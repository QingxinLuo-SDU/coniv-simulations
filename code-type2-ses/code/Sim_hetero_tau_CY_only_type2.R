
############################################################
## 0. 设置每个 R 任务使用的 CPU 线程数
## 必须放在脚本最前面
############################################################

n_threads <- 5

Sys.setenv(
  OMP_NUM_THREADS = as.character(n_threads),
  MKL_NUM_THREADS = as.character(n_threads),
  OPENBLAS_NUM_THREADS = as.character(n_threads),
  BLIS_NUM_THREADS = as.character(n_threads),
  VECLIB_MAXIMUM_THREADS = as.character(n_threads),
  NUMEXPR_NUM_THREADS = as.character(n_threads)
)

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(n_threads)
  RhpcBLASctl::omp_set_num_threads(n_threads)
}


############################################################
## 1. 数据生成
############################################################

suppressPackageStartupMessages({
  library(MASS)
})


## 1. 基础工具函数###############
## 二分类协变量生成函数：按分位数二分
Fenlei <- function(data, Clv, na.rm = TRUE) {
  
  if (Clv <= 0 || Clv >= 1) {
    stop("Clv must be between 0 and 1.")
  }
  
  cutoff <- quantile(data, probs = Clv, na.rm = na.rm, names = FALSE)
  newdata <- ifelse(data > cutoff, 1, 0)
  
  return(newdata)
}


## 安全标准化
safe_scale_vec <- function(x) {
  
  x <- as.numeric(x)
  sx <- sd(x, na.rm = TRUE)
  
  if (!is.finite(sx) || sx <= 0) {
    return(rep(0, length(x)))
  }
  
  return(as.numeric(scale(x)))
}


## 潜变量模型中根据目标发生率校准截距

calibrate_intercept <- function(eta,
                                target_prev,
                                error_sd = 1,
                                lower = -30,
                                upper = 30) {
  
  if (target_prev <= 0 || target_prev >= 1) {
    stop("target_prev must be between 0 and 1.")
  }
  
  if (error_sd <= 0) {
    stop("error_sd must be positive.")
  }
  
  eta <- as.numeric(eta)
  
  f <- function(alpha) {
    mean(pnorm((alpha + eta) / error_sd)) - target_prev
  }
  
  out <- uniroot(
    f = f,
    interval = c(lower * error_sd, upper * error_sd)
  )
  
  return(out$root)
}


## 检查 beta 长度
check_beta_length <- function(beta, expected_length, beta_name) {
  
  if (length(beta) != expected_length) {
    stop(
      paste0(
        beta_name, " length is wrong. Expected ",
        expected_length, ", but got ", length(beta), "."
      )
    )
  }
}



## 2. 固定 L 层网络结构#########


generate_layer_network <- function(parent_indices,
                                   current_layer_indices,
                                   connection_prob) {
  
  if (connection_prob < 0 || connection_prob > 1) {
    stop("connection_prob must be between 0 and 1.")
  }
  
  layer_network <- vector("list", length(current_layer_indices))
  
  for (j in seq_along(current_layer_indices)) {
    
    child_idx <- current_layer_indices[j]
    
    is_parent <- runif(length(parent_indices)) < connection_prob
    active_parent_indices <- parent_indices[is_parent]
    
    if (length(active_parent_indices) > 0) {
      
      weights <- runif(length(active_parent_indices), min = -1, max = 1)
      weights <- weights / sqrt(sum(weights^2))
      
    } else {
      
      weights <- numeric(0)
    }
    
    layer_network[[j]] <- list(
      child = child_idx,
      parents = active_parent_indices,
      weights = weights
    )
  }
  
  return(layer_network)
}


generate_layer_data_fixed_network <- function(data_matrix,
                                              layer_network,
                                              noise_sd = 0.5,
                                              standardize_output = TRUE) {
  
  if (noise_sd < 0) {
    stop("noise_sd must be non-negative.")
  }
  
  num_samples <- nrow(data_matrix)
  
  for (j in seq_along(layer_network)) {
    
    child_idx <- layer_network[[j]]$child
    active_parent_indices <- layer_network[[j]]$parents
    weights <- layer_network[[j]]$weights
    
    if (length(active_parent_indices) > 0) {
      
      parent_data <- data_matrix[, active_parent_indices, drop = FALSE]
      weighted_sum <- as.numeric(parent_data %*% weights)
      new_var <- weighted_sum + rnorm(num_samples, mean = 0, sd = noise_sd)
      
    } else {
      
      new_var <- rnorm(num_samples, mean = 0, sd = 1)
    }
    
    if (standardize_output) {
      new_var <- safe_scale_vec(new_var)
    }
    
    data_matrix[, child_idx] <- new_var
  }
  
  return(data_matrix)
}



## 3. 生成 C/U/L 协变量块########


generate_covariate_block <- function(typeC,
                                     UD = TRUE,
                                     r,
                                     N,
                                     m,
                                     k,
                                     n1,
                                     n2,
                                     n3,
                                     r_L,
                                     connection_prob,
                                     conf_L,
                                     Clv,
                                     L_network = NULL) {
  
  NC <- m + k
  total_vars <- n1 + n2 + n3
  
  ##########################################################
  ## 3.1 生成基础 C/U
  ##########################################################
  
  Sigma_CU <- matrix(rep(r, NC * NC), NC, NC)
  
  if (UD == TRUE && k > 0) {
    Sigma_CU[, (m + 1):NC] <- 0
    Sigma_CU[(m + 1):NC, ] <- 0
  }
  
  diag(Sigma_CU) <- 1
  
  mu_CU <- rep(0, NC)
  Cmat <- MASS::mvrnorm(n = N, mu = mu_CU, Sigma = Sigma_CU)
  colnames(Cmat) <- paste0("NC", 1:NC)
  
  ##########################################################
  ## 3.2 初始化 L 层数据
  ##########################################################
  
  simulated_data <- matrix(0, nrow = N, ncol = total_vars)
  colnames(simulated_data) <- c(
    paste0("L1_", 1:n1),
    paste0("L2_", 1:n2),
    paste0("L3_", 1:n3)
  )
  
  ##########################################################
  ## 3.3 生成 L1
  ##########################################################
  
  Sigma_L <- matrix(r_L, nrow = n1, ncol = n1)
  diag(Sigma_L) <- 1
  
  mu_L <- rep(0, n1)
  l1_data <- MASS::mvrnorm(n = N, mu = mu_L, Sigma = Sigma_L)
  simulated_data[, 1:n1] <- l1_data
  
  ##########################################################
  ## 3.4 固定网络生成 L2 和 L3
  ##########################################################
  
  l1_indices <- 1:n1
  l2_indices <- (n1 + 1):(n1 + n2)
  l3_indices <- (n1 + n2 + 1):total_vars
  
  if (is.null(L_network)) {
    
    L_network <- list(
      L2 = generate_layer_network(
        parent_indices = l1_indices,
        current_layer_indices = l2_indices,
        connection_prob = connection_prob
      ),
      L3 = generate_layer_network(
        parent_indices = l2_indices,
        current_layer_indices = l3_indices,
        connection_prob = connection_prob
      )
    )
  }
  
  simulated_data <- generate_layer_data_fixed_network(
    data_matrix = simulated_data,
    layer_network = L_network$L2,
    noise_sd = 0.5,
    standardize_output = TRUE
  )
  
  simulated_data <- generate_layer_data_fixed_network(
    data_matrix = simulated_data,
    layer_network = L_network$L3,
    noise_sd = 0.5,
    standardize_output = TRUE
  )
  
  df_L <- as.data.frame(simulated_data)
  
  ##########################################################
  ## 3.5 合并 C/U/L
  ##########################################################
  
  all_cov <- data.frame(Cmat, df_L)
  
  if (typeC == 2) {
    
    newc <- NULL
    
    for (ci in 1:ncol(all_cov)) {
      newc <- cbind(newc, Fenlei(all_cov[, ci], Clv))
    }
    
    all_cov <- as.data.frame(newc)
    
    all_confound <- data.frame(
      all_cov[, c(1:NC, NC + conf_L), drop = FALSE]
    )
    
  } else {
    
    all_confound <- data.frame(
      Cmat,
      df_L[, conf_L, drop = FALSE]
    )
  }
  
  colnames(all_confound) <- c(
    paste0("NC", 1:NC),
    colnames(df_L)[conf_L]
  )
  
  colnames(all_cov) <- c(
    paste0("C", 1:m),
    paste0("U", 1:k),
    colnames(df_L)
  )
  
  all_cov <- as.data.frame(all_cov)
  
  return(list(
    all_cov = all_cov,
    all_confound = all_confound,
    L_network = L_network
  ))
}



## 4. 生成单个观察性中心数据######################


Data_Generated <- function(type = 1,
                           typeC,
                           UD = TRUE,
                           r,
                           N,
                           m,
                           k,
                           n1,
                           n2,
                           n3,
                           r_L,
                           connection_prob,
                           conf_L,
                           beta_CX,
                           beta_CY,
                           beta_XY,
                           Clv,
                           L_network = NULL,
                           
                           ## 二分类潜变量模型参数
                           target_px_ob = 0.50,
                           target_py0 = 0.20,
                           conf_strength_x = 0.50,
                           conf_strength_y = 0.50,
                           error_sd_x = 1,
                           error_sd_y = 1,
                           
                           ## 是否保留潜变量，用于调试
                           keep_latent = FALSE) {
  
  NC <- m + k
  p_confound <- NC + length(conf_L)
  
  check_beta_length(beta_CX, p_confound, "beta_CX")
  check_beta_length(beta_CY, p_confound, "beta_CY")
  
  ##########################################################
  ## 4.1 生成 C/U/L
  ##########################################################
  
  cov_block <- generate_covariate_block(
    typeC = typeC,
    UD = UD,
    r = r,
    N = N,
    m = m,
    k = k,
    n1 = n1,
    n2 = n2,
    n3 = n3,
    r_L = r_L,
    connection_prob = connection_prob,
    conf_L = conf_L,
    Clv = Clv,
    L_network = L_network
  )
  
  all_cov <- cov_block$all_cov
  all_confound <- cov_block$all_confound
  
  H <- as.matrix(all_confound)
  
  ##########################################################
  ## 4.2 生成 X 和 Y
  ##########################################################
  
  if (type == 1) {
    
    ########################################################
    ## 连续 X/Y：线性结构
    ########################################################
    
    eta_x <- as.numeric(H %*% beta_CX)
    x <- eta_x + rnorm(N, mean = 0, sd = 1)
    
    eta_y <- as.numeric(H %*% beta_CY)
    y <- eta_y + beta_XY * x + rnorm(N, mean = 0, sd = 1)
    
    true_effect <- beta_XY
    
    data_total <- data.frame(x, y, all_cov)
    
    colnames(data_total)[1:(2 + NC)] <- c(
      "X",
      "Y",
      paste0("C", 1:m),
      paste0("U", 1:k)
    )
    
    diagnostics <- list(
      type = "continuous",
      true_effect = true_effect,
      mean_x = mean(x),
      mean_y = mean(y),
      sd_x = sd(x),
      sd_y = sd(y)
    )
    
  } else if (type == 2) {
    
    ########################################################
    ## 二分类 X/Y：潜变量阈值模型
    ########################################################
    
    ## X 的潜在线性部分
    eta_x_raw <- as.numeric(H %*% beta_CX)
    eta_x <- safe_scale_vec(eta_x_raw)
    eta_x <- conf_strength_x * eta_x
    
    ## 根据目标暴露率校准 alpha_x
    alpha_x <- calibrate_intercept(
      eta = eta_x,
      target_prev = target_px_ob,
      error_sd = error_sd_x
    )
    
    ## 生成观察性暴露 X
    eps_x <- rnorm(N, mean = 0, sd = error_sd_x)
    x_star <- alpha_x + eta_x + eps_x
    x <- as.integer(x_star > 0)
    
    ## Y 的潜在线性部分
    eta_y_raw <- as.numeric(H %*% beta_CY)
    eta_y <- safe_scale_vec(eta_y_raw)
    eta_y <- conf_strength_y * eta_y
    
    ## target_py0 表示 X=0 时的平均基线结局风险
    alpha_y <- calibrate_intercept(
      eta = eta_y,
      target_prev = target_py0,
      error_sd = error_sd_y
    )
    
    ## 生成观察性结局 Y
    eps_y <- rnorm(N, mean = 0, sd = error_sd_y)
    y_star <- alpha_y + eta_y + beta_XY * x + eps_y
    y <- as.integer(y_star > 0)
    
    ########################################################
    ## 计算观测风险差尺度真值
    ##
    ## 注意：
    ## beta_XY 是潜变量尺度效应；
    ## true_RD 是观测二分类结局风险差尺度效应。
    ########################################################
    
    p_y1 <- pnorm((alpha_y + eta_y + beta_XY) / error_sd_y)
    p_y0 <- pnorm((alpha_y + eta_y) / error_sd_y)
    
    true_RD <- mean(p_y1 - p_y0)
    
    data_total <- data.frame(x, y, all_cov)
    
    colnames(data_total)[1:(2 + NC)] <- c(
      "X",
      "Y",
      paste0("C", 1:m),
      paste0("U", 1:k)
    )
    
    if (keep_latent) {
      data_total$X_star <- x_star
      data_total$Y_star <- y_star
      data_total$eta_x <- eta_x
      data_total$eta_y <- eta_y
    }
    
    diagnostics <- list(
      type = "binary_latent",
      beta_XY_latent = beta_XY,
      true_RD = true_RD,
      alpha_x = alpha_x,
      alpha_y = alpha_y,
      mean_x = mean(x),
      mean_y = mean(y),
      target_px_ob = target_px_ob,
      target_py0 = target_py0,
      conf_strength_x = conf_strength_x,
      conf_strength_y = conf_strength_y,
      error_sd_x = error_sd_x,
      error_sd_y = error_sd_y
    )
    
  } else {
    
    stop("type must be 1 for continuous or 2 for binary.")
  }
  
  ##########################################################
  ## 4.3 删除不可观测 U 和可选潜变量
  ##########################################################
  
  cols_to_remove <- paste0("U", 1:k)
  
  if (!keep_latent) {
    cols_to_remove <- c(
      cols_to_remove,
      "X_star", "Y_star", "eta_x", "eta_y"
    )
  }
  
  data_total_end <- data_total[
    ,
    !names(data_total) %in% cols_to_remove,
    drop = FALSE
  ]
  
  return(list(
    data = data_total_end,
    diagnostics = diagnostics
  ))
}



## 5. 生成 RCT / do 数据######


Data_Generated_do <- function(type = 1,
                              typeC,
                              UD = TRUE,
                              r,
                              N,
                              m,
                              k,
                              n1,
                              n2,
                              n3,
                              r_L,
                              connection_prob,
                              conf_L,
                              beta_CX,
                              beta_CY,
                              beta_XY,
                              Clv,
                              p_x,
                              L_network = NULL,
                              
                              ## 二分类潜变量模型参数
                              target_py0 = 0.20,
                              conf_strength_y = 0.50,
                              error_sd_y = 1,
                              
                              ## 连续 do 暴露误差
                              error_sd_x_cont = 1,
                              
                              ## 是否保留潜变量
                              keep_latent = FALSE) {
  
  NC <- m + k
  p_confound <- NC + length(conf_L)
  
  check_beta_length(beta_CX, p_confound, "beta_CX")
  check_beta_length(beta_CY, p_confound, "beta_CY")
  
  ##########################################################
  ## 5.1 生成 C/U/L
  ##########################################################
  
  cov_block <- generate_covariate_block(
    typeC = typeC,
    UD = UD,
    r = r,
    N = N,
    m = m,
    k = k,
    n1 = n1,
    n2 = n2,
    n3 = n3,
    r_L = r_L,
    connection_prob = connection_prob,
    conf_L = conf_L,
    Clv = Clv,
    L_network = L_network
  )
  
  all_cov <- cov_block$all_cov
  all_confound <- cov_block$all_confound
  
  H <- as.matrix(all_confound)
  
  ##########################################################
  ## 5.2 处理 p_x
  ##########################################################
  
  if (length(p_x) == 1) {
    px <- rep(p_x, N)
  } else {
    px <- p_x
  }
  
  if (length(px) != N) {
    stop("p_x must be either a scalar or a vector with length N.")
  }
  
  ##########################################################
  ## 5.3 生成 X 和 Y
  ##########################################################
  
  if (type == 1) {
    
    ########################################################
    ## 连续 RCT/do：X 不由 H 决定
    ########################################################
    
    x <- px + rnorm(N, mean = 0, sd = error_sd_x_cont)
    
    eta_y <- as.numeric(H %*% beta_CY)
    y <- eta_y + beta_XY * x + rnorm(N, mean = 0, sd = 1)
    
    true_effect <- beta_XY
    
    data_total <- data.frame(x, y, all_cov)
    
    colnames(data_total)[1:(2 + NC)] <- c(
      "X",
      "Y",
      paste0("C", 1:m),
      paste0("U", 1:k)
    )
    
    diagnostics <- list(
      type = "continuous_do",
      true_effect = true_effect,
      mean_x = mean(x),
      mean_y = mean(y),
      sd_x = sd(x),
      sd_y = sd(y)
    )
    
  } else if (type == 2) {
    
    ########################################################
    ## 二分类 RCT/do：
    ## X 随机化，不由 H 决定；
    ## Y 仍由 H 的非处理通路和 X 的处理效应共同决定。
    ########################################################
    
    px <- pmin(pmax(px, 1e-6), 1 - 1e-6)
    x <- rbinom(N, size = 1, prob = px)
    
    eta_y_raw <- as.numeric(H %*% beta_CY)
    eta_y <- safe_scale_vec(eta_y_raw)
    eta_y <- conf_strength_y * eta_y
    
    alpha_y <- calibrate_intercept(
      eta = eta_y,
      target_prev = target_py0,
      error_sd = error_sd_y
    )
    
    eps_y <- rnorm(N, mean = 0, sd = error_sd_y)
    y_star <- alpha_y + eta_y + beta_XY * x + eps_y
    y <- as.integer(y_star > 0)
    
    p_y1 <- pnorm((alpha_y + eta_y + beta_XY) / error_sd_y)
    p_y0 <- pnorm((alpha_y + eta_y) / error_sd_y)
    
    true_RD <- mean(p_y1 - p_y0)
    
    data_total <- data.frame(x, y, all_cov)
    
    colnames(data_total)[1:(2 + NC)] <- c(
      "X",
      "Y",
      paste0("C", 1:m),
      paste0("U", 1:k)
    )
    
    if (keep_latent) {
      data_total$Y_star <- y_star
      data_total$eta_y <- eta_y
    }
    
    diagnostics <- list(
      type = "binary_latent_do",
      beta_XY_latent = beta_XY,
      true_RD = true_RD,
      alpha_y = alpha_y,
      mean_x = mean(x),
      mean_y = mean(y),
      p_x_mean = mean(px),
      target_py0 = target_py0,
      conf_strength_y = conf_strength_y,
      error_sd_y = error_sd_y
    )
    
  } else {
    
    stop("type must be 1 for continuous or 2 for binary.")
  }
  
  ##########################################################
  ## 5.4 删除不可观测 U 和可选潜变量
  ##########################################################
  
  cols_to_remove <- paste0("U", 1:k)
  
  if (!keep_latent) {
    cols_to_remove <- c(
      cols_to_remove,
      "X_star", "Y_star", "eta_x", "eta_y"
    )
  }
  
  data_total_end <- data_total[
    ,
    !names(data_total) %in% cols_to_remove,
    drop = FALSE
  ]
  
  return(list(
    data = data_total_end,
    diagnostics = diagnostics
  ))
}


############################################################
## Generate center-specific confounding strength
##
## Purpose:
##   Generate center-specific confounding strength for
##   heterogeneity scenarios.
##
## delta:
##   systematic shift of confounding strength
##
## tau:
##   random between-center heterogeneity
##
## rho_xy:
##   correlation between CX and CY random heterogeneity
############################################################

generate_center_confounding_strength <- function(N_center,
                                                 conf_strength_x_base = 1,
                                                 conf_strength_y_base = 1,
                                                 delta_x = 1,
                                                 delta_y = 1,
                                                 tau_x = 0,
                                                 tau_y = 0,
                                                 rho_xy = 1,
                                                 lower_bound = 0.05,
                                                 upper_bound = 5) {
  
  ##########################################################
  ## 1. Basic checks
  ##########################################################
  
  if (N_center < 1) {
    stop("N_center must be positive.")
  }
  
  if (conf_strength_x_base <= 0 || conf_strength_y_base <= 0) {
    stop("Base confounding strengths must be positive.")
  }
  
  if (delta_x <= 0 || delta_y <= 0) {
    stop("delta_x and delta_y must be positive.")
  }
  
  if (tau_x < 0 || tau_y < 0) {
    stop("tau_x and tau_y must be non-negative.")
  }
  
  if (rho_xy < -1 || rho_xy > 1) {
    stop("rho_xy must be between -1 and 1.")
  }
  
  if (lower_bound <= 0 || upper_bound <= lower_bound) {
    stop("Bounds are invalid.")
  }
  
  ##########################################################
  ## 2. No random heterogeneity
  ##########################################################
  
  if (tau_x == 0 && tau_y == 0) {
    
    mult_x <- rep(delta_x, N_center)
    mult_y <- rep(delta_y, N_center)
    
  } else {
    
    ########################################################
    ## 3. Random heterogeneity
    ##
    ## log-normal random effects:
    ##   multiplier = delta * exp(epsilon)
    ##
    ## The mean of epsilon is set to -tau^2 / 2 so that
    ## E[exp(epsilon)] is approximately 1.
    ########################################################
    
    if (tau_x > 0 && tau_y > 0) {
      
      ######################################################
      ## Special treatment for rho_xy = 1 or -1
      ## to avoid singular covariance matrix in mvrnorm.
      ######################################################
      
      if (abs(rho_xy) >= 0.999) {
        
        z <- rnorm(N_center)
        
        eps_x <- -tau_x^2 / 2 + tau_x * z
        eps_y <- -tau_y^2 / 2 + tau_y * sign(rho_xy) * z
        
      } else {
        
        Sigma <- matrix(
          c(
            tau_x^2,
            rho_xy * tau_x * tau_y,
            rho_xy * tau_x * tau_y,
            tau_y^2
          ),
          nrow = 2,
          byrow = TRUE
        )
        
        eps_mat <- MASS::mvrnorm(
          n = N_center,
          mu = c(
            -tau_x^2 / 2,
            -tau_y^2 / 2
          ),
          Sigma = Sigma
        )
        
        eps_x <- eps_mat[, 1]
        eps_y <- eps_mat[, 2]
      }
      
    } else if (tau_x > 0 && tau_y == 0) {
      
      eps_x <- rnorm(
        N_center,
        mean = -tau_x^2 / 2,
        sd = tau_x
      )
      
      eps_y <- rep(0, N_center)
      
    } else if (tau_x == 0 && tau_y > 0) {
      
      eps_x <- rep(0, N_center)
      
      eps_y <- rnorm(
        N_center,
        mean = -tau_y^2 / 2,
        sd = tau_y
      )
      
    } else {
      
      eps_x <- rep(0, N_center)
      eps_y <- rep(0, N_center)
    }
    
    mult_x <- delta_x * exp(eps_x)
    mult_y <- delta_y * exp(eps_y)
  }
  
  ##########################################################
  ## 4. Center-specific confounding strength
  ##########################################################
  
  conf_strength_x_center <- conf_strength_x_base * mult_x
  conf_strength_y_center <- conf_strength_y_base * mult_y
  
  ##########################################################
  ## 5. Bound values to avoid extreme numerical instability
  ##########################################################
  
  conf_strength_x_center <- pmin(
    pmax(conf_strength_x_center, lower_bound),
    upper_bound
  )
  
  conf_strength_y_center <- pmin(
    pmax(conf_strength_y_center, lower_bound),
    upper_bound
  )
  
  ##########################################################
  ## 6. Return
  ##########################################################
  
  return(list(
    conf_strength_x_center = conf_strength_x_center,
    conf_strength_y_center = conf_strength_y_center,
    mult_x = mult_x,
    mult_y = mult_y,
    delta_x = delta_x,
    delta_y = delta_y,
    tau_x = tau_x,
    tau_y = tau_y,
    rho_xy = rho_xy
  ))
}
############################################################
## 6. Generate all data:
##    multicenter observational data + RCT/do data
##
## Modified version:
##   1. Allows center-specific conf_strength_x
##   2. Allows center-specific conf_strength_y
############################################################

ALL_data <- function(type,
                     N_center,
                     typeC,
                     UD,
                     r,
                     Nall,
                     Nrct,
                     m,
                     k,
                     n1,
                     n2,
                     n3,
                     r_L,
                     connection_prob,
                     conf_L,
                     beta_CX,
                     beta_CY,
                     beta_XY,
                     Clv,
                     p_x,
                     
                     ## Binary latent model parameters
                     target_px_ob = 0.50,
                     target_py0 = 0.20,
                     conf_strength_x = 0.50,
                     conf_strength_y = 0.50,
                     
                     ## Center-specific parameters
                     conf_strength_x_center = NULL,
                     conf_strength_y_center = NULL,
                     conf_strength_y_rct = NULL,
                     
                     error_sd_x = 1,
                     error_sd_y = 1,
                     
                     ## Whether to keep latent variables
                     keep_latent = FALSE) {
  
  ##########################################################
  ## 0. Basic checks
  ##########################################################
  
  if (length(Nall) != N_center) {
    stop("Length of Nall must be equal to N_center.")
  }
  
  ##########################################################
  ## 0.1 Center-specific confounding strength
  ##########################################################
  
  if (is.null(conf_strength_x_center)) {
    conf_strength_x_center <- rep(conf_strength_x, N_center)
  }
  
  if (is.null(conf_strength_y_center)) {
    conf_strength_y_center <- rep(conf_strength_y, N_center)
  }
  
  if (length(conf_strength_x_center) == 1) {
    conf_strength_x_center <- rep(conf_strength_x_center, N_center)
  }
  
  if (length(conf_strength_y_center) == 1) {
    conf_strength_y_center <- rep(conf_strength_y_center, N_center)
  }
  
  if (length(conf_strength_x_center) != N_center) {
    stop("Length of conf_strength_x_center must be 1 or N_center.")
  }
  
  if (length(conf_strength_y_center) != N_center) {
    stop("Length of conf_strength_y_center must be 1 or N_center.")
  }
  
  if (any(!is.finite(conf_strength_x_center)) ||
      any(conf_strength_x_center <= 0)) {
    stop("conf_strength_x_center must contain positive finite values.")
  }
  
  if (any(!is.finite(conf_strength_y_center)) ||
      any(conf_strength_y_center <= 0)) {
    stop("conf_strength_y_center must contain positive finite values.")
  }
  
  if (is.null(conf_strength_y_rct)) {
    conf_strength_y_rct <- conf_strength_y
  }
  
  if (length(conf_strength_y_rct) != 1 ||
      !is.finite(conf_strength_y_rct) ||
      conf_strength_y_rct <= 0) {
    stop("conf_strength_y_rct must be a positive scalar.")
  }
  
  ##########################################################
  ## 1. Fix one L-layer network structure
  ##########################################################
  
  total_vars <- n1 + n2 + n3
  
  l1_indices <- 1:n1
  l2_indices <- (n1 + 1):(n1 + n2)
  l3_indices <- (n1 + n2 + 1):total_vars
  
  L_network <- list(
    L2 = generate_layer_network(
      parent_indices = l1_indices,
      current_layer_indices = l2_indices,
      connection_prob = connection_prob
    ),
    L3 = generate_layer_network(
      parent_indices = l2_indices,
      current_layer_indices = l3_indices,
      connection_prob = connection_prob
    )
  )
  
  ##########################################################
  ## 2. Generate multicenter observational data
  ##########################################################
  
  data_list <- vector("list", N_center)
  diag_list <- vector("list", N_center)
  
  for (i in 1:N_center) {
    
    gen_i <- Data_Generated(
      type = type,
      typeC = typeC,
      UD = UD,
      r = r,
      N = Nall[i],
      m = m,
      k = k,
      n1 = n1,
      n2 = n2,
      n3 = n3,
      r_L = r_L,
      connection_prob = connection_prob,
      conf_L = conf_L,
      beta_CX = beta_CX,
      beta_CY = beta_CY,
      beta_XY = beta_XY,
      Clv = Clv,
      L_network = L_network,
      target_px_ob = target_px_ob,
      target_py0 = target_py0,
      conf_strength_x = conf_strength_x_center[i],
      conf_strength_y = conf_strength_y_center[i],
      error_sd_x = error_sd_x,
      error_sd_y = error_sd_y,
      keep_latent = keep_latent
    )
    
    data_i <- gen_i$data
    data_i$centre <- i
    
    data_list[[i]] <- data_i
    
    diag_list[[i]] <- gen_i$diagnostics
    diag_list[[i]]$centre <- i
    diag_list[[i]]$N <- Nall[i]
    diag_list[[i]]$conf_strength_x_center <- conf_strength_x_center[i]
    diag_list[[i]]$conf_strength_y_center <- conf_strength_y_center[i]
  }
  
  all_data <- do.call(rbind, data_list)
  
  ##########################################################
  ## 3. Generate RCT/do data
  ##
  ## Note:
  ##   RCT/do data use conf_strength_y_rct.
  ##   This is important for delta_CY scenarios, where
  ##   observational centers can systematically deviate from
  ##   the RCT/do reference non-treatment pathway.
  ##########################################################
  
  gen_rct <- Data_Generated_do(
    type = type,
    typeC = typeC,
    UD = UD,
    r = r,
    N = Nrct,
    m = m,
    k = k,
    n1 = n1,
    n2 = n2,
    n3 = n3,
    r_L = r_L,
    connection_prob = connection_prob,
    conf_L = conf_L,
    beta_CX = beta_CX,
    beta_CY = beta_CY,
    beta_XY = beta_XY,
    Clv = Clv,
    p_x = p_x,
    L_network = L_network,
    target_py0 = target_py0,
    conf_strength_y = conf_strength_y_rct,
    error_sd_y = error_sd_y,
    keep_latent = keep_latent
  )
  
  rct_data <- gen_rct$data
  rct_diag <- gen_rct$diagnostics
  rct_diag$centre <- "RCT_do"
  rct_diag$N <- Nrct
  rct_diag$conf_strength_y_rct <- conf_strength_y_rct
  
  ##########################################################
  ## 4. Summarize truth and diagnostics
  ##
  ## Important:
  ##
  ## For binary outcome:
  ##   true_RD_ob_weighted is still calculated as
  ##   sample-size weighted average of center-specific true RD.
  ##########################################################
  
  if (type == 1) {
    
    truth <- list(
      estimand = "linear effect",
      beta_XY = beta_XY,
      true_effect = beta_XY
    )
    
  } else {
    
    center_true_RD <- sapply(diag_list, function(x) x$true_RD)
    center_N <- sapply(diag_list, function(x) x$N)
    
    weighted_true_RD_ob <- sum(center_true_RD * center_N) / sum(center_N)
    
    truth <- list(
      estimand = "observed risk difference",
      beta_XY_latent = beta_XY,
      
      true_RD_center = center_true_RD,
      true_RD_ob_weighted = weighted_true_RD_ob,
      true_RD_rct_sample = rct_diag$true_RD,
      
      target_px_ob = target_px_ob,
      target_py0 = target_py0,
      
      conf_strength_x = conf_strength_x,
      conf_strength_y = conf_strength_y,
      conf_strength_x_center = conf_strength_x_center,
      conf_strength_y_center = conf_strength_y_center,
      conf_strength_y_rct = conf_strength_y_rct,
      
      error_sd_x = error_sd_x,
      error_sd_y = error_sd_y
    )
  }
  
  diagnostics <- list(
    centers = diag_list,
    rct = rct_diag
  )
  
  ##########################################################
  ## 5. Return
  ##########################################################
  
  return(list(
    all_data = all_data,
    rct = rct_data,
    L_network = L_network,
    truth = truth,
    diagnostics = diagnostics
  ))
}




## 二分类版本：稳健线性投影变量筛选###################


## 0. 安全计算绝对相关系数##############


safe_abs_cor <- function(x, y) {
  
  cor_val <- suppressWarnings(
    cor(x, y, use = "pairwise.complete.obs")
  )
  
  if (!is.finite(cor_val) || is.na(cor_val)) {
    return(0)
  } else {
    return(abs(cor_val))
  }
}



## 1. 影响函数稳健线性回归：提取某个 term 的系数、SE、P 值##############


robust_lm_term_if <- function(y,
                              x_df,
                              term_name,
                              hc_type = "HC1",
                              min_n = 20) {
  
  ##########################################################
  ## 返回空结果
  ##########################################################
  
  empty_result <- function() {
    return(list(
      beta = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      F = NA_real_,
      IF = NULL,
      n = NA_integer_,
      valid = FALSE
    ))
  }
  
  ##########################################################
  ## 基本整理
  ##########################################################
  
  dat <- data.frame(y = y, x_df)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  
  n <- nrow(dat)
  
  if (n < min_n) {
    return(empty_result())
  }
  
  ## y 必须有变异
  if (length(unique(dat$y)) < 2) {
    return(empty_result())
  }
  
  ## 构建设计矩阵
  mm <- try(
    model.matrix(
      object = as.formula(
        paste("y ~", paste(colnames(x_df), collapse = " + "))
      ),
      data = dat
    ),
    silent = TRUE
  )
  
  if (inherits(mm, "try-error")) {
    return(empty_result())
  }
  
  if (!term_name %in% colnames(mm)) {
    return(empty_result())
  }
  
  p <- ncol(mm)
  
  if (n <= p + 2) {
    return(empty_result())
  }
  
  ## 检查秩
  qr_mm <- qr(mm)
  
  if (qr_mm$rank < p) {
    return(empty_result())
  }
  
  ##########################################################
  ## OLS 点估计
  ##########################################################
  
  beta_hat <- try(
    qr.coef(qr_mm, dat$y),
    silent = TRUE
  )
  
  if (inherits(beta_hat, "try-error")) {
    return(empty_result())
  }
  
  if (any(!is.finite(beta_hat))) {
    return(empty_result())
  }
  
  y_hat <- as.numeric(mm %*% beta_hat)
  res <- dat$y - y_hat
  
  ##########################################################
  ## 影响函数稳健方差
  ##
  ## beta_hat - beta ≈ (X'X)^(-1) sum_i x_i e_i
  ## 每个样本对 beta 的贡献：
  ## IF_i = e_i * x_i' (X'X)^(-1)
  ##########################################################
  
  XtX_inv <- try(
    solve(crossprod(mm)),
    silent = TRUE
  )
  
  if (inherits(XtX_inv, "try-error")) {
    return(empty_result())
  }
  
  IF_mat <- mm %*% XtX_inv
  IF_mat <- sweep(IF_mat, 1, res, "*")
  
  ## HC1 小样本校正
  if (hc_type == "HC1") {
    IF_mat <- IF_mat * sqrt(n / (n - p))
  }
  
  term_idx <- which(colnames(mm) == term_name)
  IF_term <- IF_mat[, term_idx]
  
  var_term <- sum(IF_term^2)
  
  if (!is.finite(var_term) || var_term <= 0) {
    return(empty_result())
  }
  
  se_term <- sqrt(var_term)
  beta_term <- beta_hat[term_idx]
  
  if (!is.finite(beta_term) || !is.finite(se_term) || se_term <= 0) {
    return(empty_result())
  }
  
  z_val <- beta_term / se_term
  p_val <- 2 * pnorm(abs(z_val), lower.tail = FALSE)
  F_val <- z_val^2
  
  return(list(
    beta = as.numeric(beta_term),
    se = as.numeric(se_term),
    z = as.numeric(z_val),
    p = as.numeric(p_val),
    F = as.numeric(F_val),
    IF = as.numeric(IF_term),
    n = n,
    valid = TRUE
  ))
}



## 2. 简单线性投影：A ~ Z##############
##    用于 X ~ Z


robust_simple_slope <- function(A,
                                Z,
                                hc_type = "HC1",
                                min_n = 20) {
  
  fit <- robust_lm_term_if(
    y = A,
    x_df = data.frame(Z = Z),
    term_name = "Z",
    hc_type = hc_type,
    min_n = min_n
  )
  
  return(fit)
}



## 3. 单中心筛选函数##############


select_binary_iv_one_center <- function(data_i,
                                        covariates_i,
                                        correlation_threshold = 0.10,
                                        selcetnum = "multi",
                                        pi_min = 0.02,
                                        F_min = 10,
                                        p_y_cutoff = 0.10,
                                        y_adjust_x = TRUE,
                                        hc_type = "HC1",
                                        min_n = 20,
                                        max_iv = Inf) {
  
  ##########################################################
  ## 空结果
  ##########################################################
  
  empty_out <- function() {
    return(list(
      selected_vars = character(0),
      score_table = data.frame()
    ))
  }
  
  if (length(covariates_i) == 0) {
    return(empty_out())
  }
  
  if (!all(c("X", "Y") %in% colnames(data_i))) {
    stop("data_i must contain X and Y.")
  }
  
  if (correlation_threshold < 0 || correlation_threshold > 1) {
    stop("correlation_threshold must be between 0 and 1.")
  }
  
  if (!(is.numeric(selcetnum) || identical(selcetnum, "multi"))) {
    stop("selcetnum must be 1, a positive integer, or 'multi'.")
  }
  
  if (is.numeric(selcetnum)) {
    if (length(selcetnum) != 1 || selcetnum < 1 || selcetnum != floor(selcetnum)) {
      stop("numeric selcetnum must be a single positive integer.")
    }
  }
  
  ##########################################################
  ## 逐个候选变量计算筛选指标
  ##########################################################
  
  score_list <- vector("list", length(covariates_i))
  
  for (j in seq_along(covariates_i)) {
    
    var_j <- covariates_i[j]
    
    if (!var_j %in% colnames(data_i)) {
      next
    }
    
    Z_j <- data_i[[var_j]]
    
    ## 变量无变异，跳过
    if (length(unique(Z_j[!is.na(Z_j)])) < 2) {
      score_list[[j]] <- data.frame(
        var = var_j,
        beta_XZ = NA_real_,
        se_XZ = NA_real_,
        F_XZ = NA_real_,
        p_XZ = NA_real_,
        beta_YZ = NA_real_,
        se_YZ = NA_real_,
        p_YZ = NA_real_,
        F_YZ = NA_real_,
        cor_XZ = NA_real_,
        pass_X = FALSE,
        pass_Y = FALSE,
        pass_all = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }
    
    ########################################################
    ## 1. 暴露关联：X ~ Z
    ########################################################
    
    fit_x <- robust_simple_slope(
      A = data_i$X,
      Z = Z_j,
      hc_type = hc_type,
      min_n = min_n
    )
    
    ########################################################
    ## 2. 结局关联：Y ~ X + Z 或 Y ~ Z
    ##
    ## 默认用 Y ~ X + Z：
    ## 目的是筛选在给定 X 后仍与 Y 相关的候选混杂工具。
    ########################################################
    
    if (y_adjust_x) {
      
      fit_y <- robust_lm_term_if(
        y = data_i$Y,
        x_df = data.frame(
          X = data_i$X,
          Z = Z_j
        ),
        term_name = "Z",
        hc_type = hc_type,
        min_n = min_n
      )
      
    } else {
      
      fit_y <- robust_simple_slope(
        A = data_i$Y,
        Z = Z_j,
        hc_type = hc_type,
        min_n = min_n
      )
    }
    
    ########################################################
    ## 3. 通过规则
    ########################################################
    
    pass_X <- isTRUE(fit_x$valid) &&
      is.finite(fit_x$beta) &&
      is.finite(fit_x$F) &&
      abs(fit_x$beta) >= pi_min &&
      fit_x$F >= F_min
    
    pass_Y <- isTRUE(fit_y$valid) &&
      is.finite(fit_y$p) &&
      fit_y$p < p_y_cutoff
    
    cor_xz <- safe_abs_cor(data_i$X, Z_j)
    
    score_list[[j]] <- data.frame(
      var = var_j,
      beta_XZ = fit_x$beta,
      se_XZ = fit_x$se,
      F_XZ = fit_x$F,
      p_XZ = fit_x$p,
      beta_YZ = fit_y$beta,
      se_YZ = fit_y$se,
      p_YZ = fit_y$p,
      F_YZ = fit_y$F,
      cor_XZ = cor_xz,
      pass_X = pass_X,
      pass_Y = pass_Y,
      pass_all = pass_X && pass_Y,
      stringsAsFactors = FALSE
    )
  }
  
  score_table <- do.call(rbind, score_list)
  
  if (is.null(score_table) || nrow(score_table) == 0) {
    return(empty_out())
  }
  
  ##########################################################
  ## 候选变量：同时通过 X 和 Y 筛选
  ##########################################################
  
  candidates <- score_table[
    score_table$pass_all == TRUE &
      is.finite(score_table$F_XZ) &
      is.finite(score_table$beta_XZ),
    ,
    drop = FALSE
  ]
  
  if (nrow(candidates) == 0) {
    return(list(
      selected_vars = character(0),
      score_table = score_table
    ))
  }
  
  ##########################################################
  ## 排序：优先 robust F_XZ，其次 |beta_XZ|
  ##########################################################
  
  candidates$abs_beta_XZ <- abs(candidates$beta_XZ)
  
  candidates <- candidates[
    order(
      candidates$F_XZ,
      candidates$abs_beta_XZ,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  ordered_vars <- candidates$var
  
  ##########################################################
  ## 去除高度相关变量
  ##
  ## 由于 candidates 已经按 F_XZ 排序，
  ## 所以遇到相关性过高时保留前面的强变量，跳过后面的弱变量。
  ##########################################################
  
  final_vars <- character(0)
  
  for (v in ordered_vars) {
    
    if (length(final_vars) == 0) {
      final_vars <- c(final_vars, v)
      next
    }
    
    cor_with_kept <- sapply(final_vars, function(v_keep) {
      safe_abs_cor(data_i[[v]], data_i[[v_keep]])
    })
    
    if (all(cor_with_kept <= correlation_threshold, na.rm = TRUE)) {
      final_vars <- c(final_vars, v)
    }
  }
  
  ##########################################################
  ## 根据 selcetnum 控制单工具或多工具
  ##
  ## 逻辑：
  ##   selcetnum = 1       ：只选最强 1 个
  ##   selcetnum = "multi" ：保留多个
  ##   selcetnum > 1       ：如果候选数 <= selcetnum，则退回最强 1 个；
  ##                         如果候选数 > selcetnum，则保留多个
  ##########################################################
  
  if (length(final_vars) > 0) {
    
    if (is.numeric(selcetnum)) {
      
      if (selcetnum == 1) {
        
        final_vars <- final_vars[1]
        
      } else if (selcetnum > 1) {
        
        if (length(final_vars) <= selcetnum) {
          final_vars <- final_vars[1]
        } else {
          final_vars <- final_vars
        }
      }
      
    } else if (identical(selcetnum, "multi")) {
      
      final_vars <- final_vars
    }
  }
  
  ##########################################################
  ## 控制最大工具数量
  ##########################################################
  
  if (is.finite(max_iv) && length(final_vars) > max_iv) {
    final_vars <- final_vars[seq_len(max_iv)]
  }
  
  return(list(
    selected_vars = final_vars,
    score_table = score_table
  ))
}



## 4. 多中心：不同中心独立筛选工具变量##############
##
## 可直接替代原 variable_selection_sorted()


variable_selection_binary_robust <- function(N_center,
                                             covariates,
                                             data,
                                             correlation_threshold = 0.10,
                                             type = 2,
                                             selcetnum = "multi",
                                             pi_min = 0.02,
                                             F_min = 10,
                                             p_y_cutoff = 0.10,
                                             y_adjust_x = TRUE,
                                             hc_type = "HC1",
                                             min_n = 20,
                                             max_iv = Inf,
                                             return_diagnostics = FALSE) {
  
  if (length(covariates) == 0) {
    stop("covariates is empty.")
  }
  
  if (!all(c("all_data", "rct") %in% names(data))) {
    stop("data must be a list containing all_data and rct.")
  }
  
  all_final_vars <- list()
  all_diagnostics <- list()
  
  for (i in 1:N_center) {
    
    data_i <- data$all_data[
      data$all_data$centre == i,
      ,
      drop = FALSE
    ]
    
    covariates_i <- covariates[
      covariates %in% colnames(data_i) &
        covariates %in% colnames(data$rct)
    ]
    
    if (length(covariates_i) == 0 || nrow(data_i) == 0) {
      
      all_final_vars[[paste0("center_", i)]] <- character(0)
      all_diagnostics[[paste0("center_", i)]] <- data.frame()
      next
    }
    
    out_i <- select_binary_iv_one_center(
      data_i = data_i,
      covariates_i = covariates_i,
      correlation_threshold = correlation_threshold,
      selcetnum = selcetnum,
      pi_min = pi_min,
      F_min = F_min,
      p_y_cutoff = p_y_cutoff,
      y_adjust_x = y_adjust_x,
      hc_type = hc_type,
      min_n = min_n,
      max_iv = max_iv
    )
    
    all_final_vars[[paste0("center_", i)]] <- out_i$selected_vars
    
    diag_i <- out_i$score_table
    
    if (nrow(diag_i) > 0) {
      diag_i$centre <- i
    }
    
    all_diagnostics[[paste0("center_", i)]] <- diag_i
  }
  
  if (return_diagnostics) {
    
    diag_all <- do.call(rbind, all_diagnostics)
    
    return(list(
      selected = all_final_vars,
      diagnostics = diag_all
    ))
    
  } else {
    
    return(all_final_vars)
  }
}



## 5. 多中心：最大中心筛选共同工具变量##############
##
## 与 variable_selection_largest_center_common() 的做法相同


variable_selection_binary_common_robust <- function(N_center,
                                                    covariates,
                                                    data,
                                                    correlation_threshold = 0.10,
                                                    type = 2,
                                                    selcetnum = "multi",
                                                    pi_min = 0.02,
                                                    F_min = 10,
                                                    p_y_cutoff = 0.10,
                                                    y_adjust_x = TRUE,
                                                    hc_type = "HC1",
                                                    min_n = 20,
                                                    max_iv = Inf,
                                                    return_diagnostics = FALSE) {
  
  if (length(covariates) == 0) {
    stop("covariates is empty.")
  }
  
  if (!all(c("all_data", "rct") %in% names(data))) {
    stop("data must be a list containing all_data and rct.")
  }
  
  ##########################################################
  ## 找到样本量最大的中心
  ##########################################################
  
  center_tab <- table(data$all_data$centre)
  
  if (length(center_tab) == 0) {
    stop("No center information found.")
  }
  
  largest_center <- as.integer(names(which.max(center_tab)))
  
  data_ref <- data$all_data[
    data$all_data$centre == largest_center,
    ,
    drop = FALSE
  ]
  
  covariates_ref <- covariates[
    covariates %in% colnames(data_ref) &
      covariates %in% colnames(data$rct)
  ]
  
  ##########################################################
  ## 在最大中心筛选
  ##########################################################
  
  out_ref <- select_binary_iv_one_center(
    data_i = data_ref,
    covariates_i = covariates_ref,
    correlation_threshold = correlation_threshold,
    selcetnum = selcetnum,
    pi_min = pi_min,
    F_min = F_min,
    p_y_cutoff = p_y_cutoff,
    y_adjust_x = y_adjust_x,
    hc_type = hc_type,
    min_n = min_n,
    max_iv = max_iv
  )
  
  common_vars <- out_ref$selected_vars
  
  ##########################################################
  ## 所有中心使用同一组工具
  ##########################################################
  
  all_final_vars <- vector("list", N_center)
  
  for (i in 1:N_center) {
    
    data_i <- data$all_data[
      data$all_data$centre == i,
      ,
      drop = FALSE
    ]
    
    common_i <- common_vars[
      common_vars %in% colnames(data_i) &
        common_vars %in% colnames(data$rct)
    ]
    
    all_final_vars[[paste0("center_", i)]] <- common_i
  }
  
  if (return_diagnostics) {
    
    diag_ref <- out_ref$score_table
    
    if (nrow(diag_ref) > 0) {
      diag_ref$screening_center <- largest_center
    }
    
    return(list(
      selected = all_final_vars,
      diagnostics = diag_ref,
      common_vars = common_vars,
      largest_center = largest_center
    ))
    
  } else {
    
    return(all_final_vars)
  }
}



############################################################
## 二分类 CON-IV 方法估计函数
##
## 适用于：
##   1. 二分类 X
##   2. 二分类 Y
##   3. 不使用 logistic
##   4. 估计观测风险差尺度效应
##   5. 方差：影响函数 robust 协方差 + delta method
##   6. 多中心合并：共享 RCT/do 协方差校正 GLS
############################################################

suppressPackageStartupMessages({
  library(MASS)
})


############################################################
## 1. 空结果函数
############################################################

empty_center_result <- function(center_id = NA,
                                reason = NA_character_) {
  
  return(list(
    valid = FALSE,
    reason = reason,
    center_id = center_id,
    
    theta = NA_real_,
    se = NA_real_,
    var = NA_real_,
    
    IV_names = character(0),
    
    pi = numeric(0),
    gamma_OB = numeric(0),
    gamma_RCT = numeric(0),
    delta = numeric(0),
    
    a_vec = numeric(0),
    b_vec = numeric(0),
    
    Sigma_pi = matrix(NA_real_, 0, 0),
    Sigma_delta = matrix(NA_real_, 0, 0),
    Sigma_delta_pi = matrix(NA_real_, 0, 0),
    
    n_OB = NA_integer_
  ))
}


############################################################
## 2. 单变量线性投影斜率 + 影响函数
##
## 估计：
##   A ~ Z
##
## 返回：
##   beta, se, p, F, IF
############################################################

slope_if_simple <- function(A,
                            Z,
                            hc_type = "HC1",
                            min_n = 20) {
  
  empty_out <- function(n_full) {
    return(list(
      beta = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      F = NA_real_,
      IF = rep(NA_real_, n_full),
      n = NA_integer_,
      valid = FALSE
    ))
  }
  
  n_full <- length(A)
  
  dat <- data.frame(A = A, Z = Z)
  ok <- complete.cases(dat)
  dat <- dat[ok, , drop = FALSE]
  
  n <- nrow(dat)
  
  if (n < min_n) {
    return(empty_out(n_full))
  }
  
  if (length(unique(dat$A)) < 2 || length(unique(dat$Z)) < 2) {
    return(empty_out(n_full))
  }
  
  Xmat <- cbind(
    "(Intercept)" = 1,
    "Z" = dat$Z
  )
  
  p_dim <- ncol(Xmat)
  
  if (qr(Xmat)$rank < p_dim) {
    return(empty_out(n_full))
  }
  
  XtX_inv <- try(
    solve(crossprod(Xmat)),
    silent = TRUE
  )
  
  if (inherits(XtX_inv, "try-error")) {
    return(empty_out(n_full))
  }
  
  beta_hat <- as.numeric(
    XtX_inv %*% crossprod(Xmat, dat$A)
  )
  
  if (any(!is.finite(beta_hat))) {
    return(empty_out(n_full))
  }
  
  fitted <- as.numeric(Xmat %*% beta_hat)
  res <- dat$A - fitted
  
  ##########################################################
  ## IF_i(beta) = e_i * x_i' (X'X)^(-1)
  ##########################################################
  
  IF_mat <- Xmat %*% XtX_inv
  IF_mat <- sweep(IF_mat, 1, res, "*")
  
  if (hc_type == "HC1") {
    IF_mat <- IF_mat * sqrt(n / (n - p_dim))
  }
  
  IF_Z <- IF_mat[, "Z"]
  
  var_beta <- sum(IF_Z^2)
  
  if (!is.finite(var_beta) || var_beta <= 0) {
    return(empty_out(n_full))
  }
  
  se_beta <- sqrt(var_beta)
  beta_Z <- beta_hat[2]
  
  z_val <- beta_Z / se_beta
  p_val <- 2 * pnorm(abs(z_val), lower.tail = FALSE)
  F_val <- z_val^2
  
  IF_full <- rep(0, n_full)
  IF_full[ok] <- IF_Z
  
  return(list(
    beta = as.numeric(beta_Z),
    se = as.numeric(se_beta),
    z = as.numeric(z_val),
    p = as.numeric(p_val),
    F = as.numeric(F_val),
    IF = as.numeric(IF_full),
    n = n,
    valid = TRUE
  ))
}


############################################################
## 3. 多个工具变量的简单斜率矩阵
##
## 对每个 Z_j 单独估计：
##   A ~ Z_j
##
## 返回：
##   beta 向量
##   IF 矩阵
##   vcov = crossprod(IF)
############################################################

slope_if_matrix <- function(A,
                            data_df,
                            iv_names,
                            hc_type = "HC1",
                            min_n = 20) {
  
  iv_names <- unique(iv_names)
  iv_names <- iv_names[
    !is.na(iv_names) &
      iv_names %in% colnames(data_df)
  ]
  
  if (length(iv_names) == 0) {
    return(list(
      beta = numeric(0),
      se = numeric(0),
      p = numeric(0),
      F = numeric(0),
      IF = matrix(NA_real_, nrow = length(A), ncol = 0),
      vcov = matrix(NA_real_, 0, 0),
      valid_names = character(0)
    ))
  }
  
  beta_vec <- numeric(0)
  se_vec <- numeric(0)
  p_vec <- numeric(0)
  F_vec <- numeric(0)
  IF_list <- list()
  
  for (iv in iv_names) {
    
    fit <- slope_if_simple(
      A = A,
      Z = data_df[[iv]],
      hc_type = hc_type,
      min_n = min_n
    )
    
    if (isTRUE(fit$valid)) {
      
      beta_vec[iv] <- fit$beta
      se_vec[iv] <- fit$se
      p_vec[iv] <- fit$p
      F_vec[iv] <- fit$F
      IF_list[[iv]] <- fit$IF
    }
  }
  
  if (length(beta_vec) == 0) {
    return(list(
      beta = numeric(0),
      se = numeric(0),
      p = numeric(0),
      F = numeric(0),
      IF = matrix(NA_real_, nrow = length(A), ncol = 0),
      vcov = matrix(NA_real_, 0, 0),
      valid_names = character(0)
    ))
  }
  
  IF_mat <- do.call(cbind, IF_list)
  colnames(IF_mat) <- names(beta_vec)
  
  vcov_mat <- crossprod(IF_mat)
  
  return(list(
    beta = beta_vec,
    se = se_vec,
    p = p_vec,
    F = F_vec,
    IF = IF_mat,
    vcov = vcov_mat,
    valid_names = names(beta_vec)
  ))
}


############################################################
## 4. RCT/do 中 Y ~ Z 的斜率协方差矩阵
##
## 所有中心共用同一个 RCT/do，因此这里的协方差
## 后续会进入多中心 V 矩阵的非对角线。
############################################################

rct_slope_vcov_if <- function(data_rct,
                              iv_names,
                              hc_type = "HC1",
                              min_n = 20) {
  
  if (!"Y" %in% colnames(data_rct)) {
    stop("data_rct must contain Y.")
  }
  
  iv_names <- unique(iv_names)
  iv_names <- iv_names[
    !is.na(iv_names) &
      iv_names %in% colnames(data_rct)
  ]
  
  if (length(iv_names) == 0) {
    stop("No valid IVs found in RCT/do data.")
  }
  
  out <- slope_if_matrix(
    A = data_rct$Y,
    data_df = data_rct,
    iv_names = iv_names,
    hc_type = hc_type,
    min_n = min_n
  )
  
  if (length(out$beta) == 0) {
    stop("No valid RCT/do slope estimates.")
  }
  
  return(list(
    beta = out$beta,
    se = out$se,
    p = out$p,
    F = out$F,
    IF = out$IF,
    vcov = out$vcov,
    valid_names = out$valid_names
  ))
}


############################################################
## 5. 权重矩阵构造函数
##
## weight_type:
##   identity    : W = I
##   diag_delta  : W = diag(1 / Var(delta_j))
##   inv_delta   : W = generalized inverse of Var(delta)
############################################################

make_weight_matrix <- function(Sigma_delta,
                               weight_type = c("diag_delta", "identity", "inv_delta"),
                               ridge = 1e-8,
                               min_var = 1e-10) {
  
  weight_type <- match.arg(weight_type)
  
  K <- nrow(Sigma_delta)
  
  if (K == 1) {
    return(matrix(1, 1, 1))
  }
  
  if (weight_type == "identity") {
    
    W <- diag(K)
    
  } else if (weight_type == "diag_delta") {
    
    v <- diag(Sigma_delta)
    v <- pmax(v, min_var)
    W <- diag(1 / v, nrow = K)
    
  } else if (weight_type == "inv_delta") {
    
    W <- MASS::ginv(
      Sigma_delta + diag(ridge, K)
    )
  }
  
  return(W)
}


############################################################
## 6. 单中心二分类 CON-IV 估计
##
## tool_mode:
##   single : 单工具 Wald ratio
##   multi  : 多工具 IVW/GMM 型估计
############################################################

estimate_binary_coniv_one_center <- function(data_OB,
                                             IV_covariates,
                                             rct_info,
                                             center_id = NA,
                                             tool_mode = c("single", "multi"),
                                             weight_type = c("diag_delta", "identity", "inv_delta"),
                                             hc_type = "HC1",
                                             min_n = 20,
                                             denom_eps = 1e-10) {
  
  tool_mode <- match.arg(tool_mode)
  weight_type <- match.arg(weight_type)
  
  n_OB <- nrow(data_OB)
  
  if (!all(c("X", "Y") %in% colnames(data_OB))) {
    return(empty_center_result(center_id, "data_OB must contain X and Y."))
  }
  
  IV_covariates <- unique(IV_covariates)
  IV_covariates <- IV_covariates[
    !is.na(IV_covariates) &
      IV_covariates %in% colnames(data_OB) &
      IV_covariates %in% names(rct_info$beta)
  ]
  
  if (length(IV_covariates) == 0) {
    return(empty_center_result(center_id, "No valid IVs for this center."))
  }
  
  if (tool_mode == "single") {
    IV_covariates <- IV_covariates[1]
  }
  
  ##########################################################
  ## 1. 观察性中心中估计：
  ##    pi       = slope(X ~ Z)
  ##    gamma_OB = slope(Y ~ Z)
  ##########################################################
  
  fit_x <- slope_if_matrix(
    A = data_OB$X,
    data_df = data_OB,
    iv_names = IV_covariates,
    hc_type = hc_type,
    min_n = min_n
  )
  
  fit_y <- slope_if_matrix(
    A = data_OB$Y,
    data_df = data_OB,
    iv_names = IV_covariates,
    hc_type = hc_type,
    min_n = min_n
  )
  
  common_iv <- IV_covariates[
    IV_covariates %in% names(fit_x$beta) &
      IV_covariates %in% names(fit_y$beta) &
      IV_covariates %in% names(rct_info$beta)
  ]
  
  if (length(common_iv) == 0) {
    return(empty_center_result(center_id, "No common valid IVs after slope estimation."))
  }
  
  if (tool_mode == "single") {
    common_iv <- common_iv[1]
  }
  
  ##########################################################
  ## 2. 提取向量和矩阵
  ##########################################################
  
  pi_vec <- fit_x$beta[common_iv]
  gamma_OB <- fit_y$beta[common_iv]
  gamma_RCT <- rct_info$beta[common_iv]
  
  delta_vec <- gamma_OB - gamma_RCT
  
  Sigma_pi <- fit_x$vcov[common_iv, common_iv, drop = FALSE]
  Sigma_g_OB <- fit_y$vcov[common_iv, common_iv, drop = FALSE]
  Sigma_g_RCT <- rct_info$vcov[common_iv, common_iv, drop = FALSE]
  
  IF_x <- fit_x$IF[, common_iv, drop = FALSE]
  IF_y <- fit_y$IF[, common_iv, drop = FALSE]
  
  ## 同一观察性中心内：
  ## Cov(delta, pi) = Cov(gamma_OB, pi)
  Sigma_delta_pi <- crossprod(IF_y, IF_x)
  
  Sigma_delta <- Sigma_g_OB + Sigma_g_RCT
  
  K <- length(common_iv)
  
  ##########################################################
  ## 3. 权重矩阵
  ##########################################################
  
  W <- make_weight_matrix(
    Sigma_delta = Sigma_delta,
    weight_type = weight_type
  )
  
  ##########################################################
  ## 4. 统一估计公式
  ##
  ## 单工具时自动退化为：
  ##   theta = delta / pi
  ##
  ## 多工具时：
  ##   theta = pi' W delta / pi' W pi
  ##########################################################
  
  pi_mat <- matrix(pi_vec, ncol = 1)
  delta_mat <- matrix(delta_vec, ncol = 1)
  
  denom <- as.numeric(t(pi_mat) %*% W %*% pi_mat)
  
  if (!is.finite(denom) || abs(denom) < denom_eps) {
    return(empty_center_result(center_id, "Denominator too small."))
  }
  
  theta_hat <- as.numeric(
    t(pi_mat) %*% W %*% delta_mat / denom
  )
  
  ##########################################################
  ## 5. Delta method
  ##
  ## a = d theta / d delta
  ## b = d theta / d pi
  ##########################################################
  
  a_vec <- as.numeric(W %*% pi_mat / denom)
  b_vec <- as.numeric(W %*% (delta_mat - 2 * theta_hat * pi_mat) / denom)
  
  names(a_vec) <- common_iv
  names(b_vec) <- common_iv
  
  var_theta <- as.numeric(
    t(a_vec) %*% Sigma_delta %*% a_vec +
      t(b_vec) %*% Sigma_pi %*% b_vec +
      2 * t(a_vec) %*% Sigma_delta_pi %*% b_vec
  )
  
  if (!is.finite(var_theta) || var_theta <= 0) {
    return(empty_center_result(center_id, "Invalid variance."))
  }
  
  se_theta <- sqrt(var_theta)
  
  return(list(
    valid = TRUE,
    reason = NA_character_,
    center_id = center_id,
    
    theta = theta_hat,
    se = se_theta,
    var = var_theta,
    
    IV_names = common_iv,
    
    pi = pi_vec,
    gamma_OB = gamma_OB,
    gamma_RCT = gamma_RCT,
    delta = delta_vec,
    
    a_vec = a_vec,
    b_vec = b_vec,
    
    Sigma_pi = Sigma_pi,
    Sigma_delta = Sigma_delta,
    Sigma_delta_pi = Sigma_delta_pi,
    
    n_OB = n_OB,
    
    weight_type = weight_type,
    tool_mode = tool_mode
  ))
}


############################################################
## 7. 多中心逐中心估计
############################################################

estimate_binary_coniv_centers <- function(data,
                                          IV_list,
                                          N_center,
                                          rct_info,
                                          tool_mode = c("single", "multi"),
                                          weight_type = c("diag_delta", "identity", "inv_delta"),
                                          hc_type = "HC1",
                                          min_n = 20) {
  
  tool_mode <- match.arg(tool_mode)
  weight_type <- match.arg(weight_type)
  
  ## 筛选结果可以是向量或 list(selected = ...)
  if (!is.null(IV_list$selected)) {
    IV_list <- IV_list$selected
  }
  
  center_results <- vector("list", N_center)
  
  for (i in 1:N_center) {
    
    data_i <- data$all_data[
      data$all_data$centre == i,
      ,
      drop = FALSE
    ]
    
    iv_i <- IV_list[[paste0("center_", i)]]
    
    center_results[[i]] <- estimate_binary_coniv_one_center(
      data_OB = data_i,
      IV_covariates = iv_i,
      rct_info = rct_info,
      center_id = i,
      tool_mode = tool_mode,
      weight_type = weight_type,
      hc_type = hc_type,
      min_n = min_n
    )
  }
  
  names(center_results) <- paste0("center_", 1:N_center)
  
  return(center_results)
}


############################################################
## 8. V 矩阵数值稳定处理
############################################################

stabilize_vcov_matrix <- function(V,
                                  ridge_start = 1e-10,
                                  max_iter = 8) {
  
  V <- as.matrix(V)
  V <- (V + t(V)) / 2
  
  ridge_added <- 0
  
  for (k in 0:max_iter) {
    
    V_try <- V + diag(ridge_added, nrow(V))
    
    ok <- try(
      chol(V_try),
      silent = TRUE
    )
    
    if (!inherits(ok, "try-error")) {
      return(list(
        V = V_try,
        ridge_added = ridge_added,
        positive_definite = TRUE
      ))
    }
    
    ridge_added <- ridge_start * 10^k
  }
  
  return(list(
    V = V + diag(ridge_added, nrow(V)),
    ridge_added = ridge_added,
    positive_definite = FALSE
  ))
}


############################################################
## 9. 多中心 GLS 合并
##
## 对角线：
##   V_cc = Var(theta_c)
##
## 非对角线：
##   V_cd = a_c' Sigma_RCT[S_c, S_d] a_d
############################################################

combine_binary_coniv_gls <- function(center_results,
                                     rct_info,
                                     ridge_start = 1e-10,
                                     use_ginv_if_needed = TRUE) {
  
  valid_idx <- which(sapply(center_results, function(x) {
    isTRUE(x$valid) &&
      is.finite(x$theta) &&
      is.finite(x$var) &&
      x$var > 0 &&
      length(x$IV_names) > 0
  }))
  
  if (length(valid_idx) == 0) {
    return(list(
      valid = FALSE,
      reason = "No valid center estimates.",
      theta_gls = NA_real_,
      se_gls = NA_real_,
      var_gls = NA_real_,
      p_value = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      V = matrix(NA_real_, 0, 0),
      center_results = center_results,
      valid_centers = integer(0)
    ))
  }
  
  valid_results <- center_results[valid_idx]
  C_valid <- length(valid_results)
  
  theta_vec <- sapply(valid_results, function(x) x$theta)
  
  V <- matrix(0, nrow = C_valid, ncol = C_valid)
  rownames(V) <- colnames(V) <- names(valid_results)
  
  ##########################################################
  ## 对角线
  ##########################################################
  
  for (i in seq_len(C_valid)) {
    V[i, i] <- valid_results[[i]]$var
  }
  
  ##########################################################
  ## 非对角线：共享 RCT/do 导致的中心间协方差
  ##########################################################
  
  if (C_valid >= 2) {
    
    for (i in 1:(C_valid - 1)) {
      for (j in (i + 1):C_valid) {
        
        iv_i <- valid_results[[i]]$IV_names
        iv_j <- valid_results[[j]]$IV_names
        
        a_i <- valid_results[[i]]$a_vec
        a_j <- valid_results[[j]]$a_vec
        
        common_ok <- all(iv_i %in% rownames(rct_info$vcov)) &&
          all(iv_j %in% colnames(rct_info$vcov))
        
        if (!common_ok) {
          
          cov_ij <- 0
          
        } else {
          
          Sigma_rct_ij <- rct_info$vcov[
            iv_i,
            iv_j,
            drop = FALSE
          ]
          
          cov_ij <- as.numeric(
            t(a_i) %*% Sigma_rct_ij %*% a_j
          )
        }
        
        V[i, j] <- cov_ij
        V[j, i] <- cov_ij
      }
    }
  }
  
  V_original <- V
  
  ##########################################################
  ## 稳定 V 矩阵
  ##########################################################
  
  stab <- stabilize_vcov_matrix(
    V = V,
    ridge_start = ridge_start
  )
  
  V_stable <- stab$V
  
  ##########################################################
  ## GLS
  ##########################################################
  
  V_inv <- try(
    solve(V_stable),
    silent = TRUE
  )
  
  if (inherits(V_inv, "try-error")) {
    
    if (use_ginv_if_needed) {
      V_inv <- MASS::ginv(V_stable)
    } else {
      return(list(
        valid = FALSE,
        reason = "V matrix cannot be inverted.",
        theta_gls = NA_real_,
        se_gls = NA_real_,
        var_gls = NA_real_,
        p_value = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        V = V_stable,
        V_original = V_original,
        center_results = center_results,
        valid_centers = valid_idx
      ))
    }
  }
  
  one <- rep(1, C_valid)
  
  denom <- as.numeric(t(one) %*% V_inv %*% one)
  
  if (!is.finite(denom) || denom <= 0) {
    return(list(
      valid = FALSE,
      reason = "Invalid GLS denominator.",
      theta_gls = NA_real_,
      se_gls = NA_real_,
      var_gls = NA_real_,
      p_value = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      V = V_stable,
      V_original = V_original,
      center_results = center_results,
      valid_centers = valid_idx
    ))
  }
  
  theta_gls <- as.numeric(
    t(one) %*% V_inv %*% theta_vec / denom
  )
  
  var_gls <- 1 / denom
  se_gls <- sqrt(var_gls)
  
  z_val <- theta_gls / se_gls
  p_val <- 2 * pnorm(abs(z_val), lower.tail = FALSE)
  
  ci_lower <- theta_gls - 1.96 * se_gls
  ci_upper <- theta_gls + 1.96 * se_gls
  
  return(list(
    valid = TRUE,
    reason = NA_character_,
    
    theta_gls = theta_gls,
    se_gls = se_gls,
    var_gls = var_gls,
    p_value = p_val,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    
    theta_centers = theta_vec,
    V = V_stable,
    V_original = V_original,
    
    ridge_added = stab$ridge_added,
    V_positive_definite = stab$positive_definite,
    
    center_results = center_results,
    valid_centers = valid_idx
  ))
}


############################################################
## 10. 已知 IV_list 时的一键运行函数
##
## 这个函数不负责变量筛选；
## 只负责：
##   1. 整理所有 IV
##   2. 构造 RCT/do 信息
##   3. 中心内估计
##   4. 多中心 GLS 合并
############################################################

run_binary_coniv_given_ivlist <- function(data,
                                          N_center,
                                          IV_list,
                                          tool_mode = c("single", "multi"),
                                          weight_type = c("diag_delta", "identity", "inv_delta"),
                                          hc_type = "HC1",
                                          min_n = 20) {
  
  tool_mode <- match.arg(tool_mode)
  weight_type <- match.arg(weight_type)
  
  if (!is.null(IV_list$selected)) {
    IV_list <- IV_list$selected
  }
  
  all_ivs <- unique(unlist(IV_list))
  all_ivs <- all_ivs[
    !is.na(all_ivs) &
      all_ivs %in% colnames(data$rct)
  ]
  
  if (length(all_ivs) == 0) {
    return(list(
      valid = FALSE,
      reason = "No IVs selected.",
      selection = IV_list
    ))
  }
  
  rct_info <- rct_slope_vcov_if(
    data_rct = data$rct,
    iv_names = all_ivs,
    hc_type = hc_type,
    min_n = min_n
  )
  
  center_results <- estimate_binary_coniv_centers(
    data = data,
    IV_list = IV_list,
    N_center = N_center,
    rct_info = rct_info,
    tool_mode = tool_mode,
    weight_type = weight_type,
    hc_type = hc_type,
    min_n = min_n
  )
  
  gls_result <- combine_binary_coniv_gls(
    center_results = center_results,
    rct_info = rct_info
  )
  
  gls_result$rct_info <- rct_info
  gls_result$selection <- IV_list
  gls_result$tool_mode <- tool_mode
  gls_result$weight_type <- weight_type
  
  return(gls_result)
}




############################################################
## 方法 1：
## Binary-CONIV-Single-Common
##
## 单工具 + 所有中心共用同一个工具
############################################################

run_binary_coniv_single_common <- function(data,
                                           N_center,
                                           covariates = NULL,
                                           correlation_threshold = 0.10,
                                           pi_min = 0.02,
                                           F_min = 10,
                                           p_y_cutoff = 0.10,
                                           y_adjust_x = TRUE,
                                           hc_type = "HC1",
                                           min_n = 20,
                                           return_selection_diagnostics = TRUE) {
  
  if (is.null(covariates)) {
    covariates <- setdiff(
      colnames(data$all_data),
      c("X", "Y", "centre")
    )
  }
  
  selection <- variable_selection_binary_common_robust(
    N_center = N_center,
    covariates = covariates,
    data = data,
    correlation_threshold = correlation_threshold,
    type = 2,
    selcetnum = 1,
    pi_min = pi_min,
    F_min = F_min,
    p_y_cutoff = p_y_cutoff,
    y_adjust_x = y_adjust_x,
    max_iv = Inf,
    hc_type = hc_type,
    min_n = min_n,
    return_diagnostics = TRUE
  )
  
  result <- run_binary_coniv_given_ivlist(
    data = data,
    N_center = N_center,
    IV_list = selection$selected,
    tool_mode = "single",
    weight_type = "identity",
    hc_type = hc_type,
    min_n = min_n
  )
  
  result$method <- "Binary-CONIV-Single-Common"
  
  if (return_selection_diagnostics) {
    result$selection_diagnostics <- selection$diagnostics
    result$common_vars <- selection$common_vars
    result$largest_center <- selection$largest_center
  }
  
  return(result)
}



############################################################
## 方法 2：
## Binary-CONIV-Single-CenterSpecific
##
## 单工具 + 每个中心单独筛选自己的工具
############################################################

run_binary_coniv_single_center_specific <- function(data,
                                                    N_center,
                                                    covariates = NULL,
                                                    correlation_threshold = 0.10,
                                                    pi_min = 0.02,
                                                    F_min = 10,
                                                    p_y_cutoff = 0.10,
                                                    y_adjust_x = TRUE,
                                                    hc_type = "HC1",
                                                    min_n = 20,
                                                    return_selection_diagnostics = TRUE) {
  
  if (is.null(covariates)) {
    covariates <- setdiff(
      colnames(data$all_data),
      c("X", "Y", "centre")
    )
  }
  
  selection <- variable_selection_binary_robust(
    N_center = N_center,
    covariates = covariates,
    data = data,
    correlation_threshold = correlation_threshold,
    type = 2,
    selcetnum = 1,
    pi_min = pi_min,
    F_min = F_min,
    p_y_cutoff = p_y_cutoff,
    y_adjust_x = y_adjust_x,
    max_iv = Inf,
    hc_type = hc_type,
    min_n = min_n,
    return_diagnostics = TRUE
  )
  
  result <- run_binary_coniv_given_ivlist(
    data = data,
    N_center = N_center,
    IV_list = selection$selected,
    tool_mode = "single",
    weight_type = "identity",
    hc_type = hc_type,
    min_n = min_n
  )
  
  result$method <- "Binary-CONIV-Single-CenterSpecific"
  
  if (return_selection_diagnostics) {
    result$selection_diagnostics <- selection$diagnostics
  }
  
  return(result)
}

############################################################
## 方法 3：
## Binary-CONIV-Multi-Common
##
## 多工具 + 所有中心共用同一组工具
############################################################

run_binary_coniv_multi_common <- function(data,
                                          N_center,
                                          covariates = NULL,
                                          correlation_threshold = 0.10,
                                          pi_min = 0.02,
                                          F_min = 10,
                                          p_y_cutoff = 0.10,
                                          y_adjust_x = TRUE,
                                          multi_iv_threshold = 5,
                                          max_iv = Inf,
                                          weight_type = "diag_delta",
                                          hc_type = "HC1",
                                          min_n = 20,
                                          return_selection_diagnostics = TRUE) {
  
  if (is.null(covariates)) {
    covariates <- setdiff(
      colnames(data$all_data),
      c("X", "Y", "centre")
    )
  }
  
  ##########################################################
  ##
  ## selcetnum = multi_iv_threshold
  ##
  ##   如果候选工具数 <= 5，自动退回最强 1 个；
  ##   如果候选工具数 >  5，保留多个工具。
  ##
  ## max_iv = Inf
  ##########################################################
  
  selection <- variable_selection_binary_common_robust(
    N_center = N_center,
    covariates = covariates,
    data = data,
    correlation_threshold = correlation_threshold,
    type = 2,
    selcetnum = multi_iv_threshold,
    pi_min = pi_min,
    F_min = F_min,
    p_y_cutoff = p_y_cutoff,
    y_adjust_x = y_adjust_x,
    max_iv = max_iv,
    hc_type = hc_type,
    min_n = min_n,
    return_diagnostics = TRUE
  )
  
  result <- run_binary_coniv_given_ivlist(
    data = data,
    N_center = N_center,
    IV_list = selection$selected,
    tool_mode = "multi",
    weight_type = weight_type,
    hc_type = hc_type,
    min_n = min_n
  )
  
  if (return_selection_diagnostics) {
    result$selection_diagnostics <- selection$diagnostics
    result$common_vars <- selection$common_vars
    result$largest_center <- selection$largest_center
  }
  
  result$selection <- selection$selected
  result$method <- "Binary-CONIV-Multi-Common"
  result$multi_iv_threshold <- multi_iv_threshold
  result$max_iv <- max_iv
  
  return(result)
}
############################################################
## 方法 4：
## Binary-CONIV-Multi-CenterSpecific
##
## 多工具 + 每个中心单独筛选自己的工具集合
############################################################

run_binary_coniv_multi_center_specific <- function(data,
                                                   N_center,
                                                   covariates = NULL,
                                                   correlation_threshold = 0.10,
                                                   pi_min = 0.02,
                                                   F_min = 10,
                                                   p_y_cutoff = 0.10,
                                                   y_adjust_x = TRUE,
                                                   multi_iv_threshold = 5,
                                                   max_iv = Inf,
                                                   weight_type = "diag_delta",
                                                   hc_type = "HC1",
                                                   min_n = 20,
                                                   return_selection_diagnostics = TRUE) {
  
  if (is.null(covariates)) {
    covariates <- setdiff(
      colnames(data$all_data),
      c("X", "Y", "centre")
    )
  }
  
  ##########################################################
  ##
  ## selcetnum = multi_iv_threshold
  ##
  ## 每个中心分别判断：
  ##   候选工具数 <= 5：该中心退回单工具
  ##   候选工具数 >  5：该中心使用多工具
  ##########################################################
  
  selection <- variable_selection_binary_robust(
    N_center = N_center,
    covariates = covariates,
    data = data,
    correlation_threshold = correlation_threshold,
    type = 2,
    selcetnum = multi_iv_threshold,
    pi_min = pi_min,
    F_min = F_min,
    p_y_cutoff = p_y_cutoff,
    y_adjust_x = y_adjust_x,
    max_iv = max_iv,
    hc_type = hc_type,
    min_n = min_n,
    return_diagnostics = TRUE
  )
  
  result <- run_binary_coniv_given_ivlist(
    data = data,
    N_center = N_center,
    IV_list = selection$selected,
    tool_mode = "multi",
    weight_type = weight_type,
    hc_type = hc_type,
    min_n = min_n
  )
  
  if (return_selection_diagnostics) {
    result$selection_diagnostics <- selection$diagnostics
  }
  
  result$selection <- selection$selected
  result$method <- "Binary-CONIV-Multi-CenterSpecific"
  result$multi_iv_threshold <- multi_iv_threshold
  result$max_iv <- max_iv
  
  return(result)
}
############################################################
## 同时运行四种二分类 CON-IV 方法
############################################################

run_binary_coniv_all_methods <- function(data,
                                         N_center,
                                         covariates = NULL,
                                         correlation_threshold = 0.10,
                                         pi_min = 0.02,
                                         F_min = 10,
                                         p_y_cutoff = 0.10,
                                         y_adjust_x = TRUE,
                                         multi_iv_threshold = 5,
                                         max_iv = Inf,
                                         weight_type_multi = "diag_delta",
                                         hc_type = "HC1",
                                         min_n = 20) {
  
  if (is.null(covariates)) {
    covariates <- setdiff(
      colnames(data$all_data),
      c("X", "Y", "centre")
    )
  }
  
  out_single_common <- run_binary_coniv_single_common(
    data = data,
    N_center = N_center,
    covariates = covariates,
    correlation_threshold = correlation_threshold,
    pi_min = pi_min,
    F_min = F_min,
    p_y_cutoff = p_y_cutoff,
    y_adjust_x = y_adjust_x,
    hc_type = hc_type,
    min_n = min_n
  )
  
  out_single_specific <- run_binary_coniv_single_center_specific(
    data = data,
    N_center = N_center,
    covariates = covariates,
    correlation_threshold = correlation_threshold,
    pi_min = pi_min,
    F_min = F_min,
    p_y_cutoff = p_y_cutoff,
    y_adjust_x = y_adjust_x,
    hc_type = hc_type,
    min_n = min_n
  )
  
  out_multi_common <- run_binary_coniv_multi_common(
    data = data,
    N_center = N_center,
    covariates = covariates,
    correlation_threshold = correlation_threshold,
    pi_min = pi_min,
    F_min = F_min,
    p_y_cutoff = p_y_cutoff,
    y_adjust_x = y_adjust_x,
    multi_iv_threshold = multi_iv_threshold,
    max_iv = max_iv,
    weight_type = weight_type_multi,
    hc_type = hc_type,
    min_n = min_n
  )
  
  out_multi_specific <- run_binary_coniv_multi_center_specific(
    data = data,
    N_center = N_center,
    covariates = covariates,
    correlation_threshold = correlation_threshold,
    pi_min = pi_min,
    F_min = F_min,
    p_y_cutoff = p_y_cutoff,
    y_adjust_x = y_adjust_x,
    multi_iv_threshold = multi_iv_threshold,
    max_iv = max_iv,
    weight_type = weight_type_multi,
    hc_type = hc_type,
    min_n = min_n
  )
  
  return(list(
    single_common = out_single_common,
    single_center_specific = out_single_specific,
    multi_common = out_multi_common,
    multi_center_specific = out_multi_specific
  ))
}

############################################################
## 提取四种方法的核心结果
############################################################

summarize_binary_coniv_results <- function(results_all,
                                           true_effect = NULL) {
  
  method_names <- names(results_all)
  
  out <- do.call(rbind, lapply(method_names, function(m) {
    
    res <- results_all[[m]]
    
    theta <- res$theta_gls
    se <- res$se_gls
    
    if (!is.null(true_effect) &&
        is.finite(theta) &&
        is.finite(se) &&
        se > 0) {
      
      bias <- theta - true_effect
      cover <- as.integer(
        res$ci_lower <= true_effect &&
          res$ci_upper >= true_effect
      )
      
    } else {
      
      bias <- NA_real_
      cover <- NA_integer_
    }
    
    data.frame(
      method = ifelse(!is.null(res$method), res$method, m),
      valid = isTRUE(res$valid),
      theta = theta,
      se = se,
      var = res$var_gls,
      p_value = res$p_value,
      ci_lower = res$ci_lower,
      ci_upper = res$ci_upper,
      bias = bias,
      cover = cover,
      n_valid_centers = length(res$valid_centers),
      ridge_added = ifelse(
        is.null(res$ridge_added),
        NA_real_,
        res$ridge_added
      ),
      stringsAsFactors = FALSE
    )
  }))
  
  rownames(out) <- NULL
  
  return(out)
}



############################################################
## 传统观察性方法变量筛选函数
##
## 筛选策略：
##
## out = "both"：
##   同时与 Y 和 X 相关，用于回归调整
##
## out = "X"：
##   只筛选与 X 相关的变量，用于 PSM 和 IPW-MSM
##
## out = "Y"：
##   只筛选与 Y 相关，一般备用
##
## type == 1：
##   连续型，使用 lm，阈值 0.05 / 候选变量数
##
## type == 2：
##   二分类，使用 logistic，阈值 0.05
############################################################

empty_covariate_list <- function(N_center) {
  out <- vector("list", N_center)
  names(out) <- paste0("center_", 1:N_center)
  for (i in 1:N_center) out[[i]] <- character(0)
  return(out)
}


Select_other <- function(data,
                         N_center,
                         type,
                         covariates,
                         out = "both") {
  
  if (!(out %in% c("both", "X", "Y"))) {
    stop("out must be 'both', 'X', or 'Y'.")
  }
  
  all_final_vars <- list()
  
  for (cen in 1:N_center) {
    
    data_i <- data$all_data[
      data$all_data$centre == cen,
      ,
      drop = FALSE
    ]
    
    covariates_i <- covariates[
      covariates %in% colnames(data_i)
    ]
    
    if (length(covariates_i) == 0) {
      all_final_vars[[paste0("center_", cen)]] <- character(0)
      next
    }
    
    p_value_Y <- rep(NA_real_, length(covariates_i))
    p_value_X <- rep(NA_real_, length(covariates_i))
    
    for (j in seq_along(covariates_i)) {
      
      var_j <- covariates_i[j]
      
      ######################################################
      ## Y 相关性模型：Y ~ X + Z_j
      ## X 相关性模型：X ~ Z_j
      ######################################################
      
      if (type == 1) {
        
        model_Y <- try(
          lm(
            reformulate(
              termlabels = c("X", var_j),
              response = "Y"
            ),
            data = data_i
          ),
          silent = TRUE
        )
        
        model_X <- try(
          lm(
            reformulate(
              termlabels = var_j,
              response = "X"
            ),
            data = data_i
          ),
          silent = TRUE
        )
        
      } else {
        
        model_Y <- try(
          glm(
            reformulate(
              termlabels = c("X", var_j),
              response = "Y"
            ),
            family = binomial(link = "logit"),
            data = data_i
          ),
          silent = TRUE
        )
        
        model_X <- try(
          glm(
            reformulate(
              termlabels = var_j,
              response = "X"
            ),
            family = binomial(link = "logit"),
            data = data_i
          ),
          silent = TRUE
        )
      }
      
      if (!inherits(model_Y, "try-error")) {
        coefs_Y <- try(summary(model_Y)$coefficients, silent = TRUE)
        
        if (!inherits(coefs_Y, "try-error") &&
            var_j %in% rownames(coefs_Y)) {
          p_value_Y[j] <- coefs_Y[var_j, 4]
        }
      }
      
      if (!inherits(model_X, "try-error")) {
        coefs_X <- try(summary(model_X)$coefficients, silent = TRUE)
        
        if (!inherits(coefs_X, "try-error") &&
            var_j %in% rownames(coefs_X)) {
          p_value_X[j] <- coefs_X[var_j, 4]
        }
      }
    }
    
    ########################################################
    ## P 值阈值：保持原策略
    ########################################################
    
    if (type == 1) {
      p_cutoff_Y <- 0.05 / length(covariates_i)
      p_cutoff_X <- 0.05 / length(covariates_i)
    } else {
      p_cutoff_Y <- 0.05
      p_cutoff_X <- 0.05
    }
    
    variable_Y <- covariates_i[
      which(!is.na(p_value_Y) & p_value_Y < p_cutoff_Y)
    ]
    
    variable_X <- covariates_i[
      which(!is.na(p_value_X) & p_value_X < p_cutoff_X)
    ]
    
    if (out == "X") {
      selected_variables <- variable_X
    } else if (out == "Y") {
      selected_variables <- variable_Y
    } else {
      selected_variables <- intersect(variable_Y, variable_X)
    }
    
    all_final_vars[[paste0("center_", cen)]] <- selected_variables
  }
  
  return(all_final_vars)
}





meta_combine_center_results <- function(center_df,
                                        method = "REML",
                                        use_t_correction = TRUE) {
  
  empty_final <- c(
    estimate = NA_real_,
    se = NA_real_,
    z_value = NA_real_,
    p_value = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    
    p_value_t = NA_real_,
    ci_lower_t = NA_real_,
    ci_upper_t = NA_real_,
    
    p_value_z = NA_real_,
    ci_lower_z = NA_real_,
    ci_upper_z = NA_real_,
    
    df_meta = NA_real_,
    
    b = NA_real_,
    p = NA_real_,
    lower = NA_real_,
    upper = NA_real_
  )
  
  if (is.null(center_df) || nrow(center_df) == 0) {
    return(list(
      final_result = empty_final,
      center_result = center_df,
      fit = NULL
    ))
  }
  
  center_df <- as.data.frame(center_df)
  
  center_df <- center_df[
    is.finite(center_df$b) &
      is.finite(center_df$se) &
      center_df$se > 0,
    ,
    drop = FALSE
  ]
  
  if (nrow(center_df) < 2) {
    return(list(
      final_result = empty_final,
      center_result = center_df,
      fit = NULL
    ))
  }
  
  fit <- try(
    metafor::rma(
      yi = b,
      sei = se,
      data = center_df,
      method = method,
      test = "knha"
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(list(
      final_result = empty_final,
      center_result = center_df,
      fit = NULL
    ))
  }
  
  estimate <- as.numeric(fit$beta)
  meta_se <- as.numeric(fit$se)
  
  if (!is.finite(estimate) || !is.finite(meta_se) || meta_se <= 0) {
    return(list(
      final_result = empty_final,
      center_result = center_df,
      fit = fit
    ))
  }
  
  z_value <- estimate / meta_se
  
  p_value_z <- 2 * pnorm(-abs(z_value))
  ci_lower_z <- estimate - 1.96 * meta_se
  ci_upper_z <- estimate + 1.96 * meta_se
  
  K <- nrow(center_df)
  df_meta <- K - 1
  
  if (use_t_correction && df_meta > 0) {
    
    crit_t <- qt(0.975, df = df_meta)
    
    p_value_t <- 2 * pt(
      abs(z_value),
      df = df_meta,
      lower.tail = FALSE
    )
    
    ci_lower_t <- estimate - crit_t * meta_se
    ci_upper_t <- estimate + crit_t * meta_se
    
    p_value <- p_value_t
    ci_lower <- ci_lower_t
    ci_upper <- ci_upper_t
    
  } else {
    
    p_value_t <- NA_real_
    ci_lower_t <- NA_real_
    ci_upper_t <- NA_real_
    
    p_value <- p_value_z
    ci_lower <- ci_lower_z
    ci_upper <- ci_upper_z
  }
  
  final_result <- c(
    estimate = estimate,
    se = meta_se,
    z_value = z_value,
    p_value = p_value,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    
    p_value_t = p_value_t,
    ci_lower_t = ci_lower_t,
    ci_upper_t = ci_upper_t,
    
    p_value_z = p_value_z,
    ci_lower_z = ci_lower_z,
    ci_upper_z = ci_upper_z,
    
    df_meta = df_meta,
    
    b = estimate,
    p = p_value,
    lower = ci_lower,
    upper = ci_upper
  )
  
  return(list(
    final_result = final_result,
    center_result = center_df,
    fit = fit,
    df_meta = df_meta,
    use_t_correction = use_t_correction
  ))
}

############################################################
## 稳健线性模型：提取 X 的 RD 和 robust SE
############################################################

robust_lm_X <- function(formula,
                        data,
                        weights = NULL,
                        min_n = 20) {
  
  empty_out <- function() {
    return(c(b = NA_real_, se = NA_real_))
  }
  
  vars <- all.vars(formula)
  
  if (!all(vars %in% colnames(data))) {
    return(empty_out())
  }
  
  dat <- data[, vars, drop = FALSE]
  
  if (is.null(weights)) {
    dat$.__w__ <- 1
  } else {
    if (length(weights) != nrow(data)) {
      return(empty_out())
    }
    dat$.__w__ <- weights
  }
  
  dat <- dat[complete.cases(dat), , drop = FALSE]
  
  if (nrow(dat) < min_n) {
    return(empty_out())
  }
  
  y <- model.response(model.frame(formula, data = dat))
  Xmat <- model.matrix(formula, data = dat)
  w <- dat$.__w__
  
  if (!"X" %in% colnames(Xmat)) {
    return(empty_out())
  }
  
  n <- length(y)
  p <- ncol(Xmat)
  
  if (n <= p + 2) {
    return(empty_out())
  }
  
  if (qr(Xmat)$rank < p) {
    return(empty_out())
  }
  
  XtWX <- crossprod(Xmat, w * Xmat)
  XtWy <- crossprod(Xmat, w * y)
  
  XtWX_inv <- try(solve(XtWX), silent = TRUE)
  
  if (inherits(XtWX_inv, "try-error")) {
    return(empty_out())
  }
  
  beta_hat <- as.numeric(XtWX_inv %*% XtWy)
  names(beta_hat) <- colnames(Xmat)
  
  fitted <- as.numeric(Xmat %*% beta_hat)
  res <- y - fitted
  
  IF_mat <- Xmat %*% XtWX_inv
  IF_mat <- sweep(IF_mat, 1, w * res, "*")
  
  ## HC1 校正
  IF_mat <- IF_mat * sqrt(n / (n - p))
  
  vcov_mat <- crossprod(IF_mat)
  
  idx <- which(colnames(Xmat) == "X")
  
  b <- beta_hat[idx]
  se <- sqrt(vcov_mat[idx, idx])
  
  if (!is.finite(b) || !is.finite(se) || se <= 0) {
    return(empty_out())
  }
  
  return(c(b = as.numeric(b), se = as.numeric(se)))
}


############################################################
## 方法一：回归调整
##
## 变量筛选：
##   covariates_both = Select_other(..., out = "both")
##
## 模型：
##   每个中心分别拟合：
##     Y ~ X + selected_covariates
##
## 效应尺度：
##   风险差 RD
############################################################

ALL_adjust_RD <- function(ALLdata,
                          N_center,
                          covariates,
                          min_n = 20) {
  
  center_result <- data.frame()
  
  for (cen in 1:N_center) {
    
    data_i <- ALLdata$all_data[
      ALLdata$all_data$centre == cen,
      ,
      drop = FALSE
    ]
    
    cov_cen <- covariates[[paste0("center_", cen)]]
    cov_cen <- cov_cen[
      cov_cen %in% colnames(data_i)
    ]
    
    ########################################################
    ## 构建模型
    ########################################################
    
    if (length(cov_cen) > 0) {
      fml <- reformulate(
        termlabels = c("X", cov_cen),
        response = "Y"
      )
    } else {
      fml <- Y ~ X
    }
    
    est <- robust_lm_X(
      formula = fml,
      data = data_i,
      weights = NULL,
      min_n = min_n
    )
    
    center_result <- rbind(
      center_result,
      data.frame(
        centre = cen,
        b = est["b"],
        se = est["se"],
        n = nrow(data_i),
        n_cov = length(cov_cen)
      )
    )
  }
  
  meta_res <- meta_combine_center_results(center_result)
  
  return(meta_res)
}


############################################################
## 计算 logit
############################################################

logit_fun <- function(p) {
  p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  return(log(p / (1 - p)))
}


############################################################
## 单中心 1:1 最近邻 PSM
##
## 不放回匹配
## caliper 默认 0.2 * SD(logit PS)
############################################################

psm_match_one_center <- function(data_i,
                                 ps,
                                 caliper_sd = 0.2) {
  
  data_i$.__ps__ <- ps
  data_i$.__logit_ps__ <- logit_fun(ps)
  data_i$.__pair_id__ <- NA_integer_
  
  treated_idx <- which(data_i$X == 1)
  control_idx <- which(data_i$X == 0)
  
  if (length(treated_idx) == 0 || length(control_idx) == 0) {
    return(data_i[0, , drop = FALSE])
  }
  
  sd_logit <- sd(data_i$.__logit_ps__, na.rm = TRUE)
  
  if (!is.finite(sd_logit) || sd_logit <= 0) {
    return(data_i[0, , drop = FALSE])
  }
  
  caliper <- caliper_sd * sd_logit
  
  available_controls <- control_idx
  matched_rows <- integer(0)
  pair_id <- 0
  
  treated_order <- treated_idx[
    order(abs(data_i$.__ps__[treated_idx] - 0.5))
  ]
  
  for (ti in treated_order) {
    
    if (length(available_controls) == 0) break
    
    dist_vec <- abs(
      data_i$.__logit_ps__[available_controls] -
        data_i$.__logit_ps__[ti]
    )
    
    best_pos <- which.min(dist_vec)
    best_dist <- dist_vec[best_pos]
    
    if (is.finite(best_dist) && best_dist <= caliper) {
      
      ci <- available_controls[best_pos]
      
      pair_id <- pair_id + 1
      
      data_i$.__pair_id__[ti] <- pair_id
      data_i$.__pair_id__[ci] <- pair_id
      
      matched_rows <- c(matched_rows, ti, ci)
      
      available_controls <- setdiff(available_controls, ci)
    }
  }
  
  if (length(matched_rows) == 0) {
    return(data_i[0, , drop = FALSE])
  }
  
  matched_data <- data_i[matched_rows, , drop = FALSE]
  
  return(matched_data)
}


############################################################
## 方法二：PSM
##
## 变量筛选：
##   covariates_X = Select_other(..., out = "X")
##
## PS 模型：
##   X ~ selected_covariates
##
## 匹配：
##   每个中心内 1:1 最近邻不放回，caliper = 0.2 SD(logit PS)
##
## 结局模型：
##   匹配样本中 Y ~ X
##
## 效应尺度：
##   风险差 RD
############################################################

ALL_PSM_RD <- function(ALLdata,
                       N_center,
                       covariates,
                       ps_trim = 0.01,
                       caliper_sd = 0.2,
                       min_n = 20) {
  
  center_result <- data.frame()
  
  for (cen in 1:N_center) {
    
    data_i <- ALLdata$all_data[
      ALLdata$all_data$centre == cen,
      ,
      drop = FALSE
    ]
    
    cov_cen <- covariates[[paste0("center_", cen)]]
    cov_cen <- cov_cen[
      cov_cen %in% colnames(data_i)
    ]
    
    ########################################################
    ## PS 模型
    ########################################################
    
    if (length(cov_cen) > 0) {
      ps_fml <- reformulate(
        termlabels = cov_cen,
        response = "X"
      )
    } else {
      ps_fml <- X ~ 1
    }
    
    vars_ps <- all.vars(ps_fml)
    dat_ps <- data_i[, unique(c("X", "Y", vars_ps)), drop = FALSE]
    dat_ps <- dat_ps[complete.cases(dat_ps), , drop = FALSE]
    
    if (nrow(dat_ps) < min_n ||
        length(unique(dat_ps$X)) < 2) {
      
      center_result <- rbind(
        center_result,
        data.frame(
          centre = cen,
          b = NA_real_,
          se = NA_real_,
          n = nrow(dat_ps),
          n_matched = 0,
          n_cov = length(cov_cen)
        )
      )
      next
    }
    
    ps_fit <- try(
      glm(
        formula = ps_fml,
        data = dat_ps,
        family = binomial(link = "logit")
      ),
      silent = TRUE
    )
    
    if (inherits(ps_fit, "try-error")) {
      center_result <- rbind(
        center_result,
        data.frame(
          centre = cen,
          b = NA_real_,
          se = NA_real_,
          n = nrow(dat_ps),
          n_matched = 0,
          n_cov = length(cov_cen)
        )
      )
      next
    }
    
    ps <- as.numeric(predict(ps_fit, type = "response"))
    ps <- pmin(pmax(ps, ps_trim), 1 - ps_trim)
    
    ########################################################
    ## 1:1 最近邻匹配
    ########################################################
    
    matched_data <- psm_match_one_center(
      data_i = dat_ps,
      ps = ps,
      caliper_sd = caliper_sd
    )
    
    if (nrow(matched_data) < min_n ||
        length(unique(matched_data$X)) < 2 ||
        length(unique(matched_data$Y)) < 2) {
      
      est <- c(b = NA_real_, se = NA_real_)
      
    } else {
      
      est <- robust_lm_X(
        formula = Y ~ X,
        data = matched_data,
        weights = NULL,
        min_n = min_n
      )
    }
    
    center_result <- rbind(
      center_result,
      data.frame(
        centre = cen,
        b = est["b"],
        se = est["se"],
        n = nrow(dat_ps),
        n_matched = nrow(matched_data),
        n_cov = length(cov_cen)
      )
    )
  }
  
  meta_res <- meta_combine_center_results(center_result)
  
  return(meta_res)
}


############################################################
## 方法三：IPW-MSM
##
## 变量筛选：
##   covariates_X = Select_other(..., out = "X")
##
## PS 模型：
##   X ~ selected_covariates
##
## 稳定化权重：
##   X=1: P(X=1) / PS
##   X=0: P(X=0) / (1-PS)
##
## MSM：
##   Y ~ X，加权稳健线性模型
##
## 效应尺度：
##   风险差 RD
############################################################

ALL_IPW_MSM_RD <- function(ALLdata,
                           N_center,
                           covariates,
                           ps_trim = 0.01,
                           weight_trunc_quantile = c(0.01, 0.99),
                           min_n = 20) {
  
  center_result <- data.frame()
  
  for (cen in 1:N_center) {
    
    data_i <- ALLdata$all_data[
      ALLdata$all_data$centre == cen,
      ,
      drop = FALSE
    ]
    
    cov_cen <- covariates[[paste0("center_", cen)]]
    cov_cen <- cov_cen[
      cov_cen %in% colnames(data_i)
    ]
    
    ########################################################
    ## PS 模型
    ########################################################
    
    if (length(cov_cen) > 0) {
      ps_fml <- reformulate(
        termlabels = cov_cen,
        response = "X"
      )
    } else {
      ps_fml <- X ~ 1
    }
    
    vars_ps <- all.vars(ps_fml)
    dat_ps <- data_i[, unique(c("X", "Y", vars_ps)), drop = FALSE]
    dat_ps <- dat_ps[complete.cases(dat_ps), , drop = FALSE]
    
    if (nrow(dat_ps) < min_n ||
        length(unique(dat_ps$X)) < 2) {
      
      center_result <- rbind(
        center_result,
        data.frame(
          centre = cen,
          b = NA_real_,
          se = NA_real_,
          n = nrow(dat_ps),
          ess = NA_real_,
          weight_mean = NA_real_,
          weight_max = NA_real_,
          n_cov = length(cov_cen)
        )
      )
      next
    }
    
    ps_fit <- try(
      glm(
        formula = ps_fml,
        data = dat_ps,
        family = binomial(link = "logit")
      ),
      silent = TRUE
    )
    
    if (inherits(ps_fit, "try-error")) {
      
      center_result <- rbind(
        center_result,
        data.frame(
          centre = cen,
          b = NA_real_,
          se = NA_real_,
          n = nrow(dat_ps),
          ess = NA_real_,
          weight_mean = NA_real_,
          weight_max = NA_real_,
          n_cov = length(cov_cen)
        )
      )
      next
    }
    
    ps <- as.numeric(predict(ps_fit, type = "response"))
    ps <- pmin(pmax(ps, ps_trim), 1 - ps_trim)
    
    ########################################################
    ## 稳定化权重
    ########################################################
    
    p_treat <- mean(dat_ps$X == 1, na.rm = TRUE)
    p_treat <- pmin(pmax(p_treat, ps_trim), 1 - ps_trim)
    
    sw <- ifelse(
      dat_ps$X == 1,
      p_treat / ps,
      (1 - p_treat) / (1 - ps)
    )
    
    ########################################################
    ## 权重截尾
    ########################################################
    
    if (!is.null(weight_trunc_quantile)) {
      q_low <- quantile(sw, weight_trunc_quantile[1], na.rm = TRUE)
      q_high <- quantile(sw, weight_trunc_quantile[2], na.rm = TRUE)
      sw <- pmin(pmax(sw, q_low), q_high)
    }
    
    ess <- sum(sw)^2 / sum(sw^2)
    
    ########################################################
    ## 加权 MSM：Y ~ X
    ########################################################
    
    est <- robust_lm_X(
      formula = Y ~ X,
      data = dat_ps,
      weights = sw,
      min_n = min_n
    )
    
    center_result <- rbind(
      center_result,
      data.frame(
        centre = cen,
        b = est["b"],
        se = est["se"],
        n = nrow(dat_ps),
        ess = ess,
        weight_mean = mean(sw),
        weight_max = max(sw),
        n_cov = length(cov_cen)
      )
    )
  }
  
  meta_res <- meta_combine_center_results(center_result)
  
  return(meta_res)
}


############################################################
## 一键运行三种传统比较方法
##
## 筛选策略：
##
## 回归调整：
##   covariates_both = Select_other(..., out = "both")
##
## PSM：
##   covariates_X = Select_other(..., out = "X")
##
## IPW-MSM：
##   covariates_X = Select_other(..., out = "X")
############################################################

run_traditional_methods_binary <- function(data,
                                           N_center,
                                           type = 2,
                                           ps_trim = 0.01,
                                           caliper_sd = 0.2,
                                           weight_trunc_quantile = c(0.01, 0.99),
                                           min_n = 20) {
  
  ##########################################################
  ## 1. 候选变量
  ##########################################################
  
  covariates_all <- setdiff(
    colnames(data$all_data),
    c("X", "Y", "centre")
  )
  
  ##########################################################
  ## 2. 传统方法变量筛选
  ##########################################################
  
  covariates_both <- try(
    Select_other(
      data = data,
      N_center = N_center,
      type = type,
      covariates = covariates_all,
      out = "both"
    ),
    silent = TRUE
  )
  
  if (inherits(covariates_both, "try-error")) {
    covariates_both <- empty_covariate_list(N_center)
  }
  
  covariates_X <- try(
    Select_other(
      data = data,
      N_center = N_center,
      type = type,
      covariates = covariates_all,
      out = "X"
    ),
    silent = TRUE
  )
  
  if (inherits(covariates_X, "try-error")) {
    covariates_X <- empty_covariate_list(N_center)
  }
  
  ##########################################################
  ## 3. 方法一：回归调整
  ##########################################################
  
  result_adjust <- try(
    ALL_adjust_RD(
      ALLdata = data,
      N_center = N_center,
      covariates = covariates_both,
      min_n = min_n
    ),
    silent = TRUE
  )
  
  if (inherits(result_adjust, "try-error")) {
    result_adjust <- list(
      final_result = NULL,
      center_result = NULL,
      fit = NULL
    )
  }
  
  ##########################################################
  ## 4. 方法二：PSM
  ##########################################################
  
  result_psm <- try(
    ALL_PSM_RD(
      ALLdata = data,
      N_center = N_center,
      covariates = covariates_X,
      ps_trim = ps_trim,
      caliper_sd = caliper_sd,
      min_n = min_n
    ),
    silent = TRUE
  )
  
  if (inherits(result_psm, "try-error")) {
    result_psm <- list(
      final_result = NULL,
      center_result = NULL,
      fit = NULL
    )
  }
  
  ##########################################################
  ## 5. 方法三：IPW-MSM
  ##########################################################
  
  result_ipw <- try(
    ALL_IPW_MSM_RD(
      ALLdata = data,
      N_center = N_center,
      covariates = covariates_X,
      ps_trim = ps_trim,
      weight_trunc_quantile = weight_trunc_quantile,
      min_n = min_n
    ),
    silent = TRUE
  )
  
  if (inherits(result_ipw, "try-error")) {
    result_ipw <- list(
      final_result = NULL,
      center_result = NULL,
      fit = NULL
    )
  }
  
  ##########################################################
  ## 6. 返回
  ##########################################################
  
  return(list(
    covariates_all = covariates_all,
    covariates_both = covariates_both,
    covariates_X = covariates_X,
    
    adjust = result_adjust,
    psm = result_psm,
    ipw_msm = result_ipw
  ))
}


############################################################
## 整理三种传统方法结果
############################################################

summarize_traditional_methods <- function(result_traditional,
                                          true_effect = NULL) {
  
  method_list <- list(
    Regression_Adjustment = result_traditional$adjust,
    Propensity_Score_Matching = result_traditional$psm,
    IPW_MSM = result_traditional$ipw_msm
  )
  
  out <- do.call(rbind, lapply(names(method_list), function(method_name) {
    
    res <- method_list[[method_name]]
    final <- extract_final_result(res)
    
    b <- final["b"]
    se <- final["se"]
    lower <- final["lower"]
    upper <- final["upper"]
    
    if (!is.null(true_effect) &&
        is.finite(b) &&
        is.finite(se)) {
      
      bias <- b - true_effect
      cover <- as.integer(lower <= true_effect & upper >= true_effect)
      
    } else {
      
      bias <- NA_real_
      cover <- NA_integer_
    }
    
    n_centers <- if (!is.null(res$center_result)) {
      sum(is.finite(res$center_result$b) &
            is.finite(res$center_result$se) &
            res$center_result$se > 0)
    } else {
      0
    }
    
    data.frame(
      method = method_name,
      b = as.numeric(b),
      se = as.numeric(se),
      p = as.numeric(final["p"]),
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      bias = as.numeric(bias),
      cover = as.integer(cover),
      n_centers = n_centers,
      stringsAsFactors = FALSE
    )
  }))
  
  rownames(out) <- NULL
  
  return(out)
}

############################################################
## 1. 安全提取 final_result
##
## 可读取的输出格式：
##   1. CON-IV GLS 输出：
##      theta_gls, se_gls, p_value, ci_lower, ci_upper
##
##   2. 完整 final_result：
##      estimate, se, z_value, p_value, ci_lower, ci_upper,
##      p_value_t, ci_lower_t, ci_upper_t,
##      p_value_z, ci_lower_z, ci_upper_z,
##      df_meta
##
##   3. 简化 final_result：
##      b, se, p, lower, upper
############################################################

extract_final_result <- function(result_obj) {
  
  empty_final <- c(
    estimate = NA_real_,
    se = NA_real_,
    z_value = NA_real_,
    p_value = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    
    p_value_t = NA_real_,
    ci_lower_t = NA_real_,
    ci_upper_t = NA_real_,
    
    p_value_z = NA_real_,
    ci_lower_z = NA_real_,
    ci_upper_z = NA_real_,
    
    df_meta = NA_real_,
    
    ## 也可读取其他列名
    b = NA_real_,
    p = NA_real_,
    lower = NA_real_,
    upper = NA_real_
  )
  
  if (inherits(result_obj, "try-error") || is.null(result_obj)) {
    return(empty_final)
  }
  
  ##########################################################
  ## 情况 1：新二分类 CON-IV GLS 输出
  ##########################################################
  
  if (is.list(result_obj) &&
      !is.null(result_obj$theta_gls)) {
    
    out <- empty_final
    
    out["estimate"] <- suppressWarnings(as.numeric(result_obj$theta_gls))
    out["se"] <- suppressWarnings(as.numeric(result_obj$se_gls))
    out["p_value"] <- suppressWarnings(as.numeric(result_obj$p_value))
    out["ci_lower"] <- suppressWarnings(as.numeric(result_obj$ci_lower))
    out["ci_upper"] <- suppressWarnings(as.numeric(result_obj$ci_upper))
    
    if (is.finite(out["estimate"]) &&
        is.finite(out["se"]) &&
        out["se"] > 0) {
      
      stat_value <- out["estimate"] / out["se"]
      
      out["z_value"] <- stat_value
      out["p_value_z"] <- 2 * pnorm(-abs(stat_value))
      out["ci_lower_z"] <- out["estimate"] - 1.96 * out["se"]
      out["ci_upper_z"] <- out["estimate"] + 1.96 * out["se"]
      
      if (!is.finite(out["p_value"])) {
        out["p_value"] <- out["p_value_z"]
      }
      
      if (!is.finite(out["ci_lower"])) {
        out["ci_lower"] <- out["ci_lower_z"]
      }
      
      if (!is.finite(out["ci_upper"])) {
        out["ci_upper"] <- out["ci_upper_z"]
      }
    }
    
    out["b"] <- out["estimate"]
    out["p"] <- out["p_value"]
    out["lower"] <- out["ci_lower"]
    out["upper"] <- out["ci_upper"]
    
    if (!is.finite(out["estimate"])) {
      return(empty_final)
    }
    
    return(out)
  }
  
  ##########################################################
  ## 情况 2：list 中包含 final_result
  ##########################################################
  
  if (is.list(result_obj) &&
      !is.null(result_obj$final_result)) {
    
    raw_out <- result_obj$final_result
    
  } else if (is.atomic(result_obj)) {
    
    raw_out <- result_obj
    
  } else {
    
    return(empty_final)
  }
  
  if (is.null(raw_out)) {
    return(empty_final)
  }
  
  ##########################################################
  ## 情况 2.1：有名字
  ##########################################################
  
  if (!is.null(names(raw_out))) {
    
    out <- empty_final
    
    raw_num <- suppressWarnings(as.numeric(raw_out))
    names(raw_num) <- names(raw_out)
    
    common_names <- intersect(
      names(out),
      names(raw_num)
    )
    
    if (length(common_names) > 0) {
      out[common_names] <- raw_num[common_names]
    }
    
    ########################################################
    ## 也可读取 b / p / lower / upper 列
    ########################################################
    
    if (!is.finite(out["estimate"]) && "b" %in% names(raw_num)) {
      out["estimate"] <- raw_num["b"]
    }
    
    if (!is.finite(out["p_value"]) && "p" %in% names(raw_num)) {
      out["p_value"] <- raw_num["p"]
    }
    
    if (!is.finite(out["ci_lower"]) && "lower" %in% names(raw_num)) {
      out["ci_lower"] <- raw_num["lower"]
    }
    
    if (!is.finite(out["ci_upper"]) && "upper" %in% names(raw_num)) {
      out["ci_upper"] <- raw_num["upper"]
    }
    
    ########################################################
    ## 补算 z 版本
    ########################################################
    
    if (is.finite(out["estimate"]) &&
        is.finite(out["se"]) &&
        out["se"] > 0) {
      
      stat_value <- out["estimate"] / out["se"]
      
      if (!is.finite(out["z_value"])) {
        out["z_value"] <- stat_value
      }
      
      if (!is.finite(out["p_value_z"])) {
        out["p_value_z"] <- 2 * pnorm(-abs(stat_value))
      }
      
      if (!is.finite(out["ci_lower_z"])) {
        out["ci_lower_z"] <- out["estimate"] - 1.96 * out["se"]
      }
      
      if (!is.finite(out["ci_upper_z"])) {
        out["ci_upper_z"] <- out["estimate"] + 1.96 * out["se"]
      }
    }
    
    if (!is.finite(out["p_value"])) {
      out["p_value"] <- out["p_value_z"]
    }
    
    if (!is.finite(out["ci_lower"])) {
      out["ci_lower"] <- out["ci_lower_z"]
    }
    
    if (!is.finite(out["ci_upper"])) {
      out["ci_upper"] <- out["ci_upper_z"]
    }
    
    out["b"] <- out["estimate"]
    out["p"] <- out["p_value"]
    out["lower"] <- out["ci_lower"]
    out["upper"] <- out["ci_upper"]
    
    if (!is.finite(out["estimate"])) {
      return(empty_final)
    }
    
    return(out)
  }
  
  ##########################################################
  ## 情况 2.2：无列名的结果
  ##
  ## 长度 >= 6：
  ##   estimate, se, z_value, p_value, ci_lower, ci_upper
  ##
  ## 长度 == 5：
  ##   estimate, se, p_value, ci_lower, ci_upper
  ##########################################################
  
  raw_num <- suppressWarnings(as.numeric(raw_out))
  
  if (length(raw_num) < 5) {
    return(empty_final)
  }
  
  out <- empty_final
  
  if (length(raw_num) >= 6) {
    
    out["estimate"] <- raw_num[1]
    out["se"] <- raw_num[2]
    out["z_value"] <- raw_num[3]
    out["p_value"] <- raw_num[4]
    out["ci_lower"] <- raw_num[5]
    out["ci_upper"] <- raw_num[6]
    
  } else {
    
    out["estimate"] <- raw_num[1]
    out["se"] <- raw_num[2]
    out["p_value"] <- raw_num[3]
    out["ci_lower"] <- raw_num[4]
    out["ci_upper"] <- raw_num[5]
  }
  
  if (is.finite(out["estimate"]) &&
      is.finite(out["se"]) &&
      out["se"] > 0) {
    
    stat_value <- out["estimate"] / out["se"]
    
    if (!is.finite(out["z_value"])) {
      out["z_value"] <- stat_value
    }
    
    out["p_value_z"] <- 2 * pnorm(-abs(stat_value))
    out["ci_lower_z"] <- out["estimate"] - 1.96 * out["se"]
    out["ci_upper_z"] <- out["estimate"] + 1.96 * out["se"]
  }
  
  out["p_value_t"] <- out["p_value"]
  out["ci_lower_t"] <- out["ci_lower"]
  out["ci_upper_t"] <- out["ci_upper"]
  
  out["b"] <- out["estimate"]
  out["p"] <- out["p_value"]
  out["lower"] <- out["ci_lower"]
  out["upper"] <- out["ci_upper"]
  
  if (!is.finite(out["estimate"])) {
    return(empty_final)
  }
  
  return(out)
}


############################################################
## 2. 安全合并 data.frame
############################################################

safe_rbind_df <- function(x, y) {
  
  if (is.null(y)) {
    return(x)
  }
  
  y <- as.data.frame(
    y,
    stringsAsFactors = FALSE
  )
  
  if (nrow(y) == 0) {
    return(x)
  }
  
  if (is.null(x)) {
    rownames(y) <- NULL
    return(y)
  }
  
  x <- as.data.frame(
    x,
    stringsAsFactors = FALSE
  )
  
  all_cols <- union(
    colnames(x),
    colnames(y)
  )
  
  for (cc in setdiff(all_cols, colnames(x))) {
    x[[cc]] <- NA
  }
  
  for (cc in setdiff(all_cols, colnames(y))) {
    y[[cc]] <- NA
  }
  
  x <- x[, all_cols, drop = FALSE]
  y <- y[, all_cols, drop = FALSE]
  
  out <- rbind(x, y)
  rownames(out) <- NULL
  
  return(out)
}


############################################################
## 3. 空变量列表
############################################################

empty_covariate_list <- function(N_center) {
  
  tmp <- vector("list", N_center)
  names(tmp) <- paste0("center_", 1:N_center)
  
  for (cen in 1:N_center) {
    tmp[[paste0("center_", cen)]] <- character(0)
  }
  
  return(tmp)
}


############################################################
## 4. 提取筛选阶段工具数量
############################################################

extract_selected_iv_count_info <- function(selected_covariates,
                                           N_center) {
  
  empty_counts <- c(
    selected_mean_iv_count = NA_real_,
    selected_median_iv_count = NA_real_,
    selected_min_iv_count = NA_real_,
    selected_max_iv_count = NA_real_
  )
  
  if (is.null(selected_covariates) ||
      inherits(selected_covariates, "try-error")) {
    return(empty_counts)
  }
  
  if (!is.null(selected_covariates$selected)) {
    selected_covariates <- selected_covariates$selected
  }
  
  iv_count <- rep(0, N_center)
  
  for (cen in 1:N_center) {
    
    tmp_iv <- selected_covariates[[paste0("center_", cen)]]
    
    if (is.null(tmp_iv)) {
      tmp_iv <- character(0)
    }
    
    iv_count[cen] <- length(tmp_iv)
  }
  
  return(c(
    selected_mean_iv_count = mean(iv_count, na.rm = TRUE),
    selected_median_iv_count = median(iv_count, na.rm = TRUE),
    selected_min_iv_count = min(iv_count, na.rm = TRUE),
    selected_max_iv_count = max(iv_count, na.rm = TRUE)
  ))
}


############################################################
## 5. 提取实际估计使用工具数量
############################################################

extract_used_iv_count_info <- function(result_obj) {
  
  empty_counts <- c(
    used_mean_iv_count = NA_real_,
    used_median_iv_count = NA_real_,
    used_min_iv_count = NA_real_,
    used_max_iv_count = NA_real_
  )
  
  if (inherits(result_obj, "try-error") || is.null(result_obj)) {
    return(empty_counts)
  }
  
  ##########################################################
  ## 新二分类 CON-IV：center_results 是 list
  ##########################################################
  
  if (!is.null(result_obj$center_results)) {
    
    iv_count <- sapply(result_obj$center_results, function(x) {
      if (isTRUE(x$valid) && !is.null(x$IV_names)) {
        return(length(x$IV_names))
      } else {
        return(NA_real_)
      }
    })
    
    iv_count <- iv_count[is.finite(iv_count)]
    
    if (length(iv_count) > 0) {
      return(c(
        used_mean_iv_count = mean(iv_count, na.rm = TRUE),
        used_median_iv_count = median(iv_count, na.rm = TRUE),
        used_min_iv_count = min(iv_count, na.rm = TRUE),
        used_max_iv_count = max(iv_count, na.rm = TRUE)
      ))
    }
  }
  
  ##########################################################
  ## center_result 中有 iv_count 时
  ##########################################################
  
  if (!is.null(result_obj$center_result)) {
    
    center_result <- as.data.frame(
      result_obj$center_result,
      stringsAsFactors = FALSE
    )
    
    if ("iv_count" %in% colnames(center_result)) {
      
      iv_count <- suppressWarnings(
        as.numeric(center_result$iv_count)
      )
      
      iv_count <- iv_count[is.finite(iv_count)]
      
      if (length(iv_count) > 0) {
        return(c(
          used_mean_iv_count = mean(iv_count, na.rm = TRUE),
          used_median_iv_count = median(iv_count, na.rm = TRUE),
          used_min_iv_count = min(iv_count, na.rm = TRUE),
          used_max_iv_count = max(iv_count, na.rm = TRUE)
        ))
      }
    }
    
    if (all(c("b", "se") %in% colnames(center_result))) {
      
      valid_center <- is.finite(center_result$b) &
        is.finite(center_result$se) &
        center_result$se > 0
      
      if (sum(valid_center, na.rm = TRUE) > 0) {
        
        used_count <- rep(1, sum(valid_center, na.rm = TRUE))
        
        return(c(
          used_mean_iv_count = mean(used_count, na.rm = TRUE),
          used_median_iv_count = median(used_count, na.rm = TRUE),
          used_min_iv_count = min(used_count, na.rm = TRUE),
          used_max_iv_count = max(used_count, na.rm = TRUE)
        ))
      }
    }
  }
  
  return(empty_counts)
}


############################################################
## 6. 提取 final_result + 工具数量
############################################################

extract_final_result_with_iv_count <- function(result_obj,
                                               selected_covariates = NULL,
                                               N_center) {
  
  final_result <- extract_final_result(result_obj)
  
  if (is.null(selected_covariates) &&
      !inherits(result_obj, "try-error") &&
      !is.null(result_obj$selection)) {
    selected_covariates <- result_obj$selection
  }
  
  selected_iv_count_info <- extract_selected_iv_count_info(
    selected_covariates = selected_covariates,
    N_center = N_center
  )
  
  used_iv_count_info <- extract_used_iv_count_info(
    result_obj = result_obj
  )
  
  return(c(
    final_result,
    selected_iv_count_info,
    used_iv_count_info
  ))
}


############################################################
## 7. 保存二分类 CON-IV 中心诊断信息
############################################################

center_results_to_df <- function(center_results,
                                 sim,
                                 method_name) {
  
  if (is.null(center_results)) {
    return(NULL)
  }
  
  out <- do.call(rbind, lapply(seq_along(center_results), function(i) {
    
    x <- center_results[[i]]
    
    data.frame(
      sim = sim,
      method = method_name,
      center = ifelse(is.null(x$center_id), i, x$center_id),
      valid = isTRUE(x$valid),
      reason = ifelse(is.null(x$reason), NA_character_, x$reason),
      theta = ifelse(is.null(x$theta), NA_real_, x$theta),
      se = ifelse(is.null(x$se), NA_real_, x$se),
      var = ifelse(is.null(x$var), NA_real_, x$var),
      n_OB = ifelse(is.null(x$n_OB), NA_integer_, x$n_OB),
      iv_count = ifelse(is.null(x$IV_names), 0, length(x$IV_names)),
      IV_names = ifelse(
        is.null(x$IV_names) || length(x$IV_names) == 0,
        NA_character_,
        paste(x$IV_names, collapse = ";")
      ),
      stringsAsFactors = FALSE
    )
  }))
  
  rownames(out) <- NULL
  
  return(out)
}


save_binary_coniv_diagnostics <- function(result_tmp,
                                          method_name,
                                          sim,
                                          all_tool_info_IV,
                                          all_center_info_IV) {
  
  if (inherits(result_tmp, "try-error") || is.null(result_tmp)) {
    return(list(
      all_tool_info_IV = all_tool_info_IV,
      all_center_info_IV = all_center_info_IV
    ))
  }
  
  ##########################################################
  ## 变量筛选诊断表
  ##########################################################
  
  if (!is.null(result_tmp$selection_diagnostics)) {
    
    tmp_tool <- as.data.frame(
      result_tmp$selection_diagnostics,
      stringsAsFactors = FALSE
    )
    
    if (nrow(tmp_tool) > 0) {
      tmp_tool$sim <- sim
      tmp_tool$method <- method_name
      
      all_tool_info_IV <- safe_rbind_df(
        all_tool_info_IV,
        tmp_tool
      )
    }
  }
  
  ##########################################################
  ## 中心估计结果
  ##########################################################
  
  if (!is.null(result_tmp$center_results)) {
    
    tmp_center <- center_results_to_df(
      center_results = result_tmp$center_results,
      sim = sim,
      method_name = method_name
    )
    
    all_center_info_IV <- safe_rbind_df(
      all_center_info_IV,
      tmp_center
    )
  }
  
  return(list(
    all_tool_info_IV = all_tool_info_IV,
    all_center_info_IV = all_center_info_IV
  ))
}

############################################################
## Summary function for causal effect simulation
##
## Used when beta_XY is not necessarily 0.
## No type I error, no power, no coverage, no rejection.
##
## Main outputs:
##   estimate_mean
##   bias
##   abs_bias
##   empirical_sd
##   mean_se
##   mse
##   rmse
##   selected / used IV counts
############################################################

resultt_FUn_binary_effect <- function(all_result,
                                      true_effect_vec) {
  
  ##########################################################
  ## safe functions
  ##########################################################
  
  safe_mean <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    mean(x)
  }
  
  safe_sd <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) <= 1) return(NA_real_)
    sd(x)
  }
  
  safe_median <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    median(x)
  }
  
  safe_min <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    min(x)
  }
  
  safe_max <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    max(x)
  }
  
  get_col_numeric <- function(df, colname) {
    if (colname %in% colnames(df)) {
      return(suppressWarnings(as.numeric(df[[colname]])))
    } else {
      return(rep(NA_real_, nrow(df)))
    }
  }
  
  ##########################################################
  ## empty output
  ##########################################################
  
  empty_out <- c(
    n_sim = NA_real_,
    n_success = NA_real_,
    success_rate = NA_real_,
    
    true_mean = NA_real_,
    true_sd = NA_real_,
    true_min = NA_real_,
    true_max = NA_real_,
    
    estimate_mean = NA_real_,
    estimate_sd = NA_real_,
    estimate_median = NA_real_,
    estimate_min = NA_real_,
    estimate_max = NA_real_,
    
    bias = NA_real_,
    abs_bias = NA_real_,
    median_bias = NA_real_,
    empirical_sd = NA_real_,
    mean_se = NA_real_,
    median_se = NA_real_,
    mse = NA_real_,
    rmse = NA_real_,
    
    selected_mean_iv_count = NA_real_,
    selected_median_iv_count = NA_real_,
    selected_min_iv_count = NA_real_,
    selected_max_iv_count = NA_real_,
    
    used_mean_iv_count = NA_real_,
    used_median_iv_count = NA_real_,
    used_min_iv_count = NA_real_,
    used_max_iv_count = NA_real_
  )
  
  if (is.null(all_result)) {
    return(empty_out)
  }
  
  df <- as.data.frame(
    all_result,
    stringsAsFactors = FALSE
  )
  
  if (nrow(df) == 0) {
    return(empty_out)
  }
  
  n_sim <- nrow(df)
  
  ##########################################################
  ## true effect vector
  ##########################################################
  
  true_effect_vec <- suppressWarnings(
    as.numeric(true_effect_vec)
  )
  
  if (length(true_effect_vec) == 1) {
    true_effect_vec <- rep(true_effect_vec, n_sim)
  }
  
  if (length(true_effect_vec) < n_sim) {
    true_effect_vec <- c(
      true_effect_vec,
      rep(NA_real_, n_sim - length(true_effect_vec))
    )
  }
  
  true_effect_vec <- true_effect_vec[seq_len(n_sim)]
  
  ##########################################################
  ## extract estimate and se
  ##
  ## Columns read:
  ##   estimate, se
  ##
  ## or:
  ##   b, se
  ##########################################################
  
  estimate <- get_col_numeric(df, "estimate")
  se <- get_col_numeric(df, "se")
  
  estimate_old <- get_col_numeric(df, "b")
  
  replace_idx <- !is.finite(estimate) & is.finite(estimate_old)
  estimate[replace_idx] <- estimate_old[replace_idx]
  
  ##########################################################
  ## valid simulations
  ##########################################################
  
  valid <- is.finite(estimate) &
    is.finite(se) &
    se > 0 &
    is.finite(true_effect_vec)
  
  n_success <- sum(valid)
  
  if (n_success == 0) {
    
    out <- empty_out
    
    out["n_sim"] <- n_sim
    out["n_success"] <- 0
    out["success_rate"] <- 0
    
    out["true_mean"] <- safe_mean(true_effect_vec)
    out["true_sd"] <- safe_sd(true_effect_vec)
    out["true_min"] <- safe_min(true_effect_vec)
    out["true_max"] <- safe_max(true_effect_vec)
    
    return(out)
  }
  
  est_v <- estimate[valid]
  se_v <- se[valid]
  true_v <- true_effect_vec[valid]
  
  bias_v <- est_v - true_v
  
  ##########################################################
  ## IV count information
  ##########################################################
  
  selected_mean_iv_count <- get_col_numeric(
    df,
    "selected_mean_iv_count"
  )
  
  selected_median_iv_count <- get_col_numeric(
    df,
    "selected_median_iv_count"
  )
  
  selected_min_iv_count <- get_col_numeric(
    df,
    "selected_min_iv_count"
  )
  
  selected_max_iv_count <- get_col_numeric(
    df,
    "selected_max_iv_count"
  )
  
  used_mean_iv_count <- get_col_numeric(
    df,
    "used_mean_iv_count"
  )
  
  used_median_iv_count <- get_col_numeric(
    df,
    "used_median_iv_count"
  )
  
  used_min_iv_count <- get_col_numeric(
    df,
    "used_min_iv_count"
  )
  
  used_max_iv_count <- get_col_numeric(
    df,
    "used_max_iv_count"
  )
  
  ##########################################################
  ## final output
  ##########################################################
  
  out <- c(
    n_sim = n_sim,
    n_success = n_success,
    success_rate = n_success / n_sim,
    
    true_mean = safe_mean(true_effect_vec),
    true_sd = safe_sd(true_effect_vec),
    true_min = safe_min(true_effect_vec),
    true_max = safe_max(true_effect_vec),
    
    estimate_mean = safe_mean(est_v),
    estimate_sd = safe_sd(est_v),
    estimate_median = safe_median(est_v),
    estimate_min = safe_min(est_v),
    estimate_max = safe_max(est_v),
    
    bias = safe_mean(bias_v),
    abs_bias = safe_mean(abs(bias_v)),
    median_bias = safe_median(bias_v),
    empirical_sd = safe_sd(est_v),
    mean_se = safe_mean(se_v),
    median_se = safe_median(se_v),
    mse = safe_mean(bias_v^2),
    rmse = sqrt(safe_mean(bias_v^2)),
    
    selected_mean_iv_count = safe_mean(selected_mean_iv_count),
    selected_median_iv_count = safe_mean(selected_median_iv_count),
    selected_min_iv_count = safe_mean(selected_min_iv_count),
    selected_max_iv_count = safe_mean(selected_max_iv_count),
    
    used_mean_iv_count = safe_mean(used_mean_iv_count),
    used_median_iv_count = safe_mean(used_median_iv_count),
    used_min_iv_count = safe_mean(used_min_iv_count),
    used_max_iv_count = safe_mean(used_max_iv_count)
  )
  
  return(out)
}


############################################################
## 9. 二分类主模拟函数
##
## 四种提出方法：
##   1. IV_single_same : 单工具 + 共用工具
##   2. IV_single_diff : 单工具 + 中心单独筛选
##   3. IV_multi_same  : 多工具 + 共用工具
##   4. IV_multi_diff  : 多工具 + 中心单独筛选
##
## 三种传统方法：
##   5. ADJ : 回归调整
##   6. PSM : 倾向得分匹配
##   7. IPW : IPW-MSM
############################################################

############################################################
## Main simulation function for binary X and binary Y
##
## Modified version:
##   1. Supports center-specific CX/CY confounding strength
##   2. Supports delta heterogeneity scenarios
##   3. Supports tau heterogeneity scenarios
############################################################

All_simulation_binary <- function(N_sim,
                                  N_center,
                                  type = 2,
                                  typeC,
                                  UD = TRUE,
                                  r,
                                  Nall,
                                  Nrct,
                                  m,
                                  k,
                                  n1,
                                  n2,
                                  n3,
                                  r_L,
                                  connection_prob,
                                  conf_L,
                                  beta_CX,
                                  beta_CY,
                                  beta_XY,
                                  Clv,
                                  p_x,
                                  
                                  ## Binary latent model parameters
                                  target_px_ob = 0.50,
                                  target_py0 = 0.20,
                                  conf_strength_x = 0.50,
                                  conf_strength_y = 0.50,
                                  
                                  ## Heterogeneity parameters
                                  hetero_delta_x = 1,
                                  hetero_delta_y = 1,
                                  hetero_tau_x = 0,
                                  hetero_tau_y = 0,
                                  hetero_rho_xy = 1,
                                  hetero_lower_bound = 0.05,
                                  hetero_upper_bound = 5,
                                  conf_strength_y_rct = NULL,
                                  
                                  error_sd_x = 1,
                                  error_sd_y = 1,
                                  
                                  ## CON-IV screening parameters
                                  correlation_threshold = 0.10,
                                  pi_min = 0.02,
                                  F_min = 10,
                                  p_y_cutoff = 0.10,
                                  y_adjust_x = TRUE,
                                  
                                  ## Multi-IV rule
                                  multi_iv_threshold = 5,
                                  max_iv = Inf,
                                  weight_type_multi = "diag_delta",
                                  
                                  ## Conventional method parameters
                                  ps_trim = 0.01,
                                  caliper_sd = 0.20,
                                  weight_trunc_quantile = c(0.01, 0.99),
                                  
                                  ## Common robust variance parameters
                                  hc_type = "HC1",
                                  min_n = 20,
                                  
                                  ## Unused parameters
                                  n_boot = NULL,
                                  trim = NULL) {
  
  ##########################################################
  ## 0. Basic checks
  ##########################################################
  
  if (type != 2) {
    stop("All_simulation_binary() is only for type = 2.")
  }
  
  if (length(Nall) != N_center) {
    stop("Length of Nall must be equal to N_center.")
  }
  
  if (multi_iv_threshold < 1) {
    stop("multi_iv_threshold must be >= 1.")
  }
  
  if (is.null(conf_strength_y_rct)) {
    conf_strength_y_rct <- conf_strength_y
  }
  
  ##########################################################
  ## 1. Storage objects
  ##########################################################
  
  all_result_IV_single_same <- NULL
  all_result_IV_single_diff <- NULL
  all_result_IV_multi_same  <- NULL
  all_result_IV_multi_diff  <- NULL
  
  all_result_ADJ <- NULL
  all_result_PSM <- NULL
  all_result_IPW <- NULL
  
  all_tool_info_IV <- NULL
  all_center_info_IV <- NULL
  all_sim_status <- NULL
  
  true_effect_vec <- rep(NA_real_, N_sim)
  
  ##########################################################
  ## 2. Local helper functions
  ##########################################################
  
  is_ok <- function(x) {
    
    if (is.null(x)) {
      return(FALSE)
    }
    
    if (!"estimate" %in% names(x)) {
      return(FALSE)
    }
    
    est <- suppressWarnings(as.numeric(x["estimate"]))
    se <- suppressWarnings(as.numeric(x["se"]))
    
    return(
      is.finite(est) &&
        is.finite(se) &&
        se > 0
    )
  }
  
  get_iv_final <- function(result_obj) {
    
    if (inherits(result_obj, "try-error") || is.null(result_obj)) {
      return(
        extract_final_result_with_iv_count(
          result_obj = NULL,
          selected_covariates = NULL,
          N_center = N_center
        )
      )
    }
    
    selected_covariates <- NULL
    
    if (!is.null(result_obj$selection)) {
      selected_covariates <- result_obj$selection
    }
    
    out <- extract_final_result_with_iv_count(
      result_obj = result_obj,
      selected_covariates = selected_covariates,
      N_center = N_center
    )
    
    return(out)
  }
  
  get_standard_final <- function(result_obj) {
    
    if (inherits(result_obj, "try-error") || is.null(result_obj)) {
      return(extract_final_result(NULL))
    }
    
    return(extract_final_result(result_obj))
  }
  
  ##########################################################
  ## 3. Simulation loop
  ##########################################################
  
  for (sim in seq_len(N_sim)) {
    
    cat("sim =", sim, "\n")
    
    ########################################################
    ## 3.1 Generate center-specific confounding strength
    ########################################################
    
    center_strength <- try(
      generate_center_confounding_strength(
        N_center = N_center,
        conf_strength_x_base = conf_strength_x,
        conf_strength_y_base = conf_strength_y,
        delta_x = hetero_delta_x,
        delta_y = hetero_delta_y,
        tau_x = hetero_tau_x,
        tau_y = hetero_tau_y,
        rho_xy = hetero_rho_xy,
        lower_bound = hetero_lower_bound,
        upper_bound = hetero_upper_bound
      ),
      silent = TRUE
    )
    
    if (inherits(center_strength, "try-error")) {
      
      empty_IV <- extract_final_result_with_iv_count(
        result_obj = NULL,
        selected_covariates = NULL,
        N_center = N_center
      )
      
      empty_standard <- extract_final_result(NULL)
      
      all_result_IV_single_same <- rbind(all_result_IV_single_same, empty_IV)
      all_result_IV_single_diff <- rbind(all_result_IV_single_diff, empty_IV)
      all_result_IV_multi_same  <- rbind(all_result_IV_multi_same,  empty_IV)
      all_result_IV_multi_diff  <- rbind(all_result_IV_multi_diff,  empty_IV)
      
      all_result_ADJ <- rbind(all_result_ADJ, empty_standard)
      all_result_PSM <- rbind(all_result_PSM, empty_standard)
      all_result_IPW <- rbind(all_result_IPW, empty_standard)
      
      all_sim_status <- rbind(
        all_sim_status,
        data.frame(
          sim = sim,
          data_success = FALSE,
          true_effect = NA_real_,
          
          IV_single_same_success = FALSE,
          IV_single_diff_success = FALSE,
          IV_multi_same_success = FALSE,
          IV_multi_diff_success = FALSE,
          
          ADJ_success = FALSE,
          PSM_success = FALSE,
          IPW_success = FALSE,
          
          conf_strength_x_center_mean = NA_real_,
          conf_strength_x_center_min = NA_real_,
          conf_strength_x_center_max = NA_real_,
          
          conf_strength_y_center_mean = NA_real_,
          conf_strength_y_center_min = NA_real_,
          conf_strength_y_center_max = NA_real_,
          
          stringsAsFactors = FALSE
        )
      )
      
      next
    }
    
    conf_strength_x_center <- center_strength$conf_strength_x_center
    conf_strength_y_center <- center_strength$conf_strength_y_center
    
    ########################################################
    ## 3.2 Generate data
    ########################################################
    
    data <- try(
      ALL_data(
        type = 2,
        N_center = N_center,
        typeC = typeC,
        UD = UD,
        r = r,
        Nall = Nall,
        Nrct = Nrct,
        m = m,
        k = k,
        n1 = n1,
        n2 = n2,
        n3 = n3,
        r_L = r_L,
        connection_prob = connection_prob,
        conf_L = conf_L,
        beta_CX = beta_CX,
        beta_CY = beta_CY,
        beta_XY = beta_XY,
        Clv = Clv,
        p_x = p_x,
        
        target_px_ob = target_px_ob,
        target_py0 = target_py0,
        conf_strength_x = conf_strength_x,
        conf_strength_y = conf_strength_y,
        
        conf_strength_x_center = conf_strength_x_center,
        conf_strength_y_center = conf_strength_y_center,
        conf_strength_y_rct = conf_strength_y_rct,
        
        error_sd_x = error_sd_x,
        error_sd_y = error_sd_y,
        keep_latent = FALSE
      ),
      silent = TRUE
    )
    
    ########################################################
    ## 3.3 If data generation failed
    ########################################################
    
    if (inherits(data, "try-error")) {
      
      empty_IV <- extract_final_result_with_iv_count(
        result_obj = NULL,
        selected_covariates = NULL,
        N_center = N_center
      )
      
      empty_standard <- extract_final_result(NULL)
      
      all_result_IV_single_same <- rbind(all_result_IV_single_same, empty_IV)
      all_result_IV_single_diff <- rbind(all_result_IV_single_diff, empty_IV)
      all_result_IV_multi_same  <- rbind(all_result_IV_multi_same,  empty_IV)
      all_result_IV_multi_diff  <- rbind(all_result_IV_multi_diff,  empty_IV)
      
      all_result_ADJ <- rbind(all_result_ADJ, empty_standard)
      all_result_PSM <- rbind(all_result_PSM, empty_standard)
      all_result_IPW <- rbind(all_result_IPW, empty_standard)
      
      all_sim_status <- rbind(
        all_sim_status,
        data.frame(
          sim = sim,
          data_success = FALSE,
          true_effect = NA_real_,
          
          IV_single_same_success = FALSE,
          IV_single_diff_success = FALSE,
          IV_multi_same_success = FALSE,
          IV_multi_diff_success = FALSE,
          
          ADJ_success = FALSE,
          PSM_success = FALSE,
          IPW_success = FALSE,
          
          conf_strength_x_center_mean = mean(conf_strength_x_center),
          conf_strength_x_center_min = min(conf_strength_x_center),
          conf_strength_x_center_max = max(conf_strength_x_center),
          
          conf_strength_y_center_mean = mean(conf_strength_y_center),
          conf_strength_y_center_min = min(conf_strength_y_center),
          conf_strength_y_center_max = max(conf_strength_y_center),
          
          stringsAsFactors = FALSE
        )
      )
      
      next
    }
    
    ########################################################
    ## 3.4 True effect for this simulation
    ##
    ## Important:
    ##   The primary truth is still true_RD_ob_weighted.
    ########################################################
    
    if (!is.null(data$truth$true_RD_ob_weighted)) {
      true_effect <- data$truth$true_RD_ob_weighted
    } else if (!is.null(data$truth$true_effect)) {
      true_effect <- data$truth$true_effect
    } else {
      true_effect <- NA_real_
    }
    
    true_effect_vec[sim] <- true_effect
    
    ########################################################
    ## 3.5 Candidate covariates
    ########################################################
    
    covariates_all <- setdiff(
      colnames(data$all_data),
      c("X", "Y", "centre")
    )
    
    ########################################################
    ## 3.6 Binary-CONIV: single IV, common IV
    ########################################################
    
    result_IV_single_same_q <- try(
      run_binary_coniv_single_common(
        data = data,
        N_center = N_center,
        covariates = covariates_all,
        correlation_threshold = correlation_threshold,
        pi_min = pi_min,
        F_min = F_min,
        p_y_cutoff = p_y_cutoff,
        y_adjust_x = y_adjust_x,
        hc_type = hc_type,
        min_n = min_n,
        return_selection_diagnostics = TRUE
      ),
      silent = TRUE
    )
    
    final_IV_single_same <- get_iv_final(result_IV_single_same_q)
    
    all_result_IV_single_same <- rbind(
      all_result_IV_single_same,
      final_IV_single_same
    )
    
    ########################################################
    ## 3.7 Binary-CONIV: single IV, center-specific IV
    ########################################################
    
    result_IV_single_diff_q <- try(
      run_binary_coniv_single_center_specific(
        data = data,
        N_center = N_center,
        covariates = covariates_all,
        correlation_threshold = correlation_threshold,
        pi_min = pi_min,
        F_min = F_min,
        p_y_cutoff = p_y_cutoff,
        y_adjust_x = y_adjust_x,
        hc_type = hc_type,
        min_n = min_n,
        return_selection_diagnostics = TRUE
      ),
      silent = TRUE
    )
    
    final_IV_single_diff <- get_iv_final(result_IV_single_diff_q)
    
    all_result_IV_single_diff <- rbind(
      all_result_IV_single_diff,
      final_IV_single_diff
    )
    
    ########################################################
    ## 3.8 Binary-CONIV: multi IV, common IV
    ########################################################
    
    result_IV_multi_same_q <- try(
      run_binary_coniv_multi_common(
        data = data,
        N_center = N_center,
        covariates = covariates_all,
        correlation_threshold = correlation_threshold,
        pi_min = pi_min,
        F_min = F_min,
        p_y_cutoff = p_y_cutoff,
        y_adjust_x = y_adjust_x,
        multi_iv_threshold = multi_iv_threshold,
        max_iv = max_iv,
        weight_type = weight_type_multi,
        hc_type = hc_type,
        min_n = min_n,
        return_selection_diagnostics = TRUE
      ),
      silent = TRUE
    )
    
    final_IV_multi_same <- get_iv_final(result_IV_multi_same_q)
    
    all_result_IV_multi_same <- rbind(
      all_result_IV_multi_same,
      final_IV_multi_same
    )
    
    ########################################################
    ## 3.9 Binary-CONIV: multi IV, center-specific IV
    ########################################################
    
    result_IV_multi_diff_q <- try(
      run_binary_coniv_multi_center_specific(
        data = data,
        N_center = N_center,
        covariates = covariates_all,
        correlation_threshold = correlation_threshold,
        pi_min = pi_min,
        F_min = F_min,
        p_y_cutoff = p_y_cutoff,
        y_adjust_x = y_adjust_x,
        multi_iv_threshold = multi_iv_threshold,
        max_iv = max_iv,
        weight_type = weight_type_multi,
        hc_type = hc_type,
        min_n = min_n,
        return_selection_diagnostics = TRUE
      ),
      silent = TRUE
    )
    
    final_IV_multi_diff <- get_iv_final(result_IV_multi_diff_q)
    
    all_result_IV_multi_diff <- rbind(
      all_result_IV_multi_diff,
      final_IV_multi_diff
    )
    
    ########################################################
    ## 3.10 Save CON-IV diagnostic information
    ########################################################
    
    method_result_list <- list(
      IV_single_same = result_IV_single_same_q,
      IV_single_diff = result_IV_single_diff_q,
      IV_multi_same  = result_IV_multi_same_q,
      IV_multi_diff  = result_IV_multi_diff_q
    )
    
    for (method_name in names(method_result_list)) {
      
      result_tmp <- method_result_list[[method_name]]
      
      if (!inherits(result_tmp, "try-error") && !is.null(result_tmp)) {
        
        tmp_save <- save_binary_coniv_diagnostics(
          result_tmp = result_tmp,
          method_name = method_name,
          sim = sim,
          all_tool_info_IV = all_tool_info_IV,
          all_center_info_IV = all_center_info_IV
        )
        
        all_tool_info_IV <- tmp_save$all_tool_info_IV
        all_center_info_IV <- tmp_save$all_center_info_IV
      }
    }
    
    ########################################################
    ## 3.11 Variable screening for conventional methods
    ##
    ## ADJ: out = "both"
    ## PSM and IPW: out = "X"
    ########################################################
    
    covariates_both <- try(
      Select_other(
        data = data,
        N_center = N_center,
        type = 2,
        covariates = covariates_all,
        out = "both"
      ),
      silent = TRUE
    )
    
    if (inherits(covariates_both, "try-error")) {
      covariates_both <- empty_covariate_list(N_center)
    }
    
    covariates_X <- try(
      Select_other(
        data = data,
        N_center = N_center,
        type = 2,
        covariates = covariates_all,
        out = "X"
      ),
      silent = TRUE
    )
    
    if (inherits(covariates_X, "try-error")) {
      covariates_X <- empty_covariate_list(N_center)
    }
    
    ########################################################
    ## 3.12 Conventional method: adjusted regression
    ########################################################
    
    result_ADJq <- try(
      ALL_adjust_RD(
        ALLdata = data,
        N_center = N_center,
        covariates = covariates_both,
        min_n = min_n
      ),
      silent = TRUE
    )
    
    final_ADJ <- get_standard_final(result_ADJq)
    
    all_result_ADJ <- rbind(
      all_result_ADJ,
      final_ADJ
    )
    
    ########################################################
    ## 3.13 Conventional method: PSM
    ########################################################
    
    result_PSMq <- try(
      ALL_PSM_RD(
        ALLdata = data,
        N_center = N_center,
        covariates = covariates_X,
        ps_trim = ps_trim,
        caliper_sd = caliper_sd,
        min_n = min_n
      ),
      silent = TRUE
    )
    
    final_PSM <- get_standard_final(result_PSMq)
    
    all_result_PSM <- rbind(
      all_result_PSM,
      final_PSM
    )
    
    ########################################################
    ## 3.14 Conventional method: IPW-MSM
    ########################################################
    
    result_IPWq <- try(
      ALL_IPW_MSM_RD(
        ALLdata = data,
        N_center = N_center,
        covariates = covariates_X,
        ps_trim = ps_trim,
        weight_trunc_quantile = weight_trunc_quantile,
        min_n = min_n
      ),
      silent = TRUE
    )
    
    final_IPW <- get_standard_final(result_IPWq)
    
    all_result_IPW <- rbind(
      all_result_IPW,
      final_IPW
    )
    
    ########################################################
    ## 3.15 Simulation status
    ########################################################
    
    all_sim_status <- rbind(
      all_sim_status,
      data.frame(
        sim = sim,
        data_success = TRUE,
        true_effect = true_effect,
        
        IV_single_same_success = is_ok(final_IV_single_same),
        IV_single_diff_success = is_ok(final_IV_single_diff),
        IV_multi_same_success  = is_ok(final_IV_multi_same),
        IV_multi_diff_success  = is_ok(final_IV_multi_diff),
        
        ADJ_success = is_ok(final_ADJ),
        PSM_success = is_ok(final_PSM),
        IPW_success = is_ok(final_IPW),
        
        conf_strength_x_center_mean = mean(conf_strength_x_center),
        conf_strength_x_center_min = min(conf_strength_x_center),
        conf_strength_x_center_max = max(conf_strength_x_center),
        
        conf_strength_y_center_mean = mean(conf_strength_y_center),
        conf_strength_y_center_min = min(conf_strength_y_center),
        conf_strength_y_center_max = max(conf_strength_y_center),
        
        hetero_delta_x = hetero_delta_x,
        hetero_delta_y = hetero_delta_y,
        hetero_tau_x = hetero_tau_x,
        hetero_tau_y = hetero_tau_y,
        hetero_rho_xy = hetero_rho_xy,
        
        stringsAsFactors = FALSE
      )
    )
  }
  
  ############################################################
  ## 4. Summarize results for causal effect simulation
  ############################################################
  
  end_result_IV_single_same <- resultt_FUn_binary_effect(
    all_result_IV_single_same,
    true_effect_vec
  )
  
  end_result_IV_single_diff <- resultt_FUn_binary_effect(
    all_result_IV_single_diff,
    true_effect_vec
  )
  
  end_result_IV_multi_same <- resultt_FUn_binary_effect(
    all_result_IV_multi_same,
    true_effect_vec
  )
  
  end_result_IV_multi_diff <- resultt_FUn_binary_effect(
    all_result_IV_multi_diff,
    true_effect_vec
  )
  
  end_result_ADJ <- resultt_FUn_binary_effect(
    all_result_ADJ,
    true_effect_vec
  )
  
  end_result_PSM <- resultt_FUn_binary_effect(
    all_result_PSM,
    true_effect_vec
  )
  
  end_result_IPW <- resultt_FUn_binary_effect(
    all_result_IPW,
    true_effect_vec
  )
  
  ##########################################################
  ## 5. Return object
  ##########################################################
  
  out <- list()
  
  out$IV_single_same <- end_result_IV_single_same
  out$IV_single_diff <- end_result_IV_single_diff
  out$IV_multi_same  <- end_result_IV_multi_same
  out$IV_multi_diff  <- end_result_IV_multi_diff
  
  out$ADJ <- end_result_ADJ
  out$PSM <- end_result_PSM
  out$IPW <- end_result_IPW
  
  out$true_effect_vec <- true_effect_vec
  out$true_mean <- mean(true_effect_vec, na.rm = TRUE)
  
  out$all_result_IV_single_same <- all_result_IV_single_same
  out$all_result_IV_single_diff <- all_result_IV_single_diff
  out$all_result_IV_multi_same  <- all_result_IV_multi_same
  out$all_result_IV_multi_diff  <- all_result_IV_multi_diff
  
  out$all_result_ADJ <- all_result_ADJ
  out$all_result_PSM <- all_result_PSM
  out$all_result_IPW <- all_result_IPW
  
  out$tool_info_IV <- all_tool_info_IV
  out$center_info_IV <- all_center_info_IV
  out$sim_status <- all_sim_status
  
  out$multi_iv_threshold <- multi_iv_threshold
  out$max_iv <- max_iv
  out$weight_type_multi <- weight_type_multi
  
  out$hetero_delta_x <- hetero_delta_x
  out$hetero_delta_y <- hetero_delta_y
  out$hetero_tau_x <- hetero_tau_x
  out$hetero_tau_y <- hetero_tau_y
  out$hetero_rho_xy <- hetero_rho_xy
  
  out$conf_strength_y_rct <- conf_strength_y_rct
  
  return(out)
}





############################################################
## Generate beta_CX and beta_CY with strong unobserved U
############################################################

generate_beta_strong_U <- function(m,
                                   k,
                                   conf_L,
                                   C_range = c(0.05, 0.20),
                                   U_range = c(1.00, 2.00),
                                   L_range = c(0.05, 0.20)) {
  
  p_confound <- m + k + length(conf_L)
  
  idx_C <- 1:m
  idx_U <- (m + 1):(m + k)
  idx_L <- (m + k + 1):p_confound
  
  beta_CX <- numeric(p_confound)
  beta_CY <- numeric(p_confound)
  
  beta_CX[idx_C] <- runif(
    length(idx_C),
    min = C_range[1],
    max = C_range[2]
  )
  
  beta_CY[idx_C] <- runif(
    length(idx_C),
    min = C_range[1],
    max = C_range[2]
  )
  
  beta_CX[idx_U] <- runif(
    length(idx_U),
    min = U_range[1],
    max = U_range[2]
  )
  
  beta_CY[idx_U] <- runif(
    length(idx_U),
    min = U_range[1],
    max = U_range[2]
  )
  
  beta_CX[idx_L] <- runif(
    length(idx_L),
    min = L_range[1],
    max = L_range[2]
  )
  
  beta_CY[idx_L] <- runif(
    length(idx_L),
    min = L_range[1],
    max = L_range[2]
  )
  
  return(list(
    beta_CX = beta_CX,
    beta_CY = beta_CY,
    p_confound = p_confound,
    idx_C = idx_C,
    idx_U = idx_U,
    idx_L = idx_L
  ))
}

############################################################
## Diagnostic 1: beta-level contribution of C, U, and L
############################################################

diagnose_beta_contribution <- function(beta_CX,
                                       beta_CY,
                                       m,
                                       k,
                                       conf_L) {
  
  p_confound <- length(beta_CX)
  
  if (length(beta_CY) != p_confound) {
    stop("beta_CX and beta_CY must have the same length.")
  }
  
  idx_C <- 1:m
  idx_U <- (m + 1):(m + k)
  idx_L <- (m + k + 1):p_confound
  
  block_names <- c("C", "U", "L")
  
  block_index <- list(
    C = idx_C,
    U = idx_U,
    L = idx_L
  )
  
  make_one_table <- function(beta, beta_name) {
    
    out <- do.call(
      rbind,
      lapply(block_names, function(block) {
        
        idx <- block_index[[block]]
        beta_block <- beta[idx]
        
        data.frame(
          beta_name = beta_name,
          block = block,
          n = length(idx),
          mean_beta = mean(beta_block),
          mean_abs_beta = mean(abs(beta_block)),
          sd_beta = sd(beta_block),
          min_beta = min(beta_block),
          max_beta = max(beta_block),
          sum_abs_beta = sum(abs(beta_block)),
          sum_beta_sq = sum(beta_block^2),
          stringsAsFactors = FALSE
        )
      })
    )
    
    out$prop_beta_sq <- out$sum_beta_sq / sum(out$sum_beta_sq)
    
    return(out)
  }
  
  table_CX <- make_one_table(beta_CX, "beta_CX")
  table_CY <- make_one_table(beta_CY, "beta_CY")
  
  table_all <- rbind(table_CX, table_CY)
  
  get_value <- function(beta_name, block, var) {
    as.numeric(
      table_all[
        table_all$beta_name == beta_name &
          table_all$block == block,
        var
      ]
    )
  }
  
  summary_vec <- c(
    beta_CX_C_n = get_value("beta_CX", "C", "n"),
    beta_CX_U_n = get_value("beta_CX", "U", "n"),
    beta_CX_L_n = get_value("beta_CX", "L", "n"),
    
    beta_CX_C_mean_abs = get_value("beta_CX", "C", "mean_abs_beta"),
    beta_CX_U_mean_abs = get_value("beta_CX", "U", "mean_abs_beta"),
    beta_CX_L_mean_abs = get_value("beta_CX", "L", "mean_abs_beta"),
    
    beta_CX_C_prop_beta_sq = get_value("beta_CX", "C", "prop_beta_sq"),
    beta_CX_U_prop_beta_sq = get_value("beta_CX", "U", "prop_beta_sq"),
    beta_CX_L_prop_beta_sq = get_value("beta_CX", "L", "prop_beta_sq"),
    
    beta_CY_C_n = get_value("beta_CY", "C", "n"),
    beta_CY_U_n = get_value("beta_CY", "U", "n"),
    beta_CY_L_n = get_value("beta_CY", "L", "n"),
    
    beta_CY_C_mean_abs = get_value("beta_CY", "C", "mean_abs_beta"),
    beta_CY_U_mean_abs = get_value("beta_CY", "U", "mean_abs_beta"),
    beta_CY_L_mean_abs = get_value("beta_CY", "L", "mean_abs_beta"),
    
    beta_CY_C_prop_beta_sq = get_value("beta_CY", "C", "prop_beta_sq"),
    beta_CY_U_prop_beta_sq = get_value("beta_CY", "U", "prop_beta_sq"),
    beta_CY_L_prop_beta_sq = get_value("beta_CY", "L", "prop_beta_sq")
  )
  
  return(list(
    table = table_all,
    summary = summary_vec
  ))
}

############################################################
## Diagnostic 2: eta-level contribution of C, U, and L
############################################################

diagnose_eta_contribution <- function(typeC,
                                      UD,
                                      r,
                                      N_diag,
                                      m,
                                      k,
                                      n1,
                                      n2,
                                      n3,
                                      r_L,
                                      connection_prob,
                                      conf_L,
                                      Clv,
                                      beta_CX,
                                      beta_CY) {
  
  cov_block <- generate_covariate_block(
    typeC = typeC,
    UD = UD,
    r = r,
    N = N_diag,
    m = m,
    k = k,
    n1 = n1,
    n2 = n2,
    n3 = n3,
    r_L = r_L,
    connection_prob = connection_prob,
    conf_L = conf_L,
    Clv = Clv,
    L_network = NULL
  )
  
  H <- as.matrix(cov_block$all_confound)
  
  p_confound <- ncol(H)
  
  idx_C <- 1:m
  idx_U <- (m + 1):(m + k)
  idx_L <- (m + k + 1):p_confound
  
  calc_one <- function(beta, beta_name) {
    
    eta_C <- as.numeric(
      H[, idx_C, drop = FALSE] %*% beta[idx_C]
    )
    
    eta_U <- as.numeric(
      H[, idx_U, drop = FALSE] %*% beta[idx_U]
    )
    
    eta_L <- as.numeric(
      H[, idx_L, drop = FALSE] %*% beta[idx_L]
    )
    
    eta_total <- eta_C + eta_U + eta_L
    
    var_C <- var(eta_C)
    var_U <- var(eta_U)
    var_L <- var(eta_L)
    
    var_sum <- var_C + var_U + var_L
    
    out <- data.frame(
      beta_name = beta_name,
      block = c("C", "U", "L"),
      sd_eta = c(
        sd(eta_C),
        sd(eta_U),
        sd(eta_L)
      ),
      var_eta = c(
        var_C,
        var_U,
        var_L
      ),
      prop_var_eta = c(
        var_C,
        var_U,
        var_L
      ) / var_sum,
      cor_with_total_eta = c(
        cor(eta_C, eta_total),
        cor(eta_U, eta_total),
        cor(eta_L, eta_total)
      ),
      stringsAsFactors = FALSE
    )
    
    return(out)
  }
  
  out_CX <- calc_one(beta_CX, "beta_CX")
  out_CY <- calc_one(beta_CY, "beta_CY")
  
  out <- rbind(out_CX, out_CY)
  
  get_value <- function(beta_name, block, var) {
    as.numeric(
      out[
        out$beta_name == beta_name &
          out$block == block,
        var
      ]
    )
  }
  
  summary_vec <- c(
    eta_CX_C_prop_var = get_value("beta_CX", "C", "prop_var_eta"),
    eta_CX_U_prop_var = get_value("beta_CX", "U", "prop_var_eta"),
    eta_CX_L_prop_var = get_value("beta_CX", "L", "prop_var_eta"),
    
    eta_CY_C_prop_var = get_value("beta_CY", "C", "prop_var_eta"),
    eta_CY_U_prop_var = get_value("beta_CY", "U", "prop_var_eta"),
    eta_CY_L_prop_var = get_value("beta_CY", "L", "prop_var_eta"),
    
    eta_CX_U_cor_total = get_value("beta_CX", "U", "cor_with_total_eta"),
    eta_CY_U_cor_total = get_value("beta_CY", "U", "cor_with_total_eta")
  )
  
  return(list(
    table = out,
    summary = summary_vec
  ))
}

############################################################
## Collect diagnostic metrics into one named vector
############################################################

collect_confounding_diagnostics <- function(beta_diag,
                                            eta_diag) {
  
  out <- c(
    beta_CX_C_prop_beta_sq = beta_diag$summary["beta_CX_C_prop_beta_sq"],
    beta_CX_U_prop_beta_sq = beta_diag$summary["beta_CX_U_prop_beta_sq"],
    beta_CX_L_prop_beta_sq = beta_diag$summary["beta_CX_L_prop_beta_sq"],
    
    beta_CY_C_prop_beta_sq = beta_diag$summary["beta_CY_C_prop_beta_sq"],
    beta_CY_U_prop_beta_sq = beta_diag$summary["beta_CY_U_prop_beta_sq"],
    beta_CY_L_prop_beta_sq = beta_diag$summary["beta_CY_L_prop_beta_sq"],
    
    beta_CX_C_mean_abs = beta_diag$summary["beta_CX_C_mean_abs"],
    beta_CX_U_mean_abs = beta_diag$summary["beta_CX_U_mean_abs"],
    beta_CX_L_mean_abs = beta_diag$summary["beta_CX_L_mean_abs"],
    
    beta_CY_C_mean_abs = beta_diag$summary["beta_CY_C_mean_abs"],
    beta_CY_U_mean_abs = beta_diag$summary["beta_CY_U_mean_abs"],
    beta_CY_L_mean_abs = beta_diag$summary["beta_CY_L_mean_abs"],
    
    eta_CX_C_prop_var = eta_diag$summary["eta_CX_C_prop_var"],
    eta_CX_U_prop_var = eta_diag$summary["eta_CX_U_prop_var"],
    eta_CX_L_prop_var = eta_diag$summary["eta_CX_L_prop_var"],
    
    eta_CY_C_prop_var = eta_diag$summary["eta_CY_C_prop_var"],
    eta_CY_U_prop_var = eta_diag$summary["eta_CY_U_prop_var"],
    eta_CY_L_prop_var = eta_diag$summary["eta_CY_L_prop_var"],
    
    eta_CX_U_cor_total = eta_diag$summary["eta_CX_U_cor_total"],
    eta_CY_U_cor_total = eta_diag$summary["eta_CY_U_cor_total"]
  )
  
  names(out) <- gsub(
    pattern = "\\..*$",
    replacement = "",
    x = names(out)
  )
  
  return(out)
}

############################################################
## 输出路径
############################################################

out_dir <- "results"
dir.create(path.expand(out_dir), recursive = TRUE, showWarnings = FALSE)


############################################################
## 1. Basic parameters
############################################################

m <- 15
k <- 15

n1 <- 5
n2 <- 5
n3 <- 5

r <- 0
UD <- TRUE
Clv <- 0.5

r_L <- 0.5
connection_prob <- 0.5


############################################################
## 2. L-layer confounder settings
############################################################

conf_L <- sample(
  1:(n1 + n2 + n3),
  n3 + 2
)

p_confound <- m + k + length(conf_L)


############################################################
## 3. Generate beta with moderate unobserved U
##
## Target:
##   U contribution around 30%-50%
############################################################

beta_obj <- generate_beta_strong_U(
  m = m,
  k = k,
  conf_L = conf_L,
  C_range = c(0.20, 0.50),
  U_range = c(0.30, 0.55),
  L_range = c(0.20, 0.50)
)

beta_CX <- beta_obj$beta_CX
beta_CY <- beta_obj$beta_CY



############################################################
## 4. Causal effect fixed
##
## r_L scenario:
##   We vary within-layer correlation of L variables only.
##
## For type I error:
##   beta_XY_fixed <- 0
##
## For power:
##   beta_XY_fixed <- 0.4
############################################################

beta_XY_fixed <- 0.4


############################################################
## 5. RCT / do data setting
############################################################

Nrct_fixed <- 10000
p_x <- 0.5


############################################################
## 6. Data type setting
############################################################

type <- 2
typeC <- 1


############################################################
## 7. Binary latent model parameters
############################################################

target_px_ob <- 0.50
target_py0 <- 0.20

conf_strength_x_fixed <- 1
conf_strength_y_fixed <- 1

error_sd_x <- 1
error_sd_y <- 1


############################################################
## 8. 模拟次数
############################################################

N_sim <- 1000


############################################################
## 9. Binary-CONIV 工具筛选参数
############################################################

correlation_threshold_fixed <- 0.10

pi_min <- 0.02
F_min <- 10
p_y_cutoff <- 0.10
y_adjust_x <- TRUE

multi_iv_threshold <- 5
weight_type_multi <- "diag_delta"
max_iv <- Inf


############################################################
## 10. 传统方法参数：PSM 和 IPW-MSM
############################################################

ps_trim <- 0.01
caliper_sd <- 0.20
weight_trunc_quantile <- c(0.01, 0.99)


############################################################
## 11. 通用稳健方差参数
############################################################

hc_type <- "HC1"
min_n <- 20


############################################################
## 12. 未使用的参数
############################################################

n_boot <- NULL
selcetnum <- NULL
trim <- NULL






############################################################
## 13. Scenario settings
##
## hetero_tau_CY_only scenario:
##   only vary between-center random heterogeneity in CY.
############################################################

scenario_type <- "hetero_tau_CY_only"

output_prefix <- "Sim_hetero_tau_CY_only_type2"

N_center_fixed <- 10
Nall_target_fixed <- 10000
Nrct_fixed <- 10000

r_L_fixed <- 0.5
connection_prob_fixed <- 0.5
correlation_threshold_fixed <- 0.10

tau_vec <- seq(0,0.5,0.05)

scenario_grid <- data.frame(
  scenario = seq_along(tau_vec),
  scenario_type = scenario_type,
  varying_parameter = "tau_CY",
  varying_value = tau_vec,
  
  r_L_current = r_L_fixed,
  connection_prob_current = connection_prob_fixed,
  correlation_threshold_current = correlation_threshold_fixed,
  
  N_center_current = N_center_fixed,
  Nall_target_current = Nall_target_fixed,
  Nrct_current = Nrct_fixed,
  
  beta_XY_current = beta_XY_fixed,
  
  conf_strength_x_current = conf_strength_x_fixed,
  conf_strength_y_current = conf_strength_y_fixed,
  
  hetero_delta_x_current = 1,
  hetero_delta_y_current = 1,
  
  hetero_tau_x_current = 0,
  hetero_tau_y_current = tau_vec,
  
  hetero_rho_xy_current = 0,
  
  stringsAsFactors = FALSE
)

all_sim_result_all <- NULL















############################################################
## 14. Loop over heterogeneity values
############################################################

for (j in seq_len(nrow(scenario_grid))) {
  
  scenario_current <- scenario_grid$scenario[j]
  scenario_type_current <- scenario_grid$scenario_type[j]
  varying_parameter_current <- scenario_grid$varying_parameter[j]
  varying_value_current <- scenario_grid$varying_value[j]
  
  r_L_current <- scenario_grid$r_L_current[j]
  connection_prob_current <- scenario_grid$connection_prob_current[j]
  correlation_threshold_current <- scenario_grid$correlation_threshold_current[j]
  
  N_center <- scenario_grid$N_center_current[j]
  Nall_target_current <- scenario_grid$Nall_target_current[j]
  Nrct_current <- scenario_grid$Nrct_current[j]
  
  beta_XY_current <- scenario_grid$beta_XY_current[j]
  
  conf_strength_x_current <- scenario_grid$conf_strength_x_current[j]
  conf_strength_y_current <- scenario_grid$conf_strength_y_current[j]
  
  hetero_delta_x_current <- scenario_grid$hetero_delta_x_current[j]
  hetero_delta_y_current <- scenario_grid$hetero_delta_y_current[j]
  
  hetero_tau_x_current <- scenario_grid$hetero_tau_x_current[j]
  hetero_tau_y_current <- scenario_grid$hetero_tau_y_current[j]
  
  hetero_rho_xy_current <- scenario_grid$hetero_rho_xy_current[j]
  
  cat(
    "Scenario =", scenario_current,
    "scenario_type =", scenario_type_current,
    "varying_parameter =", varying_parameter_current,
    "varying_value =", varying_value_current,
    "hetero_delta_x =", hetero_delta_x_current,
    "hetero_delta_y =", hetero_delta_y_current,
    "hetero_tau_x =", hetero_tau_x_current,
    "hetero_tau_y =", hetero_tau_y_current,
    "hetero_rho_xy =", hetero_rho_xy_current,
    "r_L =", r_L_current,
    "connection_prob =", connection_prob_current,
    "correlation_threshold =", correlation_threshold_current,
    "N_center =", N_center,
    "Nall_target =", Nall_target_current,
    "Nrct =", Nrct_current,
    "beta_XY_latent =", beta_XY_current,
    "type = 2",
    "\n"
  )
  
  ##########################################################
  ## Diagnostics for current confounding structure
  ##
  ## Note:
  ##   beta-level and eta-level diagnostics are still based
  ##   on the baseline beta_CX and beta_CY.
  ##
  ##   Center-level heterogeneity is recorded through
  ##   hetero_delta_x/y and hetero_tau_x/y.
  ##########################################################
  
  beta_diag <- diagnose_beta_contribution(
    beta_CX = beta_CX,
    beta_CY = beta_CY,
    m = m,
    k = k,
    conf_L = conf_L
  )
  
  eta_diag <- diagnose_eta_contribution(
    typeC = typeC,
    UD = UD,
    r = r,
    N_diag = 10000,
    m = m,
    k = k,
    n1 = n1,
    n2 = n2,
    n3 = n3,
    r_L = r_L_current,
    connection_prob = connection_prob_current,
    conf_L = conf_L,
    Clv = Clv,
    beta_CX = beta_CX,
    beta_CY = beta_CY
  )
  
  confounding_diag_vec <- collect_confounding_diagnostics(
    beta_diag = beta_diag,
    eta_diag = eta_diag
  )
  
  print(beta_diag$table)
  print(eta_diag$table)
  
  ##########################################################
  ## Generate sample size for each center
  ##########################################################
  
  Nall <- floor(
    runif(
      N_center,
      Nall_target_current * 0.8,
      Nall_target_current * 1.2
    )
  )
  
  ##########################################################
  ## Run simulation
  ##########################################################
  
  result_sim_result <- All_simulation_binary(
    N_sim = N_sim,
    N_center = N_center,
    type = type,
    typeC = typeC,
    
    UD = UD,
    r = r,
    Nall = Nall,
    Nrct = Nrct_current,
    
    m = m,
    k = k,
    n1 = n1,
    n2 = n2,
    n3 = n3,
    
    r_L = r_L_current,
    connection_prob = connection_prob_current,
    conf_L = conf_L,
    
    beta_CX = beta_CX,
    beta_CY = beta_CY,
    beta_XY = beta_XY_current,
    
    Clv = Clv,
    p_x = p_x,
    
    target_px_ob = target_px_ob,
    target_py0 = target_py0,
    conf_strength_x = conf_strength_x_current,
    conf_strength_y = conf_strength_y_current,
    
    hetero_delta_x = hetero_delta_x_current,
    hetero_delta_y = hetero_delta_y_current,
    hetero_tau_x = hetero_tau_x_current,
    hetero_tau_y = hetero_tau_y_current,
    hetero_rho_xy = hetero_rho_xy_current,
    conf_strength_y_rct = conf_strength_y_fixed,
    
    error_sd_x = error_sd_x,
    error_sd_y = error_sd_y,
    
    correlation_threshold = correlation_threshold_current,
    pi_min = pi_min,
    F_min = F_min,
    p_y_cutoff = p_y_cutoff,
    y_adjust_x = y_adjust_x,
    
    multi_iv_threshold = multi_iv_threshold,
    max_iv = max_iv,
    weight_type_multi = weight_type_multi,
    
    ps_trim = ps_trim,
    caliper_sd = caliper_sd,
    weight_trunc_quantile = weight_trunc_quantile,
    
    hc_type = hc_type,
    min_n = min_n,
    
    n_boot = n_boot,
    trim = trim
  )
  
  ##########################################################
  ## Collect one-row result
  ##########################################################
  
  result_row <- c(
    
    setNames(
      result_sim_result$IV_single_same,
      paste0("IV_single_same_", names(result_sim_result$IV_single_same))
    ),
    
    setNames(
      result_sim_result$IV_single_diff,
      paste0("IV_single_diff_", names(result_sim_result$IV_single_diff))
    ),
    
    setNames(
      result_sim_result$IV_multi_same,
      paste0("IV_multi_same_", names(result_sim_result$IV_multi_same))
    ),
    
    setNames(
      result_sim_result$IV_multi_diff,
      paste0("IV_multi_diff_", names(result_sim_result$IV_multi_diff))
    ),
    
    setNames(
      result_sim_result$ADJ,
      paste0("ADJ_", names(result_sim_result$ADJ))
    ),
    
    setNames(
      result_sim_result$PSM,
      paste0("PSM_", names(result_sim_result$PSM))
    ),
    
    setNames(
      result_sim_result$IPW,
      paste0("IPW_", names(result_sim_result$IPW))
    ),
    
    true_mean = result_sim_result$true_mean,
    true_min = min(result_sim_result$true_effect_vec, na.rm = TRUE),
    true_max = max(result_sim_result$true_effect_vec, na.rm = TRUE),
    
    scenario = scenario_current,
    scenario_type = scenario_type_current,
    
    varying_parameter = varying_parameter_current,
    varying_value = varying_value_current,
    
    beta_XY_latent = beta_XY_current,
    
    N_center = N_center,
    Nall_target = Nall_target_current,
    Nrct = Nrct_current,
    
    conf_strength_x_current = conf_strength_x_current,
    conf_strength_y_current = conf_strength_y_current,
    
    hetero_delta_x = hetero_delta_x_current,
    hetero_delta_y = hetero_delta_y_current,
    hetero_tau_x = hetero_tau_x_current,
    hetero_tau_y = hetero_tau_y_current,
    hetero_rho_xy = hetero_rho_xy_current,
    
    m = m,
    k = k,
    n1 = n1,
    n2 = n2,
    n3 = n3,
    
    confounding_diag_vec,
    
    typeC = typeC,
    type = type,
    
    p_x = p_x,
    
    Nall_mean = mean(Nall),
    Nall_min = min(Nall),
    Nall_max = max(Nall),
    Nall_total = sum(Nall),
    
    r = r,
    UD = as.numeric(UD),
    Clv = Clv,
    
    r_L = r_L_current,
    connection_prob = connection_prob_current,
    
    target_px_ob = target_px_ob,
    target_py0 = target_py0,
    
    conf_strength_x = conf_strength_x_current,
    conf_strength_y = conf_strength_y_current,
    
    error_sd_x = error_sd_x,
    error_sd_y = error_sd_y,
    
    correlation_threshold = correlation_threshold_current,
    pi_min = pi_min,
    F_min = F_min,
    p_y_cutoff = p_y_cutoff,
    y_adjust_x = as.numeric(y_adjust_x),
    
    multi_iv_threshold = multi_iv_threshold,
    max_iv = max_iv,
    weight_type_multi = weight_type_multi,
    
    ps_trim = ps_trim,
    caliper_sd = caliper_sd,
    weight_trunc_low = weight_trunc_quantile[1],
    weight_trunc_high = weight_trunc_quantile[2],
    
    N_sim = N_sim
  )
  
  result_row <- as.data.frame(
    as.list(result_row),
    stringsAsFactors = FALSE
  )
  
  all_sim_result_all <- rbind(
    all_sim_result_all,
    result_row
  )
  
  write.csv(
    all_sim_result_all,
    file = file.path(
      path.expand(out_dir),
      paste0(output_prefix, "_running.csv")
    ),
    row.names = FALSE
  )
}


############################################################
## 15. Save final result
############################################################

write.csv(
  all_sim_result_all,
  file = file.path(
    path.expand(out_dir),
    paste0(output_prefix, ".csv")
  ),
  row.names = FALSE
)