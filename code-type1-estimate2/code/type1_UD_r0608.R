
############################################################
## 0. 设置每个 R 任务使用的 CPU 线程数
## 必须放在脚本最前面
############################################################

n_threads <- 6

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
## 1. 产生数据
############################################################

########分类函数############################
Fenlei <- function(data, Clv, na.rm = TRUE) {
  
  if (Clv <= 0 || Clv >= 1) {
    stop("Clv must be between 0 and 1.")
  }
  
  cutoff <- quantile(data, probs = Clv, na.rm = na.rm, names = FALSE)
  
  newdata <- ifelse(data > cutoff, 1, 0)
  
  return(newdata)
}


## 先生成一套固定的 L 层网络结构 ###########


generate_layer_network <- function(parent_indices,
                                   current_layer_indices,
                                   connection_prob) {
  
  if (connection_prob < 0 || connection_prob > 1) {
    stop("connection_prob must be between 0 and 1.")
  }
  
  layer_network <- vector("list", length(current_layer_indices))
  
  for (j in seq_along(current_layer_indices)) {
    
    child_idx <- current_layer_indices[j]
    
    # 固定当前 child 由哪些 parent 连接而来
    is_parent <- runif(length(parent_indices)) < connection_prob
    active_parent_indices <- parent_indices[is_parent]
    
    if (length(active_parent_indices) > 0) {
      
      # 固定边权重
      weights <- runif(length(active_parent_indices), min = -1, max = 1)
      
      # 权重归一化，避免父节点越多方差越大
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



## 使用固定网络结构生成某一层数据#########################


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
      if (sd(new_var, na.rm = TRUE) > 0) {
        new_var <- as.numeric(scale(new_var))
      }
    }
    
    data_matrix[, child_idx] <- new_var
  }
  
  return(data_matrix)
}


## 生成单观察性中心数据##########

Data_Generated <- function(type = 1, typeC, UD = TRUE, r, N, m, k,
                           n1, n2, n3, r_L, connection_prob, conf_L,
                           beta_CX, beta_CY, beta_XY,
                           Clv,
                           L_network = NULL) {
  
  NC <- m + k
  
  ##########################################################
  ## 1. 生成基础 C/U 协变量
  ##########################################################
  
  Sigma_CU <- matrix(rep(r, NC * NC), NC, NC)
  
  if (UD == TRUE) {
    Sigma_CU[, (m + 1):NC] <- 0
    Sigma_CU[(m + 1):NC, ] <- 0
  }
  
  diag(Sigma_CU) <- 1
  
  mu_CU <- rep(0, NC)
  
  Cmat <- MASS::mvrnorm(n = N, mu = mu_CU, Sigma = Sigma_CU)
  colnames(Cmat) <- paste0("NC", 1:NC)
  
  ##########################################################
  ## 2. 初始化 L 层数据矩阵
  ##########################################################
  
  total_vars <- n1 + n2 + n3
  
  simulated_data <- matrix(0, nrow = N, ncol = total_vars)
  colnames(simulated_data) <- c(
    paste0("L1_", 1:n1),
    paste0("L2_", 1:n2),
    paste0("L3_", 1:n3)
  )
  
  ##########################################################
  ## 3. 生成 L1
  ##########################################################
  
  Sigma_L <- matrix(r_L, nrow = n1, ncol = n1)
  diag(Sigma_L) <- 1
  
  mu_L <- rep(0, n1)
  
  l1_data <- MASS::mvrnorm(n = N, mu = mu_L, Sigma = Sigma_L)
  simulated_data[, 1:n1] <- l1_data
  
  ##########################################################
  ## 4. 使用固定网络生成 L2 和 L3
  ##########################################################
  
  l1_indices <- 1:n1
  l2_indices <- (n1 + 1):(n1 + n2)
  l3_indices <- (n1 + n2 + 1):total_vars
  
  # 如果外部没有传入 L_network，则单独生成一套。
  # 在 ALL_data() 中会统一传入同一套 L_network。
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
  
  ##########################################################
  ## 5. 合并协变量
  ##########################################################
  
  df_L <- as.data.frame(simulated_data)
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
  
  ##########################################################
  ## 6. 生成 X 和 Y
  ##########################################################
  
  px <- as.numeric(t(c(beta_CX) %*% t(all_confound)))
  
  if (type == 1) {
    
    x <- px + rnorm(N, 0, 1)
    
    y <- as.numeric(t(c(beta_CY) %*% t(all_confound))) +
      beta_XY * x +
      rnorm(N, 0, 1)
    
    data_total <- data.frame(x, y, all_cov)
    
    colnames(data_total)[1:(2 + NC)] <- c(
      "X",
      "Y",
      paste0("C", 1:m),
      paste0("U", 1:k)
    )
    
  } else {
    
    px_prob <- pmin(pmax(px, 1e-6), 1 - 1e-6)
    x <- rbinom(N, size = 1, prob = px_prob)
    
    py <- as.numeric(t(c(beta_CY) %*% t(all_confound))) +
      beta_XY * x
    
    py_prob <- pmin(pmax(py, 1e-6), 1 - 1e-6)
    y <- rbinom(N, size = 1, prob = py_prob)
    
    data_total <- data.frame(x, y, all_cov)
    
    colnames(data_total)[1:(2 + NC)] <- c(
      "X",
      "Y",
      paste0("C", 1:m),
      paste0("U", 1:k)
    )
  }
  
  ##########################################################
  ## 7. 删除不可观测 U
  ##########################################################
  
  cols_to_remove <- paste0("U", 1:k)
  
  data_total_end <- data_total[
    ,
    -which(names(data_total) %in% cols_to_remove),
    drop = FALSE
  ]
  
  return(data_total_end)
}


## 生成 RCT / do 数据：固定 L 网络结构版########


Data_Generated_do <- function(type = 1, typeC, UD = TRUE, r, N, m, k,
                              n1, n2, n3, r_L, connection_prob, conf_L,
                              beta_CX, beta_CY, beta_XY,
                              Clv, p_x,
                              L_network = NULL) {
  
  NC <- m + k
  
  ##########################################################
  ## 1. 生成基础 C/U 协变量
  ##########################################################
  
  Sigma_CU <- matrix(rep(r, NC * NC), NC, NC)
  
  if (UD == TRUE) {
    Sigma_CU[, (m + 1):NC] <- 0
    Sigma_CU[(m + 1):NC, ] <- 0
  }
  
  diag(Sigma_CU) <- 1
  
  mu_CU <- rep(0, NC)
  
  Cmat <- MASS::mvrnorm(n = N, mu = mu_CU, Sigma = Sigma_CU)
  colnames(Cmat) <- paste0("NC", 1:NC)
  
  ##########################################################
  ## 2. 初始化 L 层数据矩阵
  ##########################################################
  
  total_vars <- n1 + n2 + n3
  
  simulated_data <- matrix(0, nrow = N, ncol = total_vars)
  colnames(simulated_data) <- c(
    paste0("L1_", 1:n1),
    paste0("L2_", 1:n2),
    paste0("L3_", 1:n3)
  )
  
  ##########################################################
  ## 3. 生成 L1
  ##########################################################
  
  Sigma_L <- matrix(r_L, nrow = n1, ncol = n1)
  diag(Sigma_L) <- 1
  
  mu_L <- rep(0, n1)
  
  l1_data <- MASS::mvrnorm(n = N, mu = mu_L, Sigma = Sigma_L)
  simulated_data[, 1:n1] <- l1_data
  
  ##########################################################
  ## 4. 使用固定网络生成 L2 和 L3
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
  
  ##########################################################
  ## 5. 合并协变量
  ##########################################################
  
  df_L <- as.data.frame(simulated_data)
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
  
  ##########################################################
  ## 6. 生成 X 和 Y
  ##########################################################
  
  # do / RCT 中 X 由 p_x 控制
  if (length(p_x) == 1) {
    px <- rep(p_x, N)
  } else {
    px <- p_x
  }
  
  if (length(px) != N) {
    stop("p_x must be either a scalar or a vector with length N.")
  }
  
  if (type == 1) {
    
    x <- px + rnorm(N, 0, 1)
    
    y <- as.numeric(t(c(beta_CY) %*% t(all_confound))) +
      beta_XY * x +
      rnorm(N, 0, 1)
    
    data_total <- data.frame(x, y, all_cov)
    
    colnames(data_total)[1:(2 + NC)] <- c(
      "X",
      "Y",
      paste0("C", 1:m),
      paste0("U", 1:k)
    )
    
  } else {
    
    px_prob <- pmin(pmax(px, 1e-6), 1 - 1e-6)
    x <- rbinom(N, size = 1, prob = px_prob)
    
    py <- as.numeric(t(c(beta_CY) %*% t(all_confound))) +
      beta_XY * x
    
    py_prob <- pmin(pmax(py, 1e-6), 1 - 1e-6)
    y <- rbinom(N, size = 1, prob = py_prob)
    
    data_total <- data.frame(x, y, all_cov)
    
    colnames(data_total)[1:(2 + NC)] <- c(
      "X",
      "Y",
      paste0("C", 1:m),
      paste0("U", 1:k)
    )
  }
  
  ##########################################################
  ## 7. 删除不可观测 U
  ##########################################################
  
  cols_to_remove <- paste0("U", 1:k)
  
  data_total_end <- data_total[
    ,
    -which(names(data_total) %in% cols_to_remove),
    drop = FALSE
  ]
  
  return(data_total_end)
}



## 总数据生成：所有中心和 RCT/do 共用同一个 L_network#######

ALL_data <- function(type, N_center, typeC, UD, r, Nall, Nrct, m, k,
                     n1, n2, n3, r_L, connection_prob, conf_L,
                     beta_CX, beta_CY, beta_XY,
                     Clv, p_x) {
  
  ##########################################################
  ## 1. 先固定一套 L 网络结构
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
  ## 2. 生成多个观察性中心
  ##########################################################
  
  data_list <- vector("list", N_center)
  
  for (i in 1:N_center) {
    
    dataI <- Data_Generated(
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
      L_network = L_network
    )
    
    dataI$centre <- i
    data_list[[i]] <- dataI
  }
  
  all_data <- do.call(rbind, data_list)
  
  ##########################################################
  ## 3. 生成 RCT / do 数据
  ##########################################################
  
  rct_data <- Data_Generated_do(
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
    L_network = L_network
  )
  
  ##########################################################
  ## 4. 返回数据和固定网络
  ##########################################################
  
  return(list(
    all_data = all_data,
    rct = rct_data,
    L_network = L_network
  ))
}



## 2 CON_LV
############################################################

###安全提取 lm 中某个变量的 P 值###########

safe_lm_pvalue <- function(formula, data, term_name) {
  
  fit <- try(lm(formula, data = data), silent = TRUE)
  
  if (inherits(fit, "try-error")) {
    return(NA_real_)
  }
  
  coefs <- try(summary(fit)$coefficients, silent = TRUE)
  
  if (inherits(coefs, "try-error")) {
    return(NA_real_)
  }
  
  # 如果该变量因为常数、共线性等原因没有进入模型，则返回 NA
  if (!term_name %in% rownames(coefs)) {
    return(NA_real_)
  }
  
  return(coefs[term_name, 4])
}

## 安全计算绝对相关系数#####


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



##  统一变量筛选函数
##
## selcetnum:
##   selcetnum = 1       ：每个中心只选最强 1 个工具变量
##   selcetnum = "multi" ：保留多个工具变量，原则与原函数一致
##   selcetnum = 3       ：可选，最多保留前 3 个


variable_selection_sorted <- function(N_center,
                                      covariates,
                                      data,
                                      correlation_threshold,
                                      type = 1,
                                      selcetnum = "multi") {
  ## 0. 基本检查
  
  if (length(covariates) == 0) {
    stop("covariates is empty.")
  }
  
  if (correlation_threshold < 0 || correlation_threshold > 1) {
    stop("correlation_threshold must be between 0 and 1.")
  }
  
  if (!(type %in% c(1, 2))) {
    stop("type must be 1 or 2.")
  }
  
  # selcetnum 可以是：
  #   1       ：始终只选最强 1 个工具
  #   >1 数字 ：多工具门槛；如果候选数 <= selcetnum，则只选最强 1 个；
  #             如果候选数 > selcetnum，则保留多个工具
  #   "multi" ：完全多工具模式，筛出多少个就保留多少个
  if (!(is.numeric(selcetnum) || selcetnum == "multi")) {
    stop("selcetnum must be 1, a positive integer, or 'multi'.")
  }
  
  if (is.numeric(selcetnum)) {
    if (length(selcetnum) != 1 || selcetnum < 1 || selcetnum != floor(selcetnum)) {
      stop("numeric selcetnum must be a single positive integer.")
    }
  }
  
  all_final_vars <- list()
  
  ## 1. 根据 type 设置 P 值阈值
  
  if (type == 1) {
    # 连续型：Bonferroni 校正
    p_cutoff <- 0.05 / length(covariates)
  } else {
    # 二分类：阈值 0.1
    p_cutoff <- 0.1
  }
  
  ## 对每个中心分别筛选变量
  for (i in 1:N_center) {
    
    ## 2.1 提取第 i 个中心的数据
    data_i <- data$all_data[
      data$all_data$centre == i,
      ,
      drop = FALSE
    ]
    
    # 只保留当前中心和 RCT 中都存在的候选变量
    covariates_i <- covariates[
      covariates %in% colnames(data_i) &
        covariates %in% colnames(data$rct)
    ]
    
    if (length(covariates_i) == 0) {
      all_final_vars[[paste0("center_", i)]] <- character(0)
      next
    }
    
    p_value_Y <- rep(NA_real_, length(covariates_i))
    p_value_X <- rep(NA_real_, length(covariates_i))
    
    ## 2.2 初筛：同时与 Y 和 X 相关
    
    for (j in seq_along(covariates_i)) {
      
      var_j <- covariates_i[j]
      
      # Y ~ X + var_j
      formula_Y <- reformulate(
        termlabels = c("X", var_j),
        response = "Y"
      )
      
      # X ~ var_j
      formula_X <- reformulate(
        termlabels = var_j,
        response = "X"
      )
      
      p_value_Y[j] <- safe_lm_pvalue(
        formula = formula_Y,
        data = data_i,
        term_name = var_j
      )
      
      p_value_X[j] <- safe_lm_pvalue(
        formula = formula_X,
        data = data_i,
        term_name = var_j
      )
    }
    
    # 与 Y 显著相关的变量
    varable1 <- covariates_i[
      which(!is.na(p_value_Y) & p_value_Y < p_cutoff)
    ]
    
    # 与 X 显著相关的变量
    varable2 <- covariates_i[
      which(!is.na(p_value_X) & p_value_X < p_cutoff)
    ]
    
    # 同时与 Y 和 X 相关的变量
    common_variables <- intersect(varable1, varable2)
    
    if (length(common_variables) == 0) {
      all_final_vars[[paste0("center_", i)]] <- character(0)
      next
    }
    
    final_vars <- common_variables
    
    ## 2.3 处理高度相关变量
    ##
    ## 如果两个变量之间相关性超过 correlation_threshold，
    ## 则保留与 X 相关性更强的那个。
    
    if (length(common_variables) > 1) {
      
      selectvar <- data_i[, common_variables, drop = FALSE]
      
      cor_matrix <- suppressWarnings(
        abs(cor(selectvar, use = "pairwise.complete.obs"))
      )
      
      vars_to_remove <- character(0)
      
      for (idx_common_var_1 in 1:(length(common_variables) - 1)) {
        for (idx_common_var_2 in (idx_common_var_1 + 1):length(common_variables)) {
          
          var1_name <- common_variables[idx_common_var_1]
          var2_name <- common_variables[idx_common_var_2]
          
          if (var1_name %in% vars_to_remove ||
              var2_name %in% vars_to_remove) {
            next
          }
          
          cor_val <- cor_matrix[var1_name, var2_name]
          
          if (!is.na(cor_val) && cor_val > correlation_threshold) {
            
            cor1_with_X <- safe_abs_cor(data_i$X, data_i[[var1_name]])
            cor2_with_X <- safe_abs_cor(data_i$X, data_i[[var2_name]])
            
            if (cor1_with_X >= cor2_with_X) {
              vars_to_remove <- c(vars_to_remove, var2_name)
            } else {
              vars_to_remove <- c(vars_to_remove, var1_name)
            }
          }
        }
      }
      
      final_vars <- setdiff(common_variables, vars_to_remove)
    }
    
    ## 2.4 按与 X 的绝对相关性从大到小排序
    if (length(final_vars) > 1) {
      
      correlations_with_X <- sapply(final_vars, function(var) {
        safe_abs_cor(data_i$X, data_i[[var]])
      })
      
      final_vars <- final_vars[
        order(correlations_with_X, decreasing = TRUE)
      ]
    }
    
    
    ## 2.5 只保留 RCT 数据中也存在的变量
    
    final_vars_end <- final_vars[
      final_vars %in% colnames(data$rct)
    ]
    
    ########################################################
    ## 2.6 根据 selcetnum 控制选择 1 个还是多个
    ##
    ## selcetnum = 1：
    ##   始终只选最强的 1 个工具变量。
    ##
    ## selcetnum = 5：
    ##   如果筛选后工具变量数 <= 5，则只选最强 1 个；
    ##   如果筛选后工具变量数 > 5，则保留多个工具变量。
    ##
    ## selcetnum = "multi"：
    ##   完全多工具模式，筛出多少个就保留多少个。
    ########################################################
    
    if (length(final_vars_end) > 0) {
      
      if (is.numeric(selcetnum)) {
        
        if (selcetnum == 1) {
          
          # 单工具模式：只选最强的 1 个
          final_vars_end <- final_vars_end[1]
          
        } else if (selcetnum > 1) {
          
          # 多工具门槛模式：
          # 如果候选工具数 <= selcetnum，则退回最强单工具；
          # 如果候选工具数 > selcetnum，则保留多个工具。
          if (length(final_vars_end) <= selcetnum) {
            final_vars_end <- final_vars_end[1]
          } else {
            final_vars_end <- final_vars_end
          }
        }
        
      } else if (selcetnum == "multi") {
        
        # 完全多工具模式：不设门槛
        final_vars_end <- final_vars_end
      }
    }
    
    all_final_vars[[paste0("center_", i)]] <- final_vars_end
  }
  
  return(all_final_vars)
}




## RCT 中多个工具变量 Y ~ Z 斜率的协方差矩阵

rct_slope_vcov <- function(data_rct, iv_names) {
  
  iv_names <- unique(iv_names)
  iv_names <- iv_names[
    !is.na(iv_names) &
      iv_names %in% colnames(data_rct)
  ]
  
  if (length(iv_names) == 0) {
    stop("No valid IVs found in RCT data.")
  }
  
  n <- nrow(data_rct)
  K <- length(iv_names)
  
  IF_mat <- matrix(NA_real_, nrow = n, ncol = K)
  colnames(IF_mat) <- iv_names
  
  beta_rct <- rep(NA_real_, K)
  names(beta_rct) <- iv_names
  
  for (j in seq_along(iv_names)) {
    
    iv_j <- iv_names[j]
    
    tmp_data <- data.frame(
      Y = data_rct$Y,
      Z = data_rct[[iv_j]]
    )
    
    fit <- try(lm(Y ~ Z, data = tmp_data), silent = TRUE)
    
    if (inherits(fit, "try-error")) {
      next
    }
    
    coefs <- summary(fit)$coefficients
    
    if (!"Z" %in% rownames(coefs)) {
      next
    }
    
    beta_rct[j] <- coefs["Z", "Estimate"]
    
    z <- tmp_data$Z
    z_centered <- z - mean(z, na.rm = TRUE)
    Sxx <- sum(z_centered^2, na.rm = TRUE)
    
    if (!is.finite(Sxx) || Sxx <= 0) {
      next
    }
    
    res <- resid(fit)
    
    # beta_hat - beta ≈ sum_i (z_i - zbar) * e_i / Sxx
    # 所以每个样本的影响贡献为：
    IF_mat[, j] <- z_centered * res / Sxx
  }
  
  valid <- is.finite(beta_rct) &
    colSums(is.finite(IF_mat)) == n
  
  beta_rct <- beta_rct[valid]
  IF_mat <- IF_mat[, valid, drop = FALSE]
  
  if (length(beta_rct) == 0) {
    stop("No valid RCT slope estimates.")
  }
  
  vcov_rct <- crossprod(IF_mat)
  
  return(list(
    beta = beta_rct,
    vcov = vcov_rct
  ))
}

## 观察性中心中 X ~ Z 与 Y ~ Z 斜率的影响函数
## 每列对应一个工具变量；crossprod() 给出斜率估计的方差和协方差

ob_slope_if <- function(data_OB, iv_names) {

  n <- nrow(data_OB)
  K <- length(iv_names)

  IF_x <- matrix(NA_real_, nrow = n, ncol = K)
  IF_y <- matrix(NA_real_, nrow = n, ncol = K)
  colnames(IF_x) <- iv_names
  colnames(IF_y) <- iv_names

  x_centered <- data_OB$X - mean(data_OB$X)
  y_centered <- data_OB$Y - mean(data_OB$Y)

  for (j in seq_len(K)) {

    z <- data_OB[[iv_names[j]]]
    z_centered <- z - mean(z)
    Sxx <- sum(z_centered^2)

    if (!is.finite(Sxx) || Sxx <= 0) {
      next
    }

    res_x <- x_centered - sum(z_centered * x_centered) / Sxx * z_centered
    res_y <- y_centered - sum(z_centered * y_centered) / Sxx * z_centered

    IF_x[, j] <- z_centered * res_x / Sxx
    IF_y[, j] <- z_centered * res_y / Sxx
  }

  return(list(IF_x = IF_x, IF_y = IF_y))
}

ratio_se_delta <- function(betaYG, betaXG, SEYG, SEXG, cov_YX = 0) {
  
  if (is.na(betaYG) || is.na(betaXG) ||
      is.na(SEYG) || is.na(SEXG) ||
      !is.finite(betaYG) || !is.finite(betaXG) ||
      !is.finite(SEYG) || !is.finite(SEXG) ||
      abs(betaXG) < 1e-8) {
    return(NA_real_)
  }
  
  var_ratio <- (SEYG^2 / betaXG^2) +
    (betaYG^2 * SEXG^2 / betaXG^4) -
    (2 * betaYG * cov_YX / betaXG^3)
  
  return(sqrt(max(var_ratio, 0)))
}

## 单中心 Wald ratio 估计：返回协方差 Meta 所需信息

############################################################
## 单中心单工具 Wald ratio 估计
## 返回共享 RCT 协方差 Meta 所需信息 + 工具诊断信息
############################################################

CON_IV_luo_top1_return <- function(data_rct_all,
                                   data_OB,
                                   IV_covariates,
                                   rct_info) {
  
  ##########################################################
  ## 内部空结果函数：保证所有失败返回列名一致
  ##########################################################
  
  empty_result <- function(IV_name = NA_character_,
                           betaXG = NA_real_,
                           SEXG = NA_real_,
                           F_XZ = NA_real_,
                           p_XZ = NA_real_,
                           n_OB = NA_integer_) {
    
    return(data.frame(
      b = NA_real_,
      se = NA_real_,
      
      betaXG = betaXG,
      SEXG = SEXG,
      F_XZ = F_XZ,
      p_XZ = p_XZ,
      
      betaY_OB = NA_real_,
      SEY_OB = NA_real_,
      
      betaY_RCT = NA_real_,
      SEY_RCT = NA_real_,
      
      betaYG = NA_real_,
      SEYG = NA_real_,
      
      IV_name = IV_name,
      n_OB = n_OB,
      
      stringsAsFactors = FALSE
    ))
  }
  
  n_OB <- nrow(data_OB)
  
  ##########################################################
  ## 0. 无工具变量时返回 NA
  ##########################################################
  
  if (length(IV_covariates) == 0) {
    return(empty_result(
      IV_name = NA_character_,
      n_OB = n_OB
    ))
  }
  
  ##########################################################
  ## 1. 单工具版本只取第一个工具变量
  ##########################################################
  
  IV_name <- IV_covariates[1]
  
  if (!IV_name %in% colnames(data_OB) ||
      !IV_name %in% colnames(data_rct_all) ||
      !IV_name %in% names(rct_info$beta)) {
    
    return(empty_result(
      IV_name = IV_name,
      n_OB = n_OB
    ))
  }
  
  ##########################################################
  ## 2. 观察性中心中估计 X ~ Z
  ##########################################################
  
  tmp_x <- data.frame(
    X = data_OB$X,
    Z = data_OB[[IV_name]]
  )
  
  modelxG <- try(
    lm(X ~ Z, data = tmp_x),
    silent = TRUE
  )
  
  if (inherits(modelxG, "try-error")) {
    return(empty_result(
      IV_name = IV_name,
      n_OB = n_OB
    ))
  }
  
  sx <- try(
    summary(modelxG)$coefficients,
    silent = TRUE
  )
  
  if (inherits(sx, "try-error") || !"Z" %in% rownames(sx)) {
    return(empty_result(
      IV_name = IV_name,
      n_OB = n_OB
    ))
  }
  
  betaXG <- sx["Z", "Estimate"]
  SEXG   <- sx["Z", "Std. Error"]
  
  if ("Pr(>|t|)" %in% colnames(sx)) {
    p_XZ <- sx["Z", "Pr(>|t|)"]
  } else {
    p_XZ <- NA_real_
  }
  
  if (!is.finite(betaXG) ||
      !is.finite(SEXG) ||
      SEXG <= 0 ||
      abs(betaXG) < 1e-8) {
    
    return(empty_result(
      IV_name = IV_name,
      betaXG = betaXG,
      SEXG = SEXG,
      p_XZ = p_XZ,
      n_OB = n_OB
    ))
  }
  
  F_XZ <- (betaXG / SEXG)^2
  
  ##########################################################
  ## 3. 观察性中心中估计 Y ~ Z
  ##########################################################
  
  tmp_y <- data.frame(
    Y = data_OB$Y,
    Z = data_OB[[IV_name]]
  )
  
  modelYG_OB <- try(
    lm(Y ~ Z, data = tmp_y),
    silent = TRUE
  )
  
  if (inherits(modelYG_OB, "try-error")) {
    return(empty_result(
      IV_name = IV_name,
      betaXG = betaXG,
      SEXG = SEXG,
      F_XZ = F_XZ,
      p_XZ = p_XZ,
      n_OB = n_OB
    ))
  }
  
  sy_ob <- try(
    summary(modelYG_OB)$coefficients,
    silent = TRUE
  )
  
  if (inherits(sy_ob, "try-error") || !"Z" %in% rownames(sy_ob)) {
    return(empty_result(
      IV_name = IV_name,
      betaXG = betaXG,
      SEXG = SEXG,
      F_XZ = F_XZ,
      p_XZ = p_XZ,
      n_OB = n_OB
    ))
  }
  
  betaY_OB <- sy_ob["Z", "Estimate"]
  SEY_OB   <- sy_ob["Z", "Std. Error"]
  
  if (!is.finite(betaY_OB) ||
      !is.finite(SEY_OB) ||
      SEY_OB <= 0) {
    
    return(empty_result(
      IV_name = IV_name,
      betaXG = betaXG,
      SEXG = SEXG,
      F_XZ = F_XZ,
      p_XZ = p_XZ,
      n_OB = n_OB
    ))
  }
  
  ##########################################################
  ## 4. RCT / do 数据中对应的 Y ~ Z 斜率和方差
  ##########################################################
  
  betaY_RCT <- rct_info$beta[IV_name]
  varY_RCT  <- rct_info$vcov[IV_name, IV_name]
  
  if (!is.finite(betaY_RCT) ||
      !is.finite(varY_RCT) ||
      varY_RCT < 0) {
    
    return(empty_result(
      IV_name = IV_name,
      betaXG = betaXG,
      SEXG = SEXG,
      F_XZ = F_XZ,
      p_XZ = p_XZ,
      n_OB = n_OB
    ))
  }
  
  SEY_RCT <- sqrt(varY_RCT)
  
  ##########################################################
  ## 5. Wald ratio
  ##
  ## betaYG = betaY_OB - betaY_RCT
  ## betaXY = betaYG / betaXG
  ##########################################################
  
  betaYG <- betaY_OB - betaY_RCT
  betaXY <- betaYG / betaXG

  ## 观察性中心两个斜率的方差及其协方差（影响函数），
  ## 与共享 RCT 斜率的方差一起进入 delta 方法
  ob_if <- ob_slope_if(data_OB, IV_name)
  varY_OB  <- sum(ob_if$IF_y[, 1]^2)
  varX_OB  <- sum(ob_if$IF_x[, 1]^2)
  covYX_OB <- sum(ob_if$IF_y[, 1] * ob_if$IF_x[, 1])

  SEYG <- sqrt(varY_OB + varY_RCT)

  SEXY <- ratio_se_delta(
    betaYG = betaYG,
    betaXG = betaXG,
    SEYG = SEYG,
    SEXG = sqrt(varX_OB),
    cov_YX = covYX_OB
  )
  
  if (!is.finite(betaXY) ||
      !is.finite(SEXY) ||
      SEXY <= 0) {
    
    return(empty_result(
      IV_name = IV_name,
      betaXG = betaXG,
      SEXG = SEXG,
      F_XZ = F_XZ,
      p_XZ = p_XZ,
      n_OB = n_OB
    ))
  }
  
  ##########################################################
  ## 6. 返回单中心结果
  ##########################################################
  
  return(data.frame(
    b = betaXY,
    se = SEXY,
    
    betaXG = betaXG,
    SEXG = SEXG,
    F_XZ = F_XZ,
    p_XZ = p_XZ,
    
    betaY_OB = betaY_OB,
    SEY_OB = SEY_OB,
    
    betaY_RCT = betaY_RCT,
    SEY_RCT = SEY_RCT,
    
    betaYG = betaYG,
    SEYG = SEYG,
    
    IV_name = IV_name,
    n_OB = n_OB,
    
    stringsAsFactors = FALSE
  ))
}




############################################################
## 单工具共享 RCT 协方差矩阵 Meta
##
## 用于：
##   ALLCONIV_cov()
##
## 方法：
##   1. 用固定效应 GLS 得到合并估计和 SE，与 meta_function_MR_cov_multi() 相同；
##   2. p 值和 95% CI 用 t 分布；
##   3. 自由度 df = 有效中心数 K - 1；
##   4. 与多工具 meta_function_MR_cov_multi() 保持一致。
############################################################

meta_function_MR_cov <- function(MR_result_all,
                                 rct_vcov,
                                 use_nearPD = TRUE,
                                 use_t_correction = TRUE) {
  
  if (!is.data.frame(MR_result_all)) {
    MR_result_all <- as.data.frame(MR_result_all)
  }
  
  required_cols <- c(
    "b",
    "se",
    "betaXG",
    "IV_name"
  )
  
  if (!all(required_cols %in% colnames(MR_result_all))) {
    stop("MR_result_all must contain b, se, betaXG, and IV_name.")
  }
  
  ##########################################################
  ## 1. 过滤无效中心
  ##########################################################
  
  MR_result_all <- MR_result_all[
    is.finite(MR_result_all$b) &
      is.finite(MR_result_all$se) &
      is.finite(MR_result_all$betaXG) &
      MR_result_all$se > 0 &
      abs(MR_result_all$betaXG) > 1e-8 &
      !is.na(MR_result_all$IV_name) &
      MR_result_all$IV_name %in% rownames(rct_vcov) &
      MR_result_all$IV_name %in% colnames(rct_vcov),
    ,
    drop = FALSE
  ]
  
  if (nrow(MR_result_all) < 2) {
    return(list(
      final_result = NULL,
      center_result = MR_result_all,
      V = NULL,
      df_meta = NA_real_,
      use_t_correction = use_t_correction
    ))
  }
  
  K <- nrow(MR_result_all)
  
  ##########################################################
  ## 2. 构造协方差矩阵 V
  ##########################################################
  
  V <- matrix(
    0,
    nrow = K,
    ncol = K
  )
  
  diag(V) <- MR_result_all$se^2
  
  for (i in 1:K) {
    for (j in 1:K) {
      
      if (i != j) {
        
        Zi <- MR_result_all$IV_name[i]
        Zj <- MR_result_all$IV_name[j]
        
        cov_rct_ij <- rct_vcov[Zi, Zj]
        
        V[i, j] <- cov_rct_ij /
          (MR_result_all$betaXG[i] * MR_result_all$betaXG[j])
      }
    }
  }
  
  ##########################################################
  ## 3. 数值保护：如果 V 非正定，做 nearPD 修正
  ##########################################################
  
  if (use_nearPD) {
    
    eig_min <- try(
      min(
        eigen(
          V,
          symmetric = TRUE,
          only.values = TRUE
        )$values
      ),
      silent = TRUE
    )
    
    if (inherits(eig_min, "try-error") ||
        !is.finite(eig_min) ||
        eig_min <= 1e-8) {
      
      if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("V is not positive definite. Please install Matrix package.")
      }
      
      V <- as.matrix(
        Matrix::nearPD(
          V,
          corr = FALSE
        )$mat
      )
    }
  }
  
  ##########################################################
  ## 4. 固定效应 GLS 合并（与 meta_function_MR_cov_multi() 相同）
  ##
  ## beta_meta = (1' V^{-1} 1)^(-1) 1' V^{-1} b
  ##########################################################

  one_vec <- rep(1, K)

  V_inv <- try(solve(V), silent = TRUE)

  if (inherits(V_inv, "try-error")) {

    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("V is singular. Please install MASS package.")
    }

    V_inv <- MASS::ginv(V)
  }

  denom <- as.numeric(
    t(one_vec) %*% V_inv %*% one_vec
  )

  if (!is.finite(denom) || denom <= 0) {
    return(list(
      final_result = NULL,
      center_result = MR_result_all,
      V = V,
      df_meta = NA_real_,
      use_t_correction = use_t_correction
    ))
  }

  estimate <- as.numeric(
    t(one_vec) %*% V_inv %*% MR_result_all$b / denom
  )

  meta_se <- sqrt(1 / denom)
  
  if (!is.finite(estimate) ||
      !is.finite(meta_se) ||
      meta_se <= 0) {
    
    return(list(
      final_result = NULL,
      center_result = MR_result_all,
      V = V,
      df_meta = NA_real_,
      use_t_correction = use_t_correction
    ))
  }
  
  stat_value <- estimate / meta_se
  
  ##########################################################
  ## 6. t 分布修正
  ##########################################################
  
  if (use_t_correction) {
    
    df_meta <- K - 1
    
    p_value <- 2 * pt(
      -abs(stat_value),
      df = df_meta
    )
    
    ci_lower <- estimate - qt(
      0.975,
      df = df_meta
    ) * meta_se
    
    ci_upper <- estimate + qt(
      0.975,
      df = df_meta
    ) * meta_se
    
  } else {
    
    df_meta <- Inf
    
    p_value <- 2 * pnorm(
      -abs(stat_value)
    )
    
    ci_lower <- estimate - 1.96 * meta_se
    ci_upper <- estimate + 1.96 * meta_se
  }
  
  ##########################################################
  ## 7. 最终结果
  ##########################################################
  
  final_result <- c(
    estimate = estimate,
    se = meta_se,
    z_value = stat_value,
    p_value = p_value,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
  
  return(list(
    final_result = final_result,
    center_result = MR_result_all,
    V = V,
    df_meta = df_meta,
    use_t_correction = use_t_correction
  ))
}





## 多工具中心内 IVW 估计
## 返回多工具协方差 Meta 所需的信息


CON_IV_luo_multi_return <- function(data_OB,
                                    IV_covariates,
                                    rct_info,
                                    center_id = NA_integer_,
                                    min_abs_betaXG = 1e-8) {
  
  empty_result <- function(IV_covariates = character(0)) {
    return(list(
      valid = FALSE,
      b = NA_real_,
      se = NA_real_,
      iv_names = character(0),
      a_vec = numeric(0),
      betaXG = numeric(0),
      betaYG = numeric(0),
      center_id = center_id,
      iv_count = 0,
      obs_y_var = NA_real_,
      x_var = NA_real_,
      rct_var = NA_real_,
      yx_cov = NA_real_
    ))
  }
  
  
  ## 0. 基本检查
  
  if (length(IV_covariates) == 0) {
    return(empty_result())
  }
  
  IV_covariates <- unique(IV_covariates)
  IV_covariates <- IV_covariates[
    !is.na(IV_covariates) &
      IV_covariates %in% colnames(data_OB) &
      IV_covariates %in% names(rct_info$beta) &
      IV_covariates %in% rownames(rct_info$vcov) &
      IV_covariates %in% colnames(rct_info$vcov)
  ]
  
  if (length(IV_covariates) == 0) {
    return(empty_result())
  }
  
  
  ## 1. 对每个工具分别估计：
  ##    X ~ Z
  ##    Y_OB ~ Z
  ##    RCT 中 Y_RCT ~ Z 的结果来自 rct_info
  
  iv_result <- data.frame(
    IV_name = character(0),
    betaXG = numeric(0),
    SEXG = numeric(0),
    betaY_OB = numeric(0),
    SEY_OB = numeric(0),
    betaY_RCT = numeric(0),
    varY_RCT = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (iv_j in IV_covariates) {
    
    
    ## X ~ Z
    
    tmp_x <- data.frame(
      X = data_OB$X,
      Z = data_OB[[iv_j]]
    )
    
    fit_x <- try(lm(X ~ Z, data = tmp_x), silent = TRUE)
    
    if (inherits(fit_x, "try-error")) {
      next
    }
    
    sx <- try(summary(fit_x)$coefficients, silent = TRUE)
    
    if (inherits(sx, "try-error") || !"Z" %in% rownames(sx)) {
      next
    }
    
    betaXG <- sx["Z", "Estimate"]
    SEXG <- sx["Z", "Std. Error"]
    
    ## Y_OB ~ Z
    
    tmp_y <- data.frame(
      Y = data_OB$Y,
      Z = data_OB[[iv_j]]
    )
    
    fit_y <- try(lm(Y ~ Z, data = tmp_y), silent = TRUE)
    
    if (inherits(fit_y, "try-error")) {
      next
    }
    
    sy <- try(summary(fit_y)$coefficients, silent = TRUE)
    
    if (inherits(sy, "try-error") || !"Z" %in% rownames(sy)) {
      next
    }
    
    betaY_OB <- sy["Z", "Estimate"]
    SEY_OB <- sy["Z", "Std. Error"]
    
    
    ## RCT 中 Y ~ Z
    
    betaY_RCT <- rct_info$beta[iv_j]
    varY_RCT <- rct_info$vcov[iv_j, iv_j]
    
    
    ## 保存
    
    iv_result <- rbind(
      iv_result,
      data.frame(
        IV_name = iv_j,
        betaXG = betaXG,
        SEXG = SEXG,
        betaY_OB = betaY_OB,
        SEY_OB = SEY_OB,
        betaY_RCT = betaY_RCT,
        varY_RCT = varY_RCT,
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (nrow(iv_result) == 0) {
    return(empty_result())
  }
  
  
  ## 2. 过滤无效工具
  
  
  keep <- is.finite(iv_result$betaXG) &
    is.finite(iv_result$SEXG) &
    is.finite(iv_result$betaY_OB) &
    is.finite(iv_result$SEY_OB) &
    is.finite(iv_result$betaY_RCT) &
    is.finite(iv_result$varY_RCT) &
    abs(iv_result$betaXG) > min_abs_betaXG &
    iv_result$SEXG > 0 &
    iv_result$SEY_OB > 0 &
    iv_result$varY_RCT >= 0
  
  iv_result <- iv_result[keep, , drop = FALSE]
  
  if (nrow(iv_result) == 0) {
    return(empty_result())
  }
  
  
  ## 3. 多工具 IVW 估计
  ##
  ## betaYG_j = betaY_OB_j - betaY_RCT_j
  ##
  ## betaXY = sum(w_j * betaXG_j * betaYG_j) /
  ##          sum(w_j * betaXG_j^2)
  
  
  betaXG <- iv_result$betaXG
  SEXG <- iv_result$SEXG
  
  betaYG <- iv_result$betaY_OB - iv_result$betaY_RCT
  
  SEYG2_diag <- iv_result$SEY_OB^2 + iv_result$varY_RCT
  
  if (any(!is.finite(SEYG2_diag)) || any(SEYG2_diag <= 0)) {
    return(empty_result())
  }
  
  w <- 1 / SEYG2_diag
  
  denom <- sum(w * betaXG^2)
  
  if (!is.finite(denom) || denom <= 0) {
    return(empty_result())
  }
  
  betaXY <- sum(w * betaXG * betaYG) / denom
  
  
  ## 4. 计算该中心估计量对 betaYG_j 的线性权重
  ##
  ## theta_i = sum_j a_ij * betaYG_ij
  ##
  ## a_j = w_j * betaXG_j / sum(w_j * betaXG_j^2)
  
  
  a_vec <- w * betaXG / denom
  names(a_vec) <- iv_result$IV_name
  
  
  ## 5. 方差
  ##
  ## obs_y_var: 观察性中心 Y~Z 斜率，含不同工具之间的协方差
  ## rct_var:   共享 RCT 中 Y~Z 斜率，含不同工具之间的协方差
  ## x_var:     观察性中心 X~Z 斜率，含不同工具之间的协方差
  ## yx_cov:    同一观察性中心中 Y~Z 与 X~Z 斜率之间的协方差

  ob_if <- ob_slope_if(data_OB, iv_result$IV_name)

  if (any(!is.finite(ob_if$IF_x)) || any(!is.finite(ob_if$IF_y))) {
    return(empty_result())
  }

  V_yy <- crossprod(ob_if$IF_y)
  V_xx <- crossprod(ob_if$IF_x)
  V_yx <- crossprod(ob_if$IF_y, ob_if$IF_x)

  obs_y_var <- as.numeric(t(a_vec) %*% V_yy %*% a_vec)

  rct_sub_vcov <- rct_info$vcov[
    iv_result$IV_name,
    iv_result$IV_name,
    drop = FALSE
  ]

  rct_var <- as.numeric(
    t(a_vec) %*% rct_sub_vcov %*% a_vec
  )

  ## betaXY 对 betaXG_j 的导数
  d_x <- w * (betaYG - 2 * betaXG * betaXY) / denom

  x_var <- as.numeric(t(d_x) %*% V_xx %*% d_x)
  yx_cov <- 2 * as.numeric(t(a_vec) %*% V_yx %*% d_x)

  total_var <- obs_y_var + rct_var + x_var + yx_cov
  
  if (!is.finite(total_var) || total_var < 0) {
    return(empty_result())
  }
  
  seXY <- sqrt(total_var)
  
  return(list(
    valid = TRUE,
    b = betaXY,
    se = seXY,
    iv_names = iv_result$IV_name,
    a_vec = a_vec,
    betaXG = betaXG,
    betaYG = betaYG,
    center_id = center_id,
    iv_count = length(iv_result$IV_name),
    obs_y_var = obs_y_var,
    x_var = x_var,
    rct_var = rct_var,
    yx_cov = yx_cov
  ))
}


############################################################
## 多工具共享 RCT 协方差矩阵 Meta
############################################################

############################################################
## 多工具共享 RCT 协方差矩阵 Meta
############################################################

meta_function_MR_cov_multi <- function(center_result_list,
                                       rct_vcov,
                                       use_nearPD = TRUE,
                                       use_t_correction = TRUE) {
  
  ##########################################################
  ## 1. 提取有效中心
  ##########################################################
  
  valid_idx <- which(sapply(center_result_list, function(x) {
    is.list(x) &&
      isTRUE(x$valid) &&
      is.finite(x$b) &&
      is.finite(x$se) &&
      x$se > 0 &&
      length(x$iv_names) > 0 &&
      length(x$a_vec) > 0
  }))
  
  if (length(valid_idx) < 2) {
    return(list(
      final_result = NULL,
      center_result = NULL,
      V = NULL,
      valid_idx = valid_idx
    ))
  }
  
  valid_results <- center_result_list[valid_idx]
  
  K <- length(valid_results)
  
  b <- sapply(valid_results, function(x) x$b)
  se <- sapply(valid_results, function(x) x$se)
  iv_count <- sapply(valid_results, function(x) x$iv_count)
  center_original <- sapply(valid_results, function(x) x$center_id)
  
  ##########################################################
  ## 2. 构造中心间协方差矩阵 V
  ##########################################################
  
  V <- matrix(0, nrow = K, ncol = K)
  
  ## 对角线：每个中心自身估计方差
  diag(V) <- se^2
  
  ## 非对角线：共享 RCT 参照误差导致的中心间协方差
  for (i in 1:K) {
    for (j in 1:K) {
      
      if (i != j) {
        
        iv_i <- valid_results[[i]]$iv_names
        iv_j <- valid_results[[j]]$iv_names
        
        a_i <- valid_results[[i]]$a_vec
        a_j <- valid_results[[j]]$a_vec
        
        common_valid <- all(iv_i %in% rownames(rct_vcov)) &&
          all(iv_j %in% colnames(rct_vcov))
        
        if (!common_valid) {
          
          V[i, j] <- 0
          
        } else {
          
          rct_sub <- rct_vcov[
            iv_i,
            iv_j,
            drop = FALSE
          ]
          
          V[i, j] <- as.numeric(
            t(a_i) %*% rct_sub %*% a_j
          )
        }
      }
    }
  }
  
  ##########################################################
  ## 3. 数值保护：如果 V 非正定，做 nearPD 修正
  ##########################################################
  
  if (use_nearPD) {
    
    eig_min <- try(
      min(eigen(V, symmetric = TRUE, only.values = TRUE)$values),
      silent = TRUE
    )
    
    if (inherits(eig_min, "try-error") ||
        !is.finite(eig_min) ||
        eig_min <= 1e-8) {
      
      if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("V is not positive definite. Please install Matrix package.")
      }
      
      V <- as.matrix(Matrix::nearPD(V, corr = FALSE)$mat)
    }
  }
  
  ##########################################################
  ## 4. 固定效应 GLS 合并
  ##
  ## beta_meta = (1' V^{-1} 1)^(-1) 1' V^{-1} b
  ##########################################################
  
  one_vec <- rep(1, K)
  
  V_inv <- try(solve(V), silent = TRUE)
  
  if (inherits(V_inv, "try-error")) {
    
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("V is singular. Please install MASS package.")
    }
    
    V_inv <- MASS::ginv(V)
  }
  
  denom <- as.numeric(
    t(one_vec) %*% V_inv %*% one_vec
  )
  
  if (!is.finite(denom) || denom <= 0) {
    return(list(
      final_result = NULL,
      center_result = NULL,
      V = V,
      valid_idx = valid_idx
    ))
  }
  
  estimate <- as.numeric(
    t(one_vec) %*% V_inv %*% b / denom
  )
  
  meta_se <- sqrt(1 / denom)
  
  if (!is.finite(meta_se) || meta_se <= 0) {
    return(list(
      final_result = NULL,
      center_result = NULL,
      V = V,
      valid_idx = valid_idx
    ))
  }
  
  ##########################################################
  ## 5. 小样本 t 分布修正
  ##
  ##   df = K - 1
  ##   p_value = 2 * pt(-abs(t), df)
  ##   CI = beta ± qt(0.975, df) * se
  ##########################################################
  
  stat_value <- estimate / meta_se
  
  if (use_t_correction) {
    
    df_meta <- K - 1
    
    p_value <- 2 * pt(
      -abs(stat_value),
      df = df_meta
    )
    
    ci_lower <- estimate - qt(
      0.975,
      df = df_meta
    ) * meta_se
    
    ci_upper <- estimate + qt(
      0.975,
      df = df_meta
    ) * meta_se
    
  } else {
    
    df_meta <- Inf
    
    p_value <- 2 * pnorm(
      -abs(stat_value)
    )
    
    ci_lower <- estimate - 1.96 * meta_se
    ci_upper <- estimate + 1.96 * meta_se
  }
  
  ##########################################################
  ## 6. 最终结果
  ##
  ## 第三个元素命名为 z_value；
  ## use_t_correction = TRUE 时它是 t 统计量。
  ##########################################################
  
  final_result <- c(
    estimate = estimate,
    se = meta_se,
    z_value = stat_value,
    p_value = p_value,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
  
  ##########################################################
  ## 7. 整理中心结果
  ##########################################################
  
  center_result <- data.frame(
    center_original = center_original,
    b = b,
    se = se,
    iv_count = iv_count,
    iv_set = sapply(valid_results, function(x) {
      paste(x$iv_names, collapse = ";")
    }),
    stringsAsFactors = FALSE
  )
  
  return(list(
    final_result = final_result,
    center_result = center_result,
    V = V,
    valid_idx = valid_idx,
    center_result_raw = valid_results,
    df_meta = df_meta,
    use_t_correction = use_t_correction
  ))
}



############################################################
## 多中心多工具 + 共享 RCT 协方差矩阵 Meta
############################################################

ALLCONIV_cov_multi <- function(ALLdata,
                               N_center,
                               type,
                               covariates,
                               n_boot = NULL) {
  
  ##########################################################
  ## 0. 当前版本先只支持连续型 type = 1
  ##########################################################
  
  if (type != 1) {
    stop("ALLCONIV_cov_multi currently supports type == 1 only.")
  }
  
  ##########################################################
  ## 1. 收集所有中心使用的所有工具变量
  ##########################################################
  
  all_selected_IVs <- character(0)
  
  for (cen in 1:N_center) {
    
    IV_cen <- covariates[[paste0("center_", cen)]]
    
    if (length(IV_cen) > 0) {
      all_selected_IVs <- c(all_selected_IVs, IV_cen)
    }
  }
  
  all_selected_IVs <- unique(all_selected_IVs)
  all_selected_IVs <- all_selected_IVs[
    !is.na(all_selected_IVs) &
      all_selected_IVs %in% colnames(ALLdata$rct)
  ]
  
  if (length(all_selected_IVs) == 0) {
    return(list(
      final_result = NULL,
      center_result = NULL,
      V = NULL
    ))
  }
  
  ##########################################################
  ## 2. 在 RCT 中估计所有工具的 Y~Z 斜率协方差矩阵
  ##########################################################
  
  rct_info <- rct_slope_vcov(
    data_rct = ALLdata$rct,
    iv_names = all_selected_IVs
  )
  
  ##########################################################
  ## 3. 每个中心做多工具 IVW 估计
  ##########################################################
  
  center_result_list <- vector("list", N_center)
  
  for (cen in 1:N_center) {
    
    IV_cen <- covariates[[paste0("center_", cen)]]
    
    data_OB_cen <- ALLdata$all_data[
      ALLdata$all_data$centre == cen,
      ,
      drop = FALSE
    ]
    
    center_result_list[[cen]] <- CON_IV_luo_multi_return(
      data_OB = data_OB_cen,
      IV_covariates = IV_cen,
      rct_info = rct_info,
      center_id = cen
    )
  }
  
  ##########################################################
  ## 4. 多工具协方差矩阵 Meta
  ##########################################################
  
  final_meta <- meta_function_MR_cov_multi(
    center_result_list = center_result_list,
    rct_vcov = rct_info$vcov
  )
  
  final_meta$rct_info <- rct_info
  final_meta$all_selected_IVs <- all_selected_IVs
  
  return(final_meta)
}




############################################################
## 1. 模拟结果汇总函数
############################################################
############################################################
## 模拟结果汇总函数
##
## 汇总：
##   estimate
##   bias
##   emp_se
##   mean_model_se
##   se_ratio
##   mse
##   coverage
##   valid_rate
##   n_valid
##
## 另外汇总工具变量个数：
##   selected_mean_iv_count
##   selected_median_iv_count
##   selected_min_iv_count
##   selected_max_iv_count
##
##   used_mean_iv_count
##   used_median_iv_count
##   used_min_iv_count
##   used_max_iv_count
############################################################

resultt_FUn <- function(result, true_effect) {
  
  empty_result <- c(
    estimate = NA_real_,
    bias = NA_real_,
    emp_se = NA_real_,
    mean_model_se = NA_real_,
    se_ratio = NA_real_,
    mse = NA_real_,
    coverage = NA_real_,
    valid_rate = 0,
    n_valid = 0,
    
    selected_mean_iv_count = NA_real_,
    selected_median_iv_count = NA_real_,
    selected_min_iv_count = NA_real_,
    selected_max_iv_count = NA_real_,
    
    used_mean_iv_count = NA_real_,
    used_median_iv_count = NA_real_,
    used_min_iv_count = NA_real_,
    used_max_iv_count = NA_real_
  )
  
  if (is.null(result)) {
    return(empty_result)
  }
  
  result <- as.data.frame(result)
  
  if (nrow(result) == 0 || ncol(result) < 2) {
    return(empty_result)
  }
  
  ##########################################################
  ## 1. 提取估计值和模型标准误
  ##########################################################
  
  est <- as.numeric(result[, 1])
  model_se <- as.numeric(result[, 2])
  
  valid <- is.finite(est)
  
  if (sum(valid, na.rm = TRUE) == 0) {
    
    ## 即使没有有效估计，也尽量保留工具筛选个数信息
    diag_names <- names(empty_result)[
      grepl("iv_count", names(empty_result))
    ]
    
    for (dn in diag_names) {
      if (dn %in% colnames(result)) {
        tmp_val <- suppressWarnings(
          mean(as.numeric(result[[dn]]), na.rm = TRUE)
        )
        if (is.finite(tmp_val)) {
          empty_result[dn] <- tmp_val
        }
      }
    }
    
    return(empty_result)
  }
  
  est_valid <- est[valid]
  
  ##########################################################
  ## 2. 估计值、Bias、经验 SE、MSE
  ##########################################################
  
  estimate <- mean(est_valid, na.rm = TRUE)
  bias <- estimate - true_effect
  emp_se <- sd(est_valid, na.rm = TRUE)
  
  mse <- mean(
    (est_valid - true_effect)^2,
    na.rm = TRUE
  )
  
  ##########################################################
  ## 3. 模型 SE 与经验 SE 比值
  ##########################################################
  
  valid_se <- valid &
    is.finite(model_se) &
    model_se > 0
  
  if (sum(valid_se, na.rm = TRUE) > 0) {
    mean_model_se <- mean(
      model_se[valid_se],
      na.rm = TRUE
    )
  } else {
    mean_model_se <- NA_real_
  }
  
  if (is.finite(emp_se) &&
      emp_se > 0 &&
      is.finite(mean_model_se)) {
    
    se_ratio <- mean_model_se / emp_se
    
  } else {
    
    se_ratio <- NA_real_
  }
  
  ##########################################################
  ## 4. Coverage
  ##########################################################
  
  if (ncol(result) >= 6) {
    
    ci_lower <- as.numeric(result[, 5])
    ci_upper <- as.numeric(result[, 6])
    
    coverage <- mean(
      is.finite(ci_lower) &
        is.finite(ci_upper) &
        ci_lower <= true_effect &
        ci_upper >= true_effect,
      na.rm = TRUE
    )
    
  } else {
    
    coverage <- NA_real_
  }
  
  ##########################################################
  ## 5. 有效率
  ##########################################################
  
  valid_rate <- mean(valid, na.rm = TRUE)
  n_valid <- sum(valid, na.rm = TRUE)
  
  ##########################################################
  ## 6. 工具变量个数汇总
  ##########################################################
  
  get_mean_col <- function(colname) {
    
    if (!colname %in% colnames(result)) {
      return(NA_real_)
    }
    
    val <- suppressWarnings(
      mean(as.numeric(result[[colname]]), na.rm = TRUE)
    )
    
    if (!is.finite(val)) {
      return(NA_real_)
    } else {
      return(val)
    }
  }
  
  selected_mean_iv_count <- get_mean_col("selected_mean_iv_count")
  selected_median_iv_count <- get_mean_col("selected_median_iv_count")
  selected_min_iv_count <- get_mean_col("selected_min_iv_count")
  selected_max_iv_count <- get_mean_col("selected_max_iv_count")
  
  used_mean_iv_count <- get_mean_col("used_mean_iv_count")
  used_median_iv_count <- get_mean_col("used_median_iv_count")
  used_min_iv_count <- get_mean_col("used_min_iv_count")
  used_max_iv_count <- get_mean_col("used_max_iv_count")
  
  ##########################################################
  ## 7. 返回
  ##########################################################
  
  return(c(
    estimate = estimate,
    bias = bias,
    emp_se = emp_se,
    mean_model_se = mean_model_se,
    se_ratio = se_ratio,
    mse = mse,
    coverage = coverage,
    valid_rate = valid_rate,
    n_valid = n_valid,
    
    selected_mean_iv_count = selected_mean_iv_count,
    selected_median_iv_count = selected_median_iv_count,
    selected_min_iv_count = selected_min_iv_count,
    selected_max_iv_count = selected_max_iv_count,
    
    used_mean_iv_count = used_mean_iv_count,
    used_median_iv_count = used_median_iv_count,
    used_min_iv_count = used_min_iv_count,
    used_max_iv_count = used_max_iv_count
  ))
}

############################################################
## 最大样本量中心筛选统一工具变量
## 所有中心使用同一组工具
##
## 要求：
##   1. 工具变量在最大样本量中心中筛选；
##   2. 工具变量必须在所有观察性中心中存在；
##   3. 工具变量必须在 RCT 数据中存在；
##   4. 可用于单工具或多工具。
############################################################

variable_selection_largest_center_common <- function(N_center,
                                                     covariates,
                                                     data,
                                                     correlation_threshold,
                                                     type = 1,
                                                     selcetnum = "multi") {
  
  ##########################################################
  ## 0. 基本检查
  ##########################################################
  
  if (length(covariates) == 0) {
    stop("covariates is empty.")
  }
  
  if (correlation_threshold < 0 || correlation_threshold > 1) {
    stop("correlation_threshold must be between 0 and 1.")
  }
  
  if (!(type %in% c(1, 2))) {
    stop("type must be 1 or 2.")
  }
  
  if (!(is.numeric(selcetnum) || selcetnum == "multi")) {
    stop("selcetnum must be 1, a positive integer, or 'multi'.")
  }
  
  if (is.numeric(selcetnum)) {
    if (length(selcetnum) != 1 ||
        selcetnum < 1 ||
        selcetnum != floor(selcetnum)) {
      stop("numeric selcetnum must be a single positive integer.")
    }
  }
  
  ##########################################################
  ## 1. 找到样本量最大的观察性中心
  ##########################################################
  
  center_sizes <- sapply(1:N_center, function(cen) {
    sum(data$all_data$centre == cen)
  })
  
  largest_center <- which.max(center_sizes)
  
  data_largest <- data$all_data[
    data$all_data$centre == largest_center,
    ,
    drop = FALSE
  ]
  
  ##########################################################
  ## 2. 只保留所有观察性中心和 RCT 中都能观测到的候选变量
  ##########################################################
  
  center_col_list <- lapply(1:N_center, function(cen) {
    colnames(
      data$all_data[
        data$all_data$centre == cen,
        ,
        drop = FALSE
      ]
    )
  })
  
  common_center_cols <- Reduce(intersect, center_col_list)
  
  covariates_common <- covariates[
    covariates %in% common_center_cols &
      covariates %in% colnames(data$rct) &
      covariates %in% colnames(data_largest)
  ]
  
  if (length(covariates_common) == 0) {
    
    all_final_vars <- list()
    
    for (cen in 1:N_center) {
      all_final_vars[[paste0("center_", cen)]] <- character(0)
    }
    
    attr(all_final_vars, "largest_center") <- largest_center
    attr(all_final_vars, "largest_center_size") <- center_sizes[largest_center]
    attr(all_final_vars, "selected_vars") <- character(0)
    attr(all_final_vars, "selection_method") <- "largest_center_common"
    
    return(all_final_vars)
  }
  
  ##########################################################
  ## 3. 设置筛选 P 值阈值
  ##########################################################
  
  if (type == 1) {
    p_cutoff <- 0.05 / length(covariates_common)
  } else {
    p_cutoff <- 0.1
  }
  
  ##########################################################
  ## 4. 在最大样本量中心中进行工具筛选
  ##########################################################
  
  p_value_Y <- rep(NA_real_, length(covariates_common))
  p_value_X <- rep(NA_real_, length(covariates_common))
  
  for (j in seq_along(covariates_common)) {
    
    var_j <- covariates_common[j]
    
    formula_Y <- reformulate(
      termlabels = c("X", var_j),
      response = "Y"
    )
    
    formula_X <- reformulate(
      termlabels = var_j,
      response = "X"
    )
    
    p_value_Y[j] <- safe_lm_pvalue(
      formula = formula_Y,
      data = data_largest,
      term_name = var_j
    )
    
    p_value_X[j] <- safe_lm_pvalue(
      formula = formula_X,
      data = data_largest,
      term_name = var_j
    )
  }
  
  variable_Y <- covariates_common[
    which(!is.na(p_value_Y) & p_value_Y < p_cutoff)
  ]
  
  variable_X <- covariates_common[
    which(!is.na(p_value_X) & p_value_X < p_cutoff)
  ]
  
  common_variables <- intersect(variable_Y, variable_X)
  
  if (length(common_variables) == 0) {
    
    selected_vars <- character(0)
    
  } else {
    
    selected_vars <- common_variables
    
    ########################################################
    ## 5. 去除高度相关的冗余工具变量
    ##    如果两个变量高度相关，保留与 X 相关性更强的变量
    ########################################################
    
    if (length(common_variables) > 1) {
      
      selectvar <- data_largest[, common_variables, drop = FALSE]
      
      cor_matrix <- suppressWarnings(
        abs(cor(selectvar, use = "pairwise.complete.obs"))
      )
      
      vars_to_remove <- character(0)
      
      for (idx_common_var_1 in 1:(length(common_variables) - 1)) {
        for (idx_common_var_2 in (idx_common_var_1 + 1):length(common_variables)) {
          
          var1_name <- common_variables[idx_common_var_1]
          var2_name <- common_variables[idx_common_var_2]
          
          if (var1_name %in% vars_to_remove ||
              var2_name %in% vars_to_remove) {
            next
          }
          
          cor_val <- cor_matrix[var1_name, var2_name]
          
          if (!is.na(cor_val) && cor_val > correlation_threshold) {
            
            cor1_with_X <- safe_abs_cor(data_largest$X, data_largest[[var1_name]])
            cor2_with_X <- safe_abs_cor(data_largest$X, data_largest[[var2_name]])
            
            if (cor1_with_X >= cor2_with_X) {
              vars_to_remove <- c(vars_to_remove, var2_name)
            } else {
              vars_to_remove <- c(vars_to_remove, var1_name)
            }
          }
        }
      }
      
      selected_vars <- setdiff(common_variables, vars_to_remove)
    }
    
    ########################################################
    ## 6. 按与 X 的绝对相关性从大到小排序
    ########################################################
    
    if (length(selected_vars) > 1) {
      
      correlations_with_X <- sapply(selected_vars, function(var) {
        safe_abs_cor(data_largest$X, data_largest[[var]])
      })
      
      selected_vars <- selected_vars[
        order(correlations_with_X, decreasing = TRUE)
      ]
    }
    
    ########################################################
    ## 7. 再次保证工具在所有中心和 RCT 中均存在
    ########################################################
    
    selected_vars <- selected_vars[
      selected_vars %in% covariates_common
    ]
    
    ########################################################
    ## 8. 根据 selcetnum 控制单工具或多工具
    ########################################################
    
    if (length(selected_vars) > 0) {
      
      if (is.numeric(selcetnum)) {
        
        if (selcetnum == 1) {
          
          selected_vars <- selected_vars[1]
          
        } else if (selcetnum > 1) {
          
          if (length(selected_vars) <= selcetnum) {
            selected_vars <- selected_vars[1]
          } else {
            selected_vars <- selected_vars
          }
        }
        
      } else if (selcetnum == "multi") {
        
        selected_vars <- selected_vars
      }
    }
  }
  
  ##########################################################
  ## 9. 所有中心使用同一组工具变量
  ##########################################################
  
  all_final_vars <- list()
  
  for (cen in 1:N_center) {
    all_final_vars[[paste0("center_", cen)]] <- selected_vars
  }
  
  attr(all_final_vars, "largest_center") <- largest_center
  attr(all_final_vars, "largest_center_size") <- center_sizes[largest_center]
  attr(all_final_vars, "selected_vars") <- selected_vars
  attr(all_final_vars, "selection_method") <- "largest_center_common"
  attr(all_final_vars, "common_candidate_vars") <- covariates_common
  
  return(all_final_vars)
}


############################################################
## 单工具 + 共享 RCT 协方差矩阵 Meta
##
## 适用于：
##   1. 单工具 + 同一工具
##   2. 单工具 + 不同中心不同工具
##
## 注意：
##   如果 covariates 中某中心有多个候选工具，
##   本函数只取第 1 个，即最强工具。
############################################################

ALLCONIV_cov <- function(ALLdata,
                         N_center,
                         type,
                         covariates,
                         n_boot = NULL) {
  
  ##########################################################
  ## 0. 当前版本只支持连续型 type = 1
  ##########################################################
  
  if (type != 1) {
    stop("ALLCONIV_cov currently supports type == 1 only.")
  }
  
  ##########################################################
  ## 1. 提取变量筛选信息
  ##########################################################
  
  largest_center <- attr(covariates, "largest_center")
  largest_center_size <- attr(covariates, "largest_center_size")
  selected_from <- attr(covariates, "selected_from")
  selected_vars <- attr(covariates, "selected_vars")
  
  if (is.null(largest_center)) {
    largest_center <- NA_integer_
  }
  
  if (is.null(largest_center_size)) {
    largest_center_size <- NA_integer_
  }
  
  if (is.null(selected_from)) {
    selected_from <- NA_character_
  }
  
  if (is.null(selected_vars)) {
    selected_vars <- character(0)
  }
  
  ##########################################################
  ## 2. 收集所有中心实际使用的工具变量
  ##    单工具版本：每个中心只取第一个工具
  ##########################################################
  
  center_iv_map <- data.frame(
    center_original = seq_len(N_center),
    IV_name = NA_character_,
    n_IV_candidates = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  all_selected_IVs <- character(0)
  
  for (cen in 1:N_center) {
    
    IV_cen <- covariates[[paste0("center_", cen)]]
    
    if (is.null(IV_cen)) {
      IV_cen <- character(0)
    }
    
    center_iv_map$n_IV_candidates[cen] <- length(IV_cen)
    
    if (length(IV_cen) > 0) {
      
      IV_top1 <- IV_cen[1]
      
      center_iv_map$IV_name[cen] <- IV_top1
      all_selected_IVs <- c(all_selected_IVs, IV_top1)
    }
  }
  
  all_selected_IVs <- unique(all_selected_IVs)
  
  all_selected_IVs <- all_selected_IVs[
    !is.na(all_selected_IVs) &
      all_selected_IVs %in% colnames(ALLdata$rct)
  ]
  
  ##########################################################
  ## 3. 如果没有有效工具变量，直接返回空结果
  ##########################################################
  
  if (length(all_selected_IVs) == 0) {
    
    iv_selection_info <- data.frame(
      largest_center = largest_center,
      largest_center_size = largest_center_size,
      selected_from = selected_from,
      
      n_selected_vars = length(selected_vars),
      selected_vars = paste(selected_vars, collapse = ";"),
      
      n_unique_IV = 0,
      unique_IVs_used = NA_character_,
      
      n_valid_center = 0,
      n_valid_unique_IV = 0,
      valid_unique_IVs_used = NA_character_,
      
      mean_F_XZ = NA_real_,
      median_F_XZ = NA_real_,
      min_F_XZ = NA_real_,
      n_weak_IV_F10 = NA_integer_,
      
      mean_iv_count = mean(ifelse(is.na(center_iv_map$IV_name), 0, 1)),
      median_iv_count = median(ifelse(is.na(center_iv_map$IV_name), 0, 1)),
      min_iv_count = min(ifelse(is.na(center_iv_map$IV_name), 0, 1)),
      max_iv_count = max(ifelse(is.na(center_iv_map$IV_name), 0, 1)),
      
      stringsAsFactors = FALSE
    )
    
    return(list(
      final_result = NULL,
      center_result = NULL,
      V = NULL,
      rct_info = NULL,
      iv_selection_info = iv_selection_info,
      center_iv_map = center_iv_map,
      all_selected_IVs = all_selected_IVs
    ))
  }
  
  ##########################################################
  ## 4. 在 RCT / do 数据中估计 Y ~ Z 斜率协方差矩阵
  ##########################################################
  
  rct_info <- try(
    rct_slope_vcov(
      data_rct = ALLdata$rct,
      iv_names = all_selected_IVs
    ),
    silent = TRUE
  )
  
  if (inherits(rct_info, "try-error")) {
    
    iv_selection_info <- data.frame(
      largest_center = largest_center,
      largest_center_size = largest_center_size,
      selected_from = selected_from,
      
      n_selected_vars = length(selected_vars),
      selected_vars = paste(selected_vars, collapse = ";"),
      
      n_unique_IV = length(all_selected_IVs),
      unique_IVs_used = paste(all_selected_IVs, collapse = ";"),
      
      n_valid_center = 0,
      n_valid_unique_IV = 0,
      valid_unique_IVs_used = NA_character_,
      
      mean_F_XZ = NA_real_,
      median_F_XZ = NA_real_,
      min_F_XZ = NA_real_,
      n_weak_IV_F10 = NA_integer_,
      
      mean_iv_count = mean(ifelse(is.na(center_iv_map$IV_name), 0, 1)),
      median_iv_count = median(ifelse(is.na(center_iv_map$IV_name), 0, 1)),
      min_iv_count = min(ifelse(is.na(center_iv_map$IV_name), 0, 1)),
      max_iv_count = max(ifelse(is.na(center_iv_map$IV_name), 0, 1)),
      
      stringsAsFactors = FALSE
    )
    
    return(list(
      final_result = NULL,
      center_result = NULL,
      V = NULL,
      rct_info = NULL,
      iv_selection_info = iv_selection_info,
      center_iv_map = center_iv_map,
      all_selected_IVs = all_selected_IVs
    ))
  }
  
  ##########################################################
  ## 5. 每个中心计算单工具 Wald ratio
  ##########################################################
  
  center_result_list <- vector("list", N_center)
  
  for (cen in 1:N_center) {
    
    IV_cen <- covariates[[paste0("center_", cen)]]
    
    if (is.null(IV_cen)) {
      IV_cen <- character(0)
    }
    
    data_OB_cen <- ALLdata$all_data[
      ALLdata$all_data$centre == cen,
      ,
      drop = FALSE
    ]
    
    center_result_list[[cen]] <- CON_IV_luo_top1_return(
      data_rct_all = ALLdata$rct,
      data_OB = data_OB_cen,
      IV_covariates = IV_cen,
      rct_info = rct_info
    )
    
    center_result_list[[cen]]$center_original <- cen
  }
  
  MR_result_all <- do.call(rbind, center_result_list)
  
  ##########################################################
  ## 6. 带共享 RCT 协方差矩阵的 Meta 合并
  ##########################################################
  
  final_meta <- meta_function_MR_cov(
    MR_result_all = MR_result_all,
    rct_vcov = rct_info$vcov
  )
  
  ##########################################################
  ## 7. 汇总工具变量诊断信息
  ##########################################################
  
  valid_center <- is.finite(MR_result_all$b) &
    is.finite(MR_result_all$se) &
    MR_result_all$se > 0
  
  if ("F_XZ" %in% colnames(MR_result_all)) {
    valid_F <- MR_result_all$F_XZ[
      is.finite(MR_result_all$F_XZ)
    ]
  } else {
    valid_F <- numeric(0)
  }
  
  IVs_used_valid <- unique(MR_result_all$IV_name[valid_center])
  IVs_used_valid <- IVs_used_valid[!is.na(IVs_used_valid)]
  
  iv_count_by_center <- ifelse(
    is.na(center_iv_map$IV_name),
    0,
    1
  )
  
  iv_selection_info <- data.frame(
    largest_center = largest_center,
    largest_center_size = largest_center_size,
    selected_from = selected_from,
    
    n_selected_vars = length(selected_vars),
    selected_vars = paste(selected_vars, collapse = ";"),
    
    n_unique_IV = length(unique(all_selected_IVs)),
    unique_IVs_used = paste(unique(all_selected_IVs), collapse = ";"),
    
    n_valid_center = sum(valid_center, na.rm = TRUE),
    n_valid_unique_IV = length(IVs_used_valid),
    valid_unique_IVs_used = paste(IVs_used_valid, collapse = ";"),
    
    mean_F_XZ = ifelse(length(valid_F) > 0, mean(valid_F, na.rm = TRUE), NA_real_),
    median_F_XZ = ifelse(length(valid_F) > 0, median(valid_F, na.rm = TRUE), NA_real_),
    min_F_XZ = ifelse(length(valid_F) > 0, min(valid_F, na.rm = TRUE), NA_real_),
    n_weak_IV_F10 = ifelse(length(valid_F) > 0, sum(valid_F < 10, na.rm = TRUE), NA_integer_),
    
    mean_iv_count = mean(iv_count_by_center, na.rm = TRUE),
    median_iv_count = median(iv_count_by_center, na.rm = TRUE),
    min_iv_count = min(iv_count_by_center, na.rm = TRUE),
    max_iv_count = max(iv_count_by_center, na.rm = TRUE),
    
    stringsAsFactors = FALSE
  )
  
  ##########################################################
  ## 8. 附加信息后返回
  ##########################################################
  
  final_meta$rct_info <- rct_info
  final_meta$iv_selection_info <- iv_selection_info
  final_meta$center_iv_map <- center_iv_map
  final_meta$all_selected_IVs <- all_selected_IVs
  
  return(final_meta)
}

##########3 other model###########


meta_function_MR <- function(MR_result_all) {
  
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Please install metafor first.")
  }
  
  if (is.null(MR_result_all)) {
    return(list(
      final_result = NULL,
      center_result = NULL,
      fit = NULL
    ))
  }
  
  if (!is.data.frame(MR_result_all)) {
    MR_result_all <- as.data.frame(MR_result_all)
  }
  
  if (ncol(MR_result_all) < 2) {
    return(list(
      final_result = NULL,
      center_result = MR_result_all,
      fit = NULL
    ))
  }
  
  colnames(MR_result_all)[1:2] <- c("b", "se")
  
  MR_result_all <- MR_result_all[
    is.finite(MR_result_all$b) &
      is.finite(MR_result_all$se) &
      MR_result_all$se > 0,
    ,
    drop = FALSE
  ]
  
  if (nrow(MR_result_all) < 2) {
    return(list(
      final_result = NULL,
      center_result = MR_result_all,
      fit = NULL
    ))
  }
  
  fit <- tryCatch({
    metafor::rma(
      yi = b,
      sei = se,
      data = MR_result_all,
      method = "REML",
      test = "knha"
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(fit)) {
    return(list(
      final_result = NULL,
      center_result = MR_result_all,
      fit = NULL
    ))
  }
  
  stat_value <- if (!is.null(fit$zval)) {
    fit$zval
  } else {
    fit$tval
  }
  
  final_result <- c(
    estimate = as.numeric(fit$beta),
    se = as.numeric(fit$se),
    z_value = as.numeric(stat_value),
    p_value = as.numeric(fit$pval),
    ci_lower = as.numeric(fit$ci.lb),
    ci_upper = as.numeric(fit$ci.ub)
  )
  
  return(list(
    final_result = final_result,
    center_result = MR_result_all,
    fit = fit
  ))
}

## 传统观察性方法变量筛选函数##############
## out = "both" : 同时与 X 和 Y 相关
## out = "X"    : 只筛与 X 相关的变量
## out = "Y"    : 只筛与 Y 相关的变量


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
    
    if (type == 1) {
      p_cutoff_Y <- 0.05 / length(covariates_i)
      p_cutoff_X <- 0.05 / length(covariates_i)
    } else {
      p_cutoff_Y <- 0.05
      p_cutoff_X <- 0.05
    }
    
    variable_Y <- covariates_i[
      !is.na(p_value_Y) & p_value_Y < p_cutoff_Y
    ]
    
    variable_X <- covariates_i[
      !is.na(p_value_X) & p_value_X < p_cutoff_X
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



## 传统观察性回归调整############
## 每个中心单独调整，然后 meta 合并############


ALL_adjust <- function(ALLdata,
                       N_center,
                       type,
                       covariates) {
  
  MR_result_all <- NULL
  
  for (cen in 1:N_center) {
    
    data_OB_i <- ALLdata$all_data[
      ALLdata$all_data$centre == cen,
      ,
      drop = FALSE
    ]
    
    cov_cen <- covariates[[paste0("center_", cen)]]
    
    if (is.null(cov_cen)) {
      cov_cen <- character(0)
    }
    
    cov_cen <- cov_cen[
      cov_cen %in% colnames(data_OB_i)
    ]
    
    terms_model <- c("X", cov_cen)
    
    formula_adj <- reformulate(
      termlabels = terms_model,
      response = "Y"
    )
    
    if (type == 1) {
      
      model_adj <- try(
        lm(formula_adj, data = data_OB_i),
        silent = TRUE
      )
      
    } else {
      
      model_adj <- try(
        glm(
          formula_adj,
          family = binomial(link = "logit"),
          data = data_OB_i
        ),
        silent = TRUE
      )
    }
    
    if (inherits(model_adj, "try-error")) {
      next
    }
    
    sum_model_adj <- try(summary(model_adj), silent = TRUE)
    
    if (inherits(sum_model_adj, "try-error")) {
      next
    }
    
    coefs <- sum_model_adj$coefficients
    
    if (!"X" %in% rownames(coefs)) {
      next
    }
    
    betaXY <- coefs["X", 1]
    SEXY   <- coefs["X", 2]
    
    if (!is.finite(betaXY) || !is.finite(SEXY) || SEXY <= 0) {
      next
    }
    
    MR_result_all <- rbind(
      MR_result_all,
      c(betaXY, SEXY)
    )
  }
  
  if (is.null(MR_result_all) || nrow(MR_result_all) < 2) {
    return(NULL)
  }
  
  end_end_result <- meta_function_MR(MR_result_all)
  
  return(end_end_result)
}


## 连续暴露 GPS adjustment###############
## 注意：这不是 PSM，而是 GPS 预测均值调整

ALL_GPS_adjust <- function(ALLdata,
                           N_center,
                           type,
                           covariates) {
  
  if (type != 1) {
    stop("ALL_GPS_adjust is designed for continuous exposure, type == 1.")
  }
  
  MR_result_all <- NULL
  
  for (cen in 1:N_center) {
    
    data_OB_i <- ALLdata$all_data[
      ALLdata$all_data$centre == cen,
      ,
      drop = FALSE
    ]
    
    cov_cen <- covariates[[paste0("center_", cen)]]
    
    if (is.null(cov_cen) || length(cov_cen) == 0) {
      next
    }
    
    cov_cen <- cov_cen[
      cov_cen %in% colnames(data_OB_i)
    ]
    
    if (length(cov_cen) == 0) {
      next
    }
    
    
    ## 1. 暴露模型：X ~ covariates
    
    formula_gps <- reformulate(
      termlabels = cov_cen,
      response = "X"
    )
    
    model_gps <- try(
      lm(formula_gps, data = data_OB_i),
      silent = TRUE
    )
    
    if (inherits(model_gps, "try-error")) {
      next
    }
    
    data_OB_i$gps_mu <- as.numeric(predict(model_gps))
    
    
    ## 2. 结局模型：Y ~ X + gps_mu
    
    model_y <- try(
      lm(Y ~ X + gps_mu, data = data_OB_i),
      silent = TRUE
    )
    
    if (inherits(model_y, "try-error")) {
      next
    }
    
    sum_model_y <- try(summary(model_y), silent = TRUE)
    
    if (inherits(sum_model_y, "try-error")) {
      next
    }
    
    coefs <- sum_model_y$coefficients
    
    if (!"X" %in% rownames(coefs)) {
      next
    }
    
    betaXY <- coefs["X", 1]
    SEXY   <- coefs["X", 2]
    
    if (!is.finite(betaXY) || !is.finite(SEXY) || SEXY <= 0) {
      next
    }
    
    MR_result_all <- rbind(
      MR_result_all,
      c(betaXY, SEXY)
    )
  }
  
  if (is.null(MR_result_all) || nrow(MR_result_all) < 2) {
    return(NULL)
  }
  
  end_end_result <- meta_function_MR(MR_result_all)
  
  return(end_end_result)
}


## 二分类暴露的 stabilized IPTW#############
##
## exposure    : 暴露变量名，例如 "X"
## denominator : 分母模型公式，例如 X ~ C1 + C2 + L1_1
## data        : 当前中心数据
## trim        : 权重截尾比例，默认 0.01；如果不想截尾，可设为 NULL 或 0
## eps         : 概率截断，避免预测概率为 0 或 1
##
## 返回：
##   与 data 行数相同的 stabilized IPTW 权重


ipw_mew_cat <- function(exposure,
                        denominator,
                        data,
                        trim = 0.01,
                        eps = 1e-6) {
  
  
  ## 0. 基本检查
  
  if (!exposure %in% colnames(data)) {
    stop("exposure is not found in data.")
  }
  
  A <- data[[exposure]]
  
  ## 当前函数只用于二分类 0/1 暴露
  valid_A <- A[!is.na(A)]
  
  if (!all(valid_A %in% c(0, 1))) {
    stop("ipw_mew_cat currently requires binary exposure coded as 0/1.")
  }
  
  n <- nrow(data)
  
  if (n == 0) {
    return(numeric(0))
  }
  
  
  ## 1. Numerator model: P(A)
  ##    稳定权重分子，只包含截距项
  
  formula_num <- as.formula(
    paste(exposure, "~ 1")
  )
  
  mod_num <- try(
    glm(
      formula = formula_num,
      family = binomial(link = "logit"),
      data = data
    ),
    silent = TRUE
  )
  
  if (inherits(mod_num, "try-error")) {
    return(rep(NA_real_, n))
  }
  
  p_num <- try(
    predict(mod_num, type = "response"),
    silent = TRUE
  )
  
  if (inherits(p_num, "try-error")) {
    return(rep(NA_real_, n))
  }
  
  p_num <- pmin(pmax(as.numeric(p_num), eps), 1 - eps)
  
  
  ## 2. Denominator model: P(A | C)
  
  mod_den <- try(
    glm(
      formula = denominator,
      family = binomial(link = "logit"),
      data = data
    ),
    silent = TRUE
  )
  
  if (inherits(mod_den, "try-error")) {
    return(rep(NA_real_, n))
  }
  
  p_den <- try(
    predict(mod_den, type = "response"),
    silent = TRUE
  )
  
  if (inherits(p_den, "try-error")) {
    return(rep(NA_real_, n))
  }
  
  p_den <- pmin(pmax(as.numeric(p_den), eps), 1 - eps)
  
  
  ## 3. 计算 stabilized IPTW
  ##
  ## A = 1: numerator = P(A=1), denominator = P(A=1 | C)
  ## A = 0: numerator = P(A=0), denominator = P(A=0 | C)
  
  
  w_num <- ifelse(A == 1, p_num, 1 - p_num)
  w_den <- ifelse(A == 1, p_den, 1 - p_den)
  
  weights <- w_num / w_den
  
  weights[!is.finite(weights)] <- NA_real_
  
  
  ## 4. 权重截尾，避免极端权重
  
  if (!is.null(trim) && trim > 0) {
    
    q <- quantile(
      weights,
      probs = c(trim, 1 - trim),
      na.rm = TRUE
    )
    
    weights <- pmin(pmax(weights, q[1]), q[2])
  }
  
  
  ## 5. 附加简单诊断信息，不影响原函数使用
  
  
  attr(weights, "diagnostics") <- list(
    mean_weight = mean(weights, na.rm = TRUE),
    sd_weight = sd(weights, na.rm = TRUE),
    min_weight = min(weights, na.rm = TRUE),
    max_weight = max(weights, na.rm = TRUE),
    n_missing_weight = sum(!is.finite(weights) | is.na(weights))
  )
  
  return(weights)
}

## 连续暴露的 stabilized generalized propensity score weight########
##
## denominator: f(X | observed covariates)
## numerator  : f(X)


ipw_mew_con <- function(data,
                        exposure = "X",
                        covariates,
                        trim = 0.01,
                        density_floor = 1e-8) {
  
  if (!exposure %in% colnames(data)) {
    stop("exposure is not found in data.")
  }
  
  n <- nrow(data)
  
  if (n == 0) {
    return(numeric(0))
  }
  
  covariates <- covariates[
    covariates %in% colnames(data)
  ]
  
  if (length(covariates) == 0) {
    
    weights <- rep(1, n)
    
    attr(weights, "diagnostics") <- list(
      mean_weight = mean(weights, na.rm = TRUE),
      sd_weight = sd(weights, na.rm = TRUE),
      min_weight = min(weights, na.rm = TRUE),
      max_weight = max(weights, na.rm = TRUE),
      n_missing_weight = sum(!is.finite(weights) | is.na(weights))
    )
    
    return(weights)
  }
  
  formula_den <- reformulate(
    termlabels = covariates,
    response = exposure
  )
  
  model_den <- try(
    lm(formula_den, data = data),
    silent = TRUE
  )
  
  if (inherits(model_den, "try-error")) {
    return(rep(NA_real_, n))
  }
  
  mu_den <- as.numeric(predict(model_den))
  resid_den <- data[[exposure]] - mu_den
  sd_den <- sqrt(mean(resid_den^2, na.rm = TRUE))
  
  if (!is.finite(sd_den) || sd_den <= 0) {
    return(rep(NA_real_, n))
  }
  
  den_density <- dnorm(
    data[[exposure]],
    mean = mu_den,
    sd = sd_den
  )
  
  mu_num <- mean(data[[exposure]], na.rm = TRUE)
  sd_num <- sd(data[[exposure]], na.rm = TRUE)
  
  if (!is.finite(sd_num) || sd_num <= 0) {
    return(rep(NA_real_, n))
  }
  
  num_density <- dnorm(
    data[[exposure]],
    mean = mu_num,
    sd = sd_num
  )
  
  den_density <- pmax(den_density, density_floor)
  num_density <- pmax(num_density, density_floor)
  
  weights <- num_density / den_density
  weights[!is.finite(weights)] <- NA_real_
  
  if (!is.null(trim) && trim > 0) {
    
    if (sum(is.finite(weights)) > 0) {
      
      q <- quantile(
        weights,
        probs = c(trim, 1 - trim),
        na.rm = TRUE
      )
      
      if (all(is.finite(q))) {
        weights <- pmin(pmax(weights, q[1]), q[2])
      }
    }
  }
  
  attr(weights, "diagnostics") <- list(
    mean_weight = mean(weights, na.rm = TRUE),
    sd_weight = sd(weights, na.rm = TRUE),
    min_weight = min(weights, na.rm = TRUE),
    max_weight = max(weights, na.rm = TRUE),
    n_missing_weight = sum(!is.finite(weights) | is.na(weights))
  )
  
  return(weights)
}


## 连续暴露 GPS-IPW MSM#######
## 每个中心单独估计权重和 MSM，然后 meta 合并


IPW_function <- function(ALLdata,
                         N_center,
                         exposure = "X",
                         covariates,
                         type,
                         trim = 0.01,
                         robust_se = TRUE) {
  
  result <- NULL
  
  for (cen in 1:N_center) {
    
    data_i <- ALLdata$all_data[
      ALLdata$all_data$centre == cen,
      ,
      drop = FALSE
    ]
    
    cov_cen <- covariates[[paste0("center_", cen)]]
    
    if (is.null(cov_cen)) {
      cov_cen <- character(0)
    }
    
    cov_cen <- cov_cen[
      cov_cen %in% colnames(data_i)
    ]
    
    
    ## 连续暴露：GPS stabilized weight
    
    
    if (type == 1) {
      
      data_i$weight <- ipw_mew_con(
        data = data_i,
        exposure = exposure,
        covariates = cov_cen,
        trim = trim
      )
      
      data_i <- data_i[
        is.finite(data_i$weight) &
          data_i$weight > 0 &
          is.finite(data_i$Y) &
          is.finite(data_i[[exposure]]),
        ,
        drop = FALSE
      ]
      
      if (nrow(data_i) < 10) {
        next
      }
      
      fit_msm <- try(
        lm(
          as.formula(paste("Y ~", exposure)),
          data = data_i,
          weights = weight
        ),
        silent = TRUE
      )
      
      if (inherits(fit_msm, "try-error")) {
        next
      }
      
      if (robust_se && requireNamespace("sandwich", quietly = TRUE)) {
        
        vc <- sandwich::vcovHC(fit_msm, type = "HC0")
        beta <- coef(fit_msm)[exposure]
        SE <- sqrt(diag(vc))[exposure]
        
      } else {
        
        sum_fit <- summary(fit_msm)
        coefs <- sum_fit$coefficients
        
        if (!exposure %in% rownames(coefs)) {
          next
        }
        
        beta <- coefs[exposure, 1]
        SE   <- coefs[exposure, 2]
      }
      
    } else {
      
      
      ## 二分类暴露：IPW
      
      
      if (length(cov_cen) == 0) {
        next
      }
      
      hh_for <- reformulate(
        termlabels = cov_cen,
        response = exposure
      )
      
      data_i$weight <- ipw_mew_cat(
        exposure = exposure,
        denominator = hh_for,
        data = data_i
      )
      
      data_i <- data_i[
        is.finite(data_i$weight) &
          data_i$weight > 0,
        ,
        drop = FALSE
      ]
      
      if (nrow(data_i) < 10) {
        next
      }
      
      fit_msm <- try(
        glm(
          as.formula(paste("Y ~", exposure)),
          data = data_i,
          family = binomial(link = "logit"),
          weights = weight
        ),
        silent = TRUE
      )
      
      if (inherits(fit_msm, "try-error")) {
        next
      }
      
      sum_fit <- summary(fit_msm)
      coefs <- sum_fit$coefficients
      
      if (!exposure %in% rownames(coefs)) {
        next
      }
      
      beta <- coefs[exposure, 1]
      SE   <- coefs[exposure, 2]
    }
    
    if (!is.finite(beta) || !is.finite(SE) || SE <= 0) {
      next
    }
    
    result <- rbind(
      result,
      c(beta, SE)
    )
  }
  
  if (is.null(result) || nrow(result) < 2) {
    return(NULL)
  }
  
  end_end_result <- meta_function_MR(result)
  
  return(end_end_result)
}





## 安全提取 final_result #######


extract_final_result <- function(result_obj) {
  
  empty_final <- c(
    estimate = NA_real_,
    se = NA_real_,
    z_value = NA_real_,
    p_value = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_
  )
  
  if (inherits(result_obj, "try-error")) {
    return(empty_final)
  }
  
  if (is.null(result_obj)) {
    return(empty_final)
  }
  
  if (is.null(result_obj$final_result)) {
    return(empty_final)
  }
  
  out <- as.numeric(result_obj$final_result)
  
  if (length(out) < 6) {
    return(empty_final)
  }
  
  out <- out[1:6]
  names(out) <- names(empty_final)
  
  if (!is.finite(out[1])) {
    return(empty_final)
  }
  
  return(out)
}




## 主模拟函数：四种提出方法 + 回归调整 + GPS调整 + GPS-IPW MSM ##########

############################################################
## 安全合并 data.frame
##
## 用途：
##   不同方法返回的 center_result / iv_selection_info 列数可能不同。
##   直接 rbind 会报错：
##     Error in rbind(deparse.level, ...) : 变量的列数不对
##
## 本函数会自动补齐缺失列，再按统一列顺序合并。
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
  
  missing_x <- setdiff(
    all_cols,
    colnames(x)
  )
  
  missing_y <- setdiff(
    all_cols,
    colnames(y)
  )
  
  if (length(missing_x) > 0) {
    for (cc in missing_x) {
      x[[cc]] <- NA
    }
  }
  
  if (length(missing_y) > 0) {
    for (cc in missing_y) {
      y[[cc]] <- NA
    }
  }
  
  x <- x[, all_cols, drop = FALSE]
  y <- y[, all_cols, drop = FALSE]
  
  out <- rbind(
    x,
    y
  )
  
  rownames(out) <- NULL
  
  return(out)
}

############################################################
## 主模拟函数：
## 四种提出方法 + 回归调整 + GPS adjustment + GPS-IPW MSM
##
## 四种提出方法：
##   1. IV_single_same : 单工具 + 不同中心同一工具
##   2. IV_single_diff : 单工具 + 不同中心不同工具
##   3. IV_multi_same  : 多工具 + 不同中心同一组工具
##   4. IV_multi_diff  : 多工具 + 不同中心不同工具
##
## 输出中增加：
##   selected_*_iv_count : 筛选阶段工具个数
##   used_*_iv_count     : 实际有效估计使用的工具个数
############################################################

All_simulation <- function(N_sim,
                           N_center,
                           type,
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
                           n_boot,
                           correlation_threshold,
                           selcetnum = 1,
                           trim = 0.01) {
  
  ##########################################################
  ## 内部小函数 1：安全合并 data.frame
  ## 用于合并不同方法返回的 center_result / iv_selection_info
  ##########################################################
  
  safe_rbind_df_local <- function(x, y) {
    
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
    
    missing_x <- setdiff(
      all_cols,
      colnames(x)
    )
    
    missing_y <- setdiff(
      all_cols,
      colnames(y)
    )
    
    if (length(missing_x) > 0) {
      for (cc in missing_x) {
        x[[cc]] <- NA
      }
    }
    
    if (length(missing_y) > 0) {
      for (cc in missing_y) {
        y[[cc]] <- NA
      }
    }
    
    x <- x[, all_cols, drop = FALSE]
    y <- y[, all_cols, drop = FALSE]
    
    out <- rbind(
      x,
      y
    )
    
    rownames(out) <- NULL
    
    return(out)
  }
  
  
  ##########################################################
  ## 内部小函数 2：生成空变量列表
  ##########################################################
  
  empty_covariate_list <- function(N_center) {
    
    tmp <- vector("list", N_center)
    
    for (cen in 1:N_center) {
      tmp[[paste0("center_", cen)]] <- character(0)
    }
    
    return(tmp)
  }
  
  
  ##########################################################
  ## 内部小函数 3：根据筛选变量列表计算“筛选阶段”的工具个数
  ##########################################################
  
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
  
  
  ##########################################################
  ## 内部小函数 4：根据估计结果计算“实际使用”的工具个数
  ##########################################################
  
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
    
    ########################################################
    ## 情况 1：多工具结果
    ## center_result 中有 iv_count
    ########################################################
    
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
    }
    
    ########################################################
    ## 情况 2：单工具结果
    ## 如果 center_result 中有有效中心，则每个有效中心使用 1 个工具
    ########################################################
    
    if (!is.null(result_obj$center_result)) {
      
      center_result <- as.data.frame(
        result_obj$center_result,
        stringsAsFactors = FALSE
      )
      
      if (all(c("b", "se") %in% colnames(center_result))) {
        
        valid_center <- is.finite(center_result$b) &
          is.finite(center_result$se) &
          center_result$se > 0
        
        if (sum(valid_center, na.rm = TRUE) > 0) {
          
          used_count <- rep(
            1,
            sum(valid_center, na.rm = TRUE)
          )
          
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
  
  
  ##########################################################
  ## 内部小函数 5：
  ## 提取 final_result + 筛选工具个数 + 实际使用工具个数
  ##########################################################
  
  extract_final_result_with_iv_count <- function(result_obj,
                                                 selected_covariates,
                                                 N_center) {
    
    if (inherits(result_obj, "try-error")) {
      final_result <- extract_final_result(NULL)
    } else {
      final_result <- extract_final_result(result_obj)
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
  
  
  ##########################################################
  ## 内部小函数 6：保存提出方法的工具信息和中心信息
  ##########################################################
  
  save_iv_diagnostics <- function(result_tmp,
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
    
    ########################################################
    ## 1. 保存工具筛选总体信息
    ########################################################
    
    if (!is.null(result_tmp$iv_selection_info)) {
      
      tmp_tool <- as.data.frame(
        result_tmp$iv_selection_info,
        stringsAsFactors = FALSE
      )
      
      tmp_tool$sim <- sim
      tmp_tool$method <- method_name
      
      all_tool_info_IV <- safe_rbind_df_local(
        all_tool_info_IV,
        tmp_tool
      )
    }
    
    ########################################################
    ## 2. 保存中心层面结果
    ## 单工具和多工具 center_result 列可能不同，必须安全合并
    ########################################################
    
    if (!is.null(result_tmp$center_result)) {
      
      tmp_center <- as.data.frame(
        result_tmp$center_result,
        stringsAsFactors = FALSE
      )
      
      if (nrow(tmp_center) > 0) {
        
        tmp_center$sim <- sim
        tmp_center$method <- method_name
        
        all_center_info_IV <- safe_rbind_df_local(
          all_center_info_IV,
          tmp_center
        )
      }
    }
    
    return(list(
      all_tool_info_IV = all_tool_info_IV,
      all_center_info_IV = all_center_info_IV
    ))
  }
  
  
  ##########################################################
  ## 1. 真实效应
  ##########################################################
  
  if (type == 1) {
    
    true_effect <- beta_XY
    
  } else {
    
    true_data <- Data_Generated_do(
      type = type,
      typeC = typeC,
      UD = UD,
      r = r,
      N = 100000,
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
      p_x = p_x
    )
    
    true_effect <-
      sum(true_data$Y == 1 & true_data$X == 1) / sum(true_data$X == 1) -
      sum(true_data$Y == 1 & true_data$X == 0) / sum(true_data$X == 0)
  }
  
  
  ##########################################################
  ## 2. 存储每次模拟结果
  ##########################################################
  
  all_result_IV_single_same <- NULL
  all_result_IV_single_diff <- NULL
  all_result_IV_multi_same <- NULL
  all_result_IV_multi_diff <- NULL
  
  all_result_ADJ <- NULL
  all_result_GPS <- NULL
  all_result_IPW <- NULL
  
  all_tool_info_IV <- NULL
  all_center_info_IV <- NULL
  all_sim_status <- NULL
  
  
  ##########################################################
  ## 3. 开始模拟
  ##########################################################
  
  for (sim in 1:N_sim) {
    
    cat("sim =", sim, "\n")
    
    
    ########################################################
    ## 3.1 生成数据
    ########################################################
    
    data <- try(
      ALL_data(
        type = type,
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
        p_x = p_x
      ),
      silent = TRUE
    )
    
    
    ########################################################
    ## 如果数据生成失败
    ########################################################
    
    if (inherits(data, "try-error")) {
      
      empty_IV <- c(
        extract_final_result(NULL),
        
        selected_mean_iv_count = NA_real_,
        selected_median_iv_count = NA_real_,
        selected_min_iv_count = NA_real_,
        selected_max_iv_count = NA_real_,
        
        used_mean_iv_count = NA_real_,
        used_median_iv_count = NA_real_,
        used_min_iv_count = NA_real_,
        used_max_iv_count = NA_real_
      )
      
      all_result_IV_single_same <- rbind(
        all_result_IV_single_same,
        empty_IV
      )
      
      all_result_IV_single_diff <- rbind(
        all_result_IV_single_diff,
        empty_IV
      )
      
      all_result_IV_multi_same <- rbind(
        all_result_IV_multi_same,
        empty_IV
      )
      
      all_result_IV_multi_diff <- rbind(
        all_result_IV_multi_diff,
        empty_IV
      )
      
      all_result_ADJ <- rbind(
        all_result_ADJ,
        extract_final_result(NULL)
      )
      
      all_result_GPS <- rbind(
        all_result_GPS,
        extract_final_result(NULL)
      )
      
      all_result_IPW <- rbind(
        all_result_IPW,
        extract_final_result(NULL)
      )
      
      all_sim_status <- rbind(
        all_sim_status,
        data.frame(
          sim = sim,
          data_success = FALSE,
          
          IV_single_same_success = FALSE,
          IV_single_diff_success = FALSE,
          IV_multi_same_success = FALSE,
          IV_multi_diff_success = FALSE,
          
          ADJ_success = FALSE,
          GPS_success = FALSE,
          IPW_success = FALSE
        )
      )
      
      next
    }
    
    
    ########################################################
    ## 3.2 候选变量
    ########################################################
    
    covariates_all <- setdiff(
      colnames(data$all_data),
      c("X", "Y", "centre")
    )
    
    
    ########################################################
    ## 3.3 传统方法变量筛选
    ########################################################
    
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
    
    
    ########################################################
    ## 3.4 提出方法一：
    ##     单工具 + 不同中心同一工具
    ########################################################
    
    select_IV_single_same <- try(
      variable_selection_largest_center_common(
        N_center = N_center,
        covariates = covariates_all,
        data = data,
        correlation_threshold = correlation_threshold,
        type = type,
        selcetnum = 1
      ),
      silent = TRUE
    )
    
    if (inherits(select_IV_single_same, "try-error")) {
      select_IV_single_same <- empty_covariate_list(N_center)
    }
    
    result_IV_single_same_q <- try(
      ALLCONIV_cov(
        ALLdata = data,
        N_center = N_center,
        type = type,
        covariates = select_IV_single_same,
        n_boot = n_boot
      ),
      silent = TRUE
    )
    
    final_IV_single_same <- extract_final_result_with_iv_count(
      result_obj = result_IV_single_same_q,
      selected_covariates = select_IV_single_same,
      N_center = N_center
    )
    
    all_result_IV_single_same <- rbind(
      all_result_IV_single_same,
      final_IV_single_same
    )
    
    
    ########################################################
    ## 3.5 提出方法二：
    ##     单工具 + 不同中心不同工具
    ########################################################
    
    select_IV_single_diff <- try(
      variable_selection_sorted(
        N_center = N_center,
        covariates = covariates_all,
        data = data,
        correlation_threshold = correlation_threshold,
        type = type,
        selcetnum = 1
      ),
      silent = TRUE
    )
    
    if (inherits(select_IV_single_diff, "try-error")) {
      select_IV_single_diff <- empty_covariate_list(N_center)
    }
    
    result_IV_single_diff_q <- try(
      ALLCONIV_cov(
        ALLdata = data,
        N_center = N_center,
        type = type,
        covariates = select_IV_single_diff,
        n_boot = n_boot
      ),
      silent = TRUE
    )
    
    final_IV_single_diff <- extract_final_result_with_iv_count(
      result_obj = result_IV_single_diff_q,
      selected_covariates = select_IV_single_diff,
      N_center = N_center
    )
    
    all_result_IV_single_diff <- rbind(
      all_result_IV_single_diff,
      final_IV_single_diff
    )
    
    
    ########################################################
    ## 3.6 提出方法三：
    ##     多工具 + 不同中心同一组工具
    ########################################################
    
    select_IV_multi_same <- try(
      variable_selection_largest_center_common(
        N_center = N_center,
        covariates = covariates_all,
        data = data,
        correlation_threshold = correlation_threshold,
        type = type,
        selcetnum = "multi"
      ),
      silent = TRUE
    )
    
    if (inherits(select_IV_multi_same, "try-error")) {
      select_IV_multi_same <- empty_covariate_list(N_center)
    }
    
    result_IV_multi_same_q <- try(
      ALLCONIV_cov_multi(
        ALLdata = data,
        N_center = N_center,
        type = type,
        covariates = select_IV_multi_same,
        n_boot = n_boot
      ),
      silent = TRUE
    )
    
    final_IV_multi_same <- extract_final_result_with_iv_count(
      result_obj = result_IV_multi_same_q,
      selected_covariates = select_IV_multi_same,
      N_center = N_center
    )
    
    all_result_IV_multi_same <- rbind(
      all_result_IV_multi_same,
      final_IV_multi_same
    )
    
    
    ########################################################
    ## 3.7 提出方法四：
    ##     多工具 + 不同中心不同工具
    ########################################################
    
    select_IV_multi_diff <- try(
      variable_selection_sorted(
        N_center = N_center,
        covariates = covariates_all,
        data = data,
        correlation_threshold = correlation_threshold,
        type = type,
        selcetnum = "multi"
      ),
      silent = TRUE
    )
    
    if (inherits(select_IV_multi_diff, "try-error")) {
      select_IV_multi_diff <- empty_covariate_list(N_center)
    }
    
    result_IV_multi_diff_q <- try(
      ALLCONIV_cov_multi(
        ALLdata = data,
        N_center = N_center,
        type = type,
        covariates = select_IV_multi_diff,
        n_boot = n_boot
      ),
      silent = TRUE
    )
    
    final_IV_multi_diff <- extract_final_result_with_iv_count(
      result_obj = result_IV_multi_diff_q,
      selected_covariates = select_IV_multi_diff,
      N_center = N_center
    )
    
    all_result_IV_multi_diff <- rbind(
      all_result_IV_multi_diff,
      final_IV_multi_diff
    )
    
    
    ########################################################
    ## 3.8 保存四种提出方法的工具信息和中心信息
    ########################################################
    
    method_result_list <- list(
      IV_single_same = result_IV_single_same_q,
      IV_single_diff = result_IV_single_diff_q,
      IV_multi_same = result_IV_multi_same_q,
      IV_multi_diff = result_IV_multi_diff_q
    )
    
    for (method_name in names(method_result_list)) {
      
      tmp_save <- save_iv_diagnostics(
        result_tmp = method_result_list[[method_name]],
        method_name = method_name,
        sim = sim,
        all_tool_info_IV = all_tool_info_IV,
        all_center_info_IV = all_center_info_IV
      )
      
      all_tool_info_IV <- tmp_save$all_tool_info_IV
      all_center_info_IV <- tmp_save$all_center_info_IV
    }
    
    
    ########################################################
    ## 3.9 方法五：传统回归调整
    ########################################################
    
    result_ADJq <- try(
      ALL_adjust(
        ALLdata = data,
        N_center = N_center,
        type = type,
        covariates = covariates_both
      ),
      silent = TRUE
    )
    
    final_ADJ <- extract_final_result(result_ADJq)
    
    all_result_ADJ <- rbind(
      all_result_ADJ,
      final_ADJ
    )
    
    
    ########################################################
    ## 3.10 方法六：GPS mean adjustment
    ########################################################
    
    result_GPSq <- try(
      ALL_GPS_adjust(
        ALLdata = data,
        N_center = N_center,
        type = type,
        covariates = covariates_X
      ),
      silent = TRUE
    )
    
    final_GPS <- extract_final_result(result_GPSq)
    
    all_result_GPS <- rbind(
      all_result_GPS,
      final_GPS
    )
    
    
    ########################################################
    ## 3.11 方法七：GPS-IPW MSM
    ########################################################
    
    result_IPWq <- try(
      IPW_function(
        ALLdata = data,
        N_center = N_center,
        exposure = "X",
        covariates = covariates_X,
        type = type,
        trim = trim
      ),
      silent = TRUE
    )
    
    final_IPW <- extract_final_result(result_IPWq)
    
    all_result_IPW <- rbind(
      all_result_IPW,
      final_IPW
    )
    
    
    ########################################################
    ## 3.12 保存本次模拟状态
    ########################################################
    
    all_sim_status <- rbind(
      all_sim_status,
      data.frame(
        sim = sim,
        data_success = TRUE,
        
        IV_single_same_success = is.finite(final_IV_single_same[1]),
        IV_single_diff_success = is.finite(final_IV_single_diff[1]),
        IV_multi_same_success = is.finite(final_IV_multi_same[1]),
        IV_multi_diff_success = is.finite(final_IV_multi_diff[1]),
        
        ADJ_success = is.finite(final_ADJ[1]),
        GPS_success = is.finite(final_GPS[1]),
        IPW_success = is.finite(final_IPW[1])
      )
    )
  }
  
  
  ##########################################################
  ## 4. 汇总模拟结果
  ##########################################################
  
  end_result_IV_single_same <- resultt_FUn(
    all_result_IV_single_same,
    true_effect
  )
  
  end_result_IV_single_diff <- resultt_FUn(
    all_result_IV_single_diff,
    true_effect
  )
  
  end_result_IV_multi_same <- resultt_FUn(
    all_result_IV_multi_same,
    true_effect
  )
  
  end_result_IV_multi_diff <- resultt_FUn(
    all_result_IV_multi_diff,
    true_effect
  )
  
  end_result_ADJ <- resultt_FUn(
    all_result_ADJ,
    true_effect
  )
  
  end_result_GPS <- resultt_FUn(
    all_result_GPS,
    true_effect
  )
  
  end_result_IPW <- resultt_FUn(
    all_result_IPW,
    true_effect
  )
  
  
  ##########################################################
  ## 5. 返回
  ##########################################################
  
  end_resulttt <- list()
  
  end_resulttt$IV_single_same <- end_result_IV_single_same
  end_resulttt$IV_single_diff <- end_result_IV_single_diff
  end_resulttt$IV_multi_same <- end_result_IV_multi_same
  end_resulttt$IV_multi_diff <- end_result_IV_multi_diff
  
  end_resulttt$ADJ <- end_result_ADJ
  end_resulttt$GPS <- end_result_GPS
  end_resulttt$IPW <- end_result_IPW
  
  end_resulttt$true <- true_effect
  
  end_resulttt$all_result_IV_single_same <- all_result_IV_single_same
  end_resulttt$all_result_IV_single_diff <- all_result_IV_single_diff
  end_resulttt$all_result_IV_multi_same <- all_result_IV_multi_same
  end_resulttt$all_result_IV_multi_diff <- all_result_IV_multi_diff
  
  end_resulttt$all_result_ADJ <- all_result_ADJ
  end_resulttt$all_result_GPS <- all_result_GPS
  end_resulttt$all_result_IPW <- all_result_IPW
  
  end_resulttt$tool_info_IV <- all_tool_info_IV
  end_resulttt$center_info_IV <- all_center_info_IV
  end_resulttt$sim_status <- all_sim_status
  
  return(end_resulttt)
}


## 参数设置############

m <- 30
k <- 5

n1 <- 5
n2 <- 5
n3 <- 15

Clv <- 0.5

r_L <- 0.5
connection_prob <- 0.5

conf_L <- sample(1:(n1 + n2 + n3), n3 + 2)

beta_CX <- runif(m + k + length(conf_L), 0, 1)
beta_CY <- runif(m + k + length(conf_L), 0, 1)
beta_XY <- 0.5

Nrct <- 10000
p_x <- 0

typeC <- 1
type <- 1

N_sim <- 1000
N_center <- 10

correlation_threshold <- 0.1
n_boot <- 1000

selcetnum <- 1
trim <- 0.01

## 这里 para_sim 表示 C/U 相关结构
## UD = TRUE, r = 0 表示 U 与 C 独立
## UD = FALSE, r > 0 表示 U 与 C 有相关
para_sim <- data.frame(
  UD = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  r = seq(0,0.6,0.06)
)


## 输出路径

out_dir <- "results"
dir.create(path.expand(out_dir), recursive = TRUE, showWarnings = FALSE)

all_sim_result_all <- NULL


############################################################
## 开始循环
############################################################

for (j in 1:nrow(para_sim)) {
  
  cat("Scenario =", j, "UD =", para_sim$UD[j], "r =", para_sim$r[j], "\n")
  
  UD <- para_sim$UD[j]
  r <- para_sim$r[j]
  
  Nall <- floor(
    runif(N_center, 10000, 20000)
  )
  
  result_sim_result <- All_simulation(
    N_sim = N_sim,
    N_center = N_center,
    type = type,
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
    n_boot = n_boot,
    correlation_threshold = correlation_threshold,
    selcetnum = selcetnum,
    trim = trim
  )
  
  result_row <- c(
    setNames(result_sim_result$IV_single_same, paste0("IV_single_same_", names(result_sim_result$IV_single_same))),
    setNames(result_sim_result$IV_single_diff, paste0("IV_single_diff_", names(result_sim_result$IV_single_diff))),
    setNames(result_sim_result$IV_multi_same, paste0("IV_multi_same_", names(result_sim_result$IV_multi_same))),
    setNames(result_sim_result$IV_multi_diff, paste0("IV_multi_diff_", names(result_sim_result$IV_multi_diff))),
    setNames(result_sim_result$ADJ, paste0("ADJ_", names(result_sim_result$ADJ))),
    setNames(result_sim_result$GPS, paste0("GPS_", names(result_sim_result$GPS))),
    setNames(result_sim_result$IPW, paste0("IPW_", names(result_sim_result$IPW))),
    
    true_effect = result_sim_result$true,
    para_sim = j,
    beta_XY = beta_XY,
    
    m = m,
    k = k,
    n1 = n1,
    n2 = n2,
    n3 = n3,
    typeC = typeC,
    type = type,
    
    N_center = N_center,
    Nrct = Nrct,
    Nall_mean = mean(Nall),
    Nall_min = min(Nall),
    Nall_max = max(Nall),
    Nall_total = sum(Nall),
    
    r = r,
    UD = as.numeric(UD),
    Clv = Clv,
    r_L = r_L,
    connection_prob = connection_prob,
    correlation_threshold = correlation_threshold,
    n_boot = n_boot,
    trim = trim
  )
  
  result_row <- as.data.frame(as.list(result_row), stringsAsFactors = FALSE)
  
  all_sim_result_all <- rbind(
    all_sim_result_all,
    result_row
  )
  
  write.csv(
    all_sim_result_all,
    file = file.path(
      path.expand(out_dir),
      "Sim_UD_r_type1_all_methods_comparison_running.csv"
    ),
    row.names = FALSE
  )
}

write.csv(
  all_sim_result_all,
  file = file.path(
    path.expand(out_dir),
    "Sim_UD_r_type1_all_methods_comparison.csv"
  ),
  row.names = FALSE
)