
############################################################
## 0. set the number of CPU threads used by each R task
## must be placed at the top of the script
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
## 1. data generation
############################################################

suppressPackageStartupMessages({
  library(MASS)
})


## 1. basic helper functions ###############
## binary covariate generator: dichotomize at quantiles
Fenlei <- function(data, Clv, na.rm = TRUE) {
  
  if (Clv <= 0 || Clv >= 1) {
    stop("Clv must be between 0 and 1.")
  }
  
  cutoff <- quantile(data, probs = Clv, na.rm = na.rm, names = FALSE)
  newdata <- ifelse(data > cutoff, 1, 0)
  
  return(newdata)
}


## safe standardization
safe_scale_vec <- function(x) {
  
  x <- as.numeric(x)
  sx <- sd(x, na.rm = TRUE)
  
  if (!is.finite(sx) || sx <= 0) {
    return(rep(0, length(x)))
  }
  
  return(as.numeric(scale(x)))
}


## calibrate the intercept to the target incidence in the latent-variable model

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


## check the length of beta
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



## 2. fix the L-layer network structure #########


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



## 3. generate C/U/L covariate blocks ########


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
  ## 3.1 generate baseline C/U
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
  ## 3.2 initialize the L-layer data
  ##########################################################
  
  simulated_data <- matrix(0, nrow = N, ncol = total_vars)
  colnames(simulated_data) <- c(
    paste0("L1_", 1:n1),
    paste0("L2_", 1:n2),
    paste0("L3_", 1:n3)
  )
  
  ##########################################################
  ## 3.3 generate L1
  ##########################################################
  
  Sigma_L <- matrix(r_L, nrow = n1, ncol = n1)
  diag(Sigma_L) <- 1
  
  mu_L <- rep(0, n1)
  l1_data <- MASS::mvrnorm(n = N, mu = mu_L, Sigma = Sigma_L)
  simulated_data[, 1:n1] <- l1_data
  
  ##########################################################
  ## 3.4 generate L2 and L3 from the fixed network
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
  ## 3.5 combine C/U/L
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



## 4. generate data for a single observational center ######################


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
                           
                           ## binary latent-variable model parameters
                           target_px_ob = 0.50,
                           target_py0 = 0.20,
                           conf_strength_x = 0.50,
                           conf_strength_y = 0.50,
                           error_sd_x = 1,
                           error_sd_y = 1,
                           
                           ## whether to keep the latent variables, for debugging
                           keep_latent = FALSE) {
  
  NC <- m + k
  p_confound <- NC + length(conf_L)
  
  check_beta_length(beta_CX, p_confound, "beta_CX")
  check_beta_length(beta_CY, p_confound, "beta_CY")
  
  ##########################################################
  ## 4.1 generate C/U/L
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
  ## 4.2 generate X and Y
  ##########################################################
  
  if (type == 1) {
    
    ########################################################
    ## continuous X/Y: linear structure
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
    ## binary X/Y: latent-threshold model
    ########################################################
    
    ## latent linear part of X
    eta_x_raw <- as.numeric(H %*% beta_CX)
    eta_x <- safe_scale_vec(eta_x_raw)
    eta_x <- conf_strength_x * eta_x
    
    ## calibrate alpha_x to the target exposure rate
    alpha_x <- calibrate_intercept(
      eta = eta_x,
      target_prev = target_px_ob,
      error_sd = error_sd_x
    )
    
    ## generate the observational exposure X
    eps_x <- rnorm(N, mean = 0, sd = error_sd_x)
    x_star <- alpha_x + eta_x + eps_x
    x <- as.integer(x_star > 0)
    
    ## latent linear part of Y
    eta_y_raw <- as.numeric(H %*% beta_CY)
    eta_y <- safe_scale_vec(eta_y_raw)
    eta_y <- conf_strength_y * eta_y
    
    ## target_py0 is the mean baseline outcome risk when X=0
    alpha_y <- calibrate_intercept(
      eta = eta_y,
      target_prev = target_py0,
      error_sd = error_sd_y
    )
    
    ## generate the observational outcome Y
    eps_y <- rnorm(N, mean = 0, sd = error_sd_y)
    y_star <- alpha_y + eta_y + beta_XY * x + eps_y
    y <- as.integer(y_star > 0)
    
    ########################################################
    ## compute the true value on the observed risk-difference scale
    ##
    ## note:
    ## beta_XY is the effect on the latent-variable scale;
    ## true_RD is the effect on the observed binary-outcome risk-difference scale.
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
  ## 4.3 drop the unobservable U and optional latent variables
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



## 5. generate RCT / do data ######


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
                              
                              ## binary latent-variable model parameters
                              target_py0 = 0.20,
                              conf_strength_y = 0.50,
                              error_sd_y = 1,
                              
                              ## continuous do exposure error
                              error_sd_x_cont = 1,
                              
                              ## whether to keep the latent variables
                              keep_latent = FALSE) {
  
  NC <- m + k
  p_confound <- NC + length(conf_L)
  
  check_beta_length(beta_CX, p_confound, "beta_CX")
  check_beta_length(beta_CY, p_confound, "beta_CY")
  
  ##########################################################
  ## 5.1 generate C/U/L
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
  ## 5.2 handle p_x
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
  ## 5.3 generate X and Y
  ##########################################################
  
  if (type == 1) {
    
    ########################################################
    ## continuous RCT/do: X is not determined by H
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
    ## binary RCT/do:
    ## X is randomized and not determined by H;
    ## Y is still determined by the non-treatment path of H and the treatment effect of X.
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
  ## 5.4 drop the unobservable U and optional latent variables
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



## 6. overall data generation: multicenter observational data + RCT/do data #################


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
                     
                     ## binary latent-variable model parameters
                     target_px_ob = 0.50,
                     target_py0 = 0.20,
                     conf_strength_x = 0.50,
                     conf_strength_y = 0.50,
                     error_sd_x = 1,
                     error_sd_y = 1,
                     
                     ## whether to keep the latent variables
                     keep_latent = FALSE) {
  
  if (length(Nall) != N_center) {
    stop("Length of Nall must be equal to N_center.")
  }
  
  ##########################################################
  ## 6.1 fix one L network structure
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
  ## 6.2 generate several observational centers
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
      conf_strength_x = conf_strength_x,
      conf_strength_y = conf_strength_y,
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
  }
  
  all_data <- do.call(rbind, data_list)
  
  ##########################################################
  ## 6.3 generate RCT / do data
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
    conf_strength_y = conf_strength_y,
    error_sd_y = error_sd_y,
    keep_latent = keep_latent
  )
  
  rct_data <- gen_rct$data
  rct_diag <- gen_rct$diagnostics
  rct_diag$centre <- "RCT_do"
  rct_diag$N <- Nrct
  
  ##########################################################
  ## 6.4 summarize true values and diagnostics
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
      error_sd_x = error_sd_x,
      error_sd_y = error_sd_y
    )
  }
  
  diagnostics <- list(
    centers = diag_list,
    rct = rct_diag
  )
  
  ##########################################################
  ## 6.5 return
  ##########################################################
  
  return(list(
    all_data = all_data,
    rct = rct_data,
    L_network = L_network,
    truth = truth,
    diagnostics = diagnostics
  ))
}





## binary version: robust linear-projection variable selection ###################


## 0. safe absolute correlation ##############


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



## 1. influence-function robust linear regression: extract coefficient, SE, and P value of a term ##############


robust_lm_term_if <- function(y,
                              x_df,
                              term_name,
                              hc_type = "HC1",
                              min_n = 20) {
  
  ##########################################################
  ## return an empty result
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
  ## basic tidying
  ##########################################################
  
  dat <- data.frame(y = y, x_df)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  
  n <- nrow(dat)
  
  if (n < min_n) {
    return(empty_result())
  }
  
  ## y must vary
  if (length(unique(dat$y)) < 2) {
    return(empty_result())
  }
  
  ## build the design matrix
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
  
  ## check the rank
  qr_mm <- qr(mm)
  
  if (qr_mm$rank < p) {
    return(empty_result())
  }
  
  ##########################################################
  ## OLS point estimate
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
  ## influence-function robust variance
  ##
  ## beta_hat - beta ≈ (X'X)^(-1) sum_i x_i e_i
  ## contribution of each observation to beta:
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
  
  ## HC1 small-sample correction
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



## 2. simple linear projection: A ~ Z ##############
##    used for X ~ Z


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



## 3. single-center selection function ##############


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
  ## empty result
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
  ## compute the selection statistics for each candidate variable
  ##########################################################
  
  score_list <- vector("list", length(covariates_i))
  
  for (j in seq_along(covariates_i)) {
    
    var_j <- covariates_i[j]
    
    if (!var_j %in% colnames(data_i)) {
      next
    }
    
    Z_j <- data_i[[var_j]]
    
    ## variable has no variation; skip
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
    ## 1. exposure association: X ~ Z
    ########################################################
    
    fit_x <- robust_simple_slope(
      A = data_i$X,
      Z = Z_j,
      hc_type = hc_type,
      min_n = min_n
    )
    
    ########################################################
    ## 2. outcome association: Y ~ X + Z or Y ~ Z
    ##
    ## Y ~ X + Z is used by default:
    ## the aim is to select candidate confounding instruments still associated with Y given X.
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
    ## 3. selection rules
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
  ## candidate variables: pass both the X and Y criteria
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
  ## ordering: robust F_XZ first, then |beta_XZ|
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
  ## remove highly correlated variables
  ##
  ## since candidates are already sorted by F_XZ,
  ## so when correlation is too high, keep the earlier, stronger variable and skip the later, weaker one.
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
  ## use selcetnum to choose single or multiple instruments
  ##
  ## logic:
  ##   selcetnum = 1       : select only the strongest one
  ##   selcetnum = "multi" : keep several
  ##   selcetnum > 1       : if the number of candidates <= selcetnum, fall back to the strongest one;
  ##                         if the number of candidates > selcetnum, keep several
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
  ## limit the maximum number of instruments
  ##########################################################
  
  if (is.finite(max_iv) && length(final_vars) > max_iv) {
    final_vars <- final_vars[seq_len(max_iv)]
  }
  
  return(list(
    selected_vars = final_vars,
    score_table = score_table
  ))
}



## 4. multicenter: each center selects its instruments independently ##############
##
## can directly replace the original variable_selection_sorted()


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



## 5. multicenter: select common instruments in the largest center ##############
##
## same approach as variable_selection_largest_center_common()


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
  ## find the center with the largest sample size
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
  ## select in the largest center
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
  ## all centers use the same instrument set
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
## binary CON-IV estimation function
##
## applies to:
##   1. binary X
##   2. binary Y
##   3. no logistic model is used
##   4. estimate the effect on the observed risk-difference scale
##   5. variance: robust influence-function covariance + delta method
##   6. multicenter combination: GLS corrected for the shared RCT/do covariance
############################################################

suppressPackageStartupMessages({
  library(MASS)
})


############################################################
## 1. empty-result function
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
## 2. single-variable linear projection slope + influence function
##
## estimation:
##   A ~ Z
##
## returns:
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
## 3. simple slope matrix for several instruments
##
## estimate each Z_j separately:
##   A ~ Z_j
##
## returns:
##   beta vector
##   IF matrix
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
## 4. covariance matrix of the Y ~ Z slopes in RCT/do
##
## all centers share the same RCT/do, so the covariance here
## these enter the off-diagonal elements of the multicenter V matrix.
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
## 5. weight-matrix constructor
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
## 6. single-center binary CON-IV estimate
##
## tool_mode:
##   single : single-instrument Wald ratio
##   multi  : multiple-instrument IVW/GMM-type estimate
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
  ## 1. estimation in observational centers:
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
  ## 2. extract vectors and matrices
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
  
  ## within the same observational center:
  ## Cov(delta, pi) = Cov(gamma_OB, pi)
  Sigma_delta_pi <- crossprod(IF_y, IF_x)
  
  Sigma_delta <- Sigma_g_OB + Sigma_g_RCT
  
  K <- length(common_iv)
  
  ##########################################################
  ## 3. weight matrix
  ##########################################################
  
  W <- make_weight_matrix(
    Sigma_delta = Sigma_delta,
    weight_type = weight_type
  )
  
  ##########################################################
  ## 4. unified estimation formula
  ##
  ## with a single instrument this reduces automatically to:
  ##   theta = delta / pi
  ##
  ## with multiple instruments:
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
## 7. multicenter center-by-center estimation
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
  
  ## selection results can be a vector or list(selected = ...)
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
## 8. numerical stabilization of the V matrix
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
## 9. multicenter GLS combination
##
## diagonal:
##   V_cc = Var(theta_c)
##
## off-diagonal:
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
  ## diagonal
  ##########################################################
  
  for (i in seq_len(C_valid)) {
    V[i, i] <- valid_results[[i]]$var
  }
  
  ##########################################################
  ## off-diagonal: between-center covariance due to the shared RCT/do
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
  ## stabilize the V matrix
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
## 10. one-call function when IV_list is known
##
## this function does not perform variable selection;
## responsible only for:
##   1. collect all IVs
##   2. construct the RCT/do information
##   3. within-center estimation
##   4. multicenter GLS combination
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
## method 1:
## Binary-CONIV-Single-Common
##
## single instrument + all centers share the same instrument
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
## method 2:
## Binary-CONIV-Single-CenterSpecific
##
## single instrument + each center selects its own instrument
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
## method 3:
## Binary-CONIV-Multi-Common
##
## multiple instruments + all centers share the same instrument set
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
  ##   if the number of candidate instruments <= 5, automatically fall back to the strongest one;
  ##   if the number of candidate instruments >  5, keep several instruments.
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
## method 4:
## Binary-CONIV-Multi-CenterSpecific
##
## multiple instruments + each center selects its own instrument set
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
  ## decide separately for each center:
  ##   number of candidate instruments <= 5: the center falls back to a single instrument
  ##   number of candidate instruments >  5: the center uses multiple instruments
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
## run the four binary CON-IV methods together
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
## extract the core results of the four methods
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
## variable selection function for conventional observational methods
##
## selection strategy:
##
## out = "both":
##   associated with both Y and X; used for regression adjustment
##
## out = "X":
##   select only variables associated with X, for PSM and IPW-MSM
##
## out = "Y":
##   select only variables associated with Y; usually a fallback
##
## type == 1:
##   continuous: lm, threshold 0.05 / number of candidates
##
## type == 2:
##   binary: logistic model, threshold 0.05
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
      ## Y association model: Y ~ X + Z_j
      ## X association model: X ~ Z_j
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
    ## P-value threshold: keep the original strategy
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
## robust linear model: extract the RD of X and its robust SE
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
  
  ## HC1 correction
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
## method 1: regression adjustment
##
## variable selection:
##   covariates_both = Select_other(..., out = "both")
##
## model:
##   fit separately in each center:
##     Y ~ X + selected_covariates
##
## effect scale:
##   risk difference RD
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
    ## build the model
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
## compute the logit
############################################################

logit_fun <- function(p) {
  p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  return(log(p / (1 - p)))
}


############################################################
## single-center 1:1 nearest-neighbor PSM
##
## matching without replacement
## caliper defaults to 0.2 * SD(logit PS)
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
## method 2: PSM
##
## variable selection:
##   covariates_X = Select_other(..., out = "X")
##
## PS model:
##   X ~ selected_covariates
##
## matching:
##   1:1 nearest-neighbor matching without replacement within each center, caliper = 0.2 SD(logit PS)
##
## outcome model:
##   Y ~ X in the matched sample
##
## effect scale:
##   risk difference RD
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
    ## PS model
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
    ## 1:1 nearest-neighbor matching
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
## method 3: IPW-MSM
##
## variable selection:
##   covariates_X = Select_other(..., out = "X")
##
## PS model:
##   X ~ selected_covariates
##
## stabilized weights:
##   X=1: P(X=1) / PS
##   X=0: P(X=0) / (1-PS)
##
## MSM:
##   Y ~ X, weighted robust linear model
##
## effect scale:
##   risk difference RD
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
    ## PS model
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
    ## stabilized weights
    ########################################################
    
    p_treat <- mean(dat_ps$X == 1, na.rm = TRUE)
    p_treat <- pmin(pmax(p_treat, ps_trim), 1 - ps_trim)
    
    sw <- ifelse(
      dat_ps$X == 1,
      p_treat / ps,
      (1 - p_treat) / (1 - ps)
    )
    
    ########################################################
    ## weight truncation
    ########################################################
    
    if (!is.null(weight_trunc_quantile)) {
      q_low <- quantile(sw, weight_trunc_quantile[1], na.rm = TRUE)
      q_high <- quantile(sw, weight_trunc_quantile[2], na.rm = TRUE)
      sw <- pmin(pmax(sw, q_low), q_high)
    }
    
    ess <- sum(sw)^2 / sum(sw^2)
    
    ########################################################
    ## weighted MSM: Y ~ X
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
## run the three conventional comparison methods in one call
##
## selection strategy:
##
## regression adjustment:
##   covariates_both = Select_other(..., out = "both")
##
## PSM:
##   covariates_X = Select_other(..., out = "X")
##
## IPW-MSM:
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
  ## 1. candidate variables
  ##########################################################
  
  covariates_all <- setdiff(
    colnames(data$all_data),
    c("X", "Y", "centre")
  )
  
  ##########################################################
  ## 2. variable selection for conventional methods
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
  ## 3. method 1: regression adjustment
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
  ## 4. method 2: PSM
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
  ## 5. method 3: IPW-MSM
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
  ## 6. return
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
## organize the results of the three conventional methods
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
## 1. safely extract final_result
##
## readable output formats:
##   1. CON-IV GLS output:
##      theta_gls, se_gls, p_value, ci_lower, ci_upper
##
##   2. full final_result:
##      estimate, se, z_value, p_value, ci_lower, ci_upper,
##      p_value_t, ci_lower_t, ci_upper_t,
##      p_value_z, ci_lower_z, ci_upper_z,
##      df_meta
##
##   3. simplified final_result:
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
    
    ## other column names can also be read
    b = NA_real_,
    p = NA_real_,
    lower = NA_real_,
    upper = NA_real_
  )
  
  if (inherits(result_obj, "try-error") || is.null(result_obj)) {
    return(empty_final)
  }
  
  ##########################################################
  ## case 1: new binary CON-IV GLS output
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
  ## case 2: list contains final_result
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
  ## case 2.1: has names
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
    ## columns b / p / lower / upper can also be read
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
    ## add the z version
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
  ## case 2.2: results without column names
  ##
  ## length >= 6:
  ##   estimate, se, z_value, p_value, ci_lower, ci_upper
  ##
  ## length == 5:
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
## 2. safely bind data.frames
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
## 3. empty variable list
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
## 4. extract the number of instruments at the selection stage
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
## 5. extract the number of instruments actually used in estimation
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
  ## new binary CON-IV: center_results is a list
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
  ## when center_result contains iv_count
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
## 6. extract final_result + numbers of instruments
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
## 7. save center diagnostics of binary CON-IV
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
  ## variable selection diagnostics table
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
  ## center estimates
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
## Summary function for type I error and power
##
## beta_XY = 0:
##   rejection = type I error
##
## beta_XY > 0:
##   rejection = power
############################################################

resultt_FUn_binary <- function(all_result,
                               true_effect_vec) {
  
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
    
    coverage = NA_real_,
    rejection = NA_real_,
    
    coverage_z = NA_real_,
    rejection_z = NA_real_,
    
    coverage_t = NA_real_,
    rejection_t = NA_real_,
    
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
  
  true_effect_vec <- suppressWarnings(as.numeric(true_effect_vec))
  
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
  ##########################################################
  
  estimate <- get_col_numeric(df, "estimate")
  se <- get_col_numeric(df, "se")
  
  estimate_old <- get_col_numeric(df, "b")
  
  replace_idx <- !is.finite(estimate) & is.finite(estimate_old)
  estimate[replace_idx] <- estimate_old[replace_idx]
  
  p_value <- get_col_numeric(df, "p_value")
  ci_lower <- get_col_numeric(df, "ci_lower")
  ci_upper <- get_col_numeric(df, "ci_upper")
  
  p_value_z <- get_col_numeric(df, "p_value_z")
  ci_lower_z <- get_col_numeric(df, "ci_lower_z")
  ci_upper_z <- get_col_numeric(df, "ci_upper_z")
  
  p_value_t <- get_col_numeric(df, "p_value_t")
  ci_lower_t <- get_col_numeric(df, "ci_lower_t")
  ci_upper_t <- get_col_numeric(df, "ci_upper_t")
  
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
  
  p_value_v <- p_value[valid]
  ci_lower_v <- ci_lower[valid]
  ci_upper_v <- ci_upper[valid]
  
  p_value_z_v <- p_value_z[valid]
  ci_lower_z_v <- ci_lower_z[valid]
  ci_upper_z_v <- ci_upper_z[valid]
  
  p_value_t_v <- p_value_t[valid]
  ci_lower_t_v <- ci_lower_t[valid]
  ci_upper_t_v <- ci_upper_t[valid]
  
  bias_v <- est_v - true_v
  
  ##########################################################
  ## coverage and rejection
  ##########################################################
  
  cover_main <- ci_lower_v <= true_v &
    ci_upper_v >= true_v
  
  reject_main <- p_value_v < 0.05
  
  cover_z <- ci_lower_z_v <= true_v &
    ci_upper_z_v >= true_v
  
  reject_z <- p_value_z_v < 0.05
  
  cover_t <- ci_lower_t_v <= true_v &
    ci_upper_t_v >= true_v
  
  reject_t <- p_value_t_v < 0.05
  
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
    
    coverage = safe_mean(as.numeric(cover_main)),
    rejection = safe_mean(as.numeric(reject_main)),
    
    coverage_z = safe_mean(as.numeric(cover_z)),
    rejection_z = safe_mean(as.numeric(reject_z)),
    
    coverage_t = safe_mean(as.numeric(cover_t)),
    rejection_t = safe_mean(as.numeric(reject_t)),
    
    selected_mean_iv_count = safe_mean(selected_mean_iv_count),
    selected_median_iv_count = safe_mean(selected_median_iv_count),
    selected_min_iv_count = safe_mean(selected_min_iv_count),
    selected_max_iv_count = safe_max(selected_max_iv_count),
    
    used_mean_iv_count = safe_mean(used_mean_iv_count),
    used_median_iv_count = safe_mean(used_median_iv_count),
    used_min_iv_count = safe_mean(used_min_iv_count),
    used_max_iv_count = safe_max(used_max_iv_count)
  )
  
  return(out)
}

############################################################
## Helper functions used in All_simulation_binary()
##
## Put this block before All_simulation_binary()
############################################################


############################################################
## 1. Empty standard result
############################################################

empty_standard_result <- function() {
  
  return(c(
    estimate = NA_real_,
    b = NA_real_,
    se = NA_real_,
    
    p_value = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    
    p_value_z = NA_real_,
    ci_lower_z = NA_real_,
    ci_upper_z = NA_real_,
    
    p_value_t = NA_real_,
    ci_lower_t = NA_real_,
    ci_upper_t = NA_real_
  ))
}


############################################################
## 2. Extract one numeric value from different object formats
############################################################

extract_numeric_value <- function(obj,
                                  candidate_names) {
  
  if (is.null(obj)) {
    return(NA_real_)
  }
  
  ##########################################################
  ## data.frame
  ##########################################################
  
  if (is.data.frame(obj)) {
    
    for (nm in candidate_names) {
      if (nm %in% colnames(obj)) {
        val <- suppressWarnings(as.numeric(obj[[nm]][1]))
        if (length(val) > 0 && is.finite(val)) {
          return(val)
        }
      }
    }
    
    return(NA_real_)
  }
  
  ##########################################################
  ## named vector
  ##########################################################
  
  if (is.atomic(obj) && !is.null(names(obj))) {
    
    for (nm in candidate_names) {
      if (nm %in% names(obj)) {
        val <- suppressWarnings(as.numeric(obj[[nm]][1]))
        if (length(val) > 0 && is.finite(val)) {
          return(val)
        }
      }
    }
    
    return(NA_real_)
  }
  
  ##########################################################
  ## list
  ##########################################################
  
  if (is.list(obj)) {
    
    for (nm in candidate_names) {
      if (nm %in% names(obj)) {
        val <- suppressWarnings(as.numeric(obj[[nm]][1]))
        if (length(val) > 0 && is.finite(val)) {
          return(val)
        }
      }
    }
    
    ## common nested result objects
    nested_names <- c(
      "final",
      "result",
      "summary",
      "estimate_result",
      "gls_result"
    )
    
    for (nested_nm in nested_names) {
      if (nested_nm %in% names(obj)) {
        val <- extract_numeric_value(
          obj = obj[[nested_nm]],
          candidate_names = candidate_names
        )
        if (is.finite(val)) {
          return(val)
        }
      }
    }
  }
  
  return(NA_real_)
}


############################################################
## Convert ADJ / PSM / IPW result into standard format
##
## Important:
##   Traditional methods return:
##     list(final_result = ..., center_result = ..., fit = ...)
##
##   Therefore we must extract result_obj$final_result.
############################################################

get_standard_final <- function(result_obj) {
  
  empty_out <- empty_standard_result()
  
  if (is.null(result_obj) || inherits(result_obj, "try-error")) {
    return(empty_out)
  }
  
  ##########################################################
  ## 1. Extract final_result from traditional method object
  ##########################################################
  
  if (is.list(result_obj) && !is.null(result_obj$final_result)) {
    final <- result_obj$final_result
  } else {
    final <- result_obj
  }
  
  if (is.null(final)) {
    return(empty_out)
  }
  
  final_num <- suppressWarnings(as.numeric(final))
  
  if (length(final_num) == 0) {
    return(empty_out)
  }
  
  if (!is.null(names(final))) {
    names(final_num) <- names(final)
  }
  
  ##########################################################
  ## 2. Extract estimate
  ##########################################################
  
  get_one <- function(candidate_names) {
    
    if (is.null(names(final_num))) {
      return(NA_real_)
    }
    
    for (nm in candidate_names) {
      if (nm %in% names(final_num)) {
        val <- suppressWarnings(as.numeric(final_num[[nm]][1]))
        if (length(val) > 0 && is.finite(val)) {
          return(val)
        }
      }
    }
    
    return(NA_real_)
  }
  
  estimate <- get_one(c(
    "estimate",
    "theta",
    "theta_gls",
    "b",
    "beta",
    "coef",
    "effect",
    "rd",
    "RD",
    "risk_difference"
  ))
  
  se <- get_one(c(
    "se",
    "se_gls",
    "std_error",
    "standard_error",
    "SE"
  ))
  
  p_value <- get_one(c(
    "p_value",
    "p",
    "P",
    "pval",
    "p_val"
  ))
  
  ci_lower <- get_one(c(
    "ci_lower",
    "lower",
    "lcl",
    "LCL",
    "lower95",
    "ci_l",
    "CI_lower"
  ))
  
  ci_upper <- get_one(c(
    "ci_upper",
    "upper",
    "ucl",
    "UCL",
    "upper95",
    "ci_u",
    "CI_upper"
  ))
  
  p_value_z <- get_one(c("p_value_z"))
  ci_lower_z <- get_one(c("ci_lower_z"))
  ci_upper_z <- get_one(c("ci_upper_z"))
  
  p_value_t <- get_one(c("p_value_t"))
  ci_lower_t <- get_one(c("ci_lower_t"))
  ci_upper_t <- get_one(c("ci_upper_t"))
  
  ##########################################################
  ## 3. If p value or CI is missing, calculate by normal CI
  ##########################################################
  
  if (is.finite(estimate) && is.finite(se) && se > 0) {
    
    z_val <- estimate / se
    
    if (!is.finite(p_value)) {
      p_value <- 2 * pnorm(abs(z_val), lower.tail = FALSE)
    }
    
    if (!is.finite(ci_lower)) {
      ci_lower <- estimate - 1.96 * se
    }
    
    if (!is.finite(ci_upper)) {
      ci_upper <- estimate + 1.96 * se
    }
    
    if (!is.finite(p_value_z)) {
      p_value_z <- 2 * pnorm(abs(z_val), lower.tail = FALSE)
    }
    
    if (!is.finite(ci_lower_z)) {
      ci_lower_z <- estimate - 1.96 * se
    }
    
    if (!is.finite(ci_upper_z)) {
      ci_upper_z <- estimate + 1.96 * se
    }
  }
  
  ##########################################################
  ## 4. Return standard result
  ##########################################################
  
  out <- c(
    estimate = estimate,
    b = estimate,
    se = se,
    
    p_value = p_value,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    
    p_value_z = p_value_z,
    ci_lower_z = ci_lower_z,
    ci_upper_z = ci_upper_z,
    
    p_value_t = p_value_t,
    ci_lower_t = ci_lower_t,
    ci_upper_t = ci_upper_t
  )
  
  return(out)
}


############################################################
## 4. Check whether one standard final result is valid
##
## This function is used in All_simulation_binary()
############################################################

is_ok <- function(x) {
  
  if (is.null(x)) {
    return(FALSE)
  }
  
  if (inherits(x, "try-error")) {
    return(FALSE)
  }
  
  x <- suppressWarnings(as.numeric(x))
  
  if (length(x) == 0) {
    return(FALSE)
  }
  
  if (!"estimate" %in% names(x)) {
    return(any(is.finite(x)))
  }
  
  return(
    is.finite(x["estimate"]) &&
      is.finite(x["se"]) &&
      x["se"] > 0
  )
}


############################################################
## 5. Empty covariate list
##
## This is used when Select_other() fails.
############################################################

empty_covariate_list <- function(N_center) {
  
  out <- vector("list", N_center)
  
  for (i in seq_len(N_center)) {
    out[[paste0("center_", i)]] <- character(0)
  }
  
  return(out)
}
############################################################
## Main simulation function for binary X and binary Y
##
## Output:
##   rejection:
##     beta_XY = 0  -> type I error
##     beta_XY > 0  -> power
##
##   coverage:
##     coverage probability of 95% CI
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
                                  
                                  ## binary latent model parameters
                                  target_px_ob = 0.50,
                                  target_py0 = 0.20,
                                  conf_strength_x = 0.50,
                                  conf_strength_y = 0.50,
                                  error_sd_x = 1,
                                  error_sd_y = 1,
                                  
                                  ## CON-IV screening parameters
                                  correlation_threshold = 0.10,
                                  pi_min = 0.02,
                                  F_min = 10,
                                  p_y_cutoff = 0.10,
                                  y_adjust_x = TRUE,
                                  multi_iv_threshold = 5,
                                  max_iv = Inf,
                                  weight_type_multi = "diag_delta",
                                  
                                  ## conventional method parameters
                                  ps_trim = 0.01,
                                  caliper_sd = 0.20,
                                  weight_trunc_quantile = c(0.01, 0.99),
                                  
                                  ## robust variance parameters
                                  hc_type = "HC1",
                                  min_n = 20,
                                  
                                  ## unused arguments
                                  n_boot = NULL,
                                  trim = NULL) {
  
  if (type != 2) {
    stop("All_simulation_binary() is designed for type = 2 binary simulations.")
  }
  
  ##########################################################
  ## 1. Store simulation results
  ##########################################################
  
  all_result_IV_single_same <- NULL
  all_result_IV_single_diff <- NULL
  all_result_IV_multi_same <- NULL
  all_result_IV_multi_diff <- NULL
  
  all_result_ADJ <- NULL
  all_result_PSM <- NULL
  all_result_IPW <- NULL
  
  all_tool_info_IV <- NULL
  all_center_info_IV <- NULL
  all_sim_status <- NULL
  
  true_effect_vec <- rep(NA_real_, N_sim)
  
  ##########################################################
  ## 2. Simulation loop
  ##########################################################
  
  for (sim in 1:N_sim) {
    
    cat("sim =", sim, "\n")
    
    ########################################################
    ## 2.1 Generate data
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
        error_sd_x = error_sd_x,
        error_sd_y = error_sd_y,
        keep_latent = FALSE
      ),
      silent = TRUE
    )
    
    ########################################################
    ## 2.2 If data generation failed
    ########################################################
    
    if (inherits(data, "try-error")) {
      
      empty_standard <- get_standard_final(NULL)
      
      empty_IV <- c(
        empty_standard,
        selected_mean_iv_count = NA_real_,
        selected_median_iv_count = NA_real_,
        selected_min_iv_count = NA_real_,
        selected_max_iv_count = NA_real_,
        used_mean_iv_count = NA_real_,
        used_median_iv_count = NA_real_,
        used_min_iv_count = NA_real_,
        used_max_iv_count = NA_real_
      )
      
      all_result_IV_single_same <- rbind(all_result_IV_single_same, empty_IV)
      all_result_IV_single_diff <- rbind(all_result_IV_single_diff, empty_IV)
      all_result_IV_multi_same <- rbind(all_result_IV_multi_same, empty_IV)
      all_result_IV_multi_diff <- rbind(all_result_IV_multi_diff, empty_IV)
      
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
          stringsAsFactors = FALSE
        )
      )
      
      next
    }
    
    ########################################################
    ## 2.3 Extract true effect
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
    ## 2.4 Candidate covariates
    ########################################################
    
    covariates_all <- setdiff(
      colnames(data$all_data),
      c("X", "Y", "centre")
    )
    
    ########################################################
    ## 2.5 IV method 1:
    ##     single IV + common IV set
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
    
    final_IV_single_same <- extract_final_result_with_iv_count(
      result_obj = result_IV_single_same_q,
      selected_covariates = if (!inherits(result_IV_single_same_q, "try-error")) {
        result_IV_single_same_q$selection
      } else {
        NULL
      },
      N_center = N_center
    )
    
    all_result_IV_single_same <- rbind(
      all_result_IV_single_same,
      final_IV_single_same
    )
    
    ########################################################
    ## 2.6 IV method 2:
    ##     single IV + center-specific IV set
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
    
    final_IV_single_diff <- extract_final_result_with_iv_count(
      result_obj = result_IV_single_diff_q,
      selected_covariates = if (!inherits(result_IV_single_diff_q, "try-error")) {
        result_IV_single_diff_q$selection
      } else {
        NULL
      },
      N_center = N_center
    )
    
    all_result_IV_single_diff <- rbind(
      all_result_IV_single_diff,
      final_IV_single_diff
    )
    
    ########################################################
    ## 2.7 IV method 3:
    ##     multiple IVs + common IV set
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
        max_iv = max_iv,
        weight_type = weight_type_multi,
        hc_type = hc_type,
        min_n = min_n,
        return_selection_diagnostics = TRUE
      ),
      silent = TRUE
    )
    
    final_IV_multi_same <- extract_final_result_with_iv_count(
      result_obj = result_IV_multi_same_q,
      selected_covariates = if (!inherits(result_IV_multi_same_q, "try-error")) {
        result_IV_multi_same_q$selection
      } else {
        NULL
      },
      N_center = N_center
    )
    
    all_result_IV_multi_same <- rbind(
      all_result_IV_multi_same,
      final_IV_multi_same
    )
    
    ########################################################
    ## 2.8 IV method 4:
    ##     multiple IVs + center-specific IV set
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
        max_iv = max_iv,
        weight_type = weight_type_multi,
        hc_type = hc_type,
        min_n = min_n,
        return_selection_diagnostics = TRUE
      ),
      silent = TRUE
    )
    
    final_IV_multi_diff <- extract_final_result_with_iv_count(
      result_obj = result_IV_multi_diff_q,
      selected_covariates = if (!inherits(result_IV_multi_diff_q, "try-error")) {
        result_IV_multi_diff_q$selection
      } else {
        NULL
      },
      N_center = N_center
    )
    
    all_result_IV_multi_diff <- rbind(
      all_result_IV_multi_diff,
      final_IV_multi_diff
    )
    
    ########################################################
    ## 2.9 Save IV diagnostics
    ########################################################
    
    method_result_list <- list(
      IV_single_same = result_IV_single_same_q,
      IV_single_diff = result_IV_single_diff_q,
      IV_multi_same = result_IV_multi_same_q,
      IV_multi_diff = result_IV_multi_diff_q
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
    ## 2.10 Variable screening for conventional methods
    ##
    ## ADJ: out = "both"
    ## PSM/IPW: out = "X"
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
    ## 2.11 Conventional method: adjusted regression
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
    ## 2.12 Conventional method: PSM
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
    ## 2.13 Conventional method: IPW-MSM
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
    ## 2.14 Simulation status
    ########################################################
    
    all_sim_status <- rbind(
      all_sim_status,
      data.frame(
        sim = sim,
        data_success = TRUE,
        true_effect = true_effect,
        
        IV_single_same_success = is_ok(final_IV_single_same),
        IV_single_diff_success = is_ok(final_IV_single_diff),
        IV_multi_same_success = is_ok(final_IV_multi_same),
        IV_multi_diff_success = is_ok(final_IV_multi_diff),
        
        ADJ_success = is_ok(final_ADJ),
        PSM_success = is_ok(final_PSM),
        IPW_success = is_ok(final_IPW),
        
        stringsAsFactors = FALSE
      )
    )
  }
  
  ##########################################################
  ## 3. Summarize results for type I error and power
  ##########################################################
  
  end_result_IV_single_same <- resultt_FUn_binary(
    all_result_IV_single_same,
    true_effect_vec
  )
  
  end_result_IV_single_diff <- resultt_FUn_binary(
    all_result_IV_single_diff,
    true_effect_vec
  )
  
  end_result_IV_multi_same <- resultt_FUn_binary(
    all_result_IV_multi_same,
    true_effect_vec
  )
  
  end_result_IV_multi_diff <- resultt_FUn_binary(
    all_result_IV_multi_diff,
    true_effect_vec
  )
  
  end_result_ADJ <- resultt_FUn_binary(
    all_result_ADJ,
    true_effect_vec
  )
  
  end_result_PSM <- resultt_FUn_binary(
    all_result_PSM,
    true_effect_vec
  )
  
  end_result_IPW <- resultt_FUn_binary(
    all_result_IPW,
    true_effect_vec
  )
  
  ##########################################################
  ## 4. Return object
  ##########################################################
  
  out <- list()
  
  out$IV_single_same <- end_result_IV_single_same
  out$IV_single_diff <- end_result_IV_single_diff
  out$IV_multi_same <- end_result_IV_multi_same
  out$IV_multi_diff <- end_result_IV_multi_diff
  
  out$ADJ <- end_result_ADJ
  out$PSM <- end_result_PSM
  out$IPW <- end_result_IPW
  
  out$true_effect_vec <- true_effect_vec
  out$true_mean <- mean(true_effect_vec, na.rm = TRUE)
  
  out$all_result_IV_single_same <- all_result_IV_single_same
  out$all_result_IV_single_diff <- all_result_IV_single_diff
  out$all_result_IV_multi_same <- all_result_IV_multi_same
  out$all_result_IV_multi_diff <- all_result_IV_multi_diff
  
  out$all_result_ADJ <- all_result_ADJ
  out$all_result_PSM <- all_result_PSM
  out$all_result_IPW <- all_result_IPW
  
  out$tool_info_IV <- all_tool_info_IV
  out$center_info_IV <- all_center_info_IV
  out$sim_status <- all_sim_status
  
  out$multi_iv_threshold <- multi_iv_threshold
  out$max_iv <- max_iv
  out$weight_type_multi <- weight_type_multi
  
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
## output path
############################################################

out_dir <- "results"
dir.create(path.expand(out_dir), recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. basic parameter settings
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
## 2. data type settings
##
## type = 2 means X/Y are binary
## typeC = 1 means C/L are continuous variables
############################################################

type <- 2
typeC <- 1


############################################################
## 3. L-layer confounder settings
############################################################

conf_L <- sample(
  1:(n1 + n2 + n3),
  n3 + 2
)

p_confound <- m + k + length(conf_L)


############################################################
## 4. confounding-path parameters: moderate unmeasured confounding
##
## Target:
##   U contribution around 40%-50%
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
## 4.1 Diagnostics for confounding strength
############################################################

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
  r_L = r_L,
  connection_prob = connection_prob,
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


############################################################
## 5. beta_XY is not fixed in the parameter section
##
## beta_XY is assigned in the scenario loop:
##   beta_XY = 0   -> type I error
##   beta_XY > 0   -> power
############################################################

causal <- NA


############################################################
## 6. RCT / do data settings
############################################################

Nrct <- 10000
p_x <- 0.5


############################################################
## 7. binary latent-variable model parameters
############################################################

target_px_ob <- 0.50
target_py0 <- 0.20

conf_strength_x <- 1
conf_strength_y <- 1

error_sd_x <- 1
error_sd_y <- 1


############################################################
## 8. number of replicates
############################################################

N_sim <- 1000


############################################################
## 9. Binary-CONIV instrument selection parameters
############################################################

correlation_threshold <- 0.10

pi_min <- 0.02
F_min <- 10
p_y_cutoff <- 0.10
y_adjust_x <- TRUE

multi_iv_threshold <- 5
weight_type_multi <- "diag_delta"
max_iv <- Inf


############################################################
## 10. conventional method parameters: PSM and IPW-MSM
############################################################

ps_trim <- 0.01
caliper_sd <- 0.20
weight_trunc_quantile <- c(0.01, 0.99)


############################################################
## 11. common robust variance parameters
############################################################

hc_type <- "HC1"
min_n <- 20


############################################################
## 12. unused parameters
############################################################

n_boot <- NULL
selcetnum <- NULL
trim <- NULL
############################################################
## 13. Scenario settings
##
## This version only varies sample size per center.
##
## Fixed:
##   1. number of centers
##   2. causal effect beta_XY
##   3. confounding strength
##   4. true causal effect calculation method
##   5. all other model parameters
##
## Varying:
##   1. sample size per center
############################################################

############################################################
## Fixed number of centers
############################################################

N_center_fixed <- 10


############################################################
## Fixed causal effect
##
## For type I error:
##   beta_XY_fixed <- 0
##
## For power:
##   beta_XY_fixed <- 0.4
##
############################################################

beta_XY_fixed <- 0

test_target <- ifelse(
  beta_XY_fixed == 0,
  "type1_error",
  "power"
)


############################################################
## Varying sample size per center
##
## N_per_center_target is the target sample size for each center.
## Actual center sample sizes are generated from:
##   80% to 120% of N_per_center_target
############################################################

N_per_center_vec <- c(
  1000,
  3000,
  5000,
  8000,
  10000,
  15000,
  20000,
  30000,
  50000
)


############################################################
## Create scenario grid
############################################################

scenario_grid <- data.frame(
  scenario = seq_along(N_per_center_vec),
  N_per_center_target = N_per_center_vec,
  stringsAsFactors = FALSE
)

scenario_grid$N_center <- N_center_fixed
scenario_grid$beta_XY_latent <- beta_XY_fixed
scenario_grid$test_target <- test_target


############################################################
## 14. Store all scenario results
############################################################

all_sim_result_all <- NULL


############################################################
## 15. Loop over sample size per center only
############################################################

for (j in seq_len(nrow(scenario_grid))) {
  
  ##########################################################
  ## Current scenario
  ##########################################################
  
  N_center <- scenario_grid$N_center[j]
  N_per_center_target <- scenario_grid$N_per_center_target[j]
  beta_XY_current <- scenario_grid$beta_XY_latent[j]
  test_target_current <- scenario_grid$test_target[j]
  
  cat(
    "Scenario =", j,
    "target =", test_target_current,
    "N_center =", N_center,
    "N_per_center_target =", N_per_center_target,
    "beta_XY_latent =", beta_XY_current,
    "type = 2",
    "\n"
  )
  
  ##########################################################
  ## Generate sample size for each center
  ##
  ## The number of centers is fixed.
  ## Each center sample size varies around the current target:
  ##   80% to 120% of N_per_center_target
  ##########################################################
  
  Nall <- floor(
    runif(
      N_center,
      N_per_center_target * 0.8,
      N_per_center_target * 1.2
    )
  )
  
  ##########################################################
  ## Run binary simulation
  ##
  ## Note:
  ##   The true causal effect calculation inside
  ##########################################################
  
  result_sim_result <- All_simulation_binary(
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
    beta_XY = beta_XY_current,
    
    Clv = Clv,
    p_x = p_x,
    
    target_px_ob = target_px_ob,
    target_py0 = target_py0,
    conf_strength_x = conf_strength_x,
    conf_strength_y = conf_strength_y,
    error_sd_x = error_sd_x,
    error_sd_y = error_sd_y,
    
    correlation_threshold = correlation_threshold,
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
    
    ########################################################
    ## Scenario information
    ########################################################
    
    scenario = j,
    scenario_type = "N_per_center",
    test_target = test_target_current,
    beta_XY_latent = beta_XY_current,
    
    ########################################################
    ## True effect information
    ##
    ########################################################
    
    true_mean = result_sim_result$true_mean,
    true_min = min(result_sim_result$true_effect_vec, na.rm = TRUE),
    true_max = max(result_sim_result$true_effect_vec, na.rm = TRUE),
    
    ########################################################
    ## Design information
    ########################################################
    
    N_center = N_center,
    N_per_center_target = N_per_center_target,
    
    Nall_mean = mean(Nall),
    Nall_min = min(Nall),
    Nall_max = max(Nall),
    Nall_total = sum(Nall),
    
    ########################################################
    ## Confounding structure
    ########################################################
    
    m = m,
    k = k,
    n1 = n1,
    n2 = n2,
    n3 = n3,
    
    confounding_diag_vec,
    
    ########################################################
    ## Data-generation parameters
    ########################################################
    
    typeC = typeC,
    type = type,
    
    Nrct = Nrct,
    p_x = p_x,
    
    r = r,
    UD = as.numeric(UD),
    Clv = Clv,
    r_L = r_L,
    connection_prob = connection_prob,
    
    target_px_ob = target_px_ob,
    target_py0 = target_py0,
    conf_strength_x = conf_strength_x,
    conf_strength_y = conf_strength_y,
    error_sd_x = error_sd_x,
    error_sd_y = error_sd_y,
    
    ########################################################
    ## Screening and estimation parameters
    ########################################################
    
    correlation_threshold = correlation_threshold,
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
  
  ##########################################################
  ## Save running result
  ##########################################################
  
  write.csv(
    all_sim_result_all,
    file = file.path(
      path.expand(out_dir),
      paste0(
        "Sim_only_Npercenter_",
        test_target,
        "_U40to50_type2_binary_running5.csv"
      )
    ),
    row.names = FALSE
  )
  
  # ##########################################################
  # ## Save detailed result for current scenario
  # ##########################################################
  # 
  # saveRDS(
  #   result_sim_result,
  #   file = file.path(
  #     path.expand(out_dir),
  #     paste0(
  #       "Detail_only_Npercenter_",
  #       test_target,
  #       "_Ncenter_",
  #       N_center,
  #       "_Nper_",
  #       N_per_center_target,
  #       "_betaXY_",
  #       beta_XY_current,
  #       "_scenario_",
  #       j,
  #       ".rds"
  #     )
  #   )
  # )
}


############################################################
## 16. Save final summary result
############################################################

write.csv(
  all_sim_result_all,
  file = file.path(
    path.expand(out_dir),
    paste0(
      "Sim_only_Npercenter_",
      test_target,
      "_U40to50_type2_binary.csv"
    )
  ),
  row.names = FALSE
)

# saveRDS(
#   all_sim_result_all,
#   file = file.path(
#     path.expand(out_dir),
#     paste0(
#       "Sim_only_Npercenter_",
#       test_target,
#       "_U40to50_type2_binary.rds"
#     )
#   )
# )