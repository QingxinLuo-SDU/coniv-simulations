
############################################################
## 0. set the number of CPU threads used by each R task
## must be placed at the top of the script
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
## 1. generate data
############################################################

######## classification function ############################
Fenlei <- function(data, Clv, na.rm = TRUE) {
  
  if (Clv <= 0 || Clv >= 1) {
    stop("Clv must be between 0 and 1.")
  }
  
  cutoff <- quantile(data, probs = Clv, na.rm = na.rm, names = FALSE)
  
  newdata <- ifelse(data > cutoff, 1, 0)
  
  return(newdata)
}


## first generate one fixed L-layer network structure ###########


generate_layer_network <- function(parent_indices,
                                   current_layer_indices,
                                   connection_prob) {
  
  if (connection_prob < 0 || connection_prob > 1) {
    stop("connection_prob must be between 0 and 1.")
  }
  
  layer_network <- vector("list", length(current_layer_indices))
  
  for (j in seq_along(current_layer_indices)) {
    
    child_idx <- current_layer_indices[j]
    
    # fix which parents each child is connected to
    is_parent <- runif(length(parent_indices)) < connection_prob
    active_parent_indices <- parent_indices[is_parent]
    
    if (length(active_parent_indices) > 0) {
      
      # fix the edge weights
      weights <- runif(length(active_parent_indices), min = -1, max = 1)
      
      # normalize the weights so that the variance does not grow with the number of parents
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



## generate one layer of data from the fixed network structure #########################


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


## generate single observational center data ##########

Data_Generated <- function(type = 1, typeC, UD = TRUE, r, N, m, k,
                           n1, n2, n3, r_L, connection_prob, conf_L,
                           beta_CX, beta_CY, beta_XY,
                           Clv,
                           L_network = NULL) {
  
  NC <- m + k
  
  ##########################################################
  ## 1. generate baseline C/U covariates
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
  ## 2. initialize the L-layer data matrix
  ##########################################################
  
  total_vars <- n1 + n2 + n3
  
  simulated_data <- matrix(0, nrow = N, ncol = total_vars)
  colnames(simulated_data) <- c(
    paste0("L1_", 1:n1),
    paste0("L2_", 1:n2),
    paste0("L3_", 1:n3)
  )
  
  ##########################################################
  ## 3. generate L1
  ##########################################################
  
  Sigma_L <- matrix(r_L, nrow = n1, ncol = n1)
  diag(Sigma_L) <- 1
  
  mu_L <- rep(0, n1)
  
  l1_data <- MASS::mvrnorm(n = N, mu = mu_L, Sigma = Sigma_L)
  simulated_data[, 1:n1] <- l1_data
  
  ##########################################################
  ## 4. generate L2 and L3 from the fixed network
  ##########################################################
  
  l1_indices <- 1:n1
  l2_indices <- (n1 + 1):(n1 + n2)
  l3_indices <- (n1 + n2 + 1):total_vars
  
  # if no L_network is supplied, generate one here.
  # ALL_data() passes the same L_network to all data sets.
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
  ## 5. combine covariates
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
  ## 6. generate X and Y
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
  ## 7. drop the unobservable U
  ##########################################################
  
  cols_to_remove <- paste0("U", 1:k)
  
  data_total_end <- data_total[
    ,
    -which(names(data_total) %in% cols_to_remove),
    drop = FALSE
  ]
  
  return(data_total_end)
}


## generate RCT / do data: fixed L network version ########


Data_Generated_do <- function(type = 1, typeC, UD = TRUE, r, N, m, k,
                              n1, n2, n3, r_L, connection_prob, conf_L,
                              beta_CX, beta_CY, beta_XY,
                              Clv, p_x,
                              L_network = NULL) {
  
  NC <- m + k
  
  ##########################################################
  ## 1. generate baseline C/U covariates
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
  ## 2. initialize the L-layer data matrix
  ##########################################################
  
  total_vars <- n1 + n2 + n3
  
  simulated_data <- matrix(0, nrow = N, ncol = total_vars)
  colnames(simulated_data) <- c(
    paste0("L1_", 1:n1),
    paste0("L2_", 1:n2),
    paste0("L3_", 1:n3)
  )
  
  ##########################################################
  ## 3. generate L1
  ##########################################################
  
  Sigma_L <- matrix(r_L, nrow = n1, ncol = n1)
  diag(Sigma_L) <- 1
  
  mu_L <- rep(0, n1)
  
  l1_data <- MASS::mvrnorm(n = N, mu = mu_L, Sigma = Sigma_L)
  simulated_data[, 1:n1] <- l1_data
  
  ##########################################################
  ## 4. generate L2 and L3 from the fixed network
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
  ## 5. combine covariates
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
  ## 6. generate X and Y
  ##########################################################
  
  # In the do / RCT data, X is controlled by p_x
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
  ## 7. drop the unobservable U
  ##########################################################
  
  cols_to_remove <- paste0("U", 1:k)
  
  data_total_end <- data_total[
    ,
    -which(names(data_total) %in% cols_to_remove),
    drop = FALSE
  ]
  
  return(data_total_end)
}



## overall data generation: all centers and RCT/do share the same L_network #######

ALL_data <- function(type, N_center, typeC, UD, r, Nall, Nrct, m, k,
                     n1, n2, n3, r_L, connection_prob, conf_L,
                     beta_CX, beta_CY, beta_XY,
                     Clv, p_x) {
  
  ##########################################################
  ## 1. first fix one L network structure
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
  ## 2. generate several observational centers
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
  ## 3. generate RCT / do data
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
  ## 4. return the data and the fixed network
  ##########################################################
  
  return(list(
    all_data = all_data,
    rct = rct_data,
    L_network = L_network
  ))
}



## 2 CON_LV
############################################################

### safely extract the P value of a variable from lm ###########

safe_lm_pvalue <- function(formula, data, term_name) {
  
  fit <- try(lm(formula, data = data), silent = TRUE)
  
  if (inherits(fit, "try-error")) {
    return(NA_real_)
  }
  
  coefs <- try(summary(fit)$coefficients, silent = TRUE)
  
  if (inherits(coefs, "try-error")) {
    return(NA_real_)
  }
  
  # return NA if the variable was dropped from the model (e.g., constant or collinear)
  if (!term_name %in% rownames(coefs)) {
    return(NA_real_)
  }
  
  return(coefs[term_name, 4])
}

## safe absolute correlation #####


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



##  unified variable selection function
##
## selcetnum:
##   selcetnum = 1       : each center selects only the strongest instrument
##   selcetnum = "multi" : keep several instruments, following the original function
##   selcetnum = 3       : optional; keep at most the first 3


variable_selection_sorted <- function(N_center,
                                      covariates,
                                      data,
                                      correlation_threshold,
                                      type = 1,
                                      selcetnum = "multi") {
  ## 0. basic checks
  
  if (length(covariates) == 0) {
    stop("covariates is empty.")
  }
  
  if (correlation_threshold < 0 || correlation_threshold > 1) {
    stop("correlation_threshold must be between 0 and 1.")
  }
  
  if (!(type %in% c(1, 2))) {
    stop("type must be 1 or 2.")
  }
  
  # selcetnum can be:
  #   1       : always select only the strongest instrument
  #   number >1 : multiple-instrument threshold; if the number of candidates <= selcetnum, select only the strongest one;
  #             if the number of candidates > selcetnum, keep several instruments
  #   "multi" : full multiple-instrument mode; keep all selected variables
  if (!(is.numeric(selcetnum) || selcetnum == "multi")) {
    stop("selcetnum must be 1, a positive integer, or 'multi'.")
  }
  
  if (is.numeric(selcetnum)) {
    if (length(selcetnum) != 1 || selcetnum < 1 || selcetnum != floor(selcetnum)) {
      stop("numeric selcetnum must be a single positive integer.")
    }
  }
  
  all_final_vars <- list()
  
  ## 1. set the P-value threshold by type
  
  if (type == 1) {
    # continuous: Bonferroni correction
    p_cutoff <- 0.05 / length(covariates)
  } else {
    # binary: threshold 0.1
    p_cutoff <- 0.1
  }
  
  ## select variables separately in each center
  for (i in 1:N_center) {
    
    ## 2.1 extract the data of center i
    data_i <- data$all_data[
      data$all_data$centre == i,
      ,
      drop = FALSE
    ]
    
    # keep only candidate variables present in both the current center and the RCT
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
    
    ## 2.2 initial selection: associated with both Y and X
    
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
    
    # variables significantly associated with Y
    varable1 <- covariates_i[
      which(!is.na(p_value_Y) & p_value_Y < p_cutoff)
    ]
    
    # variables significantly associated with X
    varable2 <- covariates_i[
      which(!is.na(p_value_X) & p_value_X < p_cutoff)
    ]
    
    # variables associated with both Y and X
    common_variables <- intersect(varable1, varable2)
    
    if (length(common_variables) == 0) {
      all_final_vars[[paste0("center_", i)]] <- character(0)
      next
    }
    
    final_vars <- common_variables
    
    ## 2.3 handle highly correlated variables
    ##
    ## if the correlation between two variables exceeds correlation_threshold,
    ## then keep the one more strongly correlated with X.
    
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
    
    ## 2.4 sort by absolute correlation with X in decreasing order
    if (length(final_vars) > 1) {
      
      correlations_with_X <- sapply(final_vars, function(var) {
        safe_abs_cor(data_i$X, data_i[[var]])
      })
      
      final_vars <- final_vars[
        order(correlations_with_X, decreasing = TRUE)
      ]
    }
    
    
    ## 2.5 keep only variables also present in the RCT data
    
    final_vars_end <- final_vars[
      final_vars %in% colnames(data$rct)
    ]
    
    ########################################################
    ## 2.6 use selcetnum to choose one or several
    ##
    ## selcetnum = 1:
    ##   always select only the strongest instrument.
    ##
    ## selcetnum = 5:
    ##   if the number of selected instruments <= 5, select only the strongest one;
    ##   if the number of selected instruments > 5, keep several instruments.
    ##
    ## selcetnum = "multi":
    ##   full multiple-instrument mode; keep all selected variables.
    ########################################################
    
    if (length(final_vars_end) > 0) {
      
      if (is.numeric(selcetnum)) {
        
        if (selcetnum == 1) {
          
          # single-instrument mode: select only the strongest one
          final_vars_end <- final_vars_end[1]
          
        } else if (selcetnum > 1) {
          
          # multiple-instrument threshold mode:
          # if the number of candidate instruments <= selcetnum, fall back to the strongest single instrument;
          # if the number of candidate instruments > selcetnum, keep several instruments.
          if (length(final_vars_end) <= selcetnum) {
            final_vars_end <- final_vars_end[1]
          } else {
            final_vars_end <- final_vars_end
          }
        }
        
      } else if (selcetnum == "multi") {
        
        # full multiple-instrument mode: no threshold
        final_vars_end <- final_vars_end
      }
    }
    
    all_final_vars[[paste0("center_", i)]] <- final_vars_end
  }
  
  return(all_final_vars)
}




## covariance matrix of the Y ~ Z slopes of several instruments in the RCT

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
    # so the influence contribution of each observation is:
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

## influence functions of the X ~ Z and Y ~ Z slopes in observational centers
## each column corresponds to one instrument; crossprod() gives the variances and covariances of the slope estimates

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

## single-center Wald ratio estimate: returns the information needed by the covariance meta-analysis

############################################################
## single-center single-instrument Wald ratio estimate
## return the information needed by the shared-RCT covariance meta-analysis + instrument diagnostics
############################################################

CON_IV_luo_top1_return <- function(data_rct_all,
                                   data_OB,
                                   IV_covariates,
                                   rct_info) {
  
  ##########################################################
  ## internal empty-result function: keep column names consistent for all failures
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
  ## 0. return NA when there is no instrument
  ##########################################################
  
  if (length(IV_covariates) == 0) {
    return(empty_result(
      IV_name = NA_character_,
      n_OB = n_OB
    ))
  }
  
  ##########################################################
  ## 1. the single-instrument version takes only the first instrument
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
  ## 2. estimate X ~ Z in observational centers
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
  ## 3. estimate Y ~ Z in observational centers
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
  ## 4. corresponding Y ~ Z slopes and variances in the RCT / do data
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

  ## variances and covariance of the two slopes in an observational center (influence functions),
  ## enters the delta method together with the variance of the shared RCT slopes
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
  ## 6. return single-center results
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
## single-instrument meta-analysis with the shared RCT covariance matrix
##
## used for:
##   ALLCONIV_cov()
##
## methods:
##   1. obtain the combined estimate and SE by fixed-effect GLS, as in meta_function_MR_cov_multi();
##   2. p values and 95% CIs use the t distribution;
##   3. degrees of freedom df = number of valid centers K - 1;
##   4. consistent with the multiple-instrument meta_function_MR_cov_multi().
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
  ## 1. filter out invalid centers
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
  ## 2. construct the covariance matrix V
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
  ## 3. numerical safeguard: apply nearPD if V is not positive definite
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
  ## 4. fixed-effect GLS combination (same as meta_function_MR_cov_multi())
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
  ## 6. t-distribution correction
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
  ## 7. final results
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





## within-center multiple-instrument IVW estimate
## return the information needed by the multiple-instrument covariance meta-analysis


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
  
  
  ## 0. basic checks
  
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
  
  
  ## 1. estimate each instrument separately:
  ##    X ~ Z
  ##    Y_OB ~ Z
  ##    the Y_RCT ~ Z results in the RCT come from rct_info
  
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
    
    
    ## Y ~ Z in the RCT
    
    betaY_RCT <- rct_info$beta[iv_j]
    varY_RCT <- rct_info$vcov[iv_j, iv_j]
    
    
    ## save
    
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
  
  
  ## 2. filter out invalid instruments
  
  
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
  
  
  ## 3. multiple-instrument IVW estimate
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
  
  
  ## 4. compute the linear weights of this center's estimator on betaYG_j
  ##
  ## theta_i = sum_j a_ij * betaYG_ij
  ##
  ## a_j = w_j * betaXG_j / sum(w_j * betaXG_j^2)
  
  
  a_vec <- w * betaXG / denom
  names(a_vec) <- iv_result$IV_name
  
  
  ## 5. variance
  ##
  ## obs_y_var: Y~Z slopes in an observational center, including covariances between instruments
  ## rct_var:   Y~Z slopes in the shared RCT, including covariances between instruments
  ## x_var:     X~Z slopes in an observational center, including covariances between instruments
  ## yx_cov:    covariance between the Y~Z and X~Z slopes within the same observational center

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

  ## derivative of betaXY with respect to betaXG_j
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
## multiple-instrument meta-analysis with the shared RCT covariance matrix
############################################################

############################################################
## multiple-instrument meta-analysis with the shared RCT covariance matrix
############################################################

meta_function_MR_cov_multi <- function(center_result_list,
                                       rct_vcov,
                                       use_nearPD = TRUE,
                                       use_t_correction = TRUE) {
  
  ##########################################################
  ## 1. extract valid centers
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
  ## 2. construct the between-center covariance matrix V
  ##########################################################
  
  V <- matrix(0, nrow = K, ncol = K)
  
  ## diagonal: variance of each center's own estimate
  diag(V) <- se^2
  
  ## off-diagonal: between-center covariance due to shared RCT reference error
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
  ## 3. numerical safeguard: apply nearPD if V is not positive definite
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
  ## 4. fixed-effect GLS combination
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
  ## 5. small-sample t-distribution correction
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
  ## 6. final results
  ##
  ## the third element is named z_value;
  ## it is a t statistic when use_t_correction = TRUE.
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
  ## 7. organize center results
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
## multicenter multiple instruments + meta-analysis with the shared RCT covariance matrix
############################################################

ALLCONIV_cov_multi <- function(ALLdata,
                               N_center,
                               type,
                               covariates,
                               n_boot = NULL) {
  
  ##########################################################
  ## 0. this version currently supports only continuous type = 1
  ##########################################################
  
  if (type != 1) {
    stop("ALLCONIV_cov_multi currently supports type == 1 only.")
  }
  
  ##########################################################
  ## 1. collect all instruments used by all centers
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
  ## 2. estimate the covariance matrix of the Y~Z slopes of all instruments in the RCT
  ##########################################################
  
  rct_info <- rct_slope_vcov(
    data_rct = ALLdata$rct,
    iv_names = all_selected_IVs
  )
  
  ##########################################################
  ## 3. multiple-instrument IVW estimate in each center
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
  ## 4. multiple-instrument covariance-matrix meta-analysis
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
## 1. simulation result summary function
############################################################
############################################################
## simulation result summary function
##
## summary:
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
## also summarize the numbers of instruments:
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
  ## 1. extract estimates and model standard errors
  ##########################################################
  
  est <- as.numeric(result[, 1])
  model_se <- as.numeric(result[, 2])
  
  valid <- is.finite(est)
  
  if (sum(valid, na.rm = TRUE) == 0) {
    
    ## keep the instrument selection counts even when there is no valid estimate
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
  ## 2. estimate, bias, empirical SE, and MSE
  ##########################################################
  
  estimate <- mean(est_valid, na.rm = TRUE)
  bias <- estimate - true_effect
  emp_se <- sd(est_valid, na.rm = TRUE)
  
  mse <- mean(
    (est_valid - true_effect)^2,
    na.rm = TRUE
  )
  
  ##########################################################
  ## 3. ratio of model SE to empirical SE
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
  ## 5. success rate
  ##########################################################
  
  valid_rate <- mean(valid, na.rm = TRUE)
  n_valid <- sum(valid, na.rm = TRUE)
  
  ##########################################################
  ## 6. instrument count summary
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
  ## 7. return
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
## select a common instrument set in the largest center
## all centers use the same instrument set
##
## requirements:
##   1. instruments are selected in the center with the largest sample size;
##   2. instruments must be present in all observational centers;
##   3. instruments must be present in the RCT data;
##   4. usable for single or multiple instruments.
############################################################

variable_selection_largest_center_common <- function(N_center,
                                                     covariates,
                                                     data,
                                                     correlation_threshold,
                                                     type = 1,
                                                     selcetnum = "multi") {
  
  ##########################################################
  ## 0. basic checks
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
  ## 1. find the observational center with the largest sample size
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
  ## 2. keep only candidate variables observed in all observational centers and the RCT
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
  ## 3. set the P-value threshold for selection
  ##########################################################
  
  if (type == 1) {
    p_cutoff <- 0.05 / length(covariates_common)
  } else {
    p_cutoff <- 0.1
  }
  
  ##########################################################
  ## 4. select instruments in the center with the largest sample size
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
    ## 5. remove highly correlated redundant instruments
    ##    if two variables are highly correlated, keep the one more strongly correlated with X
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
    ## 6. sort by absolute correlation with X in decreasing order
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
    ## 7. ensure again that instruments are present in all centers and the RCT
    ########################################################
    
    selected_vars <- selected_vars[
      selected_vars %in% covariates_common
    ]
    
    ########################################################
    ## 8. use selcetnum to choose single or multiple instruments
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
  ## 9. all centers use the same instrument set
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
## single instrument + meta-analysis with the shared RCT covariance matrix
##
## applies to:
##   1. single instrument + the same instrument
##   2. single instrument + different instruments across centers
##
## note:
##   if a center has several candidate instruments in covariates,
##   this function takes only the first one, i.e., the strongest instrument.
############################################################

ALLCONIV_cov <- function(ALLdata,
                         N_center,
                         type,
                         covariates,
                         n_boot = NULL) {
  
  ##########################################################
  ## 0. this version supports only continuous type = 1
  ##########################################################
  
  if (type != 1) {
    stop("ALLCONIV_cov currently supports type == 1 only.")
  }
  
  ##########################################################
  ## 1. extract variable selection information
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
  ## 2. collect the instruments actually used by all centers
  ##    single-instrument version: each center takes only the first instrument
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
  ## 3. if there is no valid instrument, return an empty result directly
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
  ## 4. estimate the covariance matrix of the Y ~ Z slopes in the RCT / do data
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
  ## 5. single-instrument Wald ratio in each center
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
  ## 6. meta-analysis combination with the shared RCT covariance matrix
  ##########################################################
  
  final_meta <- meta_function_MR_cov(
    MR_result_all = MR_result_all,
    rct_vcov = rct_info$vcov
  )
  
  ##########################################################
  ## 7. summarize instrument diagnostics
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
  ## 8. return with additional information
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

## variable selection function for conventional observational methods ##############
## out = "both" : associated with both X and Y
## out = "X"    : select only variables associated with X
## out = "Y"    : select only variables associated with Y


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



## conventional observational regression adjustment ############
## adjust separately in each center, then combine by meta-analysis ############


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


## continuous exposure GPS adjustment ###############
## note: this is not PSM but GPS predicted-mean adjustment

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
    
    
    ## 1. exposure model: X ~ covariates
    
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
    
    
    ## 2. outcome model: Y ~ X + gps_mu
    
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


## stabilized IPTW for a binary exposure #############
##
## exposure    : exposure variable name, e.g., "X"
## denominator : denominator model formula, e.g., X ~ C1 + C2 + L1_1
## data        : data of the current center
## trim        : weight truncation proportion, default 0.01; set to NULL or 0 for no truncation
## eps         : probability truncation to avoid predicted probabilities of 0 or 1
##
## returns:
##   stabilized IPTW weights, same length as data


ipw_mew_cat <- function(exposure,
                        denominator,
                        data,
                        trim = 0.01,
                        eps = 1e-6) {
  
  
  ## 0. basic checks
  
  if (!exposure %in% colnames(data)) {
    stop("exposure is not found in data.")
  }
  
  A <- data[[exposure]]
  
  ## this function is only for a binary 0/1 exposure
  valid_A <- A[!is.na(A)]
  
  if (!all(valid_A %in% c(0, 1))) {
    stop("ipw_mew_cat currently requires binary exposure coded as 0/1.")
  }
  
  n <- nrow(data)
  
  if (n == 0) {
    return(numeric(0))
  }
  
  
  ## 1. Numerator model: P(A)
  ##    numerator of the stabilized weight, intercept only
  
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
  
  
  ## 3. compute stabilized IPTW
  ##
  ## A = 1: numerator = P(A=1), denominator = P(A=1 | C)
  ## A = 0: numerator = P(A=0), denominator = P(A=0 | C)
  
  
  w_num <- ifelse(A == 1, p_num, 1 - p_num)
  w_den <- ifelse(A == 1, p_den, 1 - p_den)
  
  weights <- w_num / w_den
  
  weights[!is.finite(weights)] <- NA_real_
  
  
  ## 4. truncate weights to avoid extreme values
  
  if (!is.null(trim) && trim > 0) {
    
    q <- quantile(
      weights,
      probs = c(trim, 1 - trim),
      na.rm = TRUE
    )
    
    weights <- pmin(pmax(weights, q[1]), q[2])
  }
  
  
  ## 5. attach simple diagnostics; does not affect use of the original function
  
  
  attr(weights, "diagnostics") <- list(
    mean_weight = mean(weights, na.rm = TRUE),
    sd_weight = sd(weights, na.rm = TRUE),
    min_weight = min(weights, na.rm = TRUE),
    max_weight = max(weights, na.rm = TRUE),
    n_missing_weight = sum(!is.finite(weights) | is.na(weights))
  )
  
  return(weights)
}

## stabilized generalized propensity score weight for a continuous exposure ########
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


## continuous exposure GPS-IPW MSM #######
## estimate weights and the MSM separately in each center, then combine by meta-analysis


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
    
    
    ## continuous exposure: GPS stabilized weight
    
    
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
      
      
      ## binary exposure: IPW
      
      
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





## safely extract final_result #######


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




## main simulation function: four proposed methods + regression adjustment + GPS adjustment + GPS-IPW MSM ##########

############################################################
## safely bind data.frames
##
## purpose:
##   center_result / iv_selection_info returned by different methods may have different numbers of columns.
##   rbind directly would raise an error:
##     Error in rbind(deparse.level, ...) : numbers of columns of arguments do not match
##
## this function fills missing columns automatically and then binds in a common column order.
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
## main simulation function:
## four proposed methods + regression adjustment + GPS adjustment + GPS-IPW MSM
##
## four proposed methods:
##   1. IV_single_same : single instrument + the same instrument across centers
##   2. IV_single_diff : single instrument + different instruments across centers
##   3. IV_multi_same  : multiple instruments + the same instrument set across centers
##   4. IV_multi_diff  : multiple instruments + different instruments across centers
##
## additional output:
##   selected_*_iv_count : number of instruments at the selection stage
##   used_*_iv_count     : number of instruments actually used in valid estimates
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
  ## internal helper 1: safely bind data.frames
  ## used to bind center_result / iv_selection_info returned by different methods
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
  ## internal helper 2: create an empty variable list
  ##########################################################
  
  empty_covariate_list <- function(N_center) {
    
    tmp <- vector("list", N_center)
    
    for (cen in 1:N_center) {
      tmp[[paste0("center_", cen)]] <- character(0)
    }
    
    return(tmp)
  }
  
  
  ##########################################################
  ## internal helper 3: count instruments at the selection stage from the selected-variable list
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
  ## internal helper 4: count instruments actually used from the estimation results
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
    ## case 1: multiple-instrument results
    ## center_result contains iv_count
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
    ## case 2: single-instrument results
    ## if center_result has valid centers, each valid center uses one instrument
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
  ## internal helper 5:
  ## extract final_result + selected instrument count + used instrument count
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
  ## internal helper 6: save instrument and center information of the proposed methods
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
    ## 1. save overall instrument selection information
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
    ## 2. save center-level results
    ## single- and multiple-instrument center_result columns may differ and must be bound safely
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
  ## 1. true effect
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
  ## 2. store the results of each replicate
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
  ## 3. start the simulation
  ##########################################################
  
  for (sim in 1:N_sim) {
    
    cat("sim =", sim, "\n")
    
    
    ########################################################
    ## 3.1 generate data
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
    ## if data generation fails
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
    ## 3.2 candidate variables
    ########################################################
    
    covariates_all <- setdiff(
      colnames(data$all_data),
      c("X", "Y", "centre")
    )
    
    
    ########################################################
    ## 3.3 variable selection for conventional methods
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
    ## 3.4 proposed method 1:
    ##     single instrument + the same instrument across centers
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
    ## 3.5 proposed method 2:
    ##     single instrument + different instruments across centers
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
    ## 3.6 proposed method 3:
    ##     multiple instruments + the same instrument set across centers
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
    ## 3.7 proposed method 4:
    ##     multiple instruments + different instruments across centers
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
    ## 3.8 save instrument and center information of the four proposed methods
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
    ## 3.9 method 5: conventional regression adjustment
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
    ## 3.10 method 6: GPS mean adjustment
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
    ## 3.11 method 7: GPS-IPW MSM
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
    ## 3.12 save the state of this replicate
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
  ## 4. summarize simulation results
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
  ## 5. return
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


## parameter settings ############

m <- 30
k <- 5

n1 <- 5
n2 <- 5
n3 <- 15

r <- 0
UD <- TRUE
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

n_boot <- 1000

selcetnum <- 1
trim <- 0.01

## here para_sim is the correlation threshold in instrument selection
para_sim <- seq(0.1, 0.5, 0.05)


## output path

out_dir <- "results"
dir.create(path.expand(out_dir), recursive = TRUE, showWarnings = FALSE)

all_sim_result_all <- NULL


############################################################
## start the loop
############################################################

for (j in seq_along(para_sim)) {
  
  cat("Scenario =", j, "correlation_threshold =", para_sim[j], "\n")
  
  correlation_threshold <- para_sim[j]
  
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
    para_sim = para_sim[j],
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
      "Sim_correlation_threshold_type1_all_methods_comparison_running.csv"
    ),
    row.names = FALSE
  )
}

write.csv(
  all_sim_result_all,
  file = file.path(
    path.expand(out_dir),
    "Sim_correlation_threshold_type1_all_methods_comparison.csv"
  ),
  row.names = FALSE
)