extract_gridded_x <- function(x) {

  # Example 1: One irregular value
  # x <- expand.grid(V1 = 11:15, V2 = 101:105)
  # x <- rbind(x, data.frame(V1 = 16, V2 = 101))

  # Example 2: Two irregular values
  # x <- expand.grid(V1 = 11:15, V2 = 101:105)
  # x <- rbind(x, data.frame(V1 = 16, V2 = 101))
  # x <- rbind(x, data.frame(V1 = 17, V2 = 111))

  # Example 3: One duplicated value
  # For this data, the returned values should be all TRUE;
  # the duplicates will be processed in the next step of mapprox().
  # x <- expand.grid(V1 = 11:15, V2 = 101:105)
  # x <- rbind(x, data.frame(V1 = 11, V2 = 101))

  should_be_omitted <- rep(FALSE, nrow(x))
  x_now <- x

  repeat {

    # calculate the NA ratio for each data values
    data_not_exist <- tapply(x_now, as.list(x_now), function(x0) nrow(x0) == 0)
    if (!any(data_not_exist)) break
    na_ratio <- lapply(
      seq(ncol(x_now)), function(i) apply(data_not_exist, i, mean))

    # finding indices of the value with the most NAs
    maxna_index_var <- which.max(sapply(na_ratio, max))
    maxna_index_val <- which.max(na_ratio[[maxna_index_var]])

    # removing data of a value with the most NAs
    x_uniq <- lapply(x_now, unique)
    maxna_val <- x_uniq[[maxna_index_var]][maxna_index_val]
    should_be_omitted <-
      should_be_omitted | x[, maxna_index_var] == maxna_val
    x_now <- x[!should_be_omitted, , drop = FALSE]

  }

  return(!should_be_omitted)

}

detect_duplicated_x <- function(x) {

  # Example 1: Two duplicated value
  # x <- expand.grid(V1 = 11:15, V2 = 101:105)
  # x <- rbind(data.frame(V1 = 11, V2 = 101), x, data.frame(V1 = 15, V2 = 105))

  # Example 2: One variable
  # x <- data.frame(V1 = c(1, 1:10, 10))

  duplication_id <- rep(NA, nrow(x))

  # Detecting duplication
  n_data <- tapply(x, as.list(x), nrow)
  if (all(n_data == 1)) return(duplication_id)
  is_duplicated <- which(n_data > 1, arr.ind = TRUE)

  # Adding duplication ID
  for (i in seq(nrow(is_duplicated))) {
    # i <- 1

    vars_i <- mapply(
      FUN = function(name, index) as.numeric(name[index]),
      name = dimnames(n_data), index = is_duplicated[i, ], SIMPLIFY = FALSE)

    sel_i <- rep(TRUE, nrow(x))
    for (j in seq(ncol(x))) {
      # j <- 1
      sel_i <- sel_i & (x[, j] == vars_i[[j]])
    }

    duplication_id[sel_i] <- i

  }

  return(duplication_id)

}
