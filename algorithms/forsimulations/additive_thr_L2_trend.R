##Main function
## Functions for changes in the trend

additive_thr_L2_trend <- function(X, thr_const = 1.4,
                                  thr_fin = thr_const * sqrt(2 * log(nrow(X))),
                                  s = 1, e = nrow(X), points = 10, k_l = 1,
                                  k_r = 1) {
  if (!(is.matrix(X))) {
    stop("The input in `X' should be a numeric matrix, with each data sequence
        we want to investigate being a column of this matrix.")
  }
  if ((thr_const <= 0) || (points <= 0)) {
    stop("The threshold constant as well as the `points' argument that
        represents the magnitude of the expansion for the intervals should
        be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps^0.5) {
    warning("The input for `points' should be a positive integer.
        If it is a positive real number then the integer part of the
        given number is used as the value of `points'.")
  }

  points <- as.integer(points)
  l <- length(X[, 1])
  r_e_points <- seq(points, l, points)
  l_e_points <- seq(l - points + 1, 1, -points)
  chp <- 0
  if (e - s <= 1) {
    cpt <- 0
    noiseless.cpt<-0
  } else {
    pos_r <- numeric()
    CUSUM_r <- numeric()
    pos_l <- numeric()
    CUSUM_l <- numeric()
    moving_points <- IDetect::s_e_points(r_e_points, l_e_points, s, e)
    right_points <- moving_points[[1]]
    left_points <- moving_points[[2]]
    lur <- length(left_points)
    rur <- length(right_points)
    loc <- numeric()
    cpt <- 0
    noiseless.cpt <- 0
    for (i in seq_len(ncol(X))){
      if (stats::mad(diff(diff(X[,i]))/sqrt(6)) == 0){
        noiseless.cpt <- c(noiseless.cpt,all.slopechanges.are.cpts(X[,i])$cpts)
        loc <- c(loc,i)}
    }
    if (length(loc) == 0){
      loc <- ncol(X) + 1
    }
    if ((length(loc) < ncol(X)) | (loc == ncol(X) + 1)){
      if (k_r < k_l) {
        while ((chp == 0) & (k_r < min(k_l, rur))) {
          lxr <- right_points[k_r] - s + 1
          ipcr <- matrix(NA, lxr - 1, ncol(X))
          for (i in (1:ncol(X))[-loc]) {
            ipcr[, i] <- IDetect:::cumsum_lin(X[s:right_points[k_r], i])/stats::mad(diff(diff(X[,i]))/sqrt(6))
          }
          ipcr <- ipcr[, colSums(is.na(ipcr)) != nrow(ipcr), drop = FALSE]
          sipcr <- sqrt(rowSums(ipcr^2)) / sqrt(ncol(X))
          pos_r[k_r] <- which.max(sipcr) + s - 1
          CUSUM_r[k_r] <- sipcr[pos_r[k_r] - s + 1]
          if (CUSUM_r[k_r] > thr_fin) {
            chp <- pos_r[k_r]
          } else {
            k_r <- k_r + 1
          }
        }
      }
      if (k_l < k_r) {
        while ((chp == 0) & (k_l < min(k_r, lur))) {
          lxl <- e - left_points[k_l] + 1
          ipcl <- matrix(NA, lxl - 1, ncol(X))
          for (i in (1:ncol(X))[-loc]) {
            ipcl[, i] <- IDetect:::cumsum_lin(X[left_points[k_l]:e, i])/stats::mad(diff(diff(X[,i]))/sqrt(6))
          }
          ipcl <- ipcl[, colSums(is.na(ipcl)) != nrow(ipcl), drop = FALSE]
          sipcl <- sqrt(rowSums(ipcl^2)) / sqrt(ncol(X))
          pos_l[k_l] <- which.max(sipcl) + left_points[k_l] - 1
          CUSUM_l[k_l] <- sipcl[pos_l[k_l] - left_points[k_l] + 1]
          if (CUSUM_l[k_l] > thr_fin) {
            chp <- pos_l[k_l]
          } else {
            k_l <- k_l + 1
          }
        }
      }
      if (chp == 0) {
        while ((chp == 0) & (k_l <= lur) & (k_r <= rur)) {
          lxr <- right_points[k_r] - s + 1
          ipcr <- matrix(NA, lxr - 1, ncol(X))
          for (i in (1:(ncol(X)))[-loc]) {
            ipcr[, i] <- IDetect:::cumsum_lin(X[s:right_points[k_r], i])/stats::mad(diff(diff(X[,i]))/sqrt(6))
          }
          ipcr <- ipcr[, colSums(is.na(ipcr)) != nrow(ipcr), drop = FALSE]
          sipcr <- sqrt(rowSums(ipcr^2)) / sqrt(ncol(X))
          pos_r[k_r] <- which.max(sipcr) + s - 1
          CUSUM_r[k_r] <- sipcr[pos_r[k_r] - s + 1]
          if (CUSUM_r[k_r] > thr_fin) {
            chp <- pos_r[k_r]
          } else {
            lxl <- e - left_points[k_l] + 1
            ipcl <- matrix(NA, lxl - 1, ncol(X))
            for (i in (1:(ncol(X)))[-loc]) {
              ipcl[, i] <- IDetect:::cumsum_lin(X[left_points[k_l]:e, i])/stats::mad(diff(diff(X[,i]))/sqrt(6))
            }
            ipcl <- ipcl[, colSums(is.na(ipcl)) != nrow(ipcl), drop = FALSE]
            sipcl <- sqrt(rowSums(ipcl^2)) / sqrt(ncol(X))
            pos_l[k_l] <- which.max(sipcl) + left_points[k_l] - 1
            CUSUM_l[k_l] <- sipcl[pos_l[k_l] - left_points[k_l] + 1]
            if (CUSUM_l[k_l] > thr_fin) {
              chp <- pos_l[k_l]
            } else {
              k_r <- k_r + 1
              k_l <- k_l + 1
            }
          }
        }
      }
      if (chp != 0) {
        if (chp > ((e + s) / 2)) {
          r <- additive_thr_L2_trend(X, s = s, e = chp, points = points,
                                     thr_fin = thr_fin, k_r = k_r,
                                     k_l = 1)
        } else {
          r <- additive_thr_L2_trend(X, s = chp + 1, e = e,
                                     points = points, thr_fin = thr_fin,
                                     k_r = 1, k_l = max(1, k_l - 1))
        }
        cpt <- c(chp, r)
      } else {
        cpt <- chp
      }
    }
  }
  cpt <- cpt[cpt != 0]
  noiseless.cpt <- noiseless.cpt[noiseless.cpt!=0]
  return(sort(unique(c(cpt, noiseless.cpt))))}
