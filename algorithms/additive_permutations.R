Permutation_L2<-function(X,m,th,a){
  statistics <- numeric(m)
  for (s in 1:m){
    test=X[sample(1:nrow(X)),]
    ipcr <- matrix(NA, nrow(X)-1, ncol(X))
    for (i in seq_len(ncol(X))) {
      ipcr[, i] <- IDetect:::cusum_function(test[,i])/(stats::mad(diff(X[,i])/sqrt(2)) + 10^(-10))
    }
    sipcr <- sqrt(rowSums(ipcr^2)) / sqrt(ncol(X))
    statistics[s]= max(sipcr)}
  value=(length(which(statistics>th))/length(statistics))<a
  return(value)
}


Permutation_Linf<-function(X,m,th,a){
  statistics <- numeric(m)
  for (s in 1:m){
    test=X[sample(1:nrow(X)),]
    ipcl <- matrix(NA, nrow(X)-1, ncol(X))
    mpcl <- numeric()
    for (i in seq_len(ncol(X))) {
      ipcl[, i] <- IDetect:::cusum_function(test[, i])/(stats::mad(diff(X[,i])/sqrt(2))+10^(-10))
      mpcl[i] <- max(abs(ipcl[, i]))
    }
    statistics[s]=max(mpcl)}
  value=(length(which(statistics>th))/length(statistics))<a
  return(value)
}


##### define the functions
additive_thr_L2_perm <- function(X, thr_const = 1,
                                 
                                 thr_fin = thr_const * sqrt(2*log(nrow(X))),
                                 s = 1, e = nrow(X), points = 10, k_l = 1,
                                 k_r = 1,m=100,a1=0.01) {
  if (!(is.matrix(X))) {
    stop("The input in `X' should be a numeric matrix, with each data
        sequence we want to investigate being a column of this matrix.")
  }
  if ((thr_fin <= 0) || (points <= 0)) {
    stop("The threshold constant as well as the `points' argument that
        represents the magnitude of the expansion for the intervals should be
        positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps^0.5) {
    warning("The input for `points' should be a positive integer. If it is
        a positive real number then the integer part of the given number is
        used as the value of `points'.")
  }
  points <- as.integer(points)
  l <- length(X[, 1])
  r_e_points <- seq(points, l, points)
  l_e_points <- seq(l - points + 1, 1, -points)
  chp <- 0
  if (e - s <= 1) {
    cpt <- 0
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
    if (k_r < k_l) {
      while ((chp == 0) & (k_r < min(k_l, rur))) {
        lxr <- right_points[k_r] - s + 1
        ipcr <- matrix(NA, lxr - 1, ncol(X))
        for (i in seq_len(ncol(X))) {
          ipcr[, i] <- IDetect:::cusum_function(X[s:right_points[k_r], i])/(stats::mad(diff(X[,i])/sqrt(2)) + 10^(-10))
        }
        sipcr <- sqrt(rowSums(ipcr^2)) / sqrt(ncol(X))
        pos_r[k_r] <- which.max(sipcr) + s - 1
        CUSUM_r[k_r] <- sipcr[pos_r[k_r] - s + 1]
        if (Permutation_L2(X[s:right_points[k_r],],m,CUSUM_r[k_r],a=a1)){
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
        for (i in seq_len(ncol(X))) {
          ipcl[, i] <- IDetect:::cusum_function(X[left_points[k_l]:e, i])/(stats::mad(diff(X[,i])/sqrt(2)) + 10^(-10))
        }
        sipcl <- sqrt(rowSums(ipcl^2)) / sqrt(ncol(X))
        pos_l[k_l] <- which.max(sipcl) + left_points[k_l] - 1
        CUSUM_l[k_l] <- sipcl[pos_l[k_l] - left_points[k_l] + 1]
        if (Permutation_L2(X[left_points[k_l]:e,],m,CUSUM_l[k_l],a=a1)) {
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
        for (i in seq_len(ncol(X))) {
          ipcr[, i] <- IDetect:::cusum_function(X[s:right_points[k_r], i])/(stats::mad(diff(X[,i])/sqrt(2)) + 10^(-10))
        }
        sipcr <- sqrt(rowSums(ipcr^2)) / sqrt(ncol(X))
        pos_r[k_r] <- which.max(sipcr) + s - 1
        CUSUM_r[k_r] <- sipcr[pos_r[k_r] - s + 1]
        if (Permutation_L2(X[s:right_points[k_r],],m,CUSUM_r[k_r],a=a1)) {
          chp <- pos_r[k_r]
        } else {
          lxl <- e - left_points[k_l] + 1
          ipcl <- matrix(NA, lxl - 1, ncol(X))
          for (i in seq_len(ncol(X))) {
            ipcl[, i] <- IDetect:::cusum_function(X[left_points[k_l]:e, i])/(stats::mad(diff(X[,i])/sqrt(2)) + 10^(-10))
          }
          sipcl <- sqrt(rowSums(ipcl^2)) / sqrt(ncol(X))
          pos_l[k_l] <- which.max(sipcl) + left_points[k_l] - 1
          CUSUM_l[k_l] <- sipcl[pos_l[k_l] - left_points[k_l] + 1]
          if (is.na(CUSUM_l[k_l])){cat("k_r:",k_r,"k_l",k_l,"rur:",rur)}
          if (Permutation_L2(X[left_points[k_l]:e,],m,CUSUM_l[k_l],a=a1))  {
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
        r <- additive_thr_L2_perm(X, s = s, e = chp, points = points,
                                  thr_fin = thr_fin, k_r = k_r,
                                  k_l = 1,m=m,a1 = a1)
      } else {
        r <- additive_thr_L2_perm(X, s = chp + 1, e = e, points = points,
                                  thr_fin = thr_fin, k_r = 1,
                                  k_l = max(1, k_l - 1),m=m,a1 = a1)
      }
      cpt <- c(chp, r)
    } else {
      cpt <- chp
    }
  }
  cpt <- cpt[cpt != 0]
  # Remove changepoints that appear next to each other
  cpt <- cpt[!c(FALSE,diff(cpt)%in% c(1,2))]
  return(sort(cpt))
}

additive_thr_Linf_perm <- function(X, thr_const = 1,
                                   thr_fin = thr_const * sqrt(2 * log(nrow(X))),
                                   s = 1, e = nrow(X), points = 10, k_l = 1,
                                   k_r = 1,m=100,a1=0.01) {
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
    if (k_r < k_l) {
      while ((chp == 0) & (k_r < min(k_l, rur))) {
        lxr <- right_points[k_r] - s + 1
        ipcr <- matrix(NA, lxr - 1, ncol(X))
        mpcr <- numeric()
        for (i in seq_len(ncol(X))) {
          ipcr[, i] <- IDetect:::cusum_function(X[s:right_points[k_r], i])/(stats::mad(diff(X[,i])/sqrt(2))+10^(-10))
          mpcr[i] <- max(abs(ipcr[, i]))
        }
        pos.max <- which.max(mpcr)
        pos_r[k_r] <- which.max(abs(ipcr[, pos.max])) + s - 1
        CUSUM_r[k_r] <- mpcr[pos.max]
        if (Permutation_Linf(X[s:right_points[k_r],],m,CUSUM_r[k_r],a=a1)) {
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
        mpcl <- numeric()
        for (i in seq_len(ncol(X))) {
          ipcl[, i] <- IDetect:::cusum_function(X[left_points[k_l]:e, i])/(stats::mad(diff(X[,i])/sqrt(2))+10^(-10))
          mpcl[i] <- max(abs(ipcl[, i]))
        }
        pos.max <- which.max(mpcl)
        pos_l[k_l] <- which.max(abs(ipcl[, pos.max])) + left_points[k_l] - 1
        CUSUM_l[k_l] <- mpcl[pos.max]
        if (Permutation_Linf(X[left_points[k_l]:e,],m,CUSUM_l[k_l],a=a1)) {
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
        mpcr <- numeric()
        for (i in seq_len(ncol(X))) {
          ipcr[, i] <- IDetect:::cusum_function(X[s:right_points[k_r], i])/(stats::mad(diff(X[,i])/sqrt(2))+10^(-10))
          mpcr[i] <- max(abs(ipcr[, i]))
        }
        pos.max <- which.max(mpcr)
        pos_r[k_r] <- which.max(abs(ipcr[, pos.max])) + s - 1
        CUSUM_r[k_r] <- mpcr[pos.max]
        if (Permutation_Linf(X[s:right_points[k_r],],m,CUSUM_r[k_r],a=a1)) {
          chp <- pos_r[k_r]}else{
          lxl <- e - left_points[k_l] + 1
          ipcl <- matrix(NA, lxl - 1, ncol(X))
          mpcl <- numeric()
          for (i in seq_len(ncol(X))) {
            ipcl[, i] <- IDetect:::cusum_function(X[left_points[k_l]:e, i])/(stats::mad(diff(X[,i])/sqrt(2)) + 10^(-10))
            mpcl[i] <- max(abs(ipcl[, i]))
          }
          pos.max <- which.max(mpcl)
          pos_l[k_l] <- which.max(abs(ipcl[, pos.max])) + left_points[k_l] - 1
          CUSUM_l[k_l] <- mpcl[pos.max]
          if (Permutation_Linf(X[left_points[k_l]:e,],m,CUSUM_l[k_l],a=a1)){
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
        r <- additive_thr_Linf_perm(X, s = s, e = chp, points = points,
                                    thr_fin = thr_fin, k_r = k_r,
                                    k_l = 1,m=m,a1 = a1)
      } else {
        r <- additive_thr_Linf_perm(X, s = chp + 1, e = e,
                                    points = points, thr_fin = thr_fin,
                                    k_r = 1, k_l = max(1, k_l - 1),a1 = a1)
      }
      cpt <- c(chp, r)
    } else {
      cpt <- chp
    }
  }
  cpt <- cpt[cpt != 0]
  cpt <- cpt[!c(FALSE,diff(cpt)%in% c(1,2))]
  return(sort(cpt))
}
