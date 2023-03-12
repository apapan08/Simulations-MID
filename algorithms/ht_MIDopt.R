ht_MIDopt <- function(X, thr_const = 1,
                            thr_fin = thr_const * sqrt(2*log(nrow(X))),
                            s = 1, e = nrow(X), points = 10, k_l = 1,
                            k_r = 1, sc = 3) {
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
  cpt <- ht_additive_thr_Linf(X, thr_const = thr_const,
                             thr_fin = thr_const * sqrt(2*log(nrow(X))),
                             s = 1, e = nrow(X), points = points, k_l = 1,
                             k_r = 1, scal =sc)
  
  SparsityEst <- EstimatedSparsity(apply(X, 2, IDetect::normalise,sc))
  if (SparsityEst > 0.6){
    cpt <- ht_additive_thr_L2(X, thr_const = thr_const,
                             thr_fin = thr_const * sqrt(2*log(nrow(X))),
                             s = 1, e = nrow(X), points = points, k_l = 1,
                             k_r = 1, scal = sc)
  }
  return(cpt)
}
