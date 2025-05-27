###   Generic (weighted) quantile estimator
##    https://aakinshin.net/posts/weighted-quantiles/

qe.generic <- function(x, probs, cdf.gen, weights = NULL, na.rm) {
  n <- length(x)

  ##  Changes/additions based on matrixStats::weightedMad
  if (!is.null(weights)) {
    if (length(na.omit(weights)) != na.omit(n)) {
      stop(sprintf("The number of non-NA elements in arguments '%s' and '%s' does not match: %.0f != %.0f", "weights", "x", length(na.omit(weights)), na.omit(n)))
    }
  } else weights <- rep(1 / n, n)

  # Argument 'na.rm':
  na.rm <- as.logical(na.rm)
  if (is.na(na.rm)) na.rm <- FALSE

  na_value <- NA
  storage.mode(na_value) <- storage.mode(x)

  # Remove values with zero (and negative) weight. This will:
  #  1) take care of the case when all weights are zero,
  #  2) it will most likely speed up the sorting.
  tmp <- (is.na(weights) | weights > 0)
  if (!all(tmp)) {
    x <- .subset(x, tmp)
    weights <- .subset(weights, tmp)
  }
  tmp <- NULL  # Not needed anymore

  # Drop missing values?
  if (na.rm) {
    keep <- which(!is.na(x))
    x <- .subset(x, keep)
    weights <- .subset(weights, keep)
    keep <- NULL  # Not needed anymore
  } else if (anyNA(x)) {
    return(na_value)
  }

  # Are any weights Inf? Then treat them with equal weight and all others
  # with weight zero.
  tmp <- is.infinite(weights)
  if (any(tmp)) {
    keep <- tmp
    x <- .subset(x, keep)
    weights <- rep(1, times = length(x))
    keep <- NULL  # Not needed anymore
  }
  tmp <- NULL  # Not needed anymore

  # Are there any values left to calculate the weighted median of?
  # This is consistent with how stats::mad() works.
  if (length(x) == 0L) {
    return(na_value)
  } else if (length(x) == 1L) {
    zero_value <- 0
    storage.mode(zero_value) <- storage.mode(x)
    return(zero_value)
  }

  nw <- sum(weights)^2 / sum(weights^2) # Kish's effective sample size

  indexes <- order(x)
  x <- x[indexes]
  weights <- weights[indexes]

  weights <- weights / sum(weights)
  cdf.probs <- cumsum(c(0, weights))

  sapply(probs, function(p) {
    cdf <- cdf.gen(nw, p)
    q <- cdf(cdf.probs)
    w <- tail(q, -1) - head(q, -1)
    sum(w * x)
  })
}

# Weighted Harrell-Davis quantile estimator
hdQuantile <- function(x, probs, weights = NULL, na.rm=TRUE) {
  cdf.gen <- function(n, p) return(function(cdf.probs) {
    pbeta(cdf.probs, (n + 1) * p, (n + 1) * (1 - p))
  })
  qe.generic(x, probs, cdf.gen, weights, na.rm=na.rm)
}

# Weighted Type 7 quantile estimator
weighted.quantile <- function(x, probs, weights = NULL, na.rm=TRUE) {
  cdf.gen <- function(n, p) return(function(cdf.probs) {
    h <- p * (n - 1) + 1
    u <- pmax((h - 1) / n, pmin(h / n, cdf.probs))
    u * n - h + 1
  })
  qe.generic(x, probs, cdf.gen, weights, na.rm=na.rm)
}


###   Akinshin's Gamma Effect Size
##    https://aakinshin.net/posts/nonparametric-effect-size/

`gammaEffectSize` <- function(x, y, probs, weights.x = NULL, weights.y = NULL, small_n_correction = TRUE, na.rm=TRUE) {
  ###  TMRW - add in checks/changes for na.rm/weights/Kish
  (hdQuantile(y, probs = probs, weights = weights.y, na.rm=na.rm) -
   hdQuantile(x, probs = probs, weights = weights.x, na.rm=na.rm)) /
  pooled.hdMAD(x, y, weights.x, weights.y, snc = small_n_correction, rm.na = na.rm)
}


###   Utility Functions

`hdMAD` <- function(x, small_n_correction = TRUE, weights = NULL, na.rm=TRUE) {

  if (!is.logical(small_n_correction)) {
    CC <- small_n_correction
  } else CC <- 1.4826
  if (small_n_correction) {
    N <- length(x)

    if (!is.null(weights)) {
      N <- sum(weights, na.rm=na.rm)^2 / sum(weights^2, na.rm=na.rm) # Kish's effective sample size
    }

    if (N < 20000) CC <- consistency_constant(N)
  }
  CC * hdQuantile(abs(x - hdQuantile(x, probs = 0.5, weights=weights, na.rm=na.rm)), probs = 0.5, weights=weights, na.rm=na.rm)
}

pooled.hdMAD <- function(x, y, wx, wy, snc, rm.na) {
  nx <- length(x)
  ny <- length(y)

  if (!is.null(wx)) {
    if (length(na.omit(wx)) != na.omit(nx)) {
      stop(sprintf("The number of non-NA elements in arguments '%s' and '%s' does not match: %.0f != %.0f", "wx", "x", length(na.omit(wx)), na.omit(nx)))
    }
    nx <- sum(wx)^2 / sum(wx^2) # Kish's effective sample size
  }
  if (!is.null(wy)) {
    if (length(na.omit(wy)) != na.omit(ny)) {
      stop(sprintf("The number of non-NA elements in arguments '%s' and '%s' does not match: %.0f != %.0f", "wy", "y", length(na.omit(wy)), na.omit(ny)))
    }
    ny <- sum(wy)^2 / sum(wy^2) # Kish's effective sample size
  }

  sqrt(((nx - 1) * hdMAD(x, snc, wx, rm.na) ^ 2 + (ny - 1) * hdMAD(y, snc, wy, rm.na) ^ 2) / (nx + ny - 2))
}


###   Lookup Tables for MAD consistency constants with small(ish) N size
###   Taken from https://aakinshin.net/posts/unbiased-mad-hd/

small.n.lookup <- data.table::fread(
"N,A_n,C_n,N,A_n,C_n
1,NA,NA,51,0.66649,1.50039
2,0.56417,1.77250,52,0.66668,1.49998
3,0.63769,1.56816,53,0.66682,1.49966
4,0.62661,1.59589,54,0.66699,1.49926
5,0.63853,1.56611,55,0.66713,1.49895
6,0.63834,1.56656,56,0.66728,1.49863
7,0.63915,1.56458,57,0.66741,1.49833
8,0.64141,1.55908,58,0.66753,1.49805
9,0.64237,1.55675,59,0.66767,1.49774
10,0.64397,1.55288,60,0.66780,1.49746
11,0.64535,1.54955,61,0.66791,1.49720
12,0.64662,1.54651,62,0.66803,1.49694
13,0.64790,1.54346,63,0.66815,1.49667
14,0.64908,1.54064,64,0.66825,1.49644
15,0.65018,1.53803,65,0.66836,1.49621
16,0.65125,1.53552,66,0.66846,1.49597
17,0.65226,1.53313,67,0.66857,1.49574
18,0.65317,1.53101,68,0.66865,1.49555
19,0.65404,1.52896,69,0.66876,1.49531
20,0.65489,1.52698,70,0.66883,1.49514
21,0.65565,1.52520,71,0.66893,1.49493
22,0.65638,1.52351,72,0.66901,1.49475
23,0.65708,1.52190,73,0.66910,1.49456
24,0.65771,1.52043,74,0.66918,1.49437
25,0.65832,1.51902,75,0.66925,1.49422
26,0.65888,1.51772,76,0.66933,1.49402
27,0.65943,1.51647,77,0.66940,1.49387
28,0.65991,1.51536,78,0.66948,1.49370
29,0.66036,1.51433,79,0.66955,1.49354
30,0.66082,1.51328,80,0.66962,1.49339
31,0.66123,1.51233,81,0.66968,1.49325
32,0.66161,1.51146,82,0.66974,1.49312
33,0.66200,1.51057,83,0.66980,1.49298
34,0.66235,1.50977,84,0.66988,1.49281
35,0.66270,1.50899,85,0.66993,1.49270
36,0.66302,1.50824,86,0.66999,1.49257
37,0.66334,1.50753,87,0.67005,1.49244
38,0.66362,1.50688,88,0.67009,1.49233
39,0.66391,1.50623,89,0.67016,1.49219
40,0.66417,1.50563,90,0.67021,1.49207
41,0.66443,1.50504,91,0.67026,1.49196
42,0.66469,1.50447,92,0.67031,1.49185
43,0.66493,1.50393,93,0.67036,1.49174
44,0.66515,1.50341,94,0.67041,1.49161
45,0.66539,1.50289,95,0.67046,1.49152
46,0.66557,1.50246,96,0.67049,1.49144
47,0.66578,1.50200,97,0.67055,1.49131
48,0.66598,1.50155,98,0.67060,1.49121
49,0.66616,1.50115,99,0.67063,1.49114
50,0.66633,1.50076,100,0.67068,1.49102")

big.n.lookup <- data.table::fread(
"N,A_n
100,0.6707
110,0.6711
120,0.6714
130,0.6716
140,0.6719
150,0.6720
200,0.6727
250,0.6731
300,0.6733
350,0.6735
400,0.6736
450,0.6737
500,0.6738
1000,0.6742
1500,0.67425
2000,0.6743
3000,0.67435
5000,0.6744
10000,0.67445
19999,0.67449") # AVI added 3000 - 19999 to make a more gradual transition to C_n = 1.4826

big.n.lookup[, C_n := 1/A_n]

n.lookup <- rbindlist(list(small.n.lookup[, 1:3], small.n.lookup[, 4:6], big.n.lookup[-1,]))

###   Create a splinefun for interpolating values between those listed in the tables.
consistency_constant <- splinefun(n.lookup[-c(1,3:5),N], n.lookup[-c(1,3:5),C_n], method="natural")
