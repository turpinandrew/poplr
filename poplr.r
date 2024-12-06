# An efficient implementation of PoPLR
#
#> *Visual Field Progression in Glaucoma: Estimating the Overall Significance of Deterioration with Permutation Analyses of Pointwise Linear Regression (PoPLR)*
#> Neil O'Leary; Balwantray C. Chauhan; Paul H. Artes
#> *Investigative Ophthalmology & Visual Science*
#> October 2012, Vol.53, 6776-6784. https://doi.org/10.1167/iovs.12-10049


# Fri Dec  2 13:16:05 AWST 2022
#
# @param series A matrix/data.frame with one row per location in the VF
#               and one column per visit (ordered in time so that col1 is earliest
#               and the final column the most recent visit)
#               Values in the column are dB values (either raw or Total Deviation).
# @param threshold Only include PLR pvalues less than or equal to this number in S
# @param perm_count Number of permutations of series to get p-value for S
#
# @return p-value for statistic S which tracks the false positive rate of calling this series progressing. 
#
PoPLR_basic <- function(series, threshold = 0.05, perm_count = 5000) {
        # get the S statistic for a single series
    getS <- function(series, threshold = 0.05) {
        xs <- seq_len(ncol(series))
        pvals <- apply(series, 1, function(rr) {
            m <- lm(rr ~ xs)
            pt(coef(summary(m))[2, 3], m$df, lower = TRUE)  # one sided
        })
        #print(xs)
        #print(series)
        #print(pvals)
        z <- pvals <= max(threshold, min(pvals))
        sum(-log(pvals[z]))
    }

        # get the actual S for series
    S <- getS(series, threshold)

        # get S values for perm_count permutations of series
    S_perms <- sapply(1:perm_count, function(i) getS(series[, order(runif(ncol(series)))], threshold))

    #print(S)
    #print(table(round(S_perms)))
    return(1 - sum(S_perms <= S) / perm_count)
}

# test
#dbs <- round(runif(52, 0, 41))
#series <- data.frame(dbs, dbs-1, dbs-2, dbs - 3)
#colnames(series) <- NULL
#p <- PoPLR(series, perm_count = 500)
#p

# @param n Return permutations of 1:n
# @param perm_count Number of permutations to return
# @param warnings TRUE for warnings, FALSE for not.
#
# @return list of permutations of 1:n making sure that the first in the list is 1:n
#         (Note no guarantee on length of list)
#
# To find unique permutations
#     perms <- 1:n + perm_count random permutations
#     hashes[p] is list of positions of perm p in perms
#     while any hashes[p] is longer than 1
#        for pp in tail(hashes[p], -1):
#           swap two elements in pp to get p'
#           hashes[pp]--
#           hashes[p']++
#
get_perms <- function(n, perm_count = 5000, warnings = TRUE) {
        # First build a list of permutations to use
        # perms is a list of perm_count permutations of 1:n with 1,2,...n as element [[1]]
    if (perm_count >= factorial(n)) {
        if (warnings)
            warning(paste("perm_count = ", perm_count, "is greater than n!", factorial(n), "in PoPLR. Restricting to n! - 1"))
        perms <- combinat::permn(1:n)
    } else {
        if (n <= 9) {
            perms <- combinat::permn(1:n)
            perms <- perms[c(1, sample(2:length(perms), perm_count))]   # make sure to keep 1
        } else {
            perms <- c(list(1:n), lapply(1:perm_count, function(i) order(runif(n))))
            hashes <- vector("list")
            for (i_p in seq_along(perms)) {
                sig <- paste(perms[[i_p]], collapse = " ")
                hashes[[sig]] <- c(hashes[[sig]], i_p)
            }
            bad_perms <- which(unlist(lapply(hashes, length)) > 1)
            count <- 0
            while (length(bad_perms) > 0 && count < 100) {
                to_go <- unlist(lapply(hashes[bad_perms], tail, -1))  # keep the first one (it might be 1:n)
                hashes[bad_perms] <- lapply(hashes[bad_perms], head, 1)
                for (i in to_go) {
                    j <- round(runif(1, 1, n))  # swap two random elements
                    k <- round(runif(1, 1, n))
                    temp <- perms[[i]][j]
                    perms[[i]][j] <- perms[[i]][k]
                    perms[[i]][k] <- temp
                    sig <- paste(perms[[i]], collapse = " ")
                    hashes[[sig]] <- c(hashes[[sig]], i)
                }
                bad_perms <- which(unlist(lapply(hashes, length)) > 1)

                count <- count + 1
            }
            if (count == 100 && warnings)
                warning(paste("PoPLR did not find", perm_count, "unique permutations. Used", length(names(hashes))))
        }
        stopifnot(all(paste(perms[[1]], collapse = " ") == paste(1:n, collapse = " ")))
    }

    return(perms)
}

# Computation based on https://www.geo.fu-berlin.de/en/v/soga/Basics-of-statistics/Hypothesis-Tests/Inferential-Methods-in-Regression-and-Correlation/Inferences-About-the-Slope/index.html#:~:text=The%20regression%20t%2DTest%20is,(linear)%20predictor%20of%20y.
#
# ASSUMES: x will be integers or at least separated by 1/3, 1/2, etc. That is, SumXdSq > 1 ish.
# Idle thought: Could this be made faster by shuffling the x values (known to be 1:n) rather than the y values?
#
# Allows for NA values at the start of rows in the series.
#
# @param series A matrix/data.frame with one row per location in the VF
#               and one column per visit (ordered in time so that col1 is earliest
#               and the final column the most recent visit)
#               Values in the column are dB values (either raw or Total Deviation).
# @param threshold Only include PLR pvalues less than or equal to this number in S
# @param perm_count Number of permutations of series to get p-value for S
# @param warnings TRUE for warnings from PoPLR, FALSE for not.
# @param arrest EXPERIMENTAL - do not use!
#
# @return Probability that series is not progressing (kind of).
#
PoPLR <- function(series, threshold = 0.05, perm_count = 5000, warnings = TRUE, arrest = FALSE) {
        # throw out locations with 2 or less visits
    z <- apply(series, 1, function(rr) sum(!is.na(rr)) > 2)
    series <- series[z, ]

        # throw out any columns that are all NA
    z <- apply(series, 2, function(cc) any(!is.na(cc)))
    series <- series[, z]

    stopifnot(nrow(series) > 1)
    stopifnot(ncol(series) > 2)

        # get permutation vectors for largest n
    n <- ncol(series)
    perms <- get_perms(n, perm_count, warnings)
    perm_count <- length(perms) - 1  # does not include first 1:n

        # p_loc_perm[loc number, perm_number] == p value
    p_loc_perm <- matrix(NA, nrow = nrow(series), ncol = perm_count + 1)

        # Now compute all the p-values for each loc in each appropriate permutation
        # Only include perms for NA rows where no intermittent NAs
    for (loc in seq_len(nrow(series))) {
        series_loc <- series[loc, ]
        ybar <- mean(series_loc, na.rm = TRUE)

        n.loc <- sum(!is.na(series_loc))
        xs <- seq_len(n.loc)
        xbar <- mean(xs)
        xd <- xs - xbar
        sumXdSq <- sum(xd * xd)
        sqrtSumXdSq <- sqrt(sumXdSq)

        for (perm in 0:perm_count + 1) {
            ys <- series_loc[perms[[perm]]]   # this is length n >= n.loc

                # remove leading and trailing NAs and if the result
                # is not length(n.loc) skip this perm
            ii <- which(!is.na(ys))
            a <- head(ii, 1)
            b <- tail(ii, 1)
            if (b - a + 1 == n.loc) {
                ys <- ys[a:b]

                if (arrest && -1 %in% ys) {
                    ii <- n.loc
                    while (ii >= 1 && ys[[ii]] == -1)
                        ii <- ii - 1
                    if (ii > 1 && all(ys[1:ii] != -1))
                        p_loc_perm[loc, perm] <- 0.5^(2 + ii)   # {!Red}+ -> {Red}+ could exist if tt=0 with this prob
                    else
                        p_loc_perm[loc, perm] <- 1   # {.}*Red{.}* -> {Red}+
                } else {
                    yd <- ys - ybar

                    beta <- sum(xd * yd) / sumXdSq
                    alpha <- ybar - beta * xbar
                    res <- ys - alpha - beta * xs
                    sse <- sum((ys - alpha - beta * xs)^2)
                    if (sse < 1e-10)  { # often == 0 for VF analysis of a small number of points. Also assumes sqrtSumXdSq is not super small.
                        if (beta >= 0) {
                            p_loc_perm[loc, perm] <- 1   # t -> Inf  (assuming 0/0 -> Inf)
                        } else {
                            p_loc_perm[loc, perm] <- 0   # t -> -Inf
                        }
                    } else {
                        se <- sqrt(sse / (n.loc - 2)) / sqrtSumXdSq
                        t <- beta / se
                        p_loc_perm[loc, perm] <- pt(t, n.loc - 2)
            #cat(res)
            #cat(sprintf("\nloc=%2s perm=%2s n=%2.0f beta=%+4.2f sse=%5.2f se=%4.2f t=%+5.2f p=%4.2f\n",
            #      loc, perm, n.loc, beta, sse, se, t, pt(t, n.loc - 2)))
                    }
                }
            }
        }
    }

    S_values <- apply(p_loc_perm, 2, function(ps) {
        z <- ps <= max(min(ps, na.rm = TRUE), threshold)
        sum(-log(ps[z]), na.rm = TRUE)
    })

    #print("***************")
    #print(p_loc_perm)
    #print(S_values[1:7])
    #print(S_values[1])
    #print(table(round(S_values)))
    others <- tail(S_values, -1)
    pp <-  sum(others < S_values[[1]]) + sum(others == S_values[[1]])/2  # another assumption here for small discrete distributions.
    return(list(p = 1 - pp / perm_count, ss = S_values))
}

# Keep xs the same even if NAs
# This definitely alters false-pos rate as S_obs always has NAs at the start.
PoPLR3 <- function(series, threshold = 1, perm_count = 5000, warnings = TRUE, arrest = FALSE) {
        # throw out locations with 2 or less visits
    z <- apply(series, 1, function(rr) sum(!is.na(rr)) > 2)
    series <- series[z, ]

        # throw out any columns that are all NA
    z <- apply(series, 2, function(cc) any(!is.na(cc)))
    series <- series[, z]

    stopifnot(nrow(series) > 1)
    stopifnot(ncol(series) > 2)

        # get permutation vectors for largest n
    n <- ncol(series)
    perms <- get_perms(n, perm_count, warnings)
    perm_count <- length(perms) - 1  # does not include first 1:n

        # p_loc_perm[loc number, perm_number] == p value
    p_loc_perm <- matrix(NA, nrow = nrow(series), ncol = perm_count + 1)

        # setup xs
    n.loc <- n
    xs <- seq_len(n.loc)
    xbar <- mean(xs)
    xd <- xs - xbar
    sumXdSq <- sum(xd * xd)
    sqrtSumXdSq <- sqrt(sumXdSq)

    # Now compute all the p-values for each loc in each appropriate permutation
        # Only include perms for NA rows where no intermittent NAs
    for (loc in seq_len(nrow(series))) {
        series_loc <- series[loc, ]
        ybar <- mean(series_loc, na.rm = TRUE)

        for (perm in 0:perm_count + 1) {
            ys <- series_loc[perms[[perm]]]   # this is length n >= n.loc

                # remove leading and trailing NAs and if the result
                # is not length(n.loc) skip this perm
            ii <- which(!is.na(ys))
            a <- head(ii, 1)
            b <- tail(ii, 1)
            if (b - a + 1 == n.loc) {
                ys <- ys[a:b]

                if (arrest && -1 %in% ys) {
                    ii <- n.loc
                    while (ii >= 1 && ys[[ii]] == -1)
                        ii <- ii - 1
                    if (ii > 1 && all(ys[1:ii] != -1))
                        p_loc_perm[loc, perm] <- 0.5^(2 + ii)   # {!Red}+ -> {Red}+ could exist if tt=0 with this prob
                    else
                        p_loc_perm[loc, perm] <- 1   # {.}*Red{.}* -> {Red}+
                } else {
                    yd <- ys - ybar

                    beta <- sum(xd * yd, na.rm = TRUE) / sumXdSq
                    alpha <- ybar - beta * xbar
                    res <- ys - alpha - beta * xs
                    sse <- sum((ys - alpha - beta * xs)^2, na.rm = TRUE)
                    if (sse < 1e-10)  { # often == 0 for VF analysis of a small number of points. Also assumes sqrtSumXdSq is not super small.
                        if (beta >= 0) {
                            p_loc_perm[loc, perm] <- 1   # t -> Inf  (assuming 0/0 -> Inf)
                        } else {
                            p_loc_perm[loc, perm] <- 0   # t -> -Inf
                        }
                    } else {
                        se <- sqrt(sse / (n.loc - 2)) / sqrtSumXdSq
                        t <- beta / se
                        p_loc_perm[loc, perm] <- pt(t, sum(!is.na(xs)) - 2)
            #cat(res)
            #cat(sprintf("\nloc=%2s perm=%2s n=%2.0f beta=%+4.2f sse=%5.2f se=%4.2f t=%+5.2f p=%4.2f\n",
            #      loc, perm, n.loc, beta, sse, se, t, pt(t, n.loc - 2)))
                    }
                }
            }
        }
    }

    S_values <- apply(p_loc_perm, 2, function(ps) {
        z <- ps <= max(min(ps, na.rm = TRUE), threshold)
        sum(-log(ps[z]), na.rm = TRUE)
    })

    #print("***************")
    #print(p_loc_perm)
    #print(S_values[1:7])
    #print(S_values[1])
    #print(table(round(S_values)))
    others <- tail(S_values, -1)
    pp <-  sum(others < S_values[[1]]) + sum(others == S_values[[1]])/2  # another assumption here for small discrete distributions.
    return(list(p = 1 - pp / perm_count, ss = S_values))
}



#  NA, NA, ...
#  NA, NA, ...
#  NA, NA, ...
#  NA, NA, NA, ...
#  NA, NA, NA, ...
#  NA, NA, NA, ...
#  NA, NA, NA, NA, NA, ...
#  NA, NA, NA, NA, NA, ...
#  8    8   5   2   2
#
# Keep NA's at the start for each location
#' @param series A matrix/data.frame with one row per location in the VF, visit cols
#' @return List of indexes of last non-NA in a block of visits with same NA count (possibly NULL)
get_boundaries <- function(series) {
#print(series)
#print(series[seq(51, nrow(series)), , drop = FALSE])
        # num NAs for each visit
    num_NAs <- apply(series, 2, function(rr) sum(is.na(rr)))

    if (all(num_NAs == 0))
        return(NULL)

    curr_num <- num_NAs[[1]]
    block_boundaries <- c(1)
    for (i in seq(2, ncol(series))) {
#print(paste("visit=", i))
#print(block_boundaries)
        if (num_NAs[[i]] == 0)
            break

            # if last block_b column is all NA's, increment it, else add new b_b
        n <- length(block_boundaries)
        if (num_NAs[[i]] == curr_num) {
            block_boundaries[[n]] <- i
        } else {
            block_boundaries <- c(block_boundaries, i)
            curr_num <- num_NAs[[i]]
        }
    }
    return(block_boundaries)
}

# Assumes no rows with less than 2 non-NAs
# Only permutes with NA-block boundaries
get_perms4 <- function(series, perm_count, warnings = TRUE) {
    bb <- get_boundaries(series)
#print(bb)
    if (is.null(bb))
        return(get_perms(ncol(series), perm_count, warnings))

        # @param start Index of visit to start perm - 1
        # @param end   Index of visit to end perm
        # @return Perms of (start, end] = 1:(end - start) + start. Resample to get exactly perm_count perms.
        #                      Keep 1 as 1..n
    gp <- function(start, end) {
#print("------gp")
#print(paste("start=", start, "end=", end, "perm_count=", perm_count))
        p <- get_perms(end - start, perm_count - 1, warnings = FALSE)
#print("-p")
#print(p)
        l <- lapply(p, function(p) start + p)
#print("-l")
#print(l)
            # keep 1...n as perm 1
        if (length(l) < perm_count)
            l <- c(list(l[[1]]), sample(l, perm_count - 1, replace = TRUE))
#print("-l2")
#print(l)
        l
    }

    perm_blocks <- lapply(seq_along(bb), function(i_bb) {
        if (i_bb == 1)
            gp(0, bb[[1]])
        else
            gp(bb[[i_bb - 1]], bb[[i_bb]])
    })

    if (tail(bb, 1) != ncol(series))
        perm_blocks <- c(perm_blocks, list(gp(tail(bb, 1), ncol(series))))
#print(perm_blocks)
#print(unlist(lapply(perm_blocks, length)))

    # combine all perm_blocks
    res <- perm_blocks[[1]]
    for (i_pb in seq(2, length(perm_blocks)))
        res <- lapply(seq_len(perm_count), function(i) c(res[[i]], perm_blocks[[i_pb]][[i]]))

#print(res)
    return(res)
}

#' Only use permutations where NAs are at the start
#' Only use locations with at least 6 visits.
#' If arrest = TRUE, exclude any locations with <= 16 at penultimate visit.
#' ASSUMES NAs are at the start.
#
#' @param series A matrix/data.frame with one row per location in the VF
#'               and one column per visit (ordered in time so that col1 is earliest
#'               and the final column the most recent visit)
#'               Values in the column are dB values (either raw or Total Deviation).
#' @param threshold Only include PLR p-values less than or equal to this number in S
#' @param perm_count Number of permutations of series to get p-value for S
#' @param warnings TRUE for warnings from PoPLR, FALSE for not.
#' @param arrest If TRUE, discard yellow-yellow+-Orange*-Red* locations and use zest value for Green+-Yellow locations.
#
#' @return list containing
#'   * p is p-value for statistic S which tracks the false positive rate of calling this series progressing.
#'   * ss is a vector of S values for each permutation.
PoPLR4 <- function(series, threshold = 1, perm_count = 5000, warnings = TRUE, arrest = FALSE) {
        # throw out locations with 6 or less visits
    z <- apply(series, 1, function(rr) sum(!is.na(rr)) >= 6)
    series <- series[z, ]

    if (arrest) {
            # Only keep locations with 1 or fewer Yellow/Orange/Red
        z <- apply(series, 1, function(rr) sum(rr <= 16, na.rm = TRUE) <= 1)
        series <- series[z, ]
    }

        # throw out any columns that are all NA
    z <- apply(series, 2, function(cc) any(!is.na(cc)))
    series <- series[, z]

    if (nrow(series) <= 1) return(list(ss = NA, p = NA))
    if (ncol(series) <= 6) return(list(ss = NA, p = NA))

        # get permutation vectors for largest n
    n <- ncol(series)
    perms <- get_perms4(series, perm_count, warnings)
    perm_count <- length(perms) - 1  # does not include first 1:n

        # p_loc_perm[loc number, perm_number] == p value
    p_loc_perm <- matrix(NA, nrow = nrow(series), ncol = perm_count + 1)

        # Now compute all the p-values for each loc in each appropriate permutation
        # Only include perms for NA rows where no intermittent NAs
    for (loc in seq_len(nrow(series))) {
        series_loc <- series[loc, ]
        ybar <- mean(series_loc, na.rm = TRUE)

        n.loc <- sum(!is.na(series_loc))
        xs <- seq_len(n.loc)
        xbar <- mean(xs)
        xd <- xs - xbar
        sumXdSq <- sum(xd * xd)
        sqrtSumXdSq <- sqrt(sumXdSq)

        for (perm in 0:perm_count + 1) {
            ys <- series_loc[perms[[perm]]]   # this is length n >= n.loc

                # remove leading and trailing NAs and if the result
                # is not length(n.loc) skip this perm
            ii <- which(!is.na(ys))
            a <- head(ii, 1)
            b <- tail(ii, 1)
            if (b - a + 1 == n.loc) {
                ys <- ys[a:b]

                yd <- ys - ybar

                beta <- sum(xd * yd) / sumXdSq
                alpha <- ybar - beta * xbar
if (perm == 1)
print(paste("beta", beta, "alpha", alpha))
                res <- ys - alpha - beta * xs
                sse <- sum((ys - alpha - beta * xs)^2)
                if (sse < 1e-10)  { # often == 0 for VF analysis of a small number of points. Also assumes sqrtSumXdSq is not super small.
                    if (beta >= 0) {
                        p_loc_perm[loc, perm] <- 1   # t -> Inf  (assuming 0/0 -> Inf)
                    } else {
                        p_loc_perm[loc, perm] <- 0   # t -> -Inf
                    }
                } else {
                    se <- sqrt(sse / (n.loc - 2)) / sqrtSumXdSq
                    t <- beta / se
                    p_loc_perm[loc, perm] <- pt(t, n.loc - 2)
        #cat(res)
        #cat(sprintf("\nloc=%2s perm=%2s n=%2.0f beta=%+4.2f sse=%5.2f se=%4.2f t=%+5.2f p=%4.2f\n",
        #      loc, perm, n.loc, beta, sse, se, t, pt(t, n.loc - 2)))
                }
            }
        }
    }

    S_values <- apply(p_loc_perm, 2, function(ps) {
        z <- ps <= max(min(ps, na.rm = TRUE), threshold)
        sum(-log(ps[z]), na.rm = TRUE)
    })

    #print("***************")
    #print(p_loc_perm)
    #print(S_values[1:7])
    #print(S_values[1])
    #print(table(round(S_values)))
    others <- tail(S_values, -1)
    pp <-  sum(others < S_values[[1]]) + sum(others == S_values[[1]])/2  # another assumption here for small discrete distributions.
    return(list(p = 1 - pp / perm_count, ss = S_values))
}



#############
# test
#############
#dbs <- round(runif(52, 0, 41))
#series <- matrix(c(dbs, dbs - 1, dbs - 2, dbs - 3), ncol = 4)
#p <- PoPLR4(series, perm_count = 5000)
#print(p)

#for (delta in  -1 / c(0.5, 1, 2, 4, 8))
#for (prog in c(52, 25, 0)) {
#    cat(sprintf("\n delta = %6.4f num prog = %2s ", delta, prog))
#    ds <- c(rep(0, length(dbs) - prog), rep(delta, prog))
#    series <- matrix(c(dbs, dbs + ds, dbs + 2 * ds, dbs + 3 * ds), ncol = 4)
#    p <- PoPLR(series)
#    cat(sprintf("p = %s\n", p))
#}

#s <- matrix(c(
# 35.924587, 35.924587, 37.441342, 38.6314074, 38.631407, 35.445228, 38.631407,
# 11.130064,  7.920757,  7.517558,  5.2927513,  6.246617,  6.246617, 12.452755,
#  4.138881,  8.241848, 12.860583,  0.2980607,  5.292264,  6.904202,  3.293290,
# 11.871311, 12.252357, 14.636894, 10.8904141,  9.036826,  9.845077,  3.698621,
#  6.419012, 12.274877,  7.865783, 16.8246304, 12.279927,  7.865783, 10.241568,
# 11.516643, 15.330036,  4.506971, 10.2755194,  7.740855,  3.297869,  8.947669,
# 37.034851, 38.591479, 38.591479, 36.2084900, 36.505164, 35.581025, 35.417021,
# 38.508213, 38.508213, 36.375763, 34.5469293, 38.508213, 35.362040, 34.546929,
# 31.545246, 33.605193, 32.355546, 32.1473103, 30.987397, 31.545246, 32.737632,
#  6.274475, 13.429302, 11.175638,  0.2980425,  4.293997, 10.534879,  4.293997
#), ncol = 7, byrow = TRUE)
#print(paste("PoPLR_basic: ", PoPLR_basic(s)))
#print(paste("PoPLR: ", PoPLR(s)))

#s <- matrix(c(
#        NA,        NA, 37.441342, 38.6314074, 38.631407, 35.445228, 38.631407,
# 11.130064,  7.920757,  7.517558,  5.2927513,  6.246617,  6.246617, 12.452755,
#  4.138881,  8.241848, 12.860583,  0.2980607,  5.292264,  6.904202,  3.293290,
# 11.871311, 12.252357, 14.636894, 10.8904141,  9.036826,  9.845077,  3.698621,
#  6.419012, 12.274877,  7.865783, 16.8246304, 12.279927,  7.865783, 10.241568,
# 11.516643, 15.330036,  4.506971, 10.2755194,  7.740855,  3.297869,  8.947669,
# 37.034851, 38.591479, 38.591479, 36.2084900, 36.505164, 35.581025, 35.417021,
#        NA,        NA,        NA, 34.5469293, 38.508213, 35.362040, 34.546929,
#        NA,        NA,        NA,         NA,        NA,        NA, 32.737632,
#  6.274475, 13.429302, 11.175638,  0.2980425,  4.293997, 10.534879,  4.293997
#), ncol = 7, byrow = TRUE)
#for (pcount in c(2, 20, 200, 2000, 5040))
#    print(paste("PoPLR: ", PoPLR(s, perm_count = pcount)))

#series <- matrix(20, ncol = 10, nrow = 52)
#PoPLR4(series)
#series <- matrix(20, ncol = 10, nrow = 53)
#series[53, 1] <- NA
#print(PoPLR4(series, perm_count = 500))
