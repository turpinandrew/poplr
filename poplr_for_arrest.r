#
# Apply PoPLR to a series 1..n VFs measured with ARREST.
# For locations that have green at visit n, use PoPLR as normal.
# Any Red/Yellow locations in visit n do not use PoPLR but event triggers:
#     - If n is red, n > t, and n - t was not Red, p = low
#     - If n is Yellow, n > t, and n - t was Green, p = low
#   t >= 2; higher more specific (ie more baselines in effect)
#
# Fri 30 Dec 2022 07:52:09 AEDT
# Andrew Turpin

# @param series A matrix/data.frame with one row per location in the VF
#               and one column per visit (ordered in time so that col1 is earliest
#               and the final column the most recent visit)
#               Values in the column are dB values (either raw or Total Deviation).
# @param threshold Only include PLR pvalues less than or equal to this number in S
# @param perm_count Number of permutations of series to get p-value for S
# @param warnings TRUE for warnings from PoPLR, FALSE for not.
# @param numBaselines Number of visits to have before checking for event-based triggers (t)
# @param eventP p-value to return is a location has an event trigger
#
# @return p-value = prob of having this data if the series has not truly progressed at visit n
#
source('arrest6.r')
source('poplr.r')
PoPLR_ARREST <- function(series, threshold = 0.05, perm_count = 5000, warnings = TRUE,
                         numBaselines = 2) {
    n <- ncol(series)

        # Green -> Yellow/Red
        # of the green base, which are Yellow at visit n
    base <- apply(series, 1, function(rr) {
        zz <- !is.na(rr) & rr > .ArrestEnv$YELLOW
        sum(zz) >= numBaselines
    })

    z.y <- base & series[, n] == .ArrestEnv$YELLOW
    ps.y <- sapply(which(z.y), function(loc) {
        ys <- series[loc, ]
        z <- ys > .ArrestEnv$YELLOW
        mean_green <- mean(ys[z])
        1 - pr_seeing(.ArrestEnv$B, tt = mean_green, fpr = 0.15, fnr = 0.03)^2
    })

        # Yellows/Greens -> Red
    base <- apply(series, 1, function(rr) {
        zz <- !is.na(rr) & rr > .ArrestEnv$RED
        sum(zz) >= numBaselines
    })

    z.r <- base & series[, n] == .ArrestEnv$RED
    ps.r <- ifelse(any(z.r), 0, 1)

    z <- series[, n] > .ArrestEnv$YELLOW
    ps.p <- PoPLR(series[z, ], threshold = threshold, perm_count = perm_count, warnings = warnings)

    return(min(c(ps.y, ps.r, ps.p)))
}