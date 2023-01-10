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
source('arrest5.r')
source('poplr.r')
PoPLR_ARREST <- function(series, threshold = 0.05, perm_count = 5000, warnings = TRUE,
                         numBaselines = 2,
                         eventP = 0.00001) {
    n <- ncol(series)

    stopifnot(n > 2)

    z <- apply(series, 1, function(rr) sum(!is.na(rr)) > numBaselines)

        # Green -> Yellow/Red
    z.y <- z & series[, n] == .ArrestEnv$YELLOW
    if (sum(z.y) > 0 && any(series[z.y, n - 1] > .ArrestEnv$YELLOW))
        return(eventP)

        # a Orange -> Red
    z.y <- z & series[, n] == .ArrestEnv$RED
    if (sum(z.y) > 0 && any(series[z.y, n - 1] > .ArrestEnv$RED))
        return(eventP)

    z <- series[, n] > .ArrestEnv$YELLOW
    return(PoPLR(series[z, ], threshold = threshold, perm_count = perm_count, warnings = warnings))
}