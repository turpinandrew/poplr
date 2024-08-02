require(testthat)

source('poplr.r')
# two series with 52 locations: stable and steep
slopes <- c(0, -5)
baseSeries <- lapply(slopes, function(slope) 
    matrix(sapply(30 + seq(0, by = slope, length = 10), function(db) rnorm(52, db, 2)), ncol = 10)
)

#' @param series a matrix with >= 52 rows
#' @param perm A permutation of 1:ncol(series)
#' @return TRUE if all NAs in series[, perm] are at the start of rows
check_perm <- function(series, perm) {
    ss <- tail(series[, perm], -52)
    z <- apply(ss, 1, function(rr) all(diff(which(is.na(rr))) == 1))
    all(z)
}

test_that("get_boundaries_1block", {
    rr <- baseSeries[[2]][seq_len(5), ]  # take the first 5 rows
    rr[, 1:3] <- NA  # set the first 3 visits to NA
    ss <- rbind(baseSeries[[2]], rr)
    b <- get_boundaries(ss)
    expect_equal(b, c(3))
})

test_that("get_boundaries", {
    expect_null(get_boundaries(baseSeries[[1]]), label = "52 rows, slope 0")
    expect_null(get_boundaries(baseSeries[[2]]), label = "52 rows, slope 1")

    newLocations <- list(
        list(c(1, 1), c(1)),   # 1 new row, 1 NA visits
        list(c(1, 3), c(3)),   # 1 new row, 3 NA visits
        list(c(1, 10), c(10)),  # 1 new row, 10 NA visits
        list(c(5, 1), c(1)),    # 5 new rows, 1 NA visits
        list(c(1, 2), c(1, 3), c(2, 3)), # 1 new rows, 2 NA visits + 1 rows 3 NA visits
        list(c(5, 2), c(5, 3), c(2, 3)), # 5 new rows, 2 NA visits + 5 rows 3 NA visits
        list(c(1, 3), c(1, 4), c(1, 5), c(3, 4, 5)),
        list(c(2, 3), c(2, 4), c(2, 5), c(3, 4, 5)),
        list(c(5, 3), c(5, 4), c(5, 5), c(3, 4, 5)),
        list(c(1, 3), c(1, 5), c(3, 5)),
        list(c(2, 3), c(1, 5), c(3, 5)),
        list(c(1, 3), c(2, 5), c(3, 5)),
        list(c(1, 2), c(10, 4), c(2, 4)),
        list(c(1, 9), c(9)),
        list(c(1, 10), c(10))
    )

    for (i_newLocs in seq_along(newLocations)) {
        s <- baseSeries
        for (loc in head(newLocations[[i_newLocs]], -1))
            s <- lapply(s, function(base) {
                nRows <- loc[[1]]
                nVisits <- loc[[2]]
                rr <- base[seq_len(nRows), , drop = FALSE]  # take the first nRows rows
                rr[, seq_len(nVisits)] <- NA  # set the first nVisits visits to NA
                rbind(base, rr)               # append them to base
            })

        for (i_slope in seq_along(slopes)) {
            b <- get_boundaries(s[[i_slope]])

            expect_true(all(b == unlist(tail(newLocations[[i_newLocs]], 1))),
                label = paste("slope=", slopes[i_slope], "newLocations=", i_newLocs))
        }
    }
})

test_that("get_perms4", {
    rr <- baseSeries[[2]][seq_len(5), ]  # take the first 5 rows
    rr[, 1:3] <- NA  # set the first 3 visits to NA
    ss <- rbind(baseSeries[[2]], rr)
    p <- get_perms4(ss, perm_count = 20)

    expect_equal(p[[1]], 1:10, label = "first is 1..n")
})

newLocations <- list(
    list(c(1, 1), c(1)),   # 1 new row, 1 NA visits
    list(c(1, 3), c(3)),   # 1 new row, 3 NA visits
    list(c(5, 1), c(1)),    # 5 new rows, 1 NA visits
    list(c(1, 2), c(1, 3), c(2, 3)), # 1 new rows, 2 NA visits + 1 rows 3 NA visits
    list(c(5, 2), c(5, 3), c(2, 3)), # 5 new rows, 2 NA visits + 5 rows 3 NA visits
    list(c(1, 3), c(1, 4), c(1, 5), c(3, 4, 5)),
    list(c(2, 3), c(2, 4), c(2, 5), c(3, 4, 5)),
    list(c(5, 3), c(5, 4), c(5, 5), c(3, 4, 5)),
    list(c(1, 3), c(1, 5), c(3, 5)),
    list(c(2, 3), c(1, 5), c(3, 5)),
    list(c(1, 3), c(2, 5), c(3, 5)),
    list(c(1, 2), c(10, 4), c(2, 4)),
    list(c(1, 9), c(9))
)

 #' Add rows given by newLocations[[i_newLocs]] to baseSeries
 #' for each element of baseSeries.
 #'
 #' @return list of series, one for each of `slopes``
 #
 make_series <- function(i_newLocs) {
    s <- baseSeries
    for (loc in head(newLocations[[i_newLocs]], -1))
        s <- lapply(s, function(base) {
            nRows <- loc[[1]]
            nVisits <- loc[[2]]
            rr <- base[seq_len(nRows), , drop = FALSE]  # take the first nRows rows
            rr[, seq_len(nVisits)] <- NA  # set the first nVisits visits to NA
            rbind(base, rr)               # append them to base
        })
    s
 }

 for (i_newLocs in seq_along(newLocations))
    test_that(paste("newLocations=", i_newLocs), {
        s <- make_series(i_newLocs)
        n <- 20
        for (i_slope in seq_along(slopes)) {
#print(s[[i_slope]])
             p <- get_perms4(s[[i_slope]], n)
#print(p)

             expect_equal(unlist(p[[1]]), 1:10, label = "first is 1..n")
             expect_equal(length(p), min(n, factorial(ncol(s))), label = paste("has ncol! or ", n))
             expect_true(all(unlist(lapply(p, function(pp) check_perm(s[[i_slope]], pp)))), label = "check_perm")
        }
    })

for (i_newLocs in seq_along(newLocations))
    test_that(paste("PoPLR newLocations=", i_newLocs), {
        s <- make_series(i_newLocs)

        for (i_slope in seq_along(slopes)) {
            p <- PoPLR4(s[[i_slope]])$p

            expect_true(p <= 1, label = "p <= 1")
        }
    })