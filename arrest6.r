#
# Same as arrest5.r, but allows new locations between greens.
#
# Andrew Turpin
# Thu  5 Jan 2023 10:06:43 AWST
#
#require(keypress)
has_keypress_support <- function() return(FALSE)

source('graph_common.r')

############################################################
# Set up ARREST constants.
############################################################
.ArrestEnv <- new.env()

assign("FPR", 0.03,     envir=.ArrestEnv)
assign("FNR", 0.03,     envir=.ArrestEnv)
assign("MAX_PRES", 250, envir=.ArrestEnv)

assign("B",  16,        envir=.ArrestEnv)         # even and > 0. First value in the BUCKET

    # The first 4 MUST be less than B
assign("NEVER_MEASURED" , -3, envir=.ArrestEnv)
assign("RED"            , -2, envir=.ArrestEnv)  # Red at last visit, and the one before.
assign("YELLOW_RED"     , -1, envir=.ArrestEnv)  # Red at last visit (couldn't see 0 x2)
assign("YELLOW"         ,  0, envir=.ArrestEnv)  # measured in [0, B]

assign("EST_NP_NEW"     , 10, envir=.ArrestEnv)        # Estimated num presentations for an added new location
assign("EST_NP_YEL"     ,  1, envir=.ArrestEnv)        # Estimated num presentations for now yellow

assign("COLORS", function(t) {   # brewer.pal(3, 'RdYlBu')
    if (is.na(t)) return(NA)
    if (t <= -1) return("#FC8D59AA")
    if (t <= 17) return("#FFFFBFcc")
    return("#91BFDBAA")
}, envir=.ArrestEnv)

assign("STAGE_1_MISSED", 2, envir=.ArrestEnv)
assign("STAGE_2_MISSED", 2, envir=.ArrestEnv)

############################################################
# Returns a list of states, one per location to be tested.
# If there are no previous results, returns a state for each
# "connected" node in graph (as defined by get_connected_nodes()).
# If there is a previous result, returns a state for each previously
# tested location, and also extra locations if it is estimated
# that there are presentations to spare. Note this might also
# alter the graph structure as new nodes become connected.
#
# @param graph Connectivity graph to drive ARREST.
# @param prev_result Previous result from ARREST (or NULL if first test)
#                    is a data.frame with 4 columns: 'x', 'y', 'npres', 'final'
# @param makeStimHelper function(x,y) -> (function(db,n)->opistimulus)
# @param opiParams for OPI simulation: list indexed by vertex number in graph.
#
# @return list of
#       states = list of [loc (vertex number), prev_mt, makeStim, opiParams, + base_state list] lists.
#       graph  = possibly altered graph
############################################################
ARREST.start <- function(graph, prev_result, makeStimHelper, opiParams=NULL)  {
    states <- NULL

    base_state <- list(
                stage = 0,
                finished = FALSE,
                final = NA,                           # final measured value
                numPresentations = 0,                 # number of presentations so far
                stimuli = NULL,                       # vector of stims shown
                responses = NULL,                     # vector of responses (1 seen, 0 not)
                responseTimes = NULL                  # vector of response times
    )

    ten_d_2_damaged <- prev_result$x^2 +prev_result$y^2 <= 81 & prev_result$final <= .ArrestEnv$YELLOW

    if (is.null(prev_result) || all(!ten_d_2_damaged)) {
            # Create a NEVER_MEASURED state for all initial locations
            # (locations that have an edge in the start graph)
        states <- lapply(get_connected_nodes(graph), function(loc) {
            c(list(loc = loc, prev_mt = .ArrestEnv$NEVER_MEASURED), base_state)
        })
    } else {
            # Create a state for all previously measured locations (including RED)
            # and add new locations near red/yellow if we have enough presentations.

            # Create states for prev and collect mt_by_vertex_num as a side effect
        mt_by_vertex_num <- rep(NA, get_max_vertex(graph))
        states <- lapply(seq_len(nrow(prev_result)), function(i) {
            v <- find_vertex_num(graph, prev_result$x[i], prev_result$y[i])
            prev_mt <- prev_result$final[i]

            mt_by_vertex_num[v] <<- prev_mt

            s <- c(list(loc = v, prev_mt = prev_mt), base_state)
            if (prev_mt == .ArrestEnv$RED) {
                s$final <- prev_mt
                s$finished <- TRUE
            }

            return(s)
        })

            # Add locations between any location pairs if there is presentations to spare.

            # Get queue of potential vertex numbers: bucket_locs
            # and find a neighbour for each of them
        #bucket_locs <- sapply(which(prev_result$final <= .ArrestEnv$YELLOW), function(i) { states[[i]]$loc })
#sink("x")
#print("Inital")
        bucket_locs <- unlist(lapply(states, "[", "loc"), use.names = FALSE)
#print(bucket_locs)

            # neighbs[[i]] = list(gradient, u, v, list(edge vertices))
            # ASSUMES YELLOW and RED_YELLOW and RED values are used as mt
        neighbs <- NULL
        for (u in bucket_locs)
            for (possible in find_neighbour_sets6(graph, u, mt_by_vertex_num))
                if (!is.null(possible))
                    neighbs <- c(neighbs, list(possible))

        if (length(neighbs) > 0) {
            o <- order(unlist(lapply(neighbs, "[", 1)), decreasing = TRUE)
            neighbs <- neighbs[o]
        }

#print(neighbs)
#for (i in seq_along(neighbs)) {
#    if (length(neighbs[[i]]) > 0) {
#        cat(sprintf("\n%2.0f %2.0f %2.0f: <", neighbs[[i]][[1]], neighbs[[i]][[2]], neighbs[[i]][[3]]))
#        for (p in seq_along(neighbs[[i]][[4]]))
#            cat(sprintf("%2.0f ", neighbs[[i]][[4]][[p]]))
#        cat(">")
#    }
#}
#cat("\n")

            zg <- prev_result$final > .ArrestEnv$YELLOW    # green
            zy <- !zg & prev_result$final > .ArrestEnv$RED # YELLOW or YELLOW_RED
            npres_estimate <- sum(prev_result$npres[zg]) + .ArrestEnv$EST_NP_YEL * sum(zy)

            used <- NULL
            while (length(neighbs) > 0 && npres_estimate < .ArrestEnv$MAX_PRES) {
#print("=====")
#print(neighbs)
                candidates <- neighbs[[1]]          # list of "v" vertices: [gradient, u, v, list of possible new locs]
                new_loc <- head(candidates[[4]], 1)
#print(new_loc)
                candidates[[4]] <- tail(candidates[[4]], -1)  # remove new_loc from neighbs
                if (length(candidates[[4]]) == 0) {           # perhaps remove whole node too
                    neighbs[[1]] <- NULL
                } else {
                        # put updated candidates back into neighbs
                    if (length(neighbs) > 1 && neighbs[[2]][[1]] == candidates[[1]]) { # give #2 in neighbs a chance
                        neighbs[[1]] <- neighbs[[2]]
                        neighbs[[2]] <- candidates
                    } else {
                        neighbs[[1]] <- candidates
                    }
                }

                if (!new_loc %in% used) {
                    used <- c(used, new_loc)

                    s <- c(list(loc = new_loc, prev_mt = .ArrestEnv$NEVER_MEASURED), base_state)
                    states <- c(states, list(s))

                    npres_estimate <- npres_estimate + tail(.ArrestEnv$EST_NP_NEW, 1)

                    graph <- split_edge(graph, new_loc)
                }
            }# end test more while loop
        }#end if len(neighbs) > 0
#sink(NULL)

        # add in makeStim and opiParams to all states
    for (i in 1:length(states)) {
        states[[i]]$opiParams <- opiParams[[states[[i]]$loc]]
        xy <- get_vertex_xy(graph, states[[i]]$loc)
        states[[i]]$makeStim <- makeStimHelper(xy$X, xy$Y)
    }

    return(list(states = states, graph = graph))
}# ARREST.start()

############################################################
# Return TRUE if location is finished (no more presentations)
# FALSE otherwise
############################################################
ARREST.stop <- function(state) { return (state$finished) }

############################################################
# Return final estimate of threshold
############################################################
ARREST.final <- function(state) { return (state$final) }

############################################################
# Advance state by one presentation and return the new state.
# @param state as returned by ARREST.start()
# @param nextStim an opiStaticStimulus giving location of next presentation
# @return list of
#   state as returned by ARREST.start()
#   resp as returned by opiPresent
############################################################
ARREST.step <- function(state, nextStim = NULL) {
    if (ARREST.stop(state))
        stop('ARREST.step(): cannot step a stopped state')

        # Do a single presentation of x dB
    do_one_pres <- function(x) {
        params <- c(list(stim = state$makeStim(x, 0), nextStim = nextStim), state$opiParams)
        opiResp <- do.call(opiPresent, params)
        while (!is.null(opiResp$err))
            opiResp <- do.call(opiPresent, params)
        state$stimuli          <- c(state$stimuli, 0)
        state$responses        <- c(state$responses, opiResp$seen)
        state$responseTimes    <- c(state$responseTimes, opiResp$time)
        state$numPresentations <- state$numPresentations + 1

        return(opiResp)
    }

        # YELLOW_RED.
        # Stage==0 if see 0, goto YELLOW and finished
        #          if !see 0, goto stage=1
        # Stage==1 if see 0, goto YELLOW and finished
        #          if !see 0, goto RED and finished
        #
        # YELLOW.
        # Same logic but goto YELLOW_RED on !see && Stage==1
    if ((state$prev_mt == .ArrestEnv$YELLOW_RED) || (state$prev_mt == .ArrestEnv$YELLOW)) {
        opiResp <- do_one_pres(0)

        if (opiResp$seen) {
            state$finished <- TRUE
            state$final <- .ArrestEnv$YELLOW
        } else if (state$stage == 0) {
            state$stage <- 1
        } else if (state$stage == 1) {
            state$finished <- TRUE
            if (state$prev_mt == .ArrestEnv$YELLOW_RED) {
                state$final <- .ArrestEnv$RED
            } else {
                state$final <- .ArrestEnv$YELLOW_RED
            }
        }

        return(list(state = state, resp = opiResp))
    }

    # Here we are either GREEN or NEVER_MEASURED.
    # Stage == 0: ZEST.start and set stage = 1 and execute stage 1
    # Stage == 1: ZEST.step
    #             if ZEST finished
    #                if ZEST result > B
    #                   finished =TRUE, final=ZEST_final
    #                else
    #                   else stage = 2
    # Stage == 2: if B seen final=B+1, finished=TRUE
    #             else stage = 3
    # Stage == 3: if B seen final=B+1, finished=TRUE
    #             else stage = 4
    # Stage == 4: if 0 seen final=YELLOW, finished=TRUE
    #             else stage = 5
    # Stage == 5: if 0 seen final=YELLOW, finished=TRUE
    #             else final = YELLOW_RED, finished=TRUE
    if (state$stage == 0) {
        state$stage <- 1
        domain <- -5:40
        prior <- c(rep(0.01,5), 0.12, dnorm(1:20, 0, 10), 4*dnorm(21:40, 30, 3))
        prior <- prior/sum(prior)

        params <- c(list(domain=domain, prior=prior, minStimulus=0, makeStim=state$makeStim), state$opiParams)
        state$ZEST <- do.call(ZEST.start, params)
    }
    if (state$stage == 1) {
        r <- ZEST.step(state$ZEST)
        opiResp <- r$resp
        state$ZEST <- r$state
        state$stimuli          <- c(state$stimuli, tail(state$ZEST$stimuli, 1))
        state$responses        <- c(state$responses, opiResp$seen)
        state$responseTimes    <- c(state$responseTimes, opiResp$time)
        state$numPresentations <- state$numPresentations + 1

        if (ZEST.stop(state$ZEST)) {
            if (ZEST.final(state$ZEST) > .ArrestEnv$B) {
                state$final <- ZEST.final(state$ZEST)
                state$finished <- TRUE
            } else {
                state$stage <- 2
            }
        }
    } else if (state$stage == 2) {
        opiResp <- do_one_pres(.ArrestEnv$B)
        if (opiResp$seen) {
            state$finished <- TRUE
            state$final <- .ArrestEnv$B + 1
        } else {
            state$stage <- 3
        }
    } else if (state$stage == 3) {
        opiResp <- do_one_pres(.ArrestEnv$B)
        if (opiResp$seen) {
            state$finished <- TRUE
            state$final <- .ArrestEnv$B + 1
        } else {
            state$stage <- 4
        }
    } else if (state$stage == 4) {
        opiResp <- do_one_pres(0)
        if (opiResp$seen) {
            state$finished <- TRUE
            state$final <- .ArrestEnv$YELLOW
        } else {
            state$stage <- 5
        }
    } else if (state$stage == 5) {
        opiResp <- do_one_pres(0)
        if (opiResp$seen) {
            state$final <- .ArrestEnv$YELLOW
        } else {
            state$final <- .ArrestEnv$YELLOW_RED
        }
        state$finished <- TRUE
    }

    return(list(state = state, resp = opiResp))
}

###################################################################
# A whole field test. Randomly selects a location to present, allows
# for some inter-stimulus interval and also false positive catch trials.
#
# @param graph Connectivity graph to drive ARREST.
# @param prev_result Previous result from ARREST (or NULL if first test)
#                    is a data.frame with 4 columns: 'x', 'y', 'npres', 'final'
# @param makeStimHelper function(x,y) -> (function(db,n)->opiStimulus)
#
# @param opiParams for OPI simulation: list indexed by vertex number in graph.
#
# @param min_isi  Minimum inter-stimulus interval in ms.
# @param max_isi  Maximum inter-stimulus interval in ms.
#
# @param fp_check Present fp_level every fp_check presentations.
# @param fp_level dB level of false positive check.
#
# @verbose if TRUE will print information about each trial
# @progress_freq will print number of locations remaining every this presentations
# @allow_pause if TRUE, will pause on keypress.
#
# Note fp check in simulation is done with opiParams of last presentation.
#
# @return list of
#    results data.frame with 4 columns: 'x', 'y', 'npres', 'final'
#    graph   possibly altered input graph
#    fp      c(num_fp_trials_seen, num_fp_trials_presented)
#    states  The list of state per location at the end of the test
###################################################################
ARREST <- function(graph, prev_result, makeStimHelper, opiParams,
                   min_isi = 0, max_isi = 0, fp_check = 30, fp_level = 60,
                   verbose = FALSE,
                   progress_freq = NA,
                   allow_pause = FALSE) {
    if (allow_pause && !has_keypress_support()) {
        warning("Cannot use key press to pause - perhaps you are in RStudio")
        allow_pause <- FALSE
    }

    sg <- ARREST.start(graph = graph, prev_result = prev_result, makeStimHelper = makeStimHelper, opiParams = opiParams)
    states <- sg$states
    graph  <- sg$graph

    fp_presented <- fp_seen <- 0
    count <- 0
    unfinished <- which(!unlist(lapply(states, "[", "finished")))
    cur <- round(runif(1, 1, length(unfinished)))
    nxt <- round(runif(1, 1, length(unfinished)))
    while (length(unfinished) > 0) {
        if (!is.na(progress_freq) && progress_freq > 0 && count %% progress_freq == 0)
            cat(sprintf("Progress: %3.0f locations remaining of %3.0f locations (%2.0f%%)\n",
                length(unfinished), length(states), length(unfinished)/length(states)*100))

        if (allow_pause && keypress(block = FALSE) != "none") {
            if (interactive()) {
                invisible(readline(prompt = "Paused. Press <Enter> to continue..."))
            } else {
                cat("Paused. .Press <Enter> to continue...")
                invisible(readLines(file("stdin"), 1))
            }
        }

        count <- count + 1

        curr_state <- states[[unfinished[cur]]]
        next_state <- states[[unfinished[nxt]]]

        r <- ARREST.step(curr_state, nextStim = next_state$makeStim(0,0))
        curr_state <- states[[unfinished[cur]]] <- r$state

        xy <- get_vertex_xy(graph, curr_state$loc)

        if (verbose)
            cat(sprintf("x= %+5.1f y= %+5.1f db= %5.2f seen= %5s time= %6.2f type= A\n", xy$X, xy$Y,
                tail(curr_state$stimuli,1), r$resp$seen, r$resp$time))

        if (ARREST.stop(curr_state)) {
            unfinished <- unfinished[-cur]
            if (nxt > cur)
                nxt <- nxt - 1
            else if (length(unfinished) > 0 && nxt == cur)
                nxt <- round(runif(1, 1, length(unfinished)))
        }

        Sys.sleep(runif(1, min = min_isi, max = max_isi))

        if (count %% fp_check == 0) {
            fp_presented <- fp_presented + 1
            params <- c(list(stim = curr_state$makeStim(fp_level, 0), nextStim = NULL), curr_state$opiParams)
            opiResp <- do.call(opiPresent, params)
            while (!is.null(opiResp$err))
                opiResp <- do.call(opiPresent, params)

            if (opiResp$seen)
                fp_seen <- fp_seen + 1

            if (verbose)
                cat(sprintf("x= %+5.1f y= %+5.1f db= %5.2f seen= %5s time= %6.2f type= FP\n", xy$X, xy$Y,
                    fp_level, opiResp$seen, opiResp$time))
        }

        cur <- nxt
        if (length(unfinished) > 0)
            nxt <- round(runif(1, 1, length(unfinished)))
    }

    locs <- unlist(lapply(states, "[", "loc"))
    xs <- unlist(lapply(locs, function(loc) get_vertex_xy(graph, loc)$X))
    ys <- unlist(lapply(locs, function(loc) get_vertex_xy(graph, loc)$Y))

    res <- cbind(
        x = xs, y = ys,
        npres = unlist(lapply(states, "[", "numPresentations")),
        final = unlist(lapply(states, "[", "final"))
    )
    rownames(res) <- NULL
    return(list(results = as.data.frame(res), graph = graph, fp = c(fp_seen, fp_presented), states = states))
}


#######################################################################################
# Test with stables and progs
##########################################################################################
###require(OPI)
###chooseOpi("SimHenson")
###if (!is.null(opiInitialize(type = "C", cap = 6)))
###    stop("opiInitialize failed")
###
###makeStimHelper <- function(x, y) {
###    ff <- function(db, n) db + n  # dummy
###
###    body(ff) <- substitute(
###        {s <- list(x = x, y = y, level = dbTocd(db, 10000 / pi), size = 0.43, color = "white",
###                  duration = 200, responseWindow = 1500)
###         class(s) <- "opiStaticStimulus"
###         return(s)
###        }
###        , list(x = x, y = y))
###    return(ff)
###}
###
###require(parallel)

###     ########
###     # stables
###     ########
### load('graph_242_OS.Rdata')
### #load('surfaces2.Rdata')   # for xy_grid
### #load('stable_surfaces.Rdata')
### #load("stable_surfaces_MTD-6.Rdata")
### #stable_surfaces[[4]] <- stable_surface_MTD_6
### #  # just keep every 2 degrees
### #z <- xy_grid[,1] %% 2 == 1 & xy_grid[,2] %% 2 == 1
### #xy_grid <- xy_grid[z,]
### #rownames(xy_grid) <- 1:nrow(xy_grid)
### #colnames(xy_grid) <- c("X", "Y")
### #
### #stable_surfaces <- lapply(stable_surfaces, function(ss) {
### #    lapply(ss, function(m) m[z,])
### #    })
###
### load('surfaces3.Rdata')   # for xy_grid
### load('stable_surfaces4.Rdata')
###
###
### print("arrest5 on stables")
### print(Sys.time())
### res <- lapply(c(1,2,3,4,5), function(i) {
###     if (i == 2 || i == 3)
###         return (NA)
###     lapply(stable_surfaces[[i]], function(surf) {
###       mclapply(1:100, function(repi) {
###         results <- NULL
###         graph <- graph_242_2deg
###         res <- NULL
###         for (visit in 1:ncol(surf)) {
###             opiParams <- lapply(1:get_max_vertex(graph), function(v) {
###                 xy <- as.numeric(get_vertex_xy(graph_242_2deg, v))
###                 z <- which(xy_grid$X == xy[1] & xy_grid$Y == xy[2])
###                 c(tt = surf[z,visit], fpr = 0.15, fnr = 0.03)
###             })
###             rgf <- ARREST(graph = graph, prev_result = res, makeStimHelper = makeStimHelper, opiParams = opiParams)
###             res <- rgf$results
###             graph <- rgf$graph
###             results <- c(results, list(res))
###         }
###         return(results)
###       }, mc.cores = 16)
###     })
### })
### #save(res, file = arrest5.stable.2.Rdata")
### #save(res, file = arrest5.stable_type5.3.Rdata")
### #save(res, file = arrest5.stable.3.Rdata")  # btype  =  1,4,5
### #save(res, file = arrest5.stable.3.fp15.Rdata")  # btype = 1,4,5
### print(Sys.time())
###
###     ########
###     # progs
###     ########
### #load('surfaces2.Rdata')
### #z <- xy_grid[,1] %% 2 == 1 & xy_grid[,2] %% 2 == 1
### #xy_grid <- xy_grid[z,]
### #rownames(xy_grid) <- 1:nrow(xy_grid)
### #colnames(xy_grid) <- c("X", "Y")
### #surfs <- lapply(surfs, function(m) m[z,])
###
### load('surfaces3.Rdata')
### load('graph_242_OS.Rdata')
###
### print("arrest5 on progs")
### print(Sys.time())
### res <- mclapply(surfs, function(surf) {
###         lapply(1:100, function(repi) {
###             results <- NULL
###             graph <- graph_242_2deg
###             res <- NULL
###             for (visit in 1:ncol(surf)) {
###                 opiParams <- lapply(1:get_max_vertex(graph), function(v) {
###                     xy <- as.numeric(get_vertex_xy(graph_242_2deg, v))
###                     z <- which(xy_grid$X == xy[1] & xy_grid$Y == xy[2])
###                     return(c(tt=surf[z,visit], fpr=0.15, fnr=0.03))
###                 })
###                 rgf <- ARREST(graph=graph, prev_result=res, makeStimHelper=makeStimHelper, opiParams=opiParams)
###                 res <- rgf$results
###                 graph <- rgf$graph
###                 results <- c(results, list(res))
###             }
###             return(results)
###         })
### }, mc.cores=16)
### #save(res, file="arrest5.2.Rdata")
### save(res, file="arrest5.2.fp15.Rdata")
### print(Sys.time())

#####################################
#### Random Test
###source("graph_common.r")
###
###load("graph_242_102.Rdata")
###graph <- graph_242_102_OS
###
###pdf("arrest6-test.pdf", width = 10, height = 3)
###options(error = dev.off)
###
###layout(matrix(1:3, 1, 3))
###par(mar = c(3, 3, 2, 0))
###
###cols <- function(t, arrest = TRUE) {   # brewer.pal(3, "RdYlBu")
###    if (is.na(t)) return(NA)
###    if (t <= -1 || (!arrest && t <= 0)) return("#FC8D59AA")
###    if (t <= 16) return("#FFFFBFcc")
###    return("#91BFDBAA")
###}
###tts <- round(runif(1, -3, 3) + rep(30, get_max_vertex(graph)))
###seed <- sample(1:get_max_vertex(graph), 1)
###z <- (graph$xy$X - graph$xy$X[seed])^2 + (graph$xy$Y - graph$xy$Y[seed])^2 < 49
###tts[z] <- 0
###
###opiParams <- lapply(1:get_max_vertex(graph), function(v) {
###    c(tt = tts[[v]], fpr = 0.03, fnr = 0)
###})
###
###plot(graph$xy$X, graph$xy$Y, las = 1, asp = 1, pch = 19, col = sapply(tts, cols), cex = 2, main = "TRUE")
###text(graph$xy$X, graph$xy$Y, tts, cex = 0.5)
###text(graph$xy$X, graph$xy$Y, seq_along(graph$xy$Y), cex = 0.25, pos = 1)
###
###results <- NULL
###for (i in 1:20) {
###    rgf <- ARREST(graph = graph, prev_result = results, makeStimHelper = makeStimHelper, opiParams = opiParams, verbose = FALSE)
###    results <- rgf$results
###    graph <- rgf$graph
###    fp <- rgf$fp
###
###    plot(results$x, results$y, las = 1, asp = 1, pch = 19, col = sapply(results$final, cols), cex = 2, main = paste("Visit", i))
###    text(results$x, results$y, f <- round(results$final, 0), cex = 0.5)
###
###    #print(paste("Total pres=", sum(results$npres), "fp=", fp[1], "/", fp[2], "=", round(fp[1]/fp[2]*100,1), "%"))
###}
###dev.off()
###options(error = NULL)