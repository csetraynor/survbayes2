# prepare data function -----
gen_stan_dat <- function(x, status = "status", time = "time") {
  # prepare for longdat formating
  x$sample_id <- 1:nrow(x)  #create sample id
  # get unique times: only event times equivalent to Cox model
  x <- x[order(x$time), ]
  times <- x[x[["status"]], ]
  time_points <- times[order(unique(unlist(times[, "time"]))), "time"]
  dtime_points <- diff(c(0, time_points))
  x$t_observed <- x$time
  x$time <- NULL
  
  tt_dat <- plyr::ddply(x, "sample_id",
                        function(x) data.frame(time = time_points, t_id = seq_along(time_points), dtime = dtime_points, t_observed = x$t_observed )  )#create dataset with real observation
  tt_dat <- tt_dat[tt_dat$time <= tt_dat$t_observed, ]
  tt_dat$t_observed <- NULL
  tt_dat <- merge(tt_dat, x, by = "sample_id")
  
  tt_dat$status_observed <- tt_dat$status
  tt_dat$status <- 0
  tt_dat$status[tt_dat$time == tt_dat$t_observed] <- as.numeric( tt_dat$status_observed[tt_dat$time == tt_dat$t_observed] )
  tt_dat$log_dtime <- as.double( unlist( log(tt_dat$dtime) ) )
  tt_dat
}
