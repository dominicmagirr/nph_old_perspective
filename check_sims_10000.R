library(tidyverse)
library(survival)
###########################################################
####### simulate from a piece-wise exponential distribution
t_piecewise_exp <- function(n = 10, 
                            change_points = c(6, 12),
                            lambdas = c(log(2) / 9, log(2) / 9, log(2) / 9)){
  
  t_lim <- matrix(rep(c(diff(c(0, change_points)), Inf), each = n), nrow = n)
  t_sep <- do.call(cbind, purrr::map(lambdas, rexp, n = n))
  which_cells <- t(apply(t_sep < t_lim, 1, function(x){
    rep(c(T,F), c(min(which(x)), length(x) - min(which(x))))
  } ))
  rowSums(pmin(t_sep, t_lim) * which_cells)
}
##########################################################
######## simulate data according to trial parameters
sim_t_uncensored <- function(model,
                             recruitment){
  
  rec_0 <- recruitment$r_period * runif(recruitment$n_0) ^ (1 / recruitment$k)
  rec_1 <- recruitment$r_period * runif(recruitment$n_1) ^ (1 / recruitment$k)
  
  time_0 <- t_piecewise_exp(recruitment$n_0, model$change_points, model$lambdas_0)
  time_1 <- t_piecewise_exp(recruitment$n_1, model$change_points, model$lambdas_1)
  
  data.frame(time = c(time_0, time_1),
             rec = c(rec_0, rec_1),
             group = rep(c("control", "experimental"), c(recruitment$n_0, recruitment$n_1)))
  
}
#########################################################
## apply data cut off to simulated data set
apply_dco <- function(df,
                      dco = NULL,
                      events = NULL){
  
  if (is.null(dco) && is.null(events)) stop("Must specify either dco or events")
  
  df$cal_time <- df$time + df$rec
  
  if (is.null(dco)){
    dco <- sort(df$cal_time)[events]
  }
  
  df_dco <- df[df$rec < dco, ]
  df_dco$event <-  df_dco$cal_time <= dco
  df_dco$time <- pmin(df_dco$time, dco - df_dco$rec)
  df_dco$dco <- dco
  
  df_dco
  
}
##########################################################
## produce risk table
rt_row <- function(t,
                   df,
                   trt_colname,
                   time_colname,
                   event_colname){
  
  df1 = df[df[[trt_colname]] == unique(df[[trt_colname]])[1],]
  df2 = df[df[[trt_colname]] == unique(df[[trt_colname]])[2],]
  
  data.frame(time = t,
             r1 = sum(df1[[time_colname]] >= t),
             r2 = sum(df2[[time_colname]] >= t),
             e1 = sum(df1[[time_colname]] == t & df1[[event_colname]] == 1),
             e2 = sum(df2[[time_colname]] == t & df2[[event_colname]] == 1))
}
#########################################################
### apply a landmark analysis
landmark = function(df, time){
  
  ## use survival::survfit to do KM estimation by group:
  fit <- survival::survfit(Surv(time, event) ~ group, data = df)
  
  no_data <- min(fit$time[fit$n.risk == 1]) < time
  time <- time[!no_data]
  #if(any(min(fit$time[fit$n.risk == 1]) < time)) return(rep(NA, length(time)))
  
  info = summary(fit, time = time)
  
  
  ## extract survival probabilities with standard errors:
  
  info_df = data.frame(t = info$time,
                       p = info$surv,
                       se = info$std.err,
                       group = info$strata)
  
  
  
  z_stats = info_df %>% 
    group_by(t) %>%
    summarize(diff = p[group == "group=experimental"] - p[group == "group=control"],
              se = sqrt(sum(se ^ 2))) %>%
    mutate(z = diff/se)
  
  c(z_stats$z, rep(NA, sum(no_data)))
  
}
#########################################################
### Simulate 1 trial and apply all different analysis
### methods simultaneously
sim_trial <- function(dco, model, recruitment){
  
  df_uncensored <- sim_t_uncensored(model, recruitment)
  
  df_final <- apply_dco(df_uncensored, dco = dco)
  
  #########################################
  rt <- purrr::map_df(unique(df_final[["time"]]),
                      rt_row,
                      df = df_final,
                      trt_colname = "group",
                      time_colname = "time",
                      event_colname = "event")
  
  rt <- dplyr::filter(dplyr::arrange(rt, time),
                      e1 > 0 | e2 > 0)
  
  ##########################################
  #### get summary
  s_sum <- summary(survival::survfit(survival::Surv(time, event) ~ 1,
                                     data= df_final,
                                     timefix = FALSE))
  
  ### get survival probabilities for the pooled data
  s_pool <- s_sum$surv
  
  
  ##########################################
  ### get weights
  w_lr <- rep(1, length(s_pool))
  w_fh_0_1 <-  (1 - c(1, s_pool[-length(s_pool)])) 
  w_fh_0_half <- (1 - c(1, s_pool[-length(s_pool)])) ^ 0.5
  ##########################################
  if(any(s_sum$time < 12)){
    
    w_mw_12 <- pmin(1 / c(1, s_pool[-length(s_pool)]),
                    1 / s_pool[max(which(s_sum$time < 12))])
  }
  else {
    w_mw_12 <- rep(1, length(s_pool))
  }
  ###########################################
  if(any(s_sum$time < 24)){
    
    w_mw_24 <- pmin(1 / c(1, s_pool[-length(s_pool)]),
                    1 / s_pool[max(which(s_sum$time < 24))])
  }
  else {
    w_mw_24 <- rep(1, length(s_pool))
  }
  ###########################################
  get_z <- function(w, rt){
    u = sum(w * (rt$e2 - (rt$e1+rt$e2) * rt$r2 / (rt$r1+rt$r2)))
    v_u = sum(w^2 * rt$r1 * rt$r2 * (rt$e1+rt$e2) * (rt$r1+rt$r2 - (rt$e1+rt$e2)) / ((rt$r1+rt$r2)^2 * (rt$r1+rt$r2-1)), na.rm = TRUE)
    z = u / sqrt(v_u)
    z
  }
  
  list(z_lr = get_z(w_lr, rt),
       z_fh_0_1 = get_z(w_fh_0_1, rt),
       z_fh_0_half = get_z(w_fh_0_half, rt),
       z_mw_12 = get_z(w_mw_12, rt),
       z_mw_24 = get_z(w_mw_24, rt),
       z_landmark_21 = -landmark(df_final, 21),
       z_landmark_27 = -landmark(df_final, 27))
  
}
###########################################################
###########################################################
### SCENARIOS
#########################################################
recruitment = list(n_0 = 500, n_1 = 500, r_period = 12, k = 1)
##########################################################
##########################################################
model_nph = list(change_points = c(6),
                 lambdas_0 = c(log(2) / 15, log(2) / 15),
                 lambdas_1 = c(log(2) / 15, log(2) / 21))



model_weak_null = list(change_points = c(6),
                       lambdas_0 = c(log(2) / 15, log(2) / 15),
                       lambdas_1 = c(log(2) / 15, log(2) / 15))



model_strong_null = list(change_points = c(7, 27),
                         lambdas_0 = c(log(2) / 15, log(2) / 15, log(2) / 25),
                         lambdas_1 = c(log(2) / 11, log(2) / 17, log(2) / 25))


model_ph = list(change_points = c(6),
                lambdas_0 = c(log(2) / 15, log(2) / 15),
                lambdas_1 = c(log(2) / 19, log(2) / 19))



model_diminishing = list(change_points = c(9, 18),
                         lambdas_0 = c(log(2) / 15, log(2) / 15, log(2) / 15),
                         lambdas_1 = c(log(2) / 25, log(2) / 18, log(2) / 13))

#########################################################
#########################################################
## RESULTS

set.seed(246)

results_nph <- purrr::map_df(rep(36, 10000),
                             sim_trial,
                             model = model_nph, 
                             recruitment = recruitment)

results_weak_null <- purrr::map_df(rep(36, 10000),
                                   sim_trial,
                                   model = model_weak_null, 
                                   recruitment = recruitment)

results_strong_null <- purrr::map_df(rep(36, 10000),
                                     sim_trial,
                                     model = model_strong_null, 
                                     recruitment = recruitment)

results_ph <- purrr::map_df(rep(36, 10000),
                            sim_trial,
                            model = model_ph, 
                            recruitment = recruitment)

results_diminishing <- purrr::map_df(rep(36, 10000),
                                     sim_trial,
                                     model = model_diminishing, 
                                     recruitment = recruitment)


#########################################################
#########################################################
## SUMMARIZE RESULTS


colMeans(results_nph < qnorm(0.025))
colMeans(results_weak_null < qnorm(0.025))
colMeans(results_strong_null < qnorm(0.025))
colMeans(results_ph < qnorm(0.025))
colMeans(results_diminishing < qnorm(0.025))

save(results_nph,
     results_weak_null,
     results_strong_null,
     results_ph,
     results_diminishing, file = "results_10000.RData")


#########################################################
## ADDITIONAL SCENARIO (WORSE EARLY)

model_strong_null_2 = list(change_points = c(2, 6),
                           lambdas_0 = c(log(2) / 14, log(2) / 10, log(2) / 15),
                           lambdas_1 = c(log(2) / 7, log(2) / 15, log(2) / 15))

set.seed(246)

results_strong_null_2 <- purrr::map_df(rep(36, 10000),
                                       sim_trial,
                                       model = model_strong_null_2, 
                                       recruitment = recruitment)


colMeans(results_strong_null_2 < qnorm(0.025))
