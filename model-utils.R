#' Initialize forecast models
#' 
#' A function to create forecast models for each metric of interest
#' 
#' @param data a tibble containing 1 row per 6 hour interval with variables for
#' each of the metrics to be forecast
#' @param a vector containing the names of the metrics that the models will be built on
#' @param model_names a vector containing the names of the models to be tried
#' @param fourier_terms terms a list containing the results of `pick_fourier_terms()`
#' 
#' @return a list is returned with each component holding the forecast models
#'  for one of the metrics we forecast. 

library(dplyr)
library(lubridate)
library(stattleshipR)
library(forecast)
library(mgcv)
library(glue)
initialize_models <- function(data,
                              forecast_metrics = c("ctas13_arrivals", "ctas45_arrivals",
                                                   "mhesa_arrival", "admit_arrival",
                                                   "total_arrivals", "total_census", 
                                                   "ctas13_census", "ctas45_census",
                                                   "mhesa_census", "pia_total", 
                                                   "pia_ctas13", "pia_ctas45", 
                                                   "acute_total"), 
                              model_names = c("stlf", "arima",
                                              "tbats", "gam",
                                              "nnet", "theta",
                                              "snaive", "hw"),
                              fourier_terms) {
  
  # initialize model list
  all_models <- list()
  
  # for each forecast metric
  for(i in forecast_metrics) {

    # extract the features
    xreg_train <- get_features(tmp)
    
    # create a timeseries object
    ts_df <- msts(tmp[[i]], c(4,28))
    
    all_models[[i]] <- list()
    
    # fit each of the model types and store in the list
    for(j in model_names) {
      
      if(j == 'stlf') {
        
        model <- forecast::stlf(ts_df)
        
      } else if(j == "tbats") {
        
        model <- forecast::tbats(ts_df)
        
      } else if(j == "arima") {
        
        ft <- fourier_terms[[i]] %>% 
          dplyr::arrange(aicc) %>% 
          dplyr::slice(1)
        
        fourier_i <- ft$i
        fourier_j <- ft$j
        
        model <- arima(ts_df, xreg = cbind(xreg_train,
                                           fourier(ts_df, K = c(fourier_i,
                                                                fourier_j))))
      } else if(j == "gam") {
        
        model <- fit_gam(i, tmp)
        
      } else if(j == "nnet") {
        model <- forecast::nnetar(ts_df, xreg = xreg_train,
                                  decay = .8, MaxNWts= 3800)
      } else if(j == "hw") {
        
        model <- HoltWinters(ts_df)
        
      } else if(j == "snaive") {
        
        model <- snaive(ts_df)
        
      } else if(j == "theta") {
        
        model <- thetaf(ts_df)
        
      }
      
      
      all_models[[i]][[j]][["model"]] <- model
      all_models[[i]][[j]][["last_train"]] <- Sys.Date()
    }
    
  }
  
  return(all_models)
  
}


#' Get Features
#' 
#' A function that extracts the features we will use in the forecasting model
#' 
#' @param data tibble containing the variables we will potentially use
#' @param features a string vector with the names of the features we will extract
#' from the full set of data
#' 
get_features <- function(data,
                         features = c("holiday",
                                      "day_before_holiday",
                                      "day_after_holiday",
                                      "forecast_temp",
                                      "forecast_precip",
                                      "forecast_snow",
                                      "lag_heavy_snow",
                                      "city_event1",
                                      "city_event2",
                                      "sport_event1",
                                      "december_holiday",
                                      "new_year",
                                      "jan",
                                      "feb", 
                                      "mar",
                                      "apr",
                                      "may",
                                      "jun",
                                      "jul",
                                      "aug",
                                      "sep",
                                      "oct",
                                      "nov")) {
  
  xreg <- data %>% 
    dplyr::select_(features) %>% 
    as.matrix()
  return(xreg)
}


#' Function to fit our chosen GAM model
#' 
#' @param data a tibble containing the timeseries to be forecast an each
#' of the features sepecified in `get_features()``

fit_gam <- function( data, dist = "poisson") {
  
 
  
  gam_model <- gam(data ~ 
                     s(day_week, bs = "ps", k = 7, by = factor(hour)) +
                     s(hour, bs = "ps", k = 4, by = factor(weekend)) +
                     s(week_num, bs = "ps", k = 53) +
                     s(time, bs = "cs", k = 4) +
                     s(forecast_temp, bs = "cs", k = 20) +
                     forecast_precip +
                     forecast_snow +
                     forecast_snow*weekend + 
                     lag_heavy_snow + 
                     weekend + 
                     holiday +
                     day_after_holiday +
                     day_before_holiday +
                     december_holiday + 
                     new_year +
                     city_event1 +
                     city_event2 +
                   data = data,
                   select = T,
                   family = dist)
  
  return(gam_model)
  
}

#' Forecast models
#' 
#' This function will forecast each of the models over a specificed horizon
#' 
#' @param data is a tibble with each of the metrics we forecast
#' @param models is an object created with `initialize_models()`
#' @param forecast_metrics a vector with the metrics to be forecasted
#' @param model_names the names of the models contained in the models parameter
#' @param fourier_terms an object returned from `pick_fourier_terms()`
#' @param data_path the path to the data extracted 
#' @param weather_path the path to the weather forecasts
#' 
#' @return A list is returned with the forecasts and 95% confidence intervals for
#' each of the metrics and models

forecast_models <- function(data, horizon = 28,  models,
                            forecast_metrics = c("ctas13_arrivals", "ctas45_arrivals",
                                                 "mhesa_arrival", "admit_arrival",
                                                 "total_arrivals", "total_census", 
                                                 "ctas13_census", "ctas45_census",
                                                 "mhesa_census", "pia_total", 
                                                 "pia_ctas13", "pia_ctas45", 
                                                 "acute_total"), 
                            model_names = c("stlf", "arima",
                                            "tbats", "gam",
                                            "nnet", "theta",
                                            "snaive", "hw"),
                            fourier_terms,
                            data_path,
                            weather_path = "data/weather/clean_weather_forecast/") {
  
  
  
  load(glue::glue("{data_path}future_frame.Rda"))
  
  next_seven_days <- as.character(Sys.Date() + ddays(0:6))
  
  #  read in the weather forecasts ------------------------------------------
  
  fs <- list.files(weather_path, pattern = "*.csv", full.names = T)
  fs <- fs[length(fs)]
  
  weather_forecast<- readr::read_csv(i) %>% 
    mutate(file = i)
  
  weather_forecast <- weather_forecast %>% 
    filter(date %in% ymd(next_seven_days))

  weather_forecast <- prepare_weather_covs(clean_weather_forecast =weather_forecast )
  
  covs <- weather_forecast %>% 
    left_join(future_frame, by = c('date' = 'day_date'))
  
  all_forecasts <- list()
  
  for(i in forecast_metrics) {
    
    tmp <- data %>% 
      filter(!is.na(data[[i]]))
    
    xreg_train <- get_features(tmp)
    xreg_test <- get_features(covs)
    
    ts_df <- msts(tmp[[i]], c(4,28))
    
    
    all_forecasts[[i]] <- list()
    
    for(j in model_names) {
      
      if(j == 'stlf') {
        
        stream <- stlf(ts_df)
        model_mean <- round(as.numeric(forecast(stream, h = horizon)$mean))
        model_low <- round(as.numeric(forecast(stream, h = horizon)$lower[, 2]))
        model_high <- round(as.numeric(forecast(stream, h = horizon)$upper[, 2]))
      } else if(j == "tbats") {
        stream <- tbats(ts_df, model = models[[i]][[j]]$model)
        model_mean <- round(as.numeric(forecast(stream, h = horizon)$mean))
        model_low <- round(as.numeric(forecast(stream, h = horizon)$lower[, 2]))
        model_high <- round(as.numeric(forecast(stream, h = horizon)$upper[, 2]))
      } else if(j == "arima") {
        
        ft <- fourier_terms[[i]] %>% 
          dplyr::arrange(aicc) %>% 
          dplyr::slice(1)
        fourier_i <- ft$i
        fourier_j <- ft$j
        
        stream <- Arima(ts_df, model = models[[i]][[j]]$model,
                        xreg =cbind(xreg_train,
                                    fourier(ts_df, K = c(fourier_i,
                                                         fourier_j))))
        arima_test <- cbind(xreg_test, fourier(ts_df, K = c(fourier_i,
                                                            fourier_j), h = horizon))
        
        model_mean <- round(as.numeric(forecast(stream, xreg = arima_test, h = horizon)$mean))
        model_low <- round(as.numeric(forecast(stream, xreg = arima_test, h = horizon)$lower[, 2]))
        model_high <- round(as.numeric(forecast(stream, xreg = arima_test, h = horizon)$upper[, 2]))
        
      } else if(j == "gam") {
        
        model_mean <- predict_gam(models[[i]][[j]]$model, data = covs)
        model_low <- rep(0, horizon)
        model_high <- rep(0, horizon)
        
      } else if(j == "nnet") {
        stream <- nnetar(ts_df, model = models[[i]][[j]]$model,
                         xreg =cbind(xreg_train))
        model_mean <- round(as.numeric(forecast(stream, xreg = xreg_test, h = horizon)$mean))
        model_low <- rep(0, horizon)
        model_high <- rep(0, horizon)
        
      } else if(j == "hw") {
        
        stream <- HoltWinters(ts_df)
        model_mean <- round(as.numeric(forecast(stream,  h = horizon)$mean))
        model_low <- round(as.numeric(forecast(stream,  h = horizon)$lower[, 2]))
        model_high <- round(as.numeric(forecast(stream,  h = horizon)$upper[, 2]))
      } else if(j == "snaive") {
        
        stream <- snaive(ts_df)
        model_mean <- round(as.numeric(forecast(stream,  h = horizon)$mean))
        model_low <- round(as.numeric(forecast(stream,  h = horizon)$lower[, 2]))
        model_high <- round(as.numeric(forecast(stream,  h = horizon)$upper[, 2]))
        
      } else if(j == "theta") {
        
        stream <- thetaf(ts_df)
        model_mean <- round(as.numeric(forecast(stream,  h = horizon)$mean))
        model_low <- round(as.numeric(forecast(stream,  h = horizon)$lower[, 2]))
        model_high <- round(as.numeric(forecast(stream,  h = horizon)$upper[, 2]))
      }
      
      forecast_tbl <- tibble(series = i,
                             model = j,
                             day_date = covs$date,
                             index_date = covs$index_date,
                             forecast = model_mean,
                             low = model_low,
                             high = model_high)
      
      all_forecasts[[i]][[j]][["forecast"]] <- forecast_tbl
      all_forecasts[[i]][[j]][["forecast_date"]] <- Sys.Date()
      
    }
    
  }
  
  return(all_forecasts)
  
}


#' Extract weather forecasts
#' 
#' We extract weather forecasts from theweathernetwork.com 
#' 
#' @param station_id the station id to extract the forecasts
#' 
#' @return a list containing all of the forecast information found on theweathernetwork.com
#' 
#' 
get_raw_weather_forecast <- function(station_id = 'caon0696') {
  raw_weather_forecast <- jsonlite::fromJSON(glue::glue("https://www.theweathernetwork.com/api/data/{station_id}/cm?ts=1932"))
  return(raw_weather_forecast) 
}


#' Function to clean up the weather forecast data
#' 
#' a simple function to clean the list of weather forecasts
#' 
#' @param raw_weather_forecast the object returned from `get_raw_weather_forecast()`
clean_raw_forecast <- function(raw_weather_forecast) {
  forecast_data <- raw_weather_forecast$sevendays$`periods` %>% 
    tibble::as_tibble()
  
  clean_weather_forecast <- forecast_data %>% 
    mutate(date = ymd(paste0(year(Sys.Date()),"/", super_short_date))) %>% 
    select(date, rain = orig_rain, snow = orig_snow, 
           min_temp = imperial_temperatureMin,
           max_temp = imperial_temperatureMax,
           feels_like = metric_feelsLike) %>% 
    mutate(min_temp = convert_fahrenheit_to_celcius(as.numeric(min_temp)),
           max_temp = convert_fahrenheit_to_celcius(as.numeric(max_temp)),
           rain = as.numeric(rain),
           snow = as.numeric(snow))
  return(clean_weather_forecast)
}


#' function to extract sports data
#' 
#' library(stattleshipR)
set_token(Sys.getenv("stattle_token"))


# basketball data ---------------------------------------------------------
#' Function to scrape NBA game times
#' 
#' We used the stattleship API to extract NBA/NHL game times for historical 
#' exploration. It requires a contribution to get a token
#' 
#' More can be found here https://github.com/stattleship/stattleship-r
#' 
#' @param sport the sport (either basketball or hockey)
#' @param league either nba or nhl
#' @param ep the end point to extract. We always get games
#' @param team_id either nba-tor or nhl-tor
#' @param season_list the year to extract game times

extract_sport_events <- function(sport = 'basketball',
                                 league = 'nba',
                                 ep = 'games',
                                 team_id = 'nba-tor',
                                 season_list = c('nba-2015-2016', 
                                                 'nba-2016-2017',
                                                 'nba-2017-2018', 
                                                 'nba-2018-2019')) {
  
  totals <- list()
  for(i in season_list) {
    reg_season_q_body <- list(team_id= team_id, 
                              status='ended',
                              interval_type='regularseason',
                              season_id = i)
    pre_season_q_body <- list(team_id='nba-tor',
                              status='ended',
                              interval_type='preseason',
                              season_id = i)
    
    reg_season_results <- ss_get_result(sport=sport, league=league, ep=ep, query=reg_season_q_body, walk=TRUE)  
    reg_season_game_logs <- do.call('rbind', lapply(reg_season_results, function(x) x$games))  
    
    pre_season_results <- ss_get_result(sport=sport, league=league, ep=ep, query=pre_season_q_body, walk=TRUE)  
    pre_season_game_logs <- do.call('rbind', lapply(pre_season_results, function(x) x$games))  
    
    # Toronto home  venue ID is 5b32d0b9-7e6f-4751-9f8b-c4f7b2feefa3
    regular_season_game_logs <- reg_season_game_logs %>% 
      filter(venue_id == '5b32d0b9-7e6f-4751-9f8b-c4f7b2feefa3')
    pre_season_game_logs <- pre_season_game_logs %>% 
      filter(venue_id == '5b32d0b9-7e6f-4751-9f8b-c4f7b2feefa3')
    
    regular_season_game_logs <- regular_season_game_logs %>% 
      select(started_at, ended_at, duration, attendance)
    
    regular_season_game_logs$started_at <- gsub(pattern = "T", replace = " ", regular_season_game_logs$started_at)
    regular_season_game_logs$started_at <- substr(regular_season_game_logs$started_at, 1, 19)
    regular_season_game_logs$started_at <- ymd_hms(regular_season_game_logs$started_at)
    
    regular_season_game_logs$ended_at <- gsub(pattern = "T", replace = " ", regular_season_game_logs$ended_at)
    regular_season_game_logs$ended_at <- substr(regular_season_game_logs$ended_at, 1, 19)
    regular_season_game_logs$ended_at <- ymd_hms(regular_season_game_logs$ended_at)
    
    
    pre_season_game_logs <- pre_season_game_logs %>% 
      select(started_at, ended_at, duration, attendance)
    
    pre_season_game_logs$started_at <- gsub(pattern = "T", replace = " ", pre_season_game_logs$started_at)
    pre_season_game_logs$started_at <- substr(pre_season_game_logs$started_at, 1, 19)
    pre_season_game_logs$started_at <- ymd_hms(pre_season_game_logs$started_at)
    
    pre_season_game_logs$ended_at <- gsub(pattern = "T", replace = " ", pre_season_game_logs$ended_at)
    pre_season_game_logs$ended_at <- substr(pre_season_game_logs$ended_at, 1, 19)
    pre_season_game_logs$ended_at <- ymd_hms(pre_season_game_logs$ended_at)
    df <- rbind(pre_season_game_logs, regular_season_game_logs)
    
    totals[[i]] <- df
    
  }
  
  game_totals <- do.call(rbind, totals)
  
  return(game_totals)
}


