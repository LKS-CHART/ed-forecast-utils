

## ED Forecast Utility functions overview

This repo, which can be used for educational purposes, contains several functions which are used by the CHART team to explore forecast models
for our emergency department volumes. While the repo primarily contains functions
we have used to forecast ED volumes, we have included a couple of functions for
extracting weather data along with a calendar dataset we use for holidays. To understand how we calculate ED
volumes, you can see [this tutorial]

### Getting started

To use any of the functions, you can clone the repository and install the necessary 
libraries with the renv package

Loading the project for the first time:

- `renv::init()` +  select `Restore the project from the lockfile`


### Functions included

- `initialize_models` for initializing the models we typically explore
- `get_features` for extracting the variables we use in our models
- `pick_fourier_terms` for exploring different fourier terms in ARIMA models
- `fit_gam` for fitting our generalized additive model
- `forecast_models` the function we have used for extracting forecasts and 
95% confidence intervals
- `get_raw_weather_forecast` for extracting weather forecasts
- `clean_raw_forecast` for cleaning raw weather forecasts
- `extract_sport_events` for extracting sport event times (relies on an API (stattleshipR))

### Data included

- `data/calendar_variables.csv`: a dataset of all holidays
- `data/future_frame.csv`: a dataset with horizons we will forecast in the future. This is for illustrative purposes. 
While we can't include our actual data in the repo, this is the format that our data takes with 1 row per 6 hour interval. 


