# Covid19

## Creating the Shiny app

If not already installed, do:

```
devtools::install_github("terminological/uk-covid-datatools", force = TRUE)
```

Then

```
runApp('Shiny_app.R')
```

creates the app, given the current time series of R0 (R0timeseries.RData).
