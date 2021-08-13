library(tidyverse)
library(httr)
library(rvest)
library(devtools)
library(tidyverse)
library(plotly)
library(shiny)
library(shinythemes)
library(quantreg)
library(mgcv)
library(lmf)
library(dplyr)
library(shinycssloaders)
library(shinyWidgets)

#-------FUNCTIONS TO SIMULATE DATA---------------------------------------------------------------------------------------------------------------------
days_of_lactation <- 305


# Helper Functions
# INPUT:
# content of param_wilmink list
#  intercept: numeric;
#  coef_day: numeric;
#  coef_expterm: numeric;
#
# OUTPUT:
# numeric R vector length = days_of_lactation; representing Wood's Milk Curve over days_of_lactation days

signal <- function(param_wilmink = list(intercept, coef_day, coef_expterm, k), days_of_lactation = 305) {
  days <- seq(1, days_of_lactation, by = 1)
  wilmink_curve <- param_wilmink$intercept - param_wilmink$coef_day * days + param_wilmink$coef_expterm * exp(-1 * param_wilmink$k * days)
  return(wilmink_curve)
}


# INPUT ARGUMENTS:
# model: one of following strings: "WN", standing for White Noise ,"AR1", "AR(p)", an AR model of order p, "Gamma"
# param: list of necessary parameters corresponding to the model, see table below for details.
# model | param
# ---------------
# WN    | c(variance of innovations); numeric vector
# AR1   | c(variance_of_innovations, autoregression_coefficient); numeric vector
# AR(p) | c(variance_of_innovations, vector of p autoregression coefficients, order); list, with second component being a numeric vector
#
# OUTPUT:
# numeric vector of length days_of_lactation
variability_in_milk_yield <- function(VarModel, param_var = NULL, days_of_lactation = 305) {

  # check that VarModel input is appropriate
  if ((VarModel != "WN") & (VarModel != "AR1") & (VarModel != "AR(p)") & (VarModel != "Gamma")) stop("VarModel must be equal to one of: 'WN', 'AR1', 'AR(p)','Gamma'")

  if (VarModel == "WN") {
    if (is.null(param_var[1])) stop("User must input a value in the first position of the param_var vector: variance of innovations")

    noise <- arima.sim(model = list(order = c(0, 0, 0)), n = days_of_lactation, sd = sqrt(param_var[1]))
  }

  if (VarModel == "AR1") {
    if (is.null(param_var[1])) stop("User must input a value in the first position of the param_var vector: variance of innovations")
    if (is.null(param_var[2])) stop("User must input a value for autoregression coefficient the second component of the vector in the of the param_var argument")

    noise <- arima.sim(model = list(order = c(1, 0, 0), ar = param_var[2]), n = days_of_lactation, sd = sqrt(param_var[1]))
  }

  if (VarModel == "AR(p)") {
    if (is.null(param_var[[1]])) stop("User must input a value in the first position of the param_var vector: variance of innovations")
    if (length(param_var[[2]]) != param_var[[3]]) stop("User must provide vector containing p autoregression coefficients in second component of param_var argument")
    if (is.null(param_var[[3]])) stop("User must input a numeric value for autoregression param_vector (order of the AR model)")

    noise <- arima.sim(model = list(order = c(param_var[[3]], 0, 0), ar = param_var[[2]]), n = days_of_lactation, sd = sqrt(param_var[[1]]))
  }

  if (VarModel == "Gamma") {
    # param_var = c(shape, scale)

    noise <- -1 * rgamma(days_of_lactation, shape = param_var[1], scale = param_var[2])
  }

  return(noise)
}


# INPUT:
# 1st and 2nd elements of param_kickoff list
# p_ko: decimal number; rate of kickoffs as a percentage
# v_k: numeric vector of length 2 representing lower and upper bound of uniform distribution; represents range of values for possible proportion
#      of the milk that a cow would have given on a kick off day, that was collected before the milking claw was prematurely removed
# 1st and 2nd elements of param_sick list
# p_sick_each_day_of_lactation: numeric vector; probability that a cow would become ill on each day of the days_of_lactation day lactation period
# v_s: numeric vector of length 2 representing lower and upper bound of uniform distribution; range of possible values for length of illness

# OUTPUT:
# numeric vector to adjust for kickoffs and illness, length=days_of_lactation
# p_ko = 0, v_k = c(0, 2 / 3), p_sick_each_day_of_lactation = rep(1 / days_of_lactation, days_of_lactation), v_s = c(0, 0)
outliers <- function(param_kickoff = list(p_ko = 0, v_k = c(0, 2 / 3)),
                     param_sick = list(p_sick_each_day_of_lactation = rep(1 / days_of_lactation, days_of_lactation), v_s = c(0, 0)),
                     days_of_lactation = 305) {

  # extract parameters
  p_ko <- param_kickoff$p_ko
  v_k <- param_kickoff$v_k

  p_sick_each_day_of_lactation <- param_sick$p_sick_each_day_of_lactation
  v_s <- param_sick$v_s

  # kickoffs
  if (v_k[1] < 0 | v_k[2] > 1) stop("Bounds of v_k must be within [0,1]")

  total_kickoffs <- ceiling(days_of_lactation * p_ko)
  outlier_adj <- rep(1, days_of_lactation)
  kickoff_indices <- sample(1:days_of_lactation, total_kickoffs, replace = FALSE)
  outlier_adj[kickoff_indices] <- runif(total_kickoffs, v_k[1], v_k[2])

  # sick days
  if (sum(p_sick_each_day_of_lactation) != 1) stop("Sum of p_sick_each_day_of_lactation must sum to 1.")
  if (v_s[1] < 0 | v_s[2] > days_of_lactation) stop("Bounds of v_s must be within [0, days_of_lactation]")

  n_s <- runif(1, v_s[1], v_s[2])

  # ...if there are sick days
  if (n_s > 0) {
    d_s <- sample(1:days_of_lactation, 1, prob = p_sick_each_day_of_lactation)

    illness_beginning <- c(d_s, ceiling(d_s + n_s / 3))
    illness_end <- c(ceiling(d_s + n_s / 3) + 1, d_s + n_s)

    outlier_adj[illness_beginning] <- runif(length(illness_beginning), 0, 0.5)
    outlier_adj[illness_end] <- runif(length(illness_end), 0, 1)

    return(outlier_adj[1:days_of_lactation])
  }

  # ...if there are no sick days
  else {
    return(outlier_adj)
  }
}
# MAIN FUNCTION

# INPUT:
# param_wilmink: list, with below contents...
#  intercept: numeric;
#  coef_day: numeric;
#  coef_expterm: numeric;
# k
# VarModel: one of the following strings 'WN', 'AR1', 'AR(p)','Gamma'; Model to use when simulating noise
# param_var: numeric vector of parameters, contents depend on model used for variablity. see variability_in_milk_yield() for details
# param_kickoff: list, with below contents...
#   p_ko: decimal number; rate of kickoffs as a percentage
#   v_k: numeric vector of length 2 representing lower and upper bound of uniform distribution; represents range of values for possible proportion
#        of the milk that a cow would have given on a kick off day, that was collected before the milking claw was prematurely removed
# param_sick: list, with below contents...
#   p_sick_each_day_of_lactation: numeric vector; probability that a cow would become ill on each day of the days_of_lactation day lactation period
#   v_s: numeric vector of length 2 representing lower and upper bound of uniform distribution; range of possible values for length of illness
milkCow <- function(param_wilmink, VarModel = "WN", param_var = NULL, param_kickoff = list(p_ko = 0, v_k = c(0, 2 / 3)), days_of_lactation = 305,
                    param_sick = list(p_sick_each_day_of_lactation = rep(1 / days_of_lactation, days_of_lactation), v_s = c(0, 0))) {
  Signal <- signal(param_wilmink)

  Noise <- variability_in_milk_yield(VarModel = VarModel, param = param_var) # whichever of the variability functions we choose

  Outliers <- outliers(param_kickoff = param_kickoff, param_sick = param_sick)

  Curve <- (Signal + Noise) * Outliers

  ModelSpec <- list(
    "wilmink_curve_parameters" = param_wilmink,
    "VarModel" = VarModel,
    "variance_model_parameters" = param_var,
    "kickoff_outlier_parameters" = param_kickoff,
    "sick_outlier_parameters" = param_sick
  )

  Curve_Noise_Outliers_Signal_ModelSpec <- list(curve = Curve, noise = Noise, outliers = Outliers, signal = Signal, modelspec = ModelSpec)

  return(Curve_Noise_Outliers_Signal_ModelSpec)
}




# METHODS

# INPUT: milkCow() output
# OUTPUT: time series plot of simulated days_of_lactation day lactation period
plot_curve <- function(object) {
  # library(plotly)
  # library(dplyr)
  days <- seq(1, days_of_lactation, by = 1)
  milk_yield <- object$curve
  plot_ly(
    x = ~days, y = ~milk_yield, mode = "lines", type = "scatter", text = paste("Day:", days, "\n", "Milk Yield:", format(round(milk_yield, 2), nsmall = 2), "kg")
  ) %>%
    layout(
      title = paste("Daily Milk Yield Over", days_of_lactation, "Day Lactation Period"),
      xaxis = list(title = "Day of Lactation Cyle", range = c(0, days_of_lactation)),
      yaxis = list(title = "Daily Milk Yield (kg/day)", range = c(0, ceiling(max(milk_yield))))
    )
}

# INPUT: milkCow() output
# OUTPUT: time series plot of signal (Wood's Curve)
plot_signal <- function(object) {
  # library(plotly)
  # library(dplyr)
  days <- seq(1, days_of_lactation, by = 1)
  signal <- object$signal
  plot_ly(
    x = ~days, y = ~signal, mode = "lines", type = "scatter",
    text = paste("Day:", days, "\n", "Milk Yield:", format(round(signal, 2), nsmall = 2), "kg")
  ) %>%
    layout(
      title = "Daily Milk Yield Over days_of_lactation Day Lactation Period (Wood's Curve)",
      xaxis = list(title = "Day of Lactation Cyle", range = c(0, days_of_lactation)),
      yaxis = list(title = "Daily Milk Yield (kg/day)", range = c(0, ceiling(max(signal))))
    )
}




#--FUNCTION TO FIT MODEL USING RANSAC----------------------------------------------------------------------------
# input parameters


# used: https://www.r-bloggers.com/separating-the-signal-from-the-noise-robust-statistics-for-pedestrians/ as guide

# INPUT:
# data : a data frame containing milk yield as the response, specific to model
# model: one of "Wilmink", "Quantile Reg", "Spline"
# n: minimum number of data points required to fit the model
# k: maximum number of iterations allowed in the algorithm
# t: threshold value to determine when a data point fits a model; absolute value of acceptable residual from points to model
# d: number of close data points required to assert that a model fits well to data; acceptable number of inliers

# OUTPUT:
# a model
# a model
ransac <- function(model, data, n, k, t, j) {
  iterations <- 0
  bestfit <- NULL
  besterr <- 1e5

  while (iterations < k) {
    maybeinliers <- sample(nrow(data), n)

    if (model == "Woods") {
      maybemodel <- lm(log(yield) ~ log(days) + days, data = data, subset = maybeinliers) # use days_yield_reactive data table
      # a <- maybemodel$coefficients[1] %>% exp() # would Bo actually be log(a)>
      # b <- maybemodel$coefficients[2]
      # c <- -1 * maybemodel$coefficients[3] # note: 3rd coef was already negative, so use of negative sign here and then in the woods formula is redundant, using for the sake of using exact formula from references
    }

    if (model == "Wilmink") {
      maybemodel <- lm(yield ~ d + exp_term, data = data, subset = maybeinliers)
    }
    if (model == "Quantile Reg") {
      maybemodel <- rq(yield ~ poly(days, 4, raw = TRUE), tau = 0.7, data = data, subset = maybeinliers)
    }
    if (model == "Spline") {
      maybemodel <- gam(yield ~ s(days), data = data, subset = maybeinliers)
    }

    alsoinliers <- NULL

    for (rowindex in setdiff(1:nrow(data), maybeinliers)) {
      predicted <- maybemodel %>% predict(newdata = data[rowindex, ])

      if (model == "Woods") {
        predicted <- exp(predicted) # since we estimated log(yield)
      }

      actual <- data$yield[rowindex]

      if (abs(predicted - actual) < t) {
        alsoinliers <- c(alsoinliers, rowindex)
      }
    }

    # if (is.null(alsoinliers)){next}

    if (length(alsoinliers) > j) {
      if (model == "Woods") {
        bettermodel <- lm(log(yield) ~ log(days) + days, data = data, subset = c(maybeinliers, alsoinliers)) # output: exp(Bo) = a, B1= b, b2 = -c of woods
        # a <- bettermodel$coefficients[1] %>% exp() # would Bo actually be log(a)>
        # b <- bettermodel$coefficients[2]
        # c <- -1 * bettermodel$coefficients[3] # note: 3rd coef was already negative, so use of negative sign here and then in the woods formula is redundant, using for the sake of using exact formula from references
        # days <- c(maybeinliers, alsoinliers)
        # #
        predicted <- bettermodel %>%
          predict() %>%
          exp()
        actual <- data[c(maybeinliers, alsoinliers), ]$yield
        thiserr <- (predicted - actual) %>% sd()

        # thiserr <- summary(bettermodel)$sigma
      }

      if (model == "Wilmink") {
        bettermodel <- lm(yield ~ d + exp_term, data = data, subset = c(maybeinliers, alsoinliers))
        thiserr <- summary(bettermodel)$sigma
      }
      if (model == "Quantile Reg") {
        bettermodel <- rq(yield ~ poly(days, 4, raw = TRUE), tau = 0.7, data = data, subset = c(maybeinliers, alsoinliers))
        thiserr <- bettermodel$residuals %>% sd()
      }

      if (model == "Spline") {
        bettermodel <- gam(yield ~ s(days), data = data, subset = c(maybeinliers, alsoinliers))
        thiserr <- bettermodel$residuals %>% sd()
      }

      if (thiserr < besterr) {
        bestfit <- bettermodel
        besterr <- thiserr
      }
    }

    iterations <- iterations + 1
  }

  # if (model == "Woods") {
  #   #
  #   a <- bestfit$coefficients[1] %>% exp()
  #   b <- bestfit$coefficients[2]
  #   c <- -1 * bestfit$coefficients[3] # note: 3rd coef was already negative, so use of negative sign here and then in the woods formula is redundant, using for the sake of using exact formula from references
  #   days <- data$days
  #   bestfit <- a * (days^b) * exp(-c * days)
  # }
  # #

  return(bestfit)
}
#---BEGIN SHINY APP CODE -----------------------------------------------------------------------------------------------------------_____


ui <- fluidPage(
  theme = shinytheme("yeti"),

  # App Title
  titlePanel(h1("Daily Milk Yield", align = "center")),
  navbarPage(
    "Menu",
    tabPanel(
      "Plots",

      # Sidebar layout with input and output definitions
      sidebarLayout(
        position = "left",

        # Sidebar panel for inputs
        sidebarPanel(
         actionButton(inputId = "submit_input", label = "Make Plots!"),
          
          fluidRow(
          column(width = 2, materialSwitch(inputId = "RANSAC", label = "RANSAC", status = "warning")), 
          column(width = 8, offset = 2, htmlOutput("TmYmP"))
          ),
         
          
          htmlOutput("param_label"),

          #numericInput("intercept", h6("a: scale of production"), value = 30),
           sliderInput(
             inputId = "intercept",
             label = "a: scale of production",
             min = 15, max = 45,
             value =30, step = 1
           ),


          #numericInput("coef_expterm", h6("growth to peak"), value = -20),
           sliderInput(
             inputId = "coef_expterm",
             label = "b: growth to peak",
             min = -40, max = -10,
             value = -20, step = 1
           ),


          #numericInput("coef_day", h6("c: persistency"), value = -3.744),
           sliderInput("coef_day", "c: persistency",
             min = 0.01, max = .1,
             value = 0.05, step = 0.0025
           ),
          
          #numericInput("k", h6("k: day of peak"), value = 0.05)
          sliderInput("k", "k: day of peak",
                      min = 0.04, max = 0.07,
                      value = 0.05, step = 0.0025
          )
        ),




        # Main panel for displaying outputs
        mainPanel(
          tabsetPanel(
            # put it back here
            tabPanel(
              "Daily Yield",
              plotlyOutput("daily_yield"), # %>%  withSpinner(color="#0dc5c1"),
              fluidRow(htmlOutput("clicklegend"))
            ),

            tabPanel(
              "Histograms of Residuals",
              plotlyOutput("res_hist"), # %>% withSpinner(color="#0dc5c1"),
              fluidRow(htmlOutput("clicklegend2"))
            )



            # tabPanel()
          )
        )
      )
    ),

    tabPanel(
      "Optional Parameters",
      titlePanel(paste("Below are optional parameter values for the data simulation, with the default values displayed", sep = "\n") %>% h6()),
      #textOutput("Below are optional parameter values for the data simulation, with the default values displayed"),
   
      
      fluidRow(

        # Kickoff Parameters row 1 col 5
        column(
          h5("Kickoff Parameters") %>% strong(),
          # offset = 1,
          width = 3,

          sliderInput(
            inputId = "p_ko",
            label = "Rate that cow kicks off her milking claw:",
            min = 0, max = 1,
            value = 0.08, step = 0.01
          ),

          sliderInput(
            inputId = "v_k",
            label = "Range, of proportion of milk, that would have already been collected before milking claw was kicked off: ",
            min = 0, max = 1,
            value = c(0, 2 / 3), step = 0.01
          )
        ),

        # Sick Parameters, row 2 col 1
        column(
          h5("Sick Parameters")%>%strong(),
          offset = 1,
          width = 3,
          sliderInput(
            inputId = "risk_period",
            label = "Risk Period: Range of days cow is more likely to be sick ",
            min = 0, max = 305,
            value = c(0, 0), step = 1
          ),

          sliderInput(
            inputId = "p_sick",
            label = "Probability cow is sick during Risk Period (see above): ",
            min = 0, max = 1,
            value = 1 / days_of_lactation, step = .001
          ),

          sliderInput(
            inputId = "v_s",
            label = "Range of possible values for length of illness: ",
            min = 0, max = 20,
            value = c(0, 5), step = 1
          )
        ),

        # VarModel choices, row 1 col 8
        column(
          h5("Noise Parameters") %>% strong(),
          width = 4,
          offset = 1,
          radioButtons(
            inputId = "VarModel", label = "Choose one of the following models to simulate noise in the data: ",
            choices = c("WN", "AR1", "AR(p)", "Gamma"),
            selected = "WN",
            inline = TRUE
          ),


          conditionalPanel(
            condition = "input.VarModel == 'WN'",
            numericInput("var_innovations",
              h6("Variance of Innovations (must be greater than or equal to 0)"),
              value = 2
            )
          ),

          conditionalPanel(
            condition = "input.VarModel == 'AR1'",
            numericInput("var_innovations",
              h6("Variance of Innovations (must be greater than or equal to 0)"),
              value = 2
            ),
            numericInput("auto_coef1",
              h6("Autoregression Coefficient"),
              value = -0.5
            )
          ),
          conditionalPanel(
            condition = "input.VarModel == 'AR(p)'",
            numericInput("var_innovations",
                         h6("Variance of Innovations (must be greater than or equal to 0)"),
                         min = 0,
                         value = 2
            ),
            
            numericInput("order_of_model",
                         h6("Order of the model (this is the 'p' in 'Ar(p)')"),
                         value = 2
            ),
            #gsub("[[:blank:]]", "", input$auto_coef_p)  %>% str_split(",", simplify = TRUE) %>% as.numeric()
            textInput("auto_coef_p",
                      h6("Enter the p autogregression coefficients, separated by a comma."))
            
          ),

          conditionalPanel(
            condition = "input.VarModel == 'Gamma'",
            numericInput("shape",
              h6("shape"),
              value = 2
            ),

            numericInput("scale",
              h6("scale"),
              value = 2
            )
          )
        )
      )
    )
  )
)



server <- function(input, output, session) {

  # PLOTS TAB

  # ...PLOT ON DAILY YIELD SUB-TAB and PLOTS ON HISTOGRAM SUBTAB
  observeEvent(
    eventExpr = input[["submit_input"]],

    handlerExpr = {
      showModal(modalDialog(h4("Generating Plot, this may take a few seconds."), footer = NULL))

      param_kickoff <- list(p_ko = input$p_ko, v_k = input$v_k)
      
     # rp <- c(input$risk_period[1]:input$risk_period[2])
      
      p_sick_each_day_of_lactation = rep(1 / days_of_lactation, days_of_lactation)
      
      # #HELP
      #  if(input$risk_period[2] != 0){
      #  #adjust probability sick for days in risk period
      #    p_sick_each_day_of_lactation[rp)] = input$p_sick  
      #  #decrease probability sick for days outside of risk period, so p_sick_each_day_of_lactation sums to 1
      #    p_sick_each_day_of_lactation[-c(rp)] = (1 - (input$p_sick)*length(rp))/(days_of_lactation - length(rp))
      #    }
      # 
      param_sick <- list(p_sick_each_day_of_lactation = p_sick_each_day_of_lactation, v_s = c(input$v_s[1], input$v_s[2]))
      
      #param var
      if(input$VarModel == "WN"){param_var <- c(input$var_innovations)}
      if(input$VarModel =="AR1"){param_var <- c(input$var_innovations, input$auto_coef1)}
      if(input$VarModel =="AR(p)"){
        auto_coef_p <- gsub("[[:blank:]]", "", input$auto_coef_p)  %>% str_split(",", simplify = TRUE) %>% as.numeric()
        param_var <- list(input$var_innovations, auto_coef_p, input$order_of_model)
        }
      if(input$VarModel == "Gamma"){param_var <- c(input$shape, input$scale)}
      
      
      p_wilmink <- list(intercept = input$intercept, coef_day = input$coef_day, coef_expterm = input$coef_expterm, k = input$k)

      mC_object <- milkCow(param_wilmink = p_wilmink, VarModel = input$VarModel, param_var = param_var, param_kickoff = param_kickoff, param_sick = param_sick)

      sig_reactive <- mC_object$signal

      yield <- mC_object$curve
      days <- seq(1, days_of_lactation, by = 1)

      days_yield_reactive <- data.frame("days" = days, "yield" = yield) # %>% SharedData$new()


      ransac_bool <- input$RANSAC

      if (ransac_bool == FALSE) {
        # coeff of lm(): exp(Bo) = a, B1= b, b2 = -c of woods
        a <- lm.extract(log(yield) ~ log(days) + days, data = days_yield_reactive)$ajt[1] %>% exp()
        b <- lm.extract(log(yield) ~ log(days) + days, data = days_yield_reactive)$ajt[2]
        c <- lm.extract(log(yield) ~ log(days) + days, data = days_yield_reactive)$ajt[3] * -1

        woods <- a * (days^b) * exp(-1 * c * days) # replaced t with days

        poly4_curve <- rq(yield ~ poly(days, 4, raw = TRUE), tau = 0.7, data = days_yield_reactive) %>% predict()

        spline <- gam(yield ~ s(days), data = days_yield_reactive) %>% predict()

        # wilmink_df <- data.frame(yield = yield, d = days, exp_term = exp(-0.05 * days))
        # wilmink <- lm(yield ~ d + exp_term, data = wilmink_df) %>% predict()
      }

      # function(model, data, n, k, t, j)


      if (ransac_bool == TRUE) {
        woods <- ransac(model = "Woods", data = days_yield_reactive, n = 10, k = 10, t = 10, j = 30) %>%
          predict(newdata = days_yield_reactive) %>%
          exp()

        poly4_curve <- ransac(model = "Quantile Reg", data = days_yield_reactive, n = 10, k = 10, t = 10, j = 30) %>%
          predict(newdata = days_yield_reactive)

        spline <- ransac(model = "Spline", data = days_yield_reactive, n = 10, k = 10, t = 10, j = 30) %>%
          predict(newdata = days_yield_reactive)
      }


      # ...Tm Ym P values on left bar  <- will update once reactive variables are made, otherwise will just be regenerating pts
      output$TmYmP <- renderUI({
        # extract peak yield and time of peak yield from generated data points
        t_m <- which(sig_reactive == max(sig_reactive))

        y_m <- max(sig_reactive) %>%
          round(3) %>%
          format(nsmall = 3)

        p <- (100*(max(sig_reactive) - sig_reactive[305]) / max(sig_reactive)) %>%
          round(4) %>%
          format(nsmall = 3)

        Tm <- paste("Time of Peak Yield (Day): ", t_m) %>% h6()
        Ym <- paste("Peak Yield (kg): ", y_m) %>% h6()
        P <- paste("Lactation Persistency: ", p) %>% h6()
        title <- ("Descriptive Statistics") %>%
          strong() %>%
          h6()
        # woods <- ("Wilmink's Parameters") %>%
        #   strong() %>%
        #   h6()

        HTML(paste(title, Tm,  Ym, P))
      })



      output$daily_yield <- renderPlotly({
        # signal
        plot_ly(
          x = ~days,
          y = ~sig_reactive,
          mode = "lines",
          type = "scatter",
          color = I("black"),
          name = "Generating Curve (Wilmink)",
          text = paste("Day:", days, "\n", "Milk Yield:", sig_reactive, nsmall = 2, "kg")
        ) %>%
          layout(
            title = paste("Daily Milk Yield Over", days_of_lactation, "Day Lactation Period"),
            xaxis = list(title = "Days In Milk (DIM)", range = c(0, 310)),
            yaxis = list(title = "Simulated Observed Yields (kg/day)", range = c(0, ceiling(max(sig_reactive)) + 10)),
            margin = list(t = 100),
            height = 500
          ) %>%
          # points
          add_trace(
            x = ~days,
            y = ~yield,
            mode = "markers",
            type = "scatter",
            opacity = 0.5,
            color = I("blue"),
            name = "Predicted Daily \n Milk Yield",
            text = paste("Day:", days, "\n", "Milk Yield:", format(round(yield, 2), nsmall = 2), "kg")
          ) %>%
          # woods
          add_trace(
            x = ~days,
            y = ~woods,
            mode = "lines",
            type = "scatter",
            color = I("red"),
            name = "Wood Curve",
            visible = "legendonly",
            text = paste("Day:", days, "\n", "Milk Yield:", format(round(woods, 2), nsmall = 2), "kg")
          ) %>%
          # poly4
          add_trace(
            x = ~days,
            y = ~poly4_curve,
            mode = "lines",
            type = "scatter",
            color = I("green"),
            name = "Quantile Polynomial Curve",
            visible = "legendonly",
            text = paste("Day:", days, "\n", "Milk Yield:", format(round(poly4_curve, 2), nsmall = 2), "kg")
          ) %>%
          # smoothing spline
          add_lines(
            x = ~days,
            y = ~spline,
            color = I("orange"),
            name = "Smoothing Spline",
            visible = "legendonly",
            text = paste("Day:", days, "\n", "Milk Yield:", format(round(spline, 2), nsmall = 2), "kg")
          )
      })


      # ....PLOT OF HISTOGRAM OF RESIDUALS SUB-TAB
      output$res_hist <- renderPlotly({
        res_wilmink <- days_yield_reactive$yield - sig_reactive
        hist_wilmink <- plot_ly(
          x = ~res_wilmink, type = "histogram", color = I("black"), name = "Residuals from Wilmink Curve"
        ) %>%
          layout(
            annotations = list(
              text = "Histogram of Residuals from Wilmink Curve",
              xref = "paper",
              yref = "paper",
              yanchor = "bottom",
              xanchor = "center",
              align = "center",
              x = 0.5,
              y = 1,
              showarrow = FALSE
            ),
            xaxis = list(title = "Residuals"),
            yaxis = list(title = "Frequency")
          )

        res_wood <- days_yield_reactive$yield - woods

        hist_woods <- plot_ly(
          x = ~res_wood, type = "histogram", color = I("red"), name = "Residuals from Wood Curve"
        ) %>%
          layout(
            annotations = list(
              text = "Histogram of Residuals from Wood Curve",
              xref = "paper",
              yref = "paper",
              yanchor = "bottom",
              xanchor = "center",
              align = "center",
              x = 0.5,
              y = 1,
              showarrow = FALSE
            ),
            xaxis = list(title = "Residuals"),
            yaxis = list(title = "Frequency")
          )
        res_poly4 <- days_yield_reactive$yield - poly4_curve

        hist_poly4 <- plot_ly(
          x = ~res_poly4, type = "histogram", color = I("green"), name = "Residuals from Degree 4 Polynomial"
        ) %>%
          layout(
            annotations = list(
              text = "Histogram of Residuals from Degree 4 Polynomial",
              # font = f,
              xref = "paper",
              yref = "paper",
              yanchor = "bottom",
              xanchor = "center",
              align = "center",
              x = 0.5,
              y = 1,
              showarrow = FALSE
            ),
            xaxis = list(title = "Residuals"),
            yaxis = list(title = "Frequency") # ,
          )

        res_spline <- days_yield_reactive$yield - spline

        hist_spline <- plot_ly(
          x = ~res_spline, type = "histogram", color = I("orange"), name = "Residuals from Smoothing Spline"
        ) %>%
          layout(
            annotations = list(
              text = "Histogram of Residuals from Smoothing Spline",
              xref = "paper",
              yref = "paper",
              yanchor = "bottom",
              xanchor = "center",
              align = "center",
              x = 0.5,
              y = 1,
              showarrow = FALSE
            ),
            xaxis = list(title = "Residuals"),
            yaxis = list(title = "Frequency")
          )

        subplot(hist_woods, hist_wilmink, hist_poly4, hist_spline, titleX = TRUE, nrows = 2, shareY = TRUE, margin = c(0.1, 0, 0.1, 0.1)) %>%
          layout(autosize = F, width = 1100, height = 500)
      })
      removeModal()
    }
  )

  output$param_label <- renderUI({
    ("Wilmink Parameters") %>% strong() %>% h6() %>%  paste() %>% HTML()
    #HTML(paste(strong(h6("Wilmink's Parameters"))))
  }) 

  # ...TEXT UNDER DAILY YIELD PLOT
  output$clicklegend <- renderUI({
    HTML(paste("<br/><br/><br/><br/><br/>", h5("Click on key in legend to add/remove corresponding line. "), sep = "<br/>"))
  })

  # ...TEXT UNDER HISTOGRAMS
  output$clicklegend2 <- renderUI({
    HTML(paste("<br/><br/><br/><br/><br/>", h5("Click on key in legend to add/remove corresponding histogram. "), sep = "<br/>"))
  })

  # TEXT ON TAB2
  output$text <- renderText({
    "text here"
  })
}


shinyApp(ui = ui, server = server)
