#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(deSolve)
library(ggplot2)
library(ggthemes)
library(kableExtra)

# Load data
source("load_data.R")

countries <- c("Italy", "Iran", "Japan", "Korea, South", "Singapore", "Spain", "Switzerland", "US")

N <- 10000

allData <- loadData("time_series_covid19_confirmed_global.csv", "CumConfirmed") %>%
    inner_join(loadData("time_series_covid19_deaths_global.csv", "CumDeaths")) %>%
    inner_join(loadData("time_series_covid19_recovered_global.csv", "CumRecovered")) %>%
    na.omit() %>% filter(`Province/State` == "<all>" & `Country/Region` %in% countries) %>%
    group_by(date, `Country/Region`) %>% transmute(Infected = sum(CumConfirmed), Removed = sum(CumDeaths + CumRecovered)) %>%
    mutate(Infected = Infected - Removed) %>% filter(Infected != 0) %>% mutate(Infected = Infected/N, Removed = Removed/N) %>%
    transmute(Susceptible = N - Infected - Removed, Infected = Infected, Removed = Removed) %>% ungroup()

SIR <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        dS <- -beta*I*S
        dI <- beta*I*S - gamma*I
        dR <- gamma*I
        
        list(c(dS, dI, dR))
    })
}

SIR_simulation <- function(state, times, func, parameters) {
    out <- ode(y = state, times = times, func = SIR, parms = c(beta = parameters[1], gamma = parameters[2]))
    out
}

SIR_estimation <- function(true_infected, true_removed, parameters, state, times) {
    out <- SIR_simulation(state = state, times = times, func = SIR, parameters = parameters)
    l1 <- sqrt(mean((true_infected - out[, "I"])^2))
    l2 <- sqrt(mean((true_removed - out[, "R"])^2))
    alpha <- 0.5
    alpha*l1 + (1 - alpha)*l2
}

ui <- fluidPage(
    
    titlePanel("COVID19 SIR analysis"),
    
    sidebarLayout(
        
        sidebarPanel(
            
            selectInput("country", "Country:",
                        countries)),
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("Active COVID19 cases", plotOutput("plot1")),
                        tabPanel("Estimated SIR parameters", htmlOutput("estim")),
                        tabPanel("Estimated vs. Actual data plot", plotOutput("plot2")),
                        tabPanel("Estimated vs. Actual data (log-scaled)", plotOutput("plot3"))
            )
        )
    )
)

server <- function(input, output) {
    get_country <- reactive({
        input$country
    })
    
    get_df <- reactive({
        allData %>% filter(`Country/Region` == get_country()) %>% select(-`Country/Region`)
        })
    
    estimate_sir <- reactive({
        df <- get_df()
        state <- c(S = N, I = 1/N, R = 0)
        times <- 1:nrow(df)
        estim <- nlminb(rep(0, 2), objective = SIR_estimation, true_infected = df$Infected, true_removed = df$Removed,
                        state = state, times = times, lower = rep(0, 2))$par
        estim
    })
    
    output$estim <- renderText({
        estim <- estimate_sir()
        df <- data.frame(estim[1], estim[2], N*estim[1]/estim[2])
        colnames(df) <- c("beta", "gamma", "R0")
        kable(df, digits = 10) %>% kable_styling(fixed_thead = TRUE)
    })
    
    sir_simulation <- reactive({
        df <- get_df()
        estim <- estimate_sir()
        forecast <- 300
        sim <- SIR_simulation(state = c(S = N, I = 1/N, R = 0), times = 1:forecast, func = SIR, parameters = estim) %>% as_data_frame()
        sim$time <- NULL
        sim$date <- min(df$date) + 0:(forecast - 1)
        print(df %>% right_join(sim, on = "date"))
        df <- df %>% right_join(sim, on = "date") %>%
            transmute(date = date, Susceptible = Susceptible, "Susceptible Simulated" = S , Infected = Infected, "Infected Simulated" = I, Removed = Removed, "Removed Simulated" = R) %>% 
            gather(SIR, Cases, 2:7) %>% separate(SIR, into = c("SIR", "Actual/Simulated"))
        df$`Actual/Simulated` <- ifelse(is.na(df$`Actual/Simulated`), "Actual", df$`Actual/Simulated`)
        df
    })
    
    output$plot1 <- renderPlot({
        df <- get_df()
        df %>% ggplot(aes(x = date, y = Infected)) + geom_line() + theme_base() + xlab("Date") + ylab("# of Infections") + ggtitle("Number of infections per 10'000")
    })
    
    output$plot2 <- renderPlot({
        df <- sir_simulation()
        df %>% ggplot(aes(x = date, y = Cases, colour = SIR)) + geom_line(aes(linetype = `Actual/Simulated`)) + 
                         theme_base() + xlab("Date") + ylab("# of Cases") + ggtitle("Calculated SIR based on per 10'000")
    })
    
    output$plot3 <- renderPlot({
        df <- sir_simulation()
        df %>% ggplot(aes(x = date, y = Cases, colour = SIR)) + geom_line(aes(linetype = `Actual/Simulated`)) + 
            scale_y_log10() + theme_base() + xlab("Date") + ylab("# of Cases") + ggtitle("Calculated SIR based on per 10'000")
    })
}
 
shinyApp(ui = ui, server = server)