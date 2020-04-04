#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(countrycode)
library(wbstats)
library(deSolve)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(kableExtra)

# Load data
source("load_data.R")

countries <- c("Austria", "Germany", "Italy", "Iran", "Japan", "Korea, South", "Netherlands", "Portugal", "Russia", "Singapore", "Spain", 
               "Switzerland", "United Kingdom", "US")

countries_population <- data.frame(countries) %>% rowwise() %>% 
    mutate(N = wb(country = countrycode(countries, origin = "country.name", destination = "iso2c"), indicator = "SP.POP.TOTL")$value[1])
countries_population$countries <- as.character(countries_population$countries)
colnames(countries_population) <- c("Country/Region", "N")

allData <- loadData("time_series_covid19_confirmed_global.csv", "CumConfirmed") %>%
    inner_join(loadData("time_series_covid19_deaths_global.csv", "CumDeaths")) %>%
    inner_join(loadData("time_series_covid19_recovered_global.csv", "CumRecovered")) %>%
    na.omit() %>% filter(`Province/State` == "<all>" & `Country/Region` %in% countries) %>%
    group_by(date, `Country/Region`) %>% transmute(Infected = sum(CumConfirmed), Recovered = sum(CumRecovered), Deaths = sum(CumDeaths)) %>%
    mutate(Infected = Infected - Recovered - Deaths) %>% filter(Infected != 0) %>% inner_join(countries_population, on = `Country/Region`) %>%
    transmute(Susceptible = N - Infected - Recovered - Deaths, Infected = Infected, Recovered = Recovered, Deaths = Deaths, N = N) %>% ungroup()

SIR <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        dS <- -beta*I*S
        dI <- beta*I*S - gamma*I - delta*I
        dR <- gamma*I
        dD <- delta*I
        
        list(c(dS, dI, dR, dD))
    })
}

SIR_simulation <- function(state, times, func, parameters) {
    out <- ode(y = state, times = times, func = SIR, parms = c(beta = parameters[1], gamma = parameters[2], delta = parameters[3]))
    out
}

SIR_estimation <- function(true_infected, true_recovered, true_dead, parameters, state, times) {
    out <- SIR_simulation(state = state, times = times, func = SIR, parameters = parameters)
    l1 <- sqrt(mean((true_infected - out[, "I"])^2))
    l2 <- sqrt(mean((true_recovered - out[, "R"])^2))
    l3 <- sqrt(mean((true_dead - out[, "D"])^2))
    l1 + l2 + l3
}

ui <- fluidPage(
    
    titlePanel("COVID19 SIRD analysis"),
    
    sidebarLayout(
        
        sidebarPanel(
            dateRangeInput("date_range", "Estimation interval:", start = min(allData$date), end = max(allData$date), 
            min = min(allData$date), max = max(allData$date)),
            selectInput("country", "Country:", countries),
            numericInput("forecast", "Number of days to be forecasted:", value = 60, min = 1, max = 730), interval = 1),
        
        mainPanel(
            textOutput("date_range"),
            tabsetPanel(type = "tabs",
                       tabPanel("Active COVID19 cases", plotOutput("plot1")),
                       tabPanel("Active COVID19 cases across countries", plotOutput("plot2"), plotOutput("plot3")),
                       tabPanel("Estimated SIR parameters", htmlOutput("estim")),
                       tabPanel("Estimated vs. Actual data plot", plotOutput("plot4")),
                       tabPanel("Estimated vs. Actual data (log-scaled)", plotOutput("plot5")),
                       tabPanel("Simulated number of infected cases across countries", plotOutput("plot6"))
            )
        )
    )
)

server <- function(input, output) {
    get_date_range <- reactive({
        input$date_range
    })
    
    get_country <- reactive({
        input$country
    })
    
    get_forecast <- reactive({
        input$forecast
    })
    
    estimate_parameters <- reactive({
        date_range <- get_date_range()
        estim <- allData %>% filter(between(date, date_range[1], date_range[2])) %>% group_by(`Country/Region`) %>% group_modify(~{optim(rep(0, 3), fn = SIR_estimation, true_infected = .x$Infected, true_recovered = .x$Recovered,
                                                                                 true_dead = .x$Deaths, state = c(S = max(.x$N) - 1, I = 1, R = 0, D = 0), times = 1:nrow(.x))$par %>%
                tibble::enframe(name = "param", value = "value") %>% spread(param, value)}) %>% ungroup() %>% inner_join(countries_population, on = `Country/Region`)
        
        colnames(estim)[2:4] <- c("beta", "gamma", "delta")
        
        estim$R0 <- estim$N*estim$beta/(estim$gamma + estim$delta)
        estim
    })
    
    simulate_disease <- reactive({
        forecast <- get_forecast()
        estim <- estimate_parameters()
        
        sim <- estim %>% group_by(`Country/Region`) %>% group_map(function(.x, group_info){
            SIR_simulation(state = c(S = max(.x$N) - 1, I = 1, R = 0, D = 0), times = 1:forecast, func = SIR, parameters = c(.x$beta, .x$gamma, .x$delta)) %>% as.data.frame() %>%
                mutate(`Country/Region` = group_info$`Country/Region`, date = (allData %>% filter(`Country/Region` == group_info$`Country/Region`) %>%
                                                                                   select(date) %>% pull() %>% min()) + 0:(forecast - 1)) %>% select(-time)
        }) %>% bind_rows()
        
        df <- allData %>% right_join(sim, on = c("date", "Country/Region")) %>%
            transmute(date = date, `Country/Region` = `Country/Region`, Susceptible = Susceptible, "Susceptible Simulated" = S , Infected = Infected, "Infected Simulated" = I, Recovered = Recovered, "Recovered Simulated" = R,
                      Deaths = Deaths, "Deaths Simulated" = D) %>%  gather(SIR, Cases, 3:10) %>% separate(SIR, into = c("SIR", "Actual/Simulated")) %>% 
            inner_join(countries_population, on = `Country/Region`) %>%  mutate(Cases = Cases/N*1e5)
        
        df$`Actual/Simulated` <- ifelse(is.na(df$`Actual/Simulated`), "Actual", df$`Actual/Simulated`)
        df
    })
    
    get_df <- reactive({
        simulate_disease()
    })
    
    output$estim <- renderText({
        estimate_parameters() %>% filter(`Country/Region` == get_country()) %>% kable(digits = 10) %>% kable_styling(fixed_thead = TRUE)
    })
    
    output$plot1 <- renderPlot({
        get_df() %>% filter(`Country/Region` == get_country()) %>% filter(`Actual/Simulated` == "Actual" & SIR == "Infected") %>% na.omit() %>% ggplot(aes(x = date, y = Cases)) + geom_line() + 
            theme_base() + xlab("Date") + ylab("# of Infections") + ggtitle("Number of infections per 100'000")
    })
    
    output$plot2 <- renderPlot({
        get_df() %>% filter(`Actual/Simulated` == "Actual" & SIR == "Infected") %>% na.omit() %>% ggplot(aes(x = date, y = Cases, colour = `Country/Region`)) + 
            geom_line() + theme_base() + xlab("Date") + ylab("# of Infections") + ggtitle("Number of infections per 100'000")
        })
    
    output$plot3 <- renderPlot({
        get_df() %>% filter(`Actual/Simulated` == "Actual" & SIR == "Infected") %>% na.omit() %>% ggplot(aes(x = date, y = Cases, colour = `Country/Region`)) + 
            geom_line() + scale_y_log10() + theme_base() + xlab("Date") + ylab("# of Infections") + ggtitle("Number of infections per 100'000")
    })
    
    output$plot4 <- renderPlot({
        get_df() %>% filter(`Country/Region` == get_country()) %>% ggplot(aes(x = date, y = Cases, colour = SIR)) + geom_line(aes(linetype = `Actual/Simulated`)) + 
            theme_base() + xlab("Date") + ylab("# of Cases") + ggtitle("Calculated SIR per 100'000")
    })
    
    output$plot5 <- renderPlot({
        get_df() %>% filter(`Country/Region` == get_country()) %>% ggplot(aes(x = date, y = Cases, colour = SIR)) + geom_line(aes(linetype = `Actual/Simulated`)) + 
            scale_y_log10() + theme_base() + xlab("Date") + ylab("# of Cases") + ggtitle("Calculated SIR per 100'000")
    })
    
    output$plot6 <- renderPlot({
        get_df() %>% filter(SIR == "Infected" & `Actual/Simulated` == "Simulated") %>% ggplot(aes(x = date, y = Cases, colour = `Country/Region`)) + geom_line() +
            theme_base() + xlab("Date") + ylab("# of Cases") + ggtitle("Simulated # of infected cases across countries per 100'000")
    })
}
 
shinyApp(ui = ui, server = server)