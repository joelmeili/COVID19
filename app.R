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
library(gridExtra)
library(kableExtra)

# Load data
source("~/COVID19/load_data.R")

countries <- c("Italy", "Iran", "Japan", "Korea, South", "Russia", "Singapore", "Switzerland", "United Kingdom", "US")

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

state <- c(S = N - 1/N, I = 1/N, R = 0)

estim <- allData %>% group_by(`Country/Region`) %>% group_modify(~{nlminb(rep(0, 2), objective = SIR_estimation, true_infected = .x$Infected, true_removed = .x$Removed,
                                                                          state = state, times = 1:nrow(.x))$par %>%
        tibble::enframe(name = "param", value = "value") %>% spread(param, value)}) %>% ungroup()

colnames(estim)[2:3] <- c("beta", "gamma")

estim$R0 <- N*estim$beta/estim$gamma

forecast <- 300

sim <- estim %>% group_by(`Country/Region`) %>% group_map(function(.x, group_info){
    SIR_simulation(state = state, times = 1:forecast, func = SIR, parameters = c(.x$beta, .x$gamma)) %>% as.data.frame() %>%
        mutate(`Country/Region` = group_info$`Country/Region`, date = (allData %>% filter(`Country/Region` == group_info$`Country/Region`) %>%
                                                                           select(date) %>% pull() %>% min()) + 0:(forecast - 1)) %>% select(-time)
}) %>% bind_rows()


df <- allData %>% right_join(sim, on = c("date", "Country/Region")) %>%
    transmute(date = date, `Country/Region` = `Country/Region`, Susceptible = Susceptible, "Susceptible Simulated" = S , Infected = Infected, "Infected Simulated" = I, Removed = Removed, "Removed Simulated" = R) %>% 
    gather(SIR, Cases, 3:8) %>% separate(SIR, into = c("SIR", "Actual/Simulated"))

df$`Actual/Simulated` <- ifelse(is.na(df$`Actual/Simulated`), "Actual", df$`Actual/Simulated`)

ui <- fluidPage(
    
    titlePanel("COVID19 SIR analysis"),
    
    sidebarLayout(
        
        sidebarPanel(
            
            selectInput("country", "Country:",
                        countries)),
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("Active COVID19 cases", plotOutput("plot1")),
                        tabPanel("Active COVID19 cases across countries", plotOutput("plot2")),
                        tabPanel("Estimated SIR parameters", htmlOutput("estim")),
                        tabPanel("Estimated vs. Actual data plot", plotOutput("plot3")),
                        tabPanel("Estimated vs. Actual data (log-scaled)", plotOutput("plot4")),
                        tabPanel("Simulated number of infected cases across countries", plotOutput("plot5"))
            )
        )
    )
)

server <- function(input, output) {
    get_country <- reactive({
        input$country
    })
    
    get_df <- reactive({
        df %>% filter(`Country/Region` == get_country()) %>% select(-`Country/Region`)
    })
    
    output$estim <- renderText({
        estim %>% filter(`Country/Region` == get_country()) %>% kable(digits = 10) %>% kable_styling(fixed_thead = TRUE)
    })
    
    output$plot1 <- renderPlot({
        get_df() %>% filter(`Actual/Simulated` == "Actual" & SIR == "Infected") %>% na.omit() %>% ggplot(aes(x = date, y = Cases)) + geom_line() + theme_base() + xlab("Date") + ylab("# of Infections") + ggtitle("Number of infections per 10'000")
    })
    
    output$plot2 <- renderPlot({
        grid.arrange(
            df %>% filter(`Actual/Simulated` == "Actual" & SIR == "Infected") %>% na.omit() %>% ggplot(aes(x = date, y = Cases, colour = `Country/Region`)) + geom_line() + theme_base() + xlab("Date") + ylab("# of Infections") + ggtitle("Number of infections per 10'000"),
            df %>% filter(`Actual/Simulated` == "Actual" & SIR == "Infected") %>% na.omit() %>% ggplot(aes(x = date, y = Cases, colour = `Country/Region`)) + geom_line() + scale_y_log10() + theme_base() + xlab("Date") + ylab("# of Infections") + ggtitle("Number of infections per 10'000"),
            nrow = 2)
    })
    
    output$plot3 <- renderPlot({
        get_df() %>% ggplot(aes(x = date, y = Cases, colour = SIR)) + geom_line(aes(linetype = `Actual/Simulated`)) + 
            theme_base() + xlab("Date") + ylab("# of Cases") + ggtitle("Calculated SIR based on per 10'000")
    })
    
    output$plot4 <- renderPlot({
        get_df() %>% ggplot(aes(x = date, y = Cases, colour = SIR)) + geom_line(aes(linetype = `Actual/Simulated`)) + 
            scale_y_log10() + theme_base() + xlab("Date") + ylab("# of Cases") + ggtitle("Calculated SIR based on per 10'000")
    })
    
    output$plot5 <- renderPlot({
        df %>% filter(SIR == "Infected" & `Actual/Simulated` == "Simulated") %>% ggplot(aes(x = date, y = Cases, colour = `Country/Region`)) + geom_line(aes(linetype = `Actual/Simulated`)) + 
            theme_base() + xlab("Date") + ylab("# of Cases") + ggtitle("Simulated # of infected cases across countries based on per 10'000")
    })
}
 
shinyApp(ui = ui, server = server)