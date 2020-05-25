library(shiny)
library(zoo)
library(smoother)
library(lubridate)
library(tidyverse)
library(openxlsx)
library(HDInterval)
library(plotly)
library(shinydashboard)

# DT, leaflet, plotly

# r_t_range is a vector of possible values for R_t
R_T_MAX = 12
r_t_range = seq(0, R_T_MAX, length = R_T_MAX*100 + 1)

# Gamma is 1/serial interval
# https://wwwnc.cdc.gov/eid/article/26/6/20-0357_article
GAMMA = 1/4

# New cases by day
k =  c(20, 40, 55, 90)


smooth_new_cases <- function(cases){
    cases %>% 
        arrange(date) %>%
        mutate(new_cases = c(cases[1], diff(cases))) %>% 
        mutate(new_cases_smooth = round(
            smoother::smth(new_cases, window = 7, tails = TRUE)
        )) %>%
        select(date, new_cases, new_cases_smooth) # select(state, date, new_cases, new_cases_smooth)
}


compute_likelihood <- function(cases){
    likelihood <- cases %>%
        filter(new_cases_smooth > 0) %>% 
        arrange(date) %>%
        crossing(r_t = r_t_range) %>% 
        group_by(r_t) %>% 
        mutate(lambda = lag(new_cases_smooth, 1) * exp(GAMMA * (r_t - 1))) %>%
        ungroup() %>%
        mutate(likelihood_r_t = dpois(new_cases_smooth, lambda, log = TRUE)) %>% 
        filter(date > min(date))
}


compute_posterior <- function(likelihood){
    likelihood %>% 
        arrange(date) %>% 
        group_by(r_t) %>% 
        mutate(posterior = exp(
            zoo::rollapplyr(likelihood_r_t, 7, sum, partial = TRUE)
        )) %>% 
        group_by(date) %>% 
        mutate(posterior = posterior / sum(posterior, na.rm = TRUE)) %>% 
        # HACK: NaNs in the posterior create issues later on. So we remove them.
        mutate(posterior = ifelse(is.nan(posterior), 0, posterior)) %>%
        ungroup() %>%
        select(-likelihood_r_t)
}


# Estimate R_t and a 95% highest-density interval around it
estimate_rt <- function(posteriors){
    posteriors %>% 
        group_by(date) %>% # group_by(state, date) %>% 
        summarize(
            r_t_simulated = list(sample(r_t_range, 10000, replace = TRUE, prob = posterior)),
            r_t_most_likely = r_t_range[which.max(posterior)]
        ) %>% 
        mutate(
            r_t_lo = map_dbl(r_t_simulated, ~ hdi(.x)[1]),
            r_t_hi = map_dbl(r_t_simulated, ~ hdi(.x)[2])
        ) %>% 
        select(-r_t_simulated) 
}

# Define UI for application that draws a histogram
ui <- dashboardPage(
    dashboardHeader(title = "Recife - PE COVID19"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Casos", tabName = "cases", icon = icon("diagnoses")),
            menuItem("Óbitos", tabName = "deaths", icon = icon("procedures")),
            menuItem("Transmissão", tabName = "transmission", icon = icon("viruses"))
        )
    ),
    dashboardBody(tabItems(
        tabItem(tabName = "cases",
                fluidRow(
                    box(shiny::plotOutput("plot")),
                    box(shiny::plotOutput("plot_week"))
                ),
                fluidRow(
                    box(shiny::plotOutput("plot_day_pct")),
                    box(shiny::plotOutput("plot_week_pct"))
                ),
                fluidRow(
                    box(shiny::plotOutput("plot_total")),
                    box(shiny::plotOutput("plot_total_week"))
                )
        ),
        tabItem(tabName = "deaths", 
                fluidRow(
                    box(shiny::plotOutput("deaths_plot")),
                    box(shiny::plotOutput("deaths_plot_week"))
                ),
                fluidRow(
                    box(shiny::plotOutput("deaths_plot_day_pct")),
                    box(shiny::plotOutput("deaths_plot_week_pct"))
                ),
                fluidRow(
                    box(shiny::plotOutput("deaths_plot_total")),
                    box(shiny::plotOutput("deaths_plot_total_week"))
                )
        ),
        tabItem(tabName = "transmission", 
                fluidRow(
                    box(shiny::plotOutput("rt")),
                )
        )
    ))
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    data <- read.csv("https://raw.githubusercontent.com/gppeixoto/rec-covid/master/rec.csv") %>% select(c(-X)) %>% tibble() %>% 
        mutate(
            date = ymd(date),
            epi_week = lubridate::epiweek(date),
            new_cases = cases - lag(cases, 1, 0),
            new_deaths = deaths - lag(deaths, 1, 0),
            rolling_cases = smoother::smth(new_cases, window = 7, tails = TRUE),
            rolling_deaths = smoother::smth(new_deaths, window = 7, tails = TRUE),
        )
    
    output$plot <- shiny::renderPlot({
        data %>% ggplot(aes(x = date)) + theme_minimal() + 
            geom_col(aes(y = new_cases, fill = "Novos casos")) +
            geom_line(aes(y = rolling_cases, fill = "Média (7 dias)"), linetype = "dotted") +
            scale_fill_manual(name = "", values = c("Novos casos" = "orange", "Média (7 dias)" = "black")) +
            labs(x="Data", y="Novos casos", title = "Casos reportados")
    })
    
    output$plot_week <- shiny::renderPlot({
        data %>% group_by(epi_week) %>% summarise(new_cases = sum(new_cases)) %>% ggplot() + 
            geom_col(aes(epi_week, new_cases), fill="orange") + theme_minimal() +
            labs(x="Semana epidemiológica", y="Novos casos", title = "Novos casos por semana epidemiológica")
    })
    
    output$plot_week_pct <- shiny::renderPlot({
        data %>% group_by(epi_week) %>% summarise(new_cases = sum(new_cases)) %>% 
            mutate(pct_change = (new_cases / lag(new_cases, 1, 1)) - 1) %>% tail() %>% 
            ggplot() + geom_line(aes(epi_week, pct_change)) + geom_hline(aes(yintercept=0), linetype="dotted") + 
            theme_minimal() + scale_y_continuous(labels = scales::percent_format()) + 
            labs(
                x = "Semana epidemiológica", y = "Aumento percentual", 
                title = "Aumento de novos casos por semana epidemiológica", 
                subtitle = "Quantos casos a mais foram reportados que na semana anterior. Abaixo de 0% indica que a taxa de novos casos está diminuindo."
            )
    })
    
    output$plot_day_pct <- shiny::renderPlot({
        data %>% mutate(pct_change = (new_cases / lag(new_cases, 1, 1)) - 1) %>% tail(30) %>% ggplot() + 
            geom_line(aes(date, pct_change), alpha=0.2) + theme_minimal() + scale_y_continuous(labels = scales::percent_format()) + labs(
                x="Semana epidemiológica", y = "Aumento percentual", title = "Aumento relativo diário de novos casos", 
                subtitle = "Quantos casos a mais foram reportados que no dia anterior. Abaixo de 0% indica que a taxa de novos casos está diminuindo."
            ) + 
            geom_smooth(aes(date, pct_change), color="orange", fill="orange", alpha=0.3)
    })
    
    output$plot_total <- shiny::renderPlot({
        data %>% ggplot() + geom_line(aes(date, cases, color="Total acumulado"), size=1.5) + 
            geom_col(aes(date, new_cases, fill="Novos casos"), alpha=0.5) + theme_minimal() + 
            labs(y="", x="Data", title = "Total de casos reportados", subtitle = "Eixo vertical em escala linear.") + 
            scale_color_manual(values=c("Total acumulado" = "orange"), name = "") + 
            scale_fill_manual(values=c("Novos casos" = "orange"), name = "")
    })
    
    output$plot_total_week <- shiny::renderPlot({
        data %>% group_by(epi_week) %>% summarise(new_cases = sum(new_cases)) %>% mutate(total_cases = cumsum(new_cases)) %>% ggplot() + 
            geom_col(aes(epi_week, new_cases, fill = "Novos casos"), alpha=0.5) + 
            geom_line(aes(epi_week, total_cases, color="Total casos")) + theme_minimal() +
            labs(x="Semana epidemiológica", y="Novos casos", title = "Casos por semana epidemiológica", y ="", subtitle = "Eixo vertical em escala linear.") + 
            scale_fill_manual(values=c("Novos casos" = "orange"), name = "") + scale_color_manual(values=c("Total casos" = "orange"), name="")
    })
    
    ## Deaths
    ####################################################################################
    
    output$deaths_plot <- shiny::renderPlot({
        data %>% ggplot(aes(x = date)) + theme_minimal() + 
            geom_col(aes(y = new_deaths, fill = "Novos óbitos")) +
            geom_line(aes(y = rolling_deaths, fill = "Média (7 dias)"), linetype = "dashed") +
            scale_fill_manual(name = "", values = c("Novos óbitos" = "yellow4", "Média (7 dias)" = "black")) +
            labs(x="Data", y="Novos óbitos", title = "Óbitos reportados")
    })
    
    output$deaths_plot_week <- shiny::renderPlot({
        data %>% group_by(epi_week) %>% summarise(new_deaths = sum(new_deaths)) %>% ggplot() + 
            geom_col(aes(epi_week, new_deaths), fill="yellow4") + theme_minimal() +
            labs(x="Semana epidemiológica", y="Novos óbitos", title = "Novos óbitos por semana epidemiológica")
    })
    
    output$deaths_plot_week_pct <- shiny::renderPlot({
        data %>% group_by(epi_week) %>% summarise(new_deaths = sum(new_deaths)) %>% 
            mutate(pct_change = (new_deaths / lag(new_deaths, 1, 1)) - 1) %>% tail() %>% 
            ggplot() + geom_line(aes(epi_week, pct_change)) + geom_hline(aes(yintercept=0), linetype="dotted") + 
            theme_minimal() + scale_y_continuous(labels = scales::percent_format()) + 
            labs(
                x = "Semana epidemiológica", y = "Aumento percentual", 
                title = "Aumento de novos óbitos por semana epidemiológica", 
                subtitle = "Quantos óbitos a mais foram reportados que na semana anterior. Abaixo de 0% indica que a taxa de novos óbitos está diminuindo."
            )
    })
    
    output$deaths_plot_day_pct <- shiny::renderPlot({
        data %>% mutate(pct_change = (new_deaths / lag(new_deaths, 1, 1)) - 1) %>% tail(30) %>% ggplot() + 
            geom_line(aes(date, pct_change), alpha=0.2) + theme_minimal() + scale_y_continuous(labels = scales::percent_format()) + labs(
                x="Semana epidemiológica", y = "Aumento percentual", title = "Aumento relativo diário de novos óbitos", 
                subtitle = "Quantos óbitos a mais foram reportados que no dia anterior. Abaixo de 0% indica que a taxa de novos óbitos está diminuindo."
            ) + 
            geom_smooth(aes(date, pct_change), color="yellow4", fill="yellow4", alpha=0.3)
    })
    
    output$deaths_plot_total <- shiny::renderPlot({
        data %>% ggplot() + geom_line(aes(date, deaths, color="Total acumulado"), size=1.5) + 
            geom_col(aes(date, new_deaths, fill="Novos óbitos"), alpha=0.5) + theme_minimal() + 
            labs(y="", x="Data", title = "Total de óbitos reportados", subtitle = "Eixo vertical em escala linear.") + 
            scale_color_manual(values=c("Total acumulado" = "yellow4"), name = "") + 
            scale_fill_manual(values=c("Novos óbitos" = "yellow4"), name = "")
    })
    
    output$deaths_plot_total_week <- shiny::renderPlot({
        data %>% group_by(epi_week) %>% summarise(new_deaths = sum(new_deaths)) %>% mutate(total_deaths = cumsum(new_deaths)) %>% ggplot() + 
            geom_col(aes(epi_week, new_deaths, fill = "Novos óbitos"), alpha=0.5) + 
            geom_line(aes(epi_week, total_deaths, color="Total óbitos")) + theme_minimal() +
            labs(x="Semana epidemiológica", y="Novos óbitos", title = "óbitos por semana epidemiológica", y ="", subtitle = "Eixo vertical em escala linear.") + 
            scale_fill_manual(values=c("Novos óbitos" = "yellow4"), name = "") + scale_color_manual(values=c("Total óbitos" = "yellow4"), name="")
    })
    
    est <- rec %>% 
        smooth_new_cases() %>% 
        compute_likelihood() %>% 
        compute_posterior() %>% 
        estimate_rt()
    
    current_rt <- tail(est, 1)$r_t_most_likely
    
    output$rt <- shiny::renderPlot({
        est %>%
            ggplot(aes(x = date, y = r_t_most_likely)) +
            geom_point(color = "darkorange", alpha = 0.8, size = 4) +
            geom_line() +
            geom_hline(yintercept = 1, linetype = 'dashed') +
            geom_ribbon(
                aes(ymin = pmin(r_t_lo, 1), ymax = 1),
                fill = 'green4',
                alpha = 0.2
            ) +
            geom_ribbon(
                aes(ymin = 1, ymax = pmax(1, r_t_hi)),
                fill = 'red4',
                alpha = 0.2
            ) +
            labs(
                title = expression('Recife COVID19 Real time R'[t]),
                x = "Data",
                y = expression("R"[t]),
                subtitle = paste("Current expected Rt:", current_rt)
            ) + theme_minimal()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
