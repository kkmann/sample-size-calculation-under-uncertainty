library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    withMathJax(),
    # Application title
    titlePanel("Sample size calculation under uncertainty"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("prior_support",
                        "Prior support:",
                          min = -1,
                          max = 1,
                         step = 0.01,
                        ticks = FALSE,
                        value = c(-.3, .7)),
            sliderInput("prior_mean",
                        "Prior mean:",
                        min = -1,
                        max = 1,
                        step = 0.01,
                        ticks = FALSE,
                        value = .4),
            sliderInput("prior_sd",
                        "Prior standard deviation:",
                        min = 0.01,
                        max = 1,
                        step = 0.01,
                        ticks = FALSE,
                        value = .2),
            sliderInput("mcid",
                        "\\(\\theta_{MCID}\\):",
                        min = 0.0,
                        max = 1,
                        step = 0.01,
                        ticks = FALSE,
                        value = .1),
            sliderInput("alpha",
                        "\\(\\alpha\\):",
                        min = 0.01,
                        max = 0.2,
                        step = 0.01,
                        ticks = FALSE,
                        value = .025),
            sliderInput("beta",
                        "\\(\\beta\\):",
                        min = 0.05,
                        max = 0.5,
                        step = 0.01,
                        ticks = FALSE,
                        value = .2),
            sliderInput("gamma",
                        "\\(\\gamma\\):",
                        min = 0.05,
                        max = 0.95,
                        step = 0.01,
                        ticks = FALSE,
                        value = .5),
            sliderInput("lambda",
                        "\\(\\lambda\\):",
                        min = 1,
                        max = 10000,
                        step = 1,
                        ticks = FALSE,
                        value = 1500)
        ),

        mainPanel(
            includeMarkdown("description.md"),
            plotOutput("plt_operating_characteristics", height = "625px")
        )
    )
))
