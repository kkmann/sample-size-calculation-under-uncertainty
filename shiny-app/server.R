suppressPackageStartupMessages(library(tidyverse))

source("../R/priors.R")
source("../R/functions.R")
library(shiny)



shinyServer(function(input, output) {

    output$plt_operating_characteristics <- renderPlot({

                 a <- input$prior_support[1]
                 b <- input$prior_support[2]
             prior <- Normal(input$prior_mean, input$prior_sd, a, b)
        theta_null <- 0
             alpha <- input$alpha
              beta <- input$beta
        # plot the conditional/unconditional priors
        aa <- max(-1, a - 0.1)
        bb <- min(1, b + 0.1)
        # define utility optimisation
        utility <- function(n) {
            input$lambda * PoS(prior, n, qnorm(1 - alpha), input$mcid) - n
        }
        plt_prior <- tibble(
                theta = seq(aa, bb, length.out = 1000),
                `unconditional prior` = pdf(prior, theta),
                `conditional prior`   = pdf(condition(prior, lo = input$mcid), theta)
            ) %>%
            pivot_longer(-theta, names_to = "type", values_to = "PDF") %>%
            mutate(
                support = case_when(
                    theta < a & type == "unconditional prior" ~ "left",
                    theta < input$mcid & type == "conditional prior" ~ "left",
                    theta > b ~ "right",
                    TRUE ~ "support"
                )
            ) %>%
            ggplot() +
                aes(theta, PDF, linetype = type, group = interaction(support, type)) +
                geom_line() +
                scale_linetype_discrete("") +
                scale_x_continuous(expression(theta)) +
                theme_bw() +
                theme(
                    legend.position = "top"
                )

        # compute required sample sizes
        tbl_samplesizes <- tibble(
                    MCID = get_n(theta_null, input$mcid),
                    EP   = get_n_ep(theta_null, prior, mrv = input$mcid, pwr = 1 - beta, alpha = alpha),
                    MEU  =  optimize(function(n) utility(n), c(10, 1e4), maximum = TRUE)$maximum %>%
                        round,
                quantile = get_n_quantile(theta_null, prior, input$gamma, mrv = input$mcid, pwr = 1 - beta, alpha = alpha)
            ) %>%
            pivot_longer(everything(), names_to = "type", values_to = "n")

        # plot power curves
        plt_powr_curves <- full_join(
                tbl_samplesizes,
                expand_grid(
                    theta = seq(aa, bb, length.out = 1000),
                    n = tbl_samplesizes$n
                ),
                by = "n"
            ) %>%
            mutate(
                power = map2_dbl(theta, n, ~power(..1, ..2, qnorm(1 - alpha))),
                name = sprintf("%s (n = %i)", type, n)
            ) %>%
            ggplot() +
                aes(theta, power, color = name) +
                geom_line() +
                scale_color_discrete("") +
                scale_x_continuous(expression(theta), breaks = seq(-1, 1, .1)) +
                scale_y_continuous("probability to reject", breaks = seq(0, 1, .1), expand = c(0, 0)) +
                theme_bw() +
                theme(
                    panel.grid.minor = element_blank(),
                    legend.position = "top"
                )

        set.seed(42)
             n <- 1e5
        rtheta <- numeric(n)
        cprior <- condition(prior, lo = input$mcid)
        i <- 1
        while (i < n) {
            sample <- rnorm(1, mean = cprior$mu, sd = cprior$tau)
            if (between(sample, cprior$lower, cprior$upper)) {
                rtheta[i] <- sample
                i <- i + 1
            }
        }

        plt_power_cdf <- full_join(
                tbl_samplesizes,
                expand_grid(
                    rtheta = rtheta,
                    n = tbl_samplesizes$n
                ),
                by = "n"
            ) %>%
            mutate(
                rpower = map2_dbl(rtheta, n, ~power(..1, ..2, qnorm(1 - alpha))),
                name = sprintf("%s (n = %i)", type, n)
            ) %>%
            select(name, rpower) %>%
            group_by(name) %>%
            nest() %>%
            transmute(
                ecdf = map(data, ~tibble(
                    power = seq(0, 1, .01),
                    CDF   = ecdf(.$rpower)(power)
                ))
            ) %>%
            unnest(ecdf) %>%
            ggplot(aes(power, CDF, color = name)) +
            geom_line() +
            scale_color_discrete("") +
            scale_x_continuous("random power", breaks = seq(0, 1, .1)) +
            scale_y_continuous(breaks = seq(0, 1, .1)) +
            coord_cartesian(expand = FALSE) +
            theme_bw() +
            theme(
                panel.grid.minor = element_blank(),
                legend.position = "top"
            )

        legend <- cowplot::get_legend(plt_powr_curves)

        cowplot::plot_grid(
            plt_prior,
            legend,
            plt_powr_curves + theme(legend.position = "none"),
            plt_power_cdf + theme(legend.position = "none"),
            ncol = 1,
            rel_heights = c(1.15, .15, 1, 1)
        )

    })

})
