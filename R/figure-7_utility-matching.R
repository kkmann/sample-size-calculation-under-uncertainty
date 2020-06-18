options(tidyverse.quiet = TRUE)
library(tidyverse)

source("R/priors.R")
source("R/functions.R")

prior            <- Normal(0.2, 0.2, -0.3, 0.7)
# one sided maximal type one error rate
alpha            <- 0.025
# upper boundary of the null hypothesis for the location parameter
# H0: theta <= theta_null
theta_null       <- 0.0
# minimal clinically important difference (MCID)
theta_mcid       <- 0.05


utility <- function(n, lambda) {
    lambda * PoS(prior, n, qnorm(1 - alpha), theta_mcid) - n
}

get_implied_lambda <- function(expected_power) uniroot(
    function(lambda) {
        EP(
            prior,
            round(optimize(function(n) utility(n, lambda), c(10, 1e4), maximum = TRUE)$maximum),
            c = qnorm(1 - alpha),
            mrv = theta_mcid
        ) - expected_power
    },
    c(1, 1e6)
)

tbl_implied_lambda <- tibble(
        power = seq(0.6, 0.95, by = 0.01)
    ) %>%
    mutate(
        n_ep = map_dbl(
            power,
            ~get_n_ep(theta_null, prior, mrv = theta_mcid, pwr = ., alpha = alpha)
        ),
        lambda_implied = map_dbl(power, ~get_implied_lambda(.)$root)
    )

ggplot(tbl_implied_lambda) +
    aes(power, lambda_implied) +
    geom_point() +
    geom_line() +
    scale_x_continuous("expected power", breaks = seq(0.6, 0.95, by = 0.05)) +
    scale_y_continuous("implied reward [average per-patient costs]", breaks = seq(0, 20000, by = 1000)) +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank()
    )

ggsave("latex/figures/matched-reward.pdf", width = 8, height = 3.5)

optimize(function(n) utility(n, 3333), c(10, 1e4), maximum = TRUE)$maximum %>%
    round

EP(prior,
   342,
   c = qnorm(1 - alpha),
   mrv = theta_mcid
)

get_n_ep(0, prior, 0.05, pwr = 0.8, alpha = 0.025)
EP(prior, 218, qnorm(1 - alpha), 0.05)
