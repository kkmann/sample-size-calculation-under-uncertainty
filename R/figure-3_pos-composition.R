options(tidyverse.quiet = TRUE)
library(tidyverse)

source("R/priors.R")
source("R/functions.R")

# one sided maximal type one error rate
alpha            <- 0.025
# minimal power/expected power/probability of success is 1 - beta
beta             <- 0.2
# upper boundary of the null hypothesis for the location parameter
# H0: theta <= theta_null
theta_null       <- 0.0
# minimal clinically important difference (MCID)
theta_mcid       <- 0.1

PoS_all <- function(prior, n, c, null = 0, mrv = null) {
    part1 <- integrate(
        function(Delta) pdf(prior, Delta) * power(Delta, n, c), prior$lower, null
    )$value
    part2 <- integrate(
        function(Delta) pdf(prior, Delta) * power(Delta, n, c), null, mrv
    )$value
    part3 <- integrate(
        function(Delta) pdf(prior, Delta) * power(Delta, n, c), mrv, prior$upper
    )$value
    return(tibble(a = part1, b = part2, c = part3))
}

get_n_pos_prime <- function(null, prior, mrv = null, pwr = .8, alpha = .025, upper_n = 1e4) {
    f   <- function(n) PoS_all(prior, n, qnorm(1 - alpha)) %>% {.$a + .$b + .$c - pwr}
    n   <- tryCatch(
        uniroot(f, lower = 1, upper = upper_n)$root,
        error = function(e) NA_real_
    )
    return(ceiling(n))
}

get_n_pos_prime(theta_null, Normal(.3, .2, -.3, .7), mrv = theta_mcid, pwr = .8, alpha = .025)

tbl_data <- expand_grid(
        prior_mean = seq(-0.1, .2, .025),
        prior_sd   = seq(.025, .15, .025)
    ) %>%
    mutate(
        tmp = pmap(
            list(prior_mean, prior_sd),
            ~PoS_all(
                Normal(..1, ..2, -.3, .7),
                150,
                qnorm(1 - .05),
                mrv = theta_mcid)
        )
    ) %>%
    unnest(tmp) %>%
    mutate(
        id = row_number(),
        total = a + b + c,
        A = a/total,
        B = b/total,
        C = c/total
    )

ggplot(tbl_data) +
    scatterpie::geom_scatterpie(
        aes(x = prior_mean, y = prior_sd, group = id, r = .01),
        data = tbl_data, cols = c("A", "B", "C"), color = NA
    ) +
    geom_text(
        aes(x = prior_mean, y = prior_sd, label = sprintf("%.2f", total)),
        size = 2
    ) +
    theme_bw() +
    coord_equal() +
    scale_x_continuous("prior mean", breaks = seq(-.2, .3, .05)) +
    scale_y_continuous("prior standard deviation", breaks = seq(0.05, .25, .05)) +
    scale_fill_manual("",
        values = c(
            A = scales::muted("red", l = 75, c = 60),
            B = scales::muted("blue", l = 75, c = 60),
            C = scales::muted("green", l = 75, c = 60)
        )
    ) +
    theme(
        legend.position = "right",
        panel.grid      = element_blank()
    )

# save plot as pdf
dir.create("latex/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("latex/figures/pos-components.pdf", width = 7.5, height = 2.75)
