options(tidyverse.quiet = TRUE)
library(tidyverse)

source("R/priors.R")
source("R/functions.R")

# maximal sample size cut-off
nmax             <-  1000
# lowest and highest plausible values of the location parameter theta
prior_lo         <- -0.3
prior_hi         <-  0.7
# range of prior means for location parameter
prior_mean_range <- seq(prior_lo, prior_hi, by = .01)
# range of prior standard deviations for location parameter
prior_sd_range   <- seq(.01, 1, by = .01)
# one sided maximal type one error rate
alpha            <- 0.025
# minimal power/expected power/probability of success is 1 - beta
beta             <- 0.2
# upper boundary of the null hypothesis for the location parameter
# H0: theta <= theta_null
theta_null       <- 0.0
# minimal clinically important difference (MCID)
theta_mcid       <- 0.1

tbl_grid <- expand_grid(
        mu  = prior_mean_range ,
        tau = prior_sd_range
    ) %>%
    mutate(
        prior          = map2(mu, tau, ~Normal(.x, .y, lo = prior_lo, hi = prior_hi)),
        `quantile 0.5` = map_dbl(prior, ~get_n_quantile(theta_null, ., .5, pwr = 1 - beta, alpha = alpha, upper_n = nmax)),
        `quantile 0.9` = map_dbl(prior, ~get_n_quantile(theta_null, ., .9, pwr = 1 - beta, alpha = alpha, upper_n = nmax)),
        EP             = map_dbl(prior, ~get_n_ep(theta_null, ., pwr = 1 - beta, alpha = alpha, upper_n = nmax)),
        PoS            = map_dbl(prior, ~get_n_pos(theta_null, ., pwr = 1 - beta, alpha = alpha, upper_n = nmax))
    ) %>%
    pivot_longer(
        c(contains('quantile'), EP, PoS),
        values_to = 'required sample size',
        names_to  = 'criterion'
    )

tbl_grid %>%
    filter(
        `required sample size` <= nmax, # make sure maximal sample size is respected
        !is.na(`required sample size`)  # throw out instances where the maximal sample size boundary was hit
    ) %>%
    ggplot() +
        aes(mu, tau, fill = `required sample size`, z = `required sample size`) +
        geom_raster() + # geom_raster leads to some pdf viewers interpolating, do not want that!
        scale_fill_gradient(limits = c(0, 1000), low = '#FFFFFF', high = '#000000') +
        guides(
            fill = guide_colorbar("required sample size",
                                  barwidth  = grid::unit(15, "lines"),
                                  barheight = grid::unit(.5, "lines")
                )
        ) +
        coord_cartesian(expand = FALSE) +
        xlab('prior mean') +
        ylab('prior standard deviation') +
        facet_wrap(~criterion, nrow = 1) +
        theme_bw() +
        theme(
            panel.grid      = element_blank(),
            panel.spacing   = unit(1.25, 'lines'),
            legend.position = 'top'
        )

# save plot as pdf
dir.create("latex/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("latex/figures/power-constraint-comparison.pdf", width = 8, height = 3.15)
