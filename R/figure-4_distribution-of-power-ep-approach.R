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
theta_mcid       <- .1

prior <- Normal(0.3, 0.2, prior_lo, prior_hi)

plt_priors <- tibble(
        mu = .3,
        tau = .2,
        prior = list(Normal(.3, .2, lo = prior_lo, hi = prior_hi)),
        label = sprintf("mean=%.2f, sd=%.2f", mu, tau)
    ) %>%
    mutate(
        tmp = map(prior, ~tibble(
            theta = seq(prior_lo - 0.1, prior_hi + 0.1, .005),
            `unconditional prior`   = pdf(., theta) %>%
                {ifelse(theta == theta[which.min(abs(prior_lo - theta))] | (theta == theta[which.min(abs(prior_hi - theta))]), NA_real_, .)},
            `conditional prior` = pdf(condition(., lo = theta_mcid), theta) %>%
                {ifelse(theta == theta[which.min(abs(theta_mcid - theta))] | (theta == theta[which.min(abs(prior_hi - theta))]), NA_real_, .)}
        )
        )
    ) %>%
    unnest(tmp) %>%
    pivot_longer(c(`conditional prior`, `unconditional prior`)) %>%
    ggplot(aes(theta, value, linetype = label)) +
    geom_line(aes(linetype = name)) +
    scale_x_continuous(expression(theta)) +
    scale_y_continuous("PDF") +
    scale_linetype_discrete("") +
    theme_bw() +
    theme(
        legend.position  = "top",
        panel.grid.minor = element_blank()
    )

plt_priors

tbl_gamma <- expand_grid(
        gamma = c(.5, .9),
        target_power = c(0.7, 0.8)
    ) %>%
    mutate(
        n     = map2_dbl(
            gamma,
            target_power,
            ~get_n_quantile(theta_null, prior, ..1, theta_mcid, ..2, alpha)),
        label = sprintf("gamma = %.2f\n1 - beta = %.2f\nn = %i", gamma, target_power, round(n))
    )


plt_power_curves <- full_join(
        tbl_gamma,
        expand_grid(
            theta = seq(prior_lo, prior_hi, by = 0.01),
            n = tbl_gamma$n
        ),
        by = "n"
    ) %>%
    mutate(
        power = map2_dbl(theta, n, ~power(..1, ..2, qnorm(1 - alpha)))
    ) %>%
    ggplot() +
    aes(theta, power, linetype = label) +
    geom_line() +
    scale_linetype_discrete("", guide = guide_legend(nrow = 1)) +
    scale_x_continuous(expression(theta), breaks =c(0, 0.5)) +
    scale_y_continuous("probability to reject", breaks = seq(0, 1, .2)) +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.text     = element_text(size = 5),
        legend.key.size = unit(0.9,"line")
    )

plt_power_curves

set.seed(42)
tbl_samples <- tbl_gamma %>%
    mutate(
        `random probability to reject` = map(gamma, ~get_sample(prior, 1e5)),
        `random power` = map(gamma, ~get_sample(condition(prior, lo = theta_mcid), 1e5))
    ) %>%
    unnest(c(`random probability to reject`, `random power`)) %>%
    pivot_longer(
        c(`random probability to reject`, `random power`),
        names_to = "type",
        values_to = "rtheta"
    ) %>%
    mutate(
        power  = map2_dbl(rtheta, n, ~power(..1, ..2, qnorm(1 - alpha)))
    )

tbl_probs <- tbl_samples %>%
    group_by(label, type) %>%
    summarise(
        tmp = mean(power >= 1 - beta),
        power = 1 - beta
    )

plt_power_distribution <- tbl_samples %>%
    ggplot(aes(power)) +
    geom_histogram(aes(y = stat(ndensity)), bins = 25, fill = "darkgray") +
    geom_vline(aes(xintercept = 0.8), color = "black") +
    geom_text(
        aes(
            label = sprintf("%.2f", tmp),
        ),
        y     = .9,
        x     = 0.91,
        size  = 3,
        color = "black",
        data  = tbl_probs
    ) +
    scale_y_continuous("", breaks = c(), limits = c(0, 1.1)) +
    scale_x_continuous("probability to reject", breaks = seq(0, 1, .2)) +
    coord_cartesian(expand = FALSE) +
    facet_grid(type ~ label, scales = "free_y") +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank(),
        panel.spacing    = unit(1.25, "lines"),
        strip.text       = element_text(size = 6)
    )

plt_power_distribution

legend_priors <- cowplot::get_legend(plt_priors)
legend_pwr    <- cowplot::get_legend(plt_power_curves)

cowplot::plot_grid(
    cowplot::plot_grid(
        plt_priors , plt_power_curves,
        rel_widths = c(1, 1)
    ),
    plt_power_distribution,
    rel_heights = c(1, 1.5),
    ncol = 1
)

ggsave("latex/figures/power-distribution-quantile-approach.pdf", width = 7.5, height = 6)
