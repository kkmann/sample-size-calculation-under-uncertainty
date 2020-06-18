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

# redifine without rounding
get_n_quantile <- function(null, prior, prob, mrv = null, pwr = .8, alpha = .025, upper_n = 1e4) {
    qnt <- quantile(condition(prior, lo = mrv), 1 - prob)
    f   <- function(n) power(qnt, n, qnorm(1 - alpha)) - pwr
    n   <- tryCatch(
        uniroot(f, lower = 1, upper = upper_n)$root,
        error = function(e) NA_real_
    )
    return(n)
}

get_n_ep <- function(null, prior, mrv = null, pwr = .8, alpha = .025, upper_n = 1e4) {
    f   <- function(n) EP(prior, n, qnorm(1 - alpha)) - pwr
    n   <- tryCatch(
        uniroot(f, lower = 1, upper = upper_n)$root,
        error = function(e) NA_real_
    )
    return(n)
}

match_gamma <- function(n, prior) {

    uniroot(
        function(gamma) n - get_n_quantile(theta_null, prior, gamma, pwr = 1 - beta, alpha = alpha, upper_n = 10^4),
        interval = c(0.1, 0.9)
    )

}

tbl_grid <- expand_grid(
        mu  = prior_mean_range ,
        tau = prior_sd_range
    ) %>%
    mutate(
        prior         = map2(mu, tau, ~Normal(..1, ..2, lo = prior_lo, hi = prior_hi)),
        pr_relevant   = map_dbl(prior, ~(1 - cdf(., theta_mcid))),
        n_ep          = map_dbl(prior, ~get_n_ep(theta_null, ., pwr = 1 - beta, alpha = alpha, upper_n = nmax))
    ) %>%
    filter(
        !is.na(n_ep)
    ) %>%
    mutate(
        gamma_matched = map2(prior, n_ep, ~as_tibble(match_gamma(..2, ..1)[c('root', 'f.root')]))
    ) %>%
    unnest(gamma_matched) %>%
    transmute(
        mu,
        tau,
        n_ep,
        pr_relevant,
        gamma_matched = root,
        pr_reject_geq_target = root * pr_relevant
    ) %>%
    pivot_longer(c(pr_relevant, gamma_matched, pr_reject_geq_target), names_to = "probability")

print(tbl_grid)

tbl_grid %>%
    filter(probability == "pr_reject_geq_target", mu == 0.3) %>%
    arrange(desc(value))

tbl_poi <- tibble(
    mu = c(-.25, .3, .5),
    tau = c(.4, .125, .05)
)

tbl_grid %>%
    ggplot() +
        aes(mu, tau) +
        geom_raster(aes(fill = value)) +
        geom_point(data = tbl_poi,
            color = "white"
        ) +
        #geom_contour(aes(z = value), color = "white", bins = 10) +
        scale_fill_gradient(low = '#FFFFFF', high = '#000000') +
        coord_cartesian(expand = FALSE) +
        xlab('prior mean') +
        ylab('prior standard deviation') +
        facet_wrap(~probability) +
        guides(
            fill = guide_colorbar("",
                                  barwidth = grid::unit(15, "lines"),
                                  barheight = grid::unit(.5, "lines")
            )
        ) +
        theme_bw() +
        theme(
            panel.grid      = element_blank(),
            panel.spacing   = unit(1.25, 'lines'),
            legend.position = 'top'
        )

# save plot as pdf
dir.create("latex/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("latex/figures/quantile-vs-ep.pdf", width = 8, height = 3.75)

tbl_poi <- tbl_poi %>%
    mutate(
        prior = map2(mu, tau, ~Normal(..1, ..2, prior_lo, prior_hi)),
        n     = map_dbl(prior, ~get_n_ep(theta_null, ., theta_mcid, 1 - beta, alpha)),
        label = sprintf("mean=%.2f, sd=%.2f, n=%i", mu, tau, round(n))
    )

plt_priors <- tbl_poi %>%
    mutate(
        tmp = map(prior, ~tibble(
                theta = seq(prior_lo - 0.1, prior_hi + 0.1, .01),
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
        geom_line() +
        facet_wrap(~name) +
        scale_linetype_discrete("prior") +
        scale_x_continuous(expression(theta)) +
        scale_y_continuous("PDF") +
        theme_bw() +
        theme(
            legend.position  = "top",
            panel.grid.minor = element_blank(),
            strip.text       = element_text(size = 6)
        )

plt_priors


plt_power_curves <- full_join(
        tbl_poi,
        expand_grid(
            theta = seq(prior_lo, prior_hi, by = 0.01),
            n = tbl_poi$n
        ),
        by = "n"
    ) %>%
    mutate(
        power = map2_dbl(theta, n, ~power(..1, ..2, qnorm(1 - alpha)))
    ) %>%
    ggplot() +
    aes(theta, power, linetype = label) +
    geom_line() +
    scale_color_discrete("") +
    scale_x_continuous(expression(theta), breaks = c(0, 0.5)) +
    scale_y_continuous("probability to reject", breaks = seq(0, 1, .2)) +
    scale_linetype("") +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank(),
        legend.position = "top"
    )

plt_power_curves

set.seed(42)
tbl_samples <- tbl_poi %>%
    mutate(
        unconditional = map(prior, ~get_sample(., 1e5)),
        conditional   = map(prior, ~get_sample(condition(., lo = theta_mcid), 1e5))
    ) %>%
    unnest(c(conditional, unconditional)) %>%
    pivot_longer(c(conditional, unconditional), names_to = "type", values_to = "rtheta") %>%
    mutate(
        power  = map2_dbl(rtheta, n, ~power(..1, ..2, qnorm(1 - alpha)))
    )

tbl_probs <- tbl_samples %>%
    group_by(label, type) %>%
    summarise(
        tmp   = mean(power >= 1 - beta),
        power = 1 - beta
    )

plt_power_distribution <- tbl_samples %>%
    mutate(
        type = ifelse(type == "conditional",
            "random power",
            "random probability to reject"
        )
    ) %>%
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
            data  = tbl_probs %>%
                mutate(
                    type = ifelse(type == "conditional",
                                  "random power",
                                  "random probability to reject"
                    )
                )
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

legend <- cowplot::get_legend(plt_power_curves)

cowplot::plot_grid(
    legend,
    cowplot::plot_grid(
        plt_priors + theme(legend.position = "none"),
        plt_power_curves + theme(legend.position = "none"),
        rel_widths = c(2, 1.2),
        nrow = 1
    ),
    plt_power_distribution,
    rel_heights = c(.1, 1, 1.75),
    ncol = 1
)

ggsave("latex/figures/power-distribution-examples.pdf", width = 7.5, height = 6)
