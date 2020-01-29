Normal <- function(mu, tau, lo = -Inf, hi = Inf) {
    lower <- max(mu - 10, lo)
    upper <- min(mu + 10, hi)
    if (upper <= lower) stop('degenerate')
    normalizing_constant <- pnorm(upper, mu, tau) - pnorm(lower, mu, tau)
    res <- list(
        mu = mu, tau = tau,
        lower = lower, upper = upper,
        normalizing_constant = normalizing_constant
    )
    class(res) <- c('Normal', class(res))
    return(res)
}

condition <- function(dist, lo = -Inf, hi = Inf) {
    Normal(dist$mu, dist$tau, max(lo, dist$lower), min(hi, dist$upper))
}

posterior <- function(prior, estimate, n1) {
    mu    <- 1/(1/prior$tau^2 + n1/1) * (prior$mu/prior$tau^2 + n1*estimate)
    tau   <- sqrt(1/(1/prior$tau^2 + n1/1))
    Normal(mu, tau, prior$lower, prior$upper)
}

pdf <- function(dist, Delta) ifelse(
    Delta < dist$lower | Delta > dist$upper,
    0,
    dnorm(Delta, dist$mu, dist$tau) / dist$normalizing_constant
)

quantile <- function(dist, prob) {
    res <- qnorm(pnorm(dist$lower, dist$mu, dist$tau) +
        prob*(pnorm(dist$upper, dist$mu, dist$tau) - pnorm(dist$lower, dist$mu, dist$tau)), dist$mu, dist$tau)
    res
}
