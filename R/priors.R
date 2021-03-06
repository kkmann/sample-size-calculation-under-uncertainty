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

cdf <- function(dist, Delta) {
    epsilon <- (Delta - dist$mu) / dist$tau
    alpha   <- (dist$lower - Delta) / dist$tau
    beta    <- (dist$upper - Delta) / dist$tau
    z       <- pnorm(beta) - pnorm(alpha)
    return( (pnorm(epsilon) - pnorm(alpha) ) / z )
}

quantile <- function(dist, prob) {
    res <- qnorm(pnorm(dist$lower, dist$mu, dist$tau) +
        prob*(pnorm(dist$upper, dist$mu, dist$tau) - pnorm(dist$lower, dist$mu, dist$tau)), dist$mu, dist$tau)
    res
}

get_sample <- function(dist, n) {
    # sample from prior... inefficiently
    rtheta <- numeric(n)
    i <- 1
    while (i < n) {
        sample <- rnorm(1, mean = dist$mu, sd = dist$tau)
        if (between(sample, dist$lower, dist$upper)) {
            rtheta[i] <- sample
            i <- i + 1
        }
    }
    return(rtheta)
}
