# source priors.R first

power <- function(Delta, n, c) 1 - pnorm(c, mean = sqrt(n) * Delta, sd = 1)

EP <- function(prior, n, c, mrv = 0) {
    cprior <- condition(prior, mrv)
    integrate(
        function(Delta) pdf(cprior, Delta) * power(Delta, n, c), cprior$lower, cprior$upper
    )$value
}

PoS <- function(prior, n, c, mrv = 0) {
    integrate(
        function(Delta) pdf(prior, Delta) * power(Delta, n, c), mrv, prior$upper
    )$value
}

CP <- function(zn1, n1, n, c, Delta) {
    cmu <- sqrt(n)*Delta + sqrt(n1/n) * (zn1 - sqrt(n1)*Delta)
    csd <- sqrt(1 - n1/n)
    1 - pnorm(c, mean = cmu, sd = csd)
}

OCP <- function(zn1, n1, n, c, mrv = 0) CP(zn1, n1, n, c, zn1 / sqrt(n1))

CEP <- function(prior, zn1, n1, n, c, mrv = 0) {
    cprior <- condition(prior, mrv)
    cpost  <- posterior(cprior, zn1 / sqrt(n1), n1)
    integrate(
        function(Delta) pdf(cpost, Delta) * CP(zn1, n1, n, c, Delta), cpost$lower, cpost$upper
    )$value
}


get_n_ep <- function(null, prior, mrv = null, pwr = .8, alpha = .025, upper_n = 1e4) {
    f   <- function(n) EP(prior, n, qnorm(1 - alpha)) - pwr
    n   <- tryCatch(
        uniroot(f, lower = 1, upper = upper_n)$root,
        error = function(e) NA_real_
    )
    return(ceiling(n))
}

get_n_pos <- function(null, prior, mrv = null, pwr = .8, alpha = .025, upper_n = 1e4) {
    f   <- function(n) PoS(prior, n, qnorm(1 - alpha)) - pwr
    n   <- tryCatch(
        uniroot(f, lower = 1, upper = upper_n)$root,
        error = function(e) NA_real_
    )
    return(ceiling(n))
}

get_n_quantile <- function(null, prior, prob, mrv = null, pwr = .8, alpha = .025, upper_n = 1e4) {
    qnt <- quantile(condition(prior, lo = mrv), 1 - prob)
    f   <- function(n) power(qnt, n, qnorm(1 - alpha)) - pwr
    n   <- tryCatch(
        uniroot(f, lower = 1, upper = upper_n)$root,
        error = function(e) NA_real_
    )
    return(ceiling(n))
}