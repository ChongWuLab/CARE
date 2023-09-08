mr.raps.simple.robust.own <- function (b_exp, b_out, se_exp, se_out, 
                                       loss.function = c("huber", "tukey"), 
                                       k = switch(loss.function[1], huber = 1.345, tukey = 4.685), 
                                       diagnostics = FALSE) 
{
  loss.function <- match.arg(loss.function, c("huber", "tukey"))
  rho <- switch(loss.function, huber = function(r, ...) rho.huber(r, 
                                                                  k, ...), tukey = function(r, ...) rho.tukey(r, k, ...))
  delta <- integrate(function(x) x * rho(x, deriv = 1) * dnorm(x), 
                     -Inf, Inf)$value
  c1 <- integrate(function(x) rho(x, deriv = 1)^2 * dnorm(x), 
                  -Inf, Inf)$value
  c2 <- integrate(function(x) x^2 * rho(x, deriv = 1)^2 * dnorm(x), 
                  -Inf, Inf)$value - delta^2
  c3 <- integrate(function(x) x^2 * rho(x, deriv = 2) * dnorm(x), 
                  -Inf, Inf)$value
  robust.loglike <- function(beta) {
    -sum(rho((b_out - b_exp * beta)/sqrt(se_out^2 + se_exp^2 * 
                                           beta^2)))
  }
  #bound <- quantile(abs(b_out/b_exp), 0.95, na.rm = TRUE) * 
  #  2
  bound = 1
  beta.hat <- optimize(robust.loglike, bound * c(-1, 1), maximum = TRUE, 
                       tol = .Machine$double.eps^0.5)$maximum
  #while (abs(beta.hat) > 0.95 * bound) {
  #  bound <- bound * 2
  #  beta.hat <- optimize(robust.loglike, bound * c(-1, 1), 
  #                       maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum
  #}
  score.var <- c1 * sum(((b_exp^2 - se_exp^2) * se_out^2 + 
                           (b_out^2 - se_out^2) * se_exp^2 + se_exp^2 * se_out^2)/(se_out^2 + 
                                                                                     beta.hat^2 * se_exp^2)^2)
  I <- delta * sum(((b_exp^2 - se_exp^2) * se_out^2 + (b_out^2 - 
                                                         se_out^2) * se_exp^2)/(se_out^2 + beta.hat^2 * se_exp^2)^2)
  dif <- b_out - beta.hat * b_exp
  dif.var <- se_out^2 + beta.hat^2 * se_exp^2
  chi.sq.test <- sum((dif/sqrt(dif.var))^2)
  if (diagnostics) {
    std.resid <- (b_out - b_exp * beta.hat)/sqrt((se_out^2 + 
                                                    beta.hat^2 * se_exp^2))
    par(mfrow = c(1, 2))
    qqnorm(std.resid)
    abline(0, 1)
    beta.hat.loo <- rep(NA, length(b_out))
    if (length(b_out) > 100) {
      a <- quantile(abs(b_exp/se_exp), 1 - 100/length(b_out))
    }
    else {
      a <- 0
    }
    for (i in 1:length(b_out)) {
      if (abs(b_exp[i]/se_exp[i]) > a) {
        beta.hat.loo[i] <- mr.raps.simple.robust(b_exp[-i], 
                                                 b_out[-i], se_exp[-i], se_out[-i], loss.function = loss.function, 
                                                 k = k)$beta.hat
      }
    }
    plot(abs(b_exp/se_exp), beta.hat.loo)
    abline(h = beta.hat)
    print(ad.test(std.resid))
    print(shapiro.test(std.resid))
  }
  asymp.var <- solve(I) %*% score.var %*% t(solve(I))
  out <- list(beta.hat = beta.hat, beta.se = sqrt(score.var/I^2), 
              naive.se = sqrt(1/I), chi.sq.test = chi.sq.test, beta.p.value = min(1, 
                                                                                  2 * (1 - pnorm(abs(beta.hat)/sqrt(asymp.var[1, 1])))))
  if (diagnostics) {
    out$std.resid <- std.resid
    out$beta.hat.loo <- beta.hat.loo
  }
  out
}


environment(mr.raps.simple.robust.own) <- asNamespace('mr.raps')
assignInNamespace("mr.raps.simple.robust", mr.raps.simple.robust.own, ns = "mr.raps")
