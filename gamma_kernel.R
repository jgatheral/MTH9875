#######################################################
# Gamma kernel functions
#######################################################

library(gsl)
library(MittagLeffleR)


# Gamma kernel
kGamma <- function(params) {
  Vectorize(function(tau) {
    al <- params$al
    lam <- params$lam
    nu <- params$nu
    res <- (nu / gamma(al)) * tau^(al - 1) * exp(-lam * tau)
    return(res)
  })
}

# Integral \int_0^tau kGamma(s)^2 ds
K00 <- function(params) {
  Vectorize(function(tau) {
    al <- params$al
    H2 <- 2 * al - 1
    lam <- params$lam
    nu <- params$nu
    prefactor <- (nu / gamma(al))^2
    diff_gamma <- gamma(H2) - gamma_inc(H2, 2 * lam * tau)
    res <- prefactor * ifelse(lam > 0, diff_gamma / (2 * lam)^H2, tau^H2 / H2)
    return(res)
  })
}

# Integral \int_0^tau kGamma(s + j * tau)^2 ds
Kjj <- function(params) {
  Vectorize(function(j,tau) {
        kGamma2 <- function(s) {
          (kGamma(params)(s + j * tau))^2
        }
        res <- integrate(kGamma2, lower = 0, upper = tau)$value
        return(res)
      },vectorize.args = c("j","tau"))
}

# Integral \int_0^tau kGamma(s + tau)^2 ds
K11 <- function(params)function(tau){Kjj(params)(1,tau)}

# Integral \int_0^tau kGamma(s) ds
K0 <- function(params) {
  Vectorize(function(tau) {
    al <- params$al
    lam <- params$lam
    nu <- params$nu
    prefactor <- (nu / gamma(al)) / (lam^al)
    bkt <- gamma(al) - gamma_inc(al, lam * tau)
    res2 <- (nu / gamma(al)) / al * tau^(al)
    res <- ifelse(lam > 0, prefactor * bkt, res2)
    return(res)
  })
}

# Integral \int_0^tau kGamma(s) kGamma(s + tau) ds
K01 <- function(params)function(t){
  
  kp <- kGamma(params)
  integr <- function(s){kp(s)*kp(s+t)}
  res <- integrate(integr, lower=0,upper=t)$value
  return(res)
  
}

# Integral \int_0^tau kGamma(s + j * tau) ds
Kj <- function(params)function(j)function(dt){
  
  al <- params$al
  H <- al - 1/2
  lam <- params$lam
  nu <- params$nu
  
  prefactor <- nu/(lam^al*gamma(al))
  bkt <- gamma_inc(al,j*lam*dt)- gamma_inc(al,(j+1)*lam*dt)
  res2 <- nu/gamma(al+1)*dt^al*((j+1)^al-j^al)
  res <- ifelse(lam>0,prefactor * bkt,res2)
  return(res)
  
}

# Integral \int_0^tau kGamma(s + tau) ds
K1 <- function(params)Kj(params)(1)

# Resolvent kernel of kGamma^2
bigK <- function(params) {
  function(tau) {
    al <- params$al
    H2 <- 2 * al - 1
    lam <- params$lam
    nu <- params$nu
    nuHat2 <- nu^2 * gamma(H2) / gamma(al)^2
    res <- nuHat2 * exp(-2 * lam * tau) * tau^(H2 - 1) * mlf(nuHat2 * tau^H2, H2, H2)
    return(res)
  }
}



# Integral \int_0^tau bigK(s) ds
bigK0 <- function(params) {
  Vectorize(function(tau) {
    bigKp <- bigK(params)
    integ <- function(s) {
      bigKp(s)
    }
    res <- ifelse(tau > 0, integrate(integ, lower = 0, upper = tau)$value, 0)
    return(res)
  })
}
