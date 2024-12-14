library(extraDistr)
library(rjags)
library(coda)
library(ggplot2)
library(dplyr)

##############################################
# Slice Sampling Algorithm
##############################################

# Doubling Procedure
doubling <- function( f, x0, y, w, p ){
  U <- runif( 1 )
  L <- x0 - w * U
  R <- L + w
  K <- p
  while ( K > 0 && ( y < f(L) || y< f(R) ) ) {
    V <- runif( 1, 0, 1 )
    if( V < 1/2 ){
      L <- L - ( R - L )
    }else{
      R <- R + ( R - L )
    }
    K <- K-1
  }
  return( c(L, R) )
}

# Test Acceptance Procedure
test.acceptance <- function( f, x0, x1, y, w, L, R ){
  L.hat <- L
  R.hat <- R
  D = FALSE
  while ( R.hat  - L.hat > 1.1 * w ) {
    M <- ( L.hat + R.hat ) / 2
    if ( ( x0 < M && x1 >= M )  || ( x0 >= M && x1 < M ) ){
      D <- TRUE
    }
    if ( x1 < M ){
      R.hat <- M
    }else{
      L.hat <- M
    }
    if ( D &&  y >= f( L.hat ) && y >= f( R.hat ) ){
      return( FALSE )
    }
  }
  return( TRUE )
}

# Shrinkage Procedure
shrinkage <- function( f, x0, L, R, y, w ){
  L.bar <- L
  R.bar <- R
  while(TRUE) {
    U <- runif(1)
    x1 <- L.bar + U * ( R.bar - L.bar )
    if ( y < f( x1 ) && test.acceptance( f, x0, x1, y, w, L.bar, R.bar )) {
      return( x1 )
    }
    if ( x1 < x0 ){
      L.bar <- x1
    }else{
      R.bar <- x1
    }
  }
}

# Slice Sampler
slice.sampling <- function( f, x0, w, p, n.samples ){
  samples <- numeric(n.samples)
  samples[1] <- x0  
  for ( i in 2:n.samples ){
    y <- runif( 1, 0, f( x0 ) )  
    L_R <- doubling( f, x0, y, w, p )
    L <- L_R[1]
    R <- L_R[2]
    x1 <- shrinkage( f, x0, L, R, y, w )
    samples[i] <- x1
    x0 <- x1  
  }
  return ( samples )
}




##############################################
# Metric Functions
##############################################

# L2 Distance
l2.distance<- function(samples, theoretical.cdf, l.limit, u.limit){
  empirical.cdf <- ecdf(samples)
  x.vals <- seq(l.limit, u.limit, length.out = 1000)
  empirical.values <- empirical.cdf(x.vals)
  theoretical.values <- theoretical.cdf(x.vals)
  squared.diff <- (empirical.values - theoretical.values)^2
  L2.distance <- sum(squared.diff) * (x.vals[2] - x.vals[1])
  return(L2.distance)
}

# CDF Plot Function
cdf.plot <- function(samples, theoretical.cdf, l.limit, u.limit){
  empirical.cdf <- ecdf(samples)
  x.vals <- seq(l.limit, u.limit, length.out = 1000)
  empirical.values <- empirical.cdf(x.vals)
  theoretical.values <- theoretical.cdf(x.vals)
  plot.data <- data.frame(
    x = rep(x.vals, 2),
    CDF = c(theoretical.values, empirical.values),
    Type = rep(c("Theoretical CDF", "Empirical CDF"),each = length(x.vals))
  )
  p<-ggplot(plot.data, aes(x = x, y = CDF, color = Type)) +
    geom_line(linewidth = 1.2) +
    labs(
      x = "x",
      y = "CDF",
      color = "CDF Type"
    ) +
    scale_color_manual(values = c("black", "red")) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
  return(plot = p)
}

# PDF Plot Function
pdf.plot <- function(sample.density, true.density, x.vals){
  true.density.data <- data.frame(
    x = x.vals, 
    density = true.density,
    Type = "True Density"
  )
  sample.density.data <- data.frame(
    x = sample.density$x, 
    density = sample.density$y,
    Type = "Slice Sample Density"
  )
  hist.data <- data.frame(x = samples)
  p <- ggplot() +
    geom_histogram(
      data = hist.data, aes(x = x, y = ..density..), 
      bins = 30, fill = "gray", color = "black", alpha = 0.6
    ) +
    geom_line(
      data = true.density.data, aes(x = x, y = density, color = Type), 
      size = 1
    ) +
    geom_line(
      data = sample.density.data, aes(x = x, y = density, color = Type), 
      linetype = "dashed", size = 1
    ) +
    labs(
      x = "x",
      y = "Density",
      color = "Legend"
    ) +
    scale_color_manual(values = c("red", "black")) +
    theme_bw()+
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )+
    coord_cartesian(xlim = c(min(x.vals), max(x.vals)))
  return(plot = p)
}

# ACF Plot
acf.plot <- function(samples){
  acf.result <- acf(samples, plot = FALSE)
  acf.data <- data.frame(
    Lag = acf.result$lag[,1,1],         
    Autocorrelation = acf.result$acf[,1,1]
  )
  p<-ggplot(acf.data, aes(x = Lag, y = Autocorrelation)) +
    geom_bar(stat = "identity", fill = "gray", color = "black", width = 0.1) +
    geom_hline(yintercept = c(-1.96 / sqrt(length(samples)), 
                              1.96 / sqrt(length(samples))), 
               linetype = "dashed", color = "red") +
    theme_bw() +
    labs(x = "Lag", y = "Autocorrelation")
  return(plot=p)
}

# Posterior Summary and ESS
post.ess<-function(samples){
  mcmc.samples <- as.mcmc(samples)
  s.summary<-summary(mcmc.samples)
  s.size<-effectiveSize(mcmc.samples)
  return(list(post = s.summary, ess = s.size))
}




##############################################
# Sampling
##############################################

# Sampling From Beta(2,6)
alpha <- 2
beta <- 6
target.density <- function(x) {
  if (x < 0 || x > 1) return(0) 
  dbeta(x, shape1 = alpha, shape2 = beta)
}
theoretical.cdf <- function(x) pbeta(x, shape1 = alpha, shape2 = beta)
x0 <- 0.5
w <- 0.1
p <- 10
n.samples <- 5000
samples = slice.sampling(target.density, x0, w, p, n.samples)
mean.samples <- mean(samples)
sd.samples <- sd(samples)
#statistics and l2 distance
l.limit <- 0
u.limit <- 1
l2.samples <- l2.distance(samples, theoretical.cdf, 0, 1)
cat("l2: ", l2.samples, "\n")
#cdf plot
cdf.samples <- cdf.plot(samples, theoretical.cdf, 0, 1)
print(cdf.samples)
#pdf plot
x.vals <- seq(l.limit, u.limit, length.out = 1000)
true.density <- dbeta(x.vals, shape1 = 2, shape2 = 6)
sample.density <- density(samples, from = l.limit, to = u.limit)
pdf.samples<-pdf.plot(sample.density, true.density, x.vals)
print(pdf.samples)
#acf plot
acf.samples<-acf.plot(samples)
print(acf.samples)
#posterior summary ess
post.ess.samples<-post.ess(samples)
print(post.ess.samples)
#theoretical statistics
cat("t mean", (alpha/(alpha+beta)),"\n")
cat("t sd", sqrt( ( alpha*beta ) / ( (alpha+beta)^2 * (alpha + beta +1) ) ),"\n")




# Sampling From Wald(3,2)
mu<-3
lambda <- 2
target.density <- function(x) {
  if (x <= 0) return(0)
  sqrt(lambda / (2 * pi * x^3)) * exp(-lambda * (x - mu)^2 / (2 * mu^2 * x))
}
theoretical.cdf <- function(x) pwald(x, mu, lambda)
x0 <- 3
w <- 0.1
p <- 10
n.samples <- 5000
samples = slice.sampling(target.density, x0, w, p, n.samples)
mean.samples <- mean(samples)
sd.samples <- sd(samples)
#statistics and l2 distance
l.limit <- 0
u.limit <- 15
l2.samples <- l2.distance(samples, theoretical.cdf, 0, 15)
cat("l2: ", l2.samples, "\n")
#cdf plot
cdf.samples <- cdf.plot(samples, theoretical.cdf, 0, 15)
print(cdf.samples)
#pdf plot
x.vals <- seq(l.limit, u.limit, length.out = 1000)
true.density <- dwald(x.vals, mu = 3, lambda = 2)
samples <- samples[samples >= l.limit & samples <= u.limit]
sample.density <- density(samples, from = l.limit, to = u.limit)
pdf.samples<-pdf.plot(sample.density, true.density, x.vals)
print(pdf.samples)
#acf plot
acf.samples<-acf.plot(samples)
print(acf.samples)
#posterior summary ess
post.ess.samples<-post.ess(samples)
print(post.ess.samples)
#theoretical statistics
cat("t mean: ", mu, "\n")
cat("t sd: ", sqrt(mu^3/lambda), "\n")




# Sampling From Sinu-Exp
target.density <- function(x) {
  ifelse(x <= 0, 0, abs(sin(x)) * exp(-x / 4))
}
normalizing.constant <- integrate(
  function(x) abs(sin(x)) * exp(-x / 4), 
  lower = 0, upper = Inf
)$value
theoretical.pdf <- function(x) {
  ifelse(x <= 0, 0, (abs(sin(x)) * exp(-x / 4)) / normalizing.constant)
}
theoretical.cdf <- function(x) {
  sapply(x, function(xi) {
    if (xi <= 0) {
      0
    } else {
      integrate(theoretical.pdf, lower = 0, upper = xi)$value
    }
  })
}
x0 <- 0.5
w <- 0.1
p <- 10
n.samples <- 5000
samples = slice.sampling(target.density, x0, w, p, n.samples)
mean.samples <- mean(samples)
sd.samples <- sd(samples)
#statistics and l2 distance
l.limit <- 0
u.limit <- 15
l2.samples <- l2.distance(samples, theoretical.cdf, 0, 15)
cat("l2: ", l2.samples, "\n")
#cdf plot
cdf.samples <- cdf.plot(samples, theoretical.cdf, 0, 15)
print(cdf.samples)
#pdf plot
x.vals <- seq(l.limit, u.limit, length.out = 1000)
true.density <- theoretical.pdf(x.vals)
samples <- samples[samples >= l.limit & samples <= u.limit]
sample.density <- density(samples, from = l.limit, to = u.limit)
pdf.samples<-pdf.plot(sample.density, true.density, x.vals)
print(pdf.samples)
#acf plot
acf.samples<-acf.plot(samples)
print(acf.samples)
#posterior summary ess
post.ess.samples<-post.ess(samples)
print(post.ess.samples)
#theoretical statistics
t.mean <- integrate(
  function(x) x * theoretical.pdf(x),
  lower = 0,
  upper = Inf
)$value
second.moment <- integrate(
  function(x) x^2 * theoretical.pdf(x),
  lower = 0,
  upper = Inf,
  subdivisions = 1000
)$value
t.variance <- second.moment - (theoretical.mean)^2
t.sd <- sqrt(theoretical.variance)
cat("t mean: ", t.mean, "\n")
cat("t sd: ", t.sd, "\n")




# Sampling From Bimodal-GM
w1 <- 0.6
mu1 <- 0
sigma1 <- 1
w2 <- 0.4
mu2 <- 5
sigma2 <- 1.5
gaussian.pdf <- function(x, mu, sigma) {
  (1 / (sqrt(2 * pi) * sigma)) * exp(-0.5 * ((x - mu) / sigma)^2)
}
target.density<- function(x) {
  w1 * gaussian.pdf(x, mu1, sigma1) + w2 * gaussian.pdf(x, mu2, sigma2)
}
gaussian.cdf <- function(x, mu, sigma) {
  pnorm(x, mean = mu, sd = sigma)
}
theoretical.cdf <- function(x) {
  w1 * gaussian.cdf(x, mu1, sigma1) + w2 * gaussian.cdf(x, mu2, sigma2)
}
x0 <- 0.5
w <- 0.1
p <- 10
n.samples <- 5000
samples = slice.sampling(target.density, x0, w, p, n.samples)
mean.samples <- mean(samples)
sd.samples <- sd(samples)
#statistics and l2 distance
l.limit <- -10
u.limit <- 15
l2.samples <- l2.distance(samples, theoretical.cdf, -10, 15)
cat("l2: ", l2.samples, "\n")
#cdf plot
cdf.samples <- cdf.plot(samples, theoretical.cdf, -10, 15)
print(cdf.samples)
#pdf plot
x.vals <- seq(l.limit, u.limit, length.out = 1000)
true.density <- target.density(x.vals)
samples <- samples[samples >= l.limit & samples <= u.limit]
sample.density <- density(samples, from = l.limit, to = u.limit)
pdf.samples<-pdf.plot(sample.density, true.density, x.vals)
print(pdf.samples)
#acf plot
acf.samples<-acf.plot(samples)
print(acf.samples)
#posterior summary ess
post.ess.samples<-post.ess(samples)
print(post.ess.samples)
#theoretical statistics
t.mean <- w1 * mu1 + w2 * mu2
t.variance <- w1 * (sigma1^2 + mu1^2) + w2 * (sigma2^2 + mu2^2) - gmm.mean^2
t.sd <- sqrt(gmm.variance)
cat("t mean: ", t.mean, "\n")
cat("t sd: ", t.sd, "\n")






