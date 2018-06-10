## Simulation 

library(tidyverse); library(rjags); library(jagsUI)

# Simulation of the imputing process ------------------------------------------

# 1. Construct the completed dataset
n <- 1000
x2 <- seq(from = 0, to = 1, length = 1000)
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10
f <- f2(x2)
e <- rnorm(n, 0, 2)
y <- f + e
dat1 <- data.frame(y = y, x = x2, f = f)
rm(list=ls()[!str_detect(ls(), 'dat|n')])
ggplot(dat1, aes(y=y, x=x)) + geom_point()

# 2. Remove some points to get missing dataset
missrate <- 0.05
na.index <- sort(sample(1:n, missrate*n, replace = FALSE))
dat2 <- dat1 %>% mutate(ymiss = replace(y, na.index, NA), 
                        fmiss = replace(f, na.index, NA),
                        IsNa  = as.factor(ifelse(is.na(ymiss), 1, 0)))
ggplot(dat2) + 
  geom_point(aes(y = y, x = x, col = IsNa)) + 
  geom_line(aes(y = f, x = x), size = 1) + 
  scale_color_manual(values = c('grey60', 'red')) + 
  theme_bw()

# 3-1. Imputation - random walk
s1 <- sd(dat2$fmiss, na.rm=T)
s2 <- sd(-diff(dat2$fmiss), na.rm = TRUE)  

dat.use <- list(x = dat2$fmiss, n = nrow(dat2), sigma = s2) 
initial <- NULL

model <- jags.model(file     = 'similation_model_1.txt',
                    data     = dat.use, 
                    inits    = initial,
                    n.chains = 3, 
                    n.adapt  = 1000)
update(model, n.iter = 2000)
sample <- coda.samples(model, 
                       variable.names = paste('x[', na.index, ']', sep=''), 
                       n.iter = 2000)

res.m <- summary(sample)

impute.median <- res.m$quantiles[,3]
truevalue <- dat2$f[na.index]
summary(truevalue - impute.median)

# 3-2. Imputation - measurement error
getmu <- function(v) {
  n <- length(v)
  v1 <- v2 <- v
  for (i in 2:n) { v1[i] <- ifelse(is.na(v1[i]), v1[i-1], v1[i]) }
  for (i in (n-1):1) { v2[i] <- ifelse(is.na(v2[i]), v2[i+1], v2[i]) }
  vnew <- apply(data.frame(v1, v2), 1, mean, na.rm=T)
  return(vnew)
}
s2 <- sd(-diff(dat2$fmiss), na.rm = TRUE) 
mu <- getmu(dat2$fmiss)
dat.use <- list(x = dat2$fmiss, n = nrow(dat2), mu = mu, sigma = s2) 
initial <- NULL

model <- jags.model(file     = 'similation_model_2.txt',
                    data     = dat.use, 
                    inits    = initial,
                    n.chains = 3, 
                    n.adapt  = 1000)
update(model, n.iter = 2000)
sample <- coda.samples(model, 
                       variable.names = paste('x[', na.index, ']', sep=''), 
                       n.iter = 2000)

res.m <- summary(sample)

impute.median <- res.m$quantiles[,3]
truevalue <- dat2$f[na.index]
summary(truevalue - impute.median)

# Simulation of the estimation process -------------------------------------

## Generating data
n <- 1000; k <- 5; alpha <- 0.1; beta <- 0.03; 
g1 <- rep(1/k, k)
g2 <- c(0.05, 0.1, 0.1, 0.5, 0.25)

b1 <- beta * g1
b2 <- beta * g2

xvct <- sin(1.2*seq(0, 50, length = n*k)) * 10 + rnorm(n*k, mean = 50, sd = 10)
xmat <- matrix(xvct, nrow = n, byrow = TRUE)

xmean1 <- as.vector(xmat %*% g1)
xmean2 <- as.vector(xmat %*% g2)

g <- g2
xmean <- xmean2

theta <- beta * xmean
mu <- exp(theta)
y  <- rpois(n, mu)

## Explore the simulated dataset
plot(xmean, theta)
plot(xmean, log(y))

## Generalized linear model
dat1 <- data.frame(y, xmean, xmat)
glm(y ~ xmean, data = dat1, family = poisson(link = log))
glm(y ~ X1 + X2 + X3 + X4 + X5, data = dat1, poisson(link = log))

## Bayesian model
hmat <- matrix(1:(n*k), nrow = n, ncol = k, byrow = TRUE)
nhour <- n*k

dat.use <- list(y = y, x1 = xvct, hmat = hmat, nhour = nhour, N = n, K = k) 
params  <- c('want')

out <- jagsUI::jags(data = dat.use, inits = NULL, parameters.to.save = params,
                    model.file = "simulation_model_3_f1p1.txt",
                    n.chains = 3, n.adapt = 1000, n.iter = 10000, n.burnin = 2000, n.thin = 5,
                    parallel = TRUE, n.cores = 3, verbose = T)
summary(out)
out$Rhat
out$summary
out$samples

## Compare the true weighted coeffecients with the estimated values
want.true <- beta * g
out$mean$want
out$q50$want

## GET the posterior distribution of weights & compare with the true gamma (g1 or g2)
postsample <- out$sims.list$want %>% as.tibble()
colnames(postsample) <- paste('w', 1:k, sep='')

weights <- t(apply(postsample, 1, prop.table))%>% as.tibble()
colnames(weights) <- paste('w', 1:k, sep='')

# the mean of posterior distribution of weights (gamma)
colMeans(weights)
# the median of posterior distribution of weights (gamma)
apply(weights, 2, median)
# posterior density of weights (gamma)
psts <- gather(weights, key = wname, value = weights)
trgm <- tibble(wname = paste('w', 1:5, sep=''), g)

ggplot(data = psts) + 
  geom_density(aes(x = weights), size = 0.5, col = 'darkblue', fill = 'lightblue') + 
  geom_vline(data = trgm, aes(xintercept = g), color="blue", linetype="dashed", size=1) +
  facet_grid(wname~., scales = 'free_y') + 
  scale_y_continuous(name = '') + 
  scale_x_continuous(name = 'gamma', limits = c(-0.01, 0.51), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) + 
  theme_bw()