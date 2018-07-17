## Simulation 
## 2018-07-15

setwd('E:/GitHub/Paper-relative-coefficients-of-hourly-PM2.5')

# Scenario 1-1 --------------------------------------------------------------
# Exposures are independent, coefficients are independent.

## generate data
mean24 <- apply(part1[,str_which(colnames(part1), 'Hour')], 2, function(x) {mean(x, na.rm=TRUE)})
sd24   <- apply(part1[,str_which(colnames(part1), 'Hour')], 2, function(x) {sd(x, na.rm=TRUE)})

beta1 <- rnorm(24, 0.5, 0.1)

bstar <- 0.03
w <- rep(1/24, 24) + 0.002 * c(-((24/2):1), 1:(24/2))
beta2 <- bstar * w

beta3 <- beta2 * 10
beta4 <- beta2 * 100

generate.data <- function(n, beta, mean24, sd24) {
  xmat <- matrix(NA, nrow = n, ncol = 24)
  for (k in 1:24) {
    xmat[,k] <- rnorm(n = n, mean = mean24[k], sd = sd24[k] - 40)
  }
  meanf <- xmat %*% beta
  y <- meanf + rnorm(n, 0.1, 0.1)
  # mu <- exp(meanf)
  # y  <- rpois(n, mu)
  dat <- data.frame(y, xmat)
  return(dat)
}

dat1 <- generate.data(n = 400, beta = beta1, mean24, sd24)
dat2 <- generate.data(n = 400, beta = beta2, mean24, sd24)
dat3 <- generate.data(n = 400, beta = beta3, mean24, sd24)
dat4 <- generate.data(n = 400, beta = beta4, mean24, sd24)

dat2.1200 <- generate.data(n = 1200, beta = beta2, mean24, sd24)
dat4.1200 <- generate.data(n = 1200, beta = beta4, mean24, sd24)

## Frequentist model

xnam <- paste0("X", 1:24)
fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"), "-1"))

m1 <- glm(fmla, family = gaussian, dat1)
m2 <- glm(fmla, family = gaussian, dat2)
m3 <- glm(fmla, family = gaussian, dat3)
m4 <- glm(fmla, family = gaussian, dat4)

m2.1200 <- glm(fmla, family = gaussian, dat2.1200)

round(m1$coefficients, 5)
round(m2$coefficients, 5)
round(m3$coefficients, 5)
round(m4$coefficients, 5)


## Bayesian model
dat <- dat4.1200
n <- nrow(dat)
hmat <- matrix(1:(n*24), nrow = n, ncol = 24, byrow = TRUE)
xvct <- dat[,-1] %>%
  mutate(date = 1:nrow(.)) %>% 
  gather(key = hour, value = pm, -date) %>% 
  arrange(date) %>% 
  select(pm) %>% 
  as.vector

dat.use <- list(y = dat$y, xvct = xvct, hmat = hmat, N = n, K = 24) 
params  <- c('beta')

out <- jagsUI::jags(data = dat.use, inits = NULL, parameters.to.save = params,
                    model.file = "simulation_model_3_f2p5.txt",
                    n.chains = 3, n.adapt = 1000, n.iter = 10000, n.burnin = 2000, n.thin = 5,
                    parallel = TRUE, n.cores = 3, verbose = T)
out$mean
out$q50


# Scenario 1-2 ------------------------------------------------------------
beta1 <- rnorm(24, 0.002, 0.0001)
generate.data <- function(n, beta, mean24, sd24, family = 'gaussian') {
  xmat <- matrix(NA, nrow = n, ncol = 24)
  for (k in 1:24) {
    xmat[,k] <- rnorm(n = n, mean = mean24[k], sd = sd24[k] - 40)
  }
  meanf <- xmat %*% beta
  if (family == 'gaussian')
    y <- meanf + rnorm(n, 0.1, 0.1)
  if (family == 'poisson')
    y  <- rpois(n, exp(meanf))
  dat <- data.frame(y, xmat)
  return(dat)
}

dat1 <- generate.data(400, beta1, mean24, sd24, 'poisson')
dat1.1200 <- generate.data(1200, beta1, mean24, sd24, 'poisson')

## Frequentist model

xnam <- paste0("X", 1:24)
fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"), "-1"))

m1 <- glm(fmla, family = poisson, dat1)
m1.1200 <- glm(fmla, family = poisson, dat1.1200)

round(m1$coefficients, 5)
round(m1.1200$coefficients, 5)



## Bayesian model
dat <- dat4.1200
n <- nrow(dat)
hmat <- matrix(1:(n*24), nrow = n, ncol = 24, byrow = TRUE)
xvct <- dat[,-1] %>%
  mutate(date = 1:nrow(.)) %>% 
  gather(key = hour, value = pm, -date) %>% 
  arrange(date) %>% 
  select(pm) %>% 
  as.vector

dat.use <- list(y = dat$y, xvct = xvct, hmat = hmat, N = n, K = 24) 
params  <- c('beta')

out <- jagsUI::jags(data = dat.use, inits = NULL, parameters.to.save = params,
                    model.file = "simulation_model_3_f2p5.txt",
                    n.chains = 3, n.adapt = 1000, n.iter = 10000, n.burnin = 2000, n.thin = 5,
                    parallel = TRUE, n.cores = 3, verbose = T)
out$mean
out$q50

# Scenario 2-1 --------------------------------------------------------------
beta1 <- rnorm(24, 0.06, 0.01)
beta2 <- rnorm(24, 0.6, 0.01)

xmat <- part1 %>% 
  filter(city.x == 5101) %>% 
  select(str_which(colnames(.), 'Hour')) %>% 
  drop_na()
meanf <- as.matrix(xmat) %*% beta1
y <- meanf + rnorm(length(meanf), 1, 0.02)

dat <- data.frame(y, xmat)

xnam <- paste0("Hour", 0:23)
fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"), "-1"))

m1 <- glm(fmla, family = gaussian, dat)
m1$coefficients

## Bayesian model
n <- nrow(dat)
hmat <- matrix(1:(n*24), nrow = n, ncol = 24, byrow = TRUE)
xvct <- dat[,-1] %>%
  mutate(date = 1:nrow(.)) %>% 
  gather(key = hour, value = pm, -date) %>% 
  arrange(date) %>% 
  select(pm) %>% 
  as.vector

dat.use <- list(y = dat$y, xvct = xvct, hmat = hmat, N = n, K = 24) 
params  <- c('beta')

out <- jagsUI::jags(data = dat.use, inits = NULL, parameters.to.save = params,
                    model.file = "simulation_model_3_f2p5.txt",
                    n.chains = 3, n.adapt = 1000, n.iter = 10000, n.burnin = 2000, n.thin = 5,
                    parallel = TRUE, n.cores = 3, verbose = T)
out$mean$beta
out$q50$beta


# Scenario 2-2 ------------------------------------------------------------
beta <- sin(seq(0.1, 3, 0.10))[4:27] * 0.02
plot(beta)

hpm <- part1 %>% 
  filter(city.x == 5101) %>% 
  select(str_which(colnames(.), 'Hour')) %>% 
  drop_na()
xmat <- rbind(hpm, hpm, hpm) + rnorm(nrow(hpm)*3, 2, 1)
plot(apply(xmat, 1, mean))

meanf <- as.matrix(xmat) %*% beta
y <- meanf + rnorm(length(meanf), 2, 0.05)

dat <- data.frame(y, xmat)
bstar <- 0.005
gamma <- beta/bstar

## Bayesian model
n <- nrow(dat)
hmat <- matrix(1:(n*24), nrow = n, ncol = 24, byrow = TRUE)
xvct <- dat[,-1] %>%
  mutate(date = 1:n) %>% 
  gather(key = hour, value = pm, -date) %>% 
  arrange(date) %>% 
  select(pm) %>% 
  as.vector

dat.use <- list(y = dat$y, xvct = xvct, hmat = hmat, N = n, K = 24) 
params  <- c('beta', 'bstar', 'gamma')

out <- jagsUI::jags(data = dat.use, inits = NULL, parameters.to.save = params,
                    model.file = "simulation_model_3_f2p6.txt",
                    n.chains = 3, n.adapt = 1000, n.iter = 10000, n.burnin = 2000, n.thin = 5,
                    parallel = TRUE, n.cores = 3, verbose = T)

beta.post <- out$mean$beta
plot(beta.post/sum(beta.post))

## constrained
out2 <- jagsUI::jags(data = dat.use, inits = NULL, parameters.to.save = params,
                    model.file = "simulation_model_3_f2p7.txt",
                    n.chains = 3, n.adapt = 1000, n.iter = 10000, n.burnin = 2000, n.thin = 5,
                    parallel = TRUE, n.cores = 3, verbose = T)

beta.post <- out2$mean$beta
plot(beta.post/sum(beta.post))
