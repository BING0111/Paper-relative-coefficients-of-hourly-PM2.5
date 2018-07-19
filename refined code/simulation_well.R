## Simulation 
## 2018-07-15

setwd('E:/GitHub/Paper-relative-coefficients-of-hourly-PM2.5')
load("xmat.RData")

mean24 <- apply(xmat, 2, function(x) {mean(x, na.rm=TRUE)})
sd24   <- apply(xmat, 2, function(x) {sd(x, na.rm=TRUE) - 40})

bstar <- 0.03
w1 <- rnorm(24, 0.55, 0.13) # independent
w2 <- sin(seq(0.1, 3, 0.10))[4:27] # constrained
beta1 <- round(x = bstar * w1, digits = 4)
beta2 <- round(x = bstar * w2, digits = 4)

# function - generate data  -----------------------------------------------

generate.data <- 
  function(Xstr = 'independent', N = NULL, H = NULL, Xmean = NULL, Xsd = NULL, xmat = NULL, 
           bstr = 'independent', beta, 
           family = 'gaussian') {
    
    #### two kinds of xmat ####
    
    if (Xstr == 'independent') {
      
      if (is.null(N) || is.null(H) || is.null(Xmean) || is.null(Xsd))
        stop('N, H, Xmean, or Xsd is missing in this scenario.')
      
      ## generate xmat from normal distribution.
      xmat <- matrix(NA, nrow = N, ncol = H)
      for (i in 1:H) {
        xmat[,i] <- rnorm(n = N, mean = Xmean[i], sd = Xsd[i])
      }
    }
    
    if (Xstr == 'dependent') {
      
      if (is.null(xmat))
        stop('xmat is missing in this scenario')
        
      ## assign true air quality monitoring data as xmat.
      xmat <- as.matrix(xmat)
      H <- ncol(xmat)
      
      if (is.null(N)) { 
        N <- nrow(xmat) 
      } else {
        nx <- nrow(xmat)
        divisor <- floor(N/nx)
        if ( divisor== 0 ) { xmat <- xmat[1:(N-nx),] }
        if ( divisor > 0 ) { 
          xmat <- xmat[c(rep(1:nx, divisor), 1:(N - nx*divisor)),] 
          xmat <- xmat + rnorm(nrow(xmat)*ncol(xmat), 2, 1)}
        }
    }
    
    
    #### two kinds of beta structure ####
    
    # if (bstr == 'independent') {
    #   beta1 <- rnorm(24, 0.06, 0.01)
    #   beta2 <- rnorm(24, 0.6, 0.01)
    # }
    # 
    # if (bstr == 'constrained') {
    #   beta <- sin(seq(0.1, 3, 0.10))[4:27] * 0.02
    # }
    
    
    #### mean function ####
    
    meanf <- xmat %*% beta
    
    
    #### two kinds of family & link function ####
    
    if (family == 'gaussian')
      y <- meanf + rnorm(N, 0.1, 0.1)
    if (family == 'poisson')
      y  <- rpois(N, exp(meanf))
    
    
    #### the data we need ####
    
    dat <- data.frame(y, meanf, xmat)
    return(dat)
  }

# Scenario 1-1-1 --------------------------------------------------------------

## Exposures are independent; 
## coefficients are independent;
## family of y is gaussian.

dat111.400 <- generate.data(Xstr = 'independent', N = 400, H = 24, Xmean = mean24, Xsd = sd24,
                            bstr = 'independent', beta = beta1,
                            family = 'gaussian')

# Scenario 1-1-2 ----------------------------------------------------------

## Exposures are independent; 
## coefficients are independent;
## family of y is poisson.

dat112.400 <- generate.data(Xstr = 'independent', N = 400, H = 24, Xmean = mean24, Xsd = sd24,
                            bstr = 'independent', beta = beta1,
                            family = 'poisson')

# Scenario 1-2-1 ----------------------------------------------------------

## Exposures are independent; 
## coefficients are constrained;
## family of y is gaussian

dat121.400 <- generate.data(Xstr = 'independent', N = 400, H = 24, Xmean = mean24, Xsd = sd24,
                            bstr = 'constrained', beta = beta2,
                            family = 'gaussian')

# Scenario 1-2-2 ----------------------------------------------------------

## Exposures are independent; 
## coefficients are constrained;
## family of y is poisson.


dat122.400 <- generate.data(Xstr = 'independent', N = 400, H = 24, Xmean = mean24, Xsd = sd24,
                            bstr = 'constrained', beta = beta,
                            family = 'poisson')

# Scenario 2-1-1 ----------------------------------------------------------

## Exposures are dependent; 
## coefficients are independent;
## family of y is gaussian.

dat211.400 <- generate.data(Xstr = 'dependent', N = 400, H = 24, Xmean = mean24, Xsd = sd24,
                            bstr = 'independent', beta = beta,
                            family = 'gaussian')

# Scenario 2-1-2 ----------------------------------------------------------

## Exposures are dependent; 
## coefficients are independent;
## family of y is poisson.

dat212.400 <- generate.data(Xstr = 'dependent', N = 400, H = 24, Xmean = mean24, Xsd = sd24,
                            bstr = 'independent', beta = beta,
                            family = 'poisson')

# Scenario 2-2-1 ----------------------------------------------------------

## Exposures are dependent; 
## coefficients are constrained;
## family of y is gaussian.

dat221 <- generate.data(Xstr = 'dependent', xmat = xmat1,
                        bstr = 'constrained', beta = beta2,
                        family = 'gaussian')

# Scenario 2-2-2 ----------------------------------------------------------

## Exposures are dependent; 
## coefficients are constrained;
## family of y is poisson

dat222 <- generate.data(Xstr = 'dependent', xmat = xmat1,
                        bstr = 'constrained', beta = beta2 * 0.05,
                        family = 'poisson')


# Analyze -----------------------------------------------------------------

## Frequentist model
xnam <- paste0("Hour", 0:23)
fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"), "-1"))

m221 <- glm(fmla, family = gaussian, dat221)

summary(m221)
beta.hat <- m221$coefficients %>% round(., 5) %>% unname

coef.res <- data.frame(b = c(beta2, beta.hat), 
                       type = factor(rep(c('beta', 'beta.hat'), each=24)), 
                       hour = rep(1:24, 2))

ggplot(coef.res, aes(x = hour, y = b)) + 
  geom_point(aes(col = type, group = type), size=2) +
  ggtitle('True X, smoothed beta, normal family with Identity link') + 
  theme_bw()


## Bayesian model
dat <- dat221
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
                    model.file = "refined code/simulation_normal_2prior_RandomWalk_0.txt",
                    n.chains = 3, n.adapt = 1000, n.iter = 10000, n.burnin = 2000, n.thin = 5,
                    parallel = TRUE, n.cores = 3, verbose = T)
out$mean$beta
out$q50$beta
summary(out)


coef.res2 <- data.frame(b = c(beta2, out$q50$beta), 
                       type = factor(rep(c('beta', 'beta.hat'), each=24)), 
                       hour = rep(1:24, 2))

ggplot(coef.res2, aes(x = hour, y = b)) + 
  geom_point(aes(col = type, group = type), size=2) +
  ggtitle('data: true X, smoothed beta, normal family with identity link \nmodel: normal_2prior_RandomWalk_0') + 
  theme(plot.title = element_text(vjust = 0.1)) + 
  theme_bw()
