library(bayesplot)

### Reference:
### https://cran.r-project.org/web/packages/bayesplot/vignettes/plotting-mcmc-draws.html#
###

load('E:/GitHub/Paper-relative-coefficients-of-hourly-PM2.5/visualization/res0805.RData')
post.tmp <- as.array(MultipleCity.m2)
post.samples <- aperm(post.tmp, c(1,3,2))
dim(post.samples)

# 1. Posterior uncertainty intervals
color_scheme_set("red")
mcmc_intervals(post.samples, 
               pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", 
                        "beta[5]", "beta[6]", "beta[7]", "beta[8]", 
                        "beta[9]", "beta[10]", "beta[11]", "beta[12]"))

mcmc_areas(
  post.samples, 
  pars = c("beta[1]", "beta[6]", "beta[12]", "beta[18]"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

# 2. Univariate marginal posterior distributions
color_scheme_set("green")
mcmc_hist(post.samples, pars = c("beta[9]", "beta[12]"))

color_scheme_set("brightblue")
mcmc_hist_by_chain(post.samples, pars = c("beta[9]", "beta[12]"))

color_scheme_set("purple")
mcmc_dens(post.samples, pars = c("beta[9]", "beta[12]"))

mcmc_dens_overlay(post.samples, pars = c("beta[9]", "beta[12]"))

color_scheme_set("teal")
mcmc_violin(post.samples, pars = c("beta[9]", "beta[12]"), probs = c(0.1, 0.5, 0.9))


# 3. Bivariate plots
color_scheme_set("gray")
mcmc_scatter(post.samples, pars = c("beta[9]", "beta[12]"), 
             size = 1.5, alpha = 0.5)

mcmc_hex(post.samples, pars = c("beta[9]", "beta[12]"))

color_scheme_set("pink")
mcmc_pairs(post.samples, pars = c("beta[9]", "beta[12]", "beta[20]"),
           off_diag_args = list(size = 1.5))


# 4. Trace plots
color_scheme_set("blue")
mcmc_trace(post.samples, pars = c("beta[9]", "beta[12]", "beta[20]"))

color_scheme_set("mix-blue-red")
mcmc_trace(post.samples, pars = c("beta[9]", "beta[12]", "beta[20]"), 
           facet_args = list(ncol = 1, strip.position = "left"))

mcmc_trace_highlight(post.samples, pars = c("beta[9]"), highlight = 3)


### Reference:
### https://cran.r-project.org/web/packages/bayesplot/vignettes/visual-mcmc-diagnostics.html
### 

# 1. Diagnostics for the No-U-Turn Sampler
# 2. General MCMC diagnostics
# 2.1 Rhat: potential scale reduction statistic
rhats <- gelman.diag(MultipleCity.m2)$psrf[,1]
color_scheme_set("brightblue") 
mcmc_rhat(rhats) + yaxis_text(hjust = 1)

# 2.2 Effective sample size
Neff <- effectiveSize(MultipleCity.m2)
ratios <- Neff/(4000 * 3)
mcmc_neff(ratios, size = 2)

# 2.3 Autocorrelation
mcmc_acf(MultipleCity.m2, pars = "beta[1]", lags = 10)
mcmc_acf(MultipleCity.m2, pars = "beta[9]", lags = 10)


## Reference
## https://cran.r-project.org/web/packages/bayesplot/vignettes/graphical-ppcs.html
## 

## TO LEARN
