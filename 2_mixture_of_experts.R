
# Mixture of experts







# Packages ----------------------------------------------------------------

library(rstan)
library(bayesplot)



# Read data ---------------------------------------------------------------

dat = read.csv("data/Transcript_m6A_15.csv")


# log2FoldChange = observed effect size
# lfc_SE = standard error
# mean_rel_location = moderator

# remove genes without m6A, i.e. position is NA
dat = dat[!is.na(dat$mean_rel_location) & !is.na(dat$log2FoldChange), ]


# only genes that lost all sites
pdf("output/overview_sites.pdf")
plot(n_sites_lost ~ n_sites_all, dat, 
     col = rgb(0.7*(dat$n_sites_lost >= dat$n_sites_all*0.8), 0, 0, 0.05),
     xlab = "number of sites",
     ylab = "number of sites lost",
     pch = 16)
abline(0, 1)
dev.off()
# dat = dat[dat$n_sites_lost == dat$n_sites_all, ] # 100% loss
dat = dat[dat$n_sites_lost >= dat$n_sites_all*0.8, ] # 80% loss

plot(mean_rel_location ~ n_sites_lost, dat)
plot(mean_rel_location_all ~ n_sites_lost, dat)



# scale predictor
# subtract 0.5
# divide by sd
sd_loc = sd(dat$mean_rel_location)
dat$scaled_mean_rel_location = (dat$mean_rel_location-0.5)/sd_loc
hist(dat$scaled_mean_rel_location, breaks = 50)
plot(log2FoldChange ~ mean_rel_location, dat)
plot(log2FoldChange_shrunk ~ mean_rel_location, dat)



# Prepare stan data -------------------------------------------------------

stan_dat <- list(
  N  = nrow(dat),
  # y  = as.numeric(dat$log2FoldChange),
  # se = as.numeric(dat$lfc_SE),
  y  = as.numeric(dat$log2FoldChange_shrunk),
  se = as.numeric(dat$lfc_SE_shrunk),
  x  = as.numeric(dat$scaled_mean_rel_location),
  prior_only = 0
)


# Test model --------------------------------------------------------------

# Compile the model
stan_model <- stan_model(file = "model.stan")

# Run prior predictive checks
fit_prior <- sampling(
  stan_model,
  data = within(stan_dat, prior_only <- 1),
  chains = 3, iter = 500, cores = 3
)
fit_prior





# Run model with rstan ----------------------------------------------------


# Fit the model
iter = 6000
fit <- sampling(
  stan_model,
  data = stan_dat,
  seed = 1,
  chains = 4,
  cores = 4,
  iter = iter,
  warmup = iter/2,
  thin = 6,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 10
  )
)
save(fit, file = "output/fit.Rdata")



# Check convergence -------------------------------------------------------

sumfit = summary(fit)
sumfit$summary[, c("n_eff", "Rhat")]

traceplot(fit)
# traceplot(fit, pars = c("b0", "z_b[1]", "z_b[2]", "tau_b"))
traceplot(fit, pars = c("b[1]", "b[2]"))

# parameters
summary(fit)

stan_dens(fit, pars = c("gamma10", "gamma11", "gamma30", "gamma31"))
stan_dens(fit, pars = c("sigma", "nu"))

stan_dens(fit, pars = c("gamma10", "gamma11", "gamma30", "gamma31"), separate_chains = T)

stan_scat(fit, pars = c("gamma10", "gamma11"))
stan_scat(fit, pars = c("z_mu[1]", "b[1]"))
stan_scat(fit, pars = c("gamma31", "b[2]"))
stan_scat(fit, pars = c("gamma11", "b[1]"))




# Posterior samples and x_grid --------------------------------------------


post <- rstan::extract(fit)
x_grid <- seq(min(dat$scaled_mean_rel_location), 
              max(dat$scaled_mean_rel_location), 
              length.out = 100)



# Plot probabilites -------------------------------------------------------

gamma10 <- post$gamma10
gamma11 <- post$gamma11
gamma30 <- post$gamma30
gamma31 <- post$gamma31


# initialize matrix to store probabilities
p_up <- matrix(NA, nrow = length(gamma10), ncol = length(x_grid))
p_null <- matrix(NA, nrow = length(gamma10), ncol = length(x_grid))
p_down <- matrix(NA, nrow = length(gamma10), ncol = length(x_grid))

for (i in seq_along(x_grid)) {
  xg <- x_grid[i]
  
  logit_p <- cbind(
    gamma10 + gamma11 * xg,  # down
    0,                        # null
    gamma30 + gamma31 * xg    # up
  )
  
  log_softmax_matrix <- function(M) {
    M_centered <- M - apply(M, 1, max)
    M_centered - log(rowSums(exp(M_centered)))
  }
  
  # log_softmax over columns
  log_p <- log_softmax_matrix(logit_p)  # pseudo function: do row-wise log_softmax
  p_down[, i] <- exp(log_p[,1])
  p_null[, i] <- exp(log_p[,2])
  p_up[, i]   <- exp(log_p[,3])
}

pdf("output/modelresults.pdf", 
    width = 10)
par(mfrow = c(1, 2))
plot(0.5 + sd_loc*x_grid,
     apply(p_up, 2, median),
     type = "l", lwd = 2,
     ylim = c(0, 1),
     ylab = "Probability",
     xlab = "m6A position",
     col = "darkred")
polygon(
  0.5 + sd_loc*c(x_grid, rev(x_grid)),
  c(apply(p_up, 2, quantile, 0.1),
    rev(apply(p_up, 2, quantile, 0.9))),
  col = rgb(1, 0.2, 0.2, 0.1),
  border = NA
)
lines(0.5 + sd_loc*x_grid,
     apply(p_down, 2, median),
     type = "l", lwd = 2,
     col = "darkblue")
polygon(
  0.5 + sd_loc*c(x_grid, rev(x_grid)),
  c(apply(p_down, 2, quantile, 0.1),
    rev(apply(p_down, 2, quantile, 0.9))),
  col = rgb(0.1, 0.2, 1, 0.2),
  border = NA
)
lines(0.5 + sd_loc*x_grid,
      apply(p_null, 2, median),
      type = "l", lwd = 2,
      col = "grey20")
polygon(
  0.5 + sd_loc*c(x_grid, rev(x_grid)),
  c(apply(p_null, 2, quantile, 0.1),
    rev(apply(p_null, 2, quantile, 0.9))),
  col = rgb(0, 0, 0, 0.2),
  border = NA
)

# extract percentages
extract_percentage = function(x) round(100*quantile(x, p = c(0.1, 0.5, 0.9)), 2)

# across range of m6A positions
extract_percentage(p_null[, ])
extract_percentage(p_down[, ])
extract_percentage(p_up[, ])

# at position 0.5 -> for paper
extract_percentage(p_null[, 50])
extract_percentage(p_down[, 50])
extract_percentage(p_up[, 50])



# Write component probabilities to file -----------------------------------


# Extract values
vals_null <- extract_percentage(p_null[, 50])
vals_down <- extract_percentage(p_down[, 50])
vals_up   <- extract_percentage(p_up[, 50])

# Helper to format one vector
format_line <- function(name, vals) {
  paste0(
    name, ": ",
    paste0(names(vals), " = ", sprintf("%.2f", vals), collapse = " | ")
  )
}

# Create lines
lines <- c(
  format_line("Null", vals_null),
  format_line("Down", vals_down),
  format_line("Up", vals_up)
)

# Write to file
writeLines(lines, "output/component_probabilities.txt")





# Plot effects ------------------------------------------------------------

min_effect <- 0.5

# conditional effect
mu_down_grid <- sapply(x_grid, function(xg) {
  - min_effect - log1p(exp(post$mu0 + post$tau_mu * post$z_mu[,1] + post$b[,1] * xg))
})

mu_up_grid <- sapply(x_grid, function(xg) {
  min_effect + log1p(exp(post$mu0 + post$tau_mu * post$z_mu[,2] + post$b[,2] * xg))
})


plot(0.5 + sd_loc*x_grid,
     apply(mu_up_grid, 2, median),
     type = "l", lwd = 2,
     ylim = c(-3, 3),
     # ylim = range(dat$log2FoldChange_shrunk),
     col = "darkred",
     ylab = "log2FC",
     xlab = "m6A position")
polygon(
  0.5 + sd_loc*c(x_grid, rev(x_grid)),
  c(apply(mu_up_grid, 2, quantile, 0.1),
    rev(apply(mu_up_grid, 2, quantile, 0.9))),
  col = rgb(1, 0.2, 0.2, 0.1),
  border = NA)

polygon(
  0.5 + sd_loc*c(x_grid, rev(x_grid)),
  c(apply(mu_down_grid, 2, quantile, 0.1),
    rev(apply(mu_down_grid, 2, quantile, 0.9))),
  col = rgb(0.1, 0.2, 1, 0.2),
  border = NA)  
lines(0.5 + sd_loc*x_grid,
      apply(mu_down_grid, 2, median),
      col = "darkblue",
      lwd = 2)

polygon(
  0.5 + sd_loc*c(x_grid, rev(x_grid)),
  c(rep(quantile(post$Intercept_mu2, 0.1), length(x_grid)),
    rep(quantile(post$Intercept_mu2, 0.9), length(x_grid))),
  col = rgb(0, 0, 0, 0.2),
  border = NA)
abline(h = median(post$Intercept_mu2), 
       col = "darkgray", 
       lwd = 2)

w <- 1 / dat$lfc_SE
w_scaled <- w / max(w)
points(log2FoldChange_shrunk ~ mean_rel_location, dat,
       pch = 16, 
       # cex = 0.5, 
       cex = 0.2 + 1.2 * w_scaled,   # controlled range
       col = rgb(0, 0, 0, 0.1)
       # col = rgb(0, 0, 0, 0.01/dat$lfc_SE)
)

abline(h = c(-min_effect, min_effect), lty = 2, lwd = 2, col = "red")
dev.off()





# Posterior predictive checks ---------------------------------------------

# extract
# y_obs <- dat$log2FoldChange
y_obs <- dat$log2FoldChange_shrunk
x_obs <- dat$mean_rel_location
y_rep <- post$y_rep   # iterations x N
z_rep <- post$z_rep   # iterations x N, if kept

dim(y_rep)
length(y_obs)
length(x_obs)

# use a subset of posterior draws for plotting
set.seed(1)
draws_plot <- sample(seq_len(nrow(y_rep)), min(500, nrow(y_rep)))

pdf("output/PPC_nufitted.pdf")

# density overlay
ppc_dens_overlay(y_obs, y_rep[draws_plot, ])

# summary stats
ppc_stat(y_obs, y_rep, stat = "median")
ppc_stat(y_obs, y_rep, stat = "sd")

# tail behavior
ppc_stat(y_obs, y_rep, stat = function(z) mean(abs(z) > 0.5))
ppc_stat(y_obs, y_rep, stat = function(z) mean(abs(z) > 1))
ppc_stat(y_obs, y_rep, stat = function(z) mean(abs(z) > 2))

# directional mass
ppc_stat(y_obs, y_rep, stat = function(z) mean(z < -0.5))
ppc_stat(y_obs, y_rep, stat = function(z) mean(abs(z) <= 0.5))
ppc_stat(y_obs, y_rep, stat = function(z) mean(z > 0.5))


dev.off()






