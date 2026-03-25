
data {
  int<lower=1> N;
  vector[N] y;                       // log2FC
  vector<lower=0>[N] se;             // standard errors
  vector[N] x;                       // m6A position
  int<lower=0, upper=1> prior_only;  // should the likelihood be ignored?
}


parameters {

  // ---- hierarchical magnitudes (log scale) ----
  real mu0;
  vector[2] z_mu;
  real<lower=0> tau_mu;

  // ---- structural null ----
  real Intercept_mu2;

  // ---- slopes (log scale) ----
  vector[2] b;

  // ---- noise ----
  real<lower=0> sigma;
  
  // ---- t degrees of freedom ----
  real<lower=2> nu;

  // ---- mixture-of-experts gating ----
  real gamma10;   // intercept: comp 1 vs null
  real gamma11;   // slope on x

  real gamma30;   // intercept: comp 3 vs null
  real gamma31;   // slope on x
}

transformed parameters {

  real lprior = 0;


  // ---- priors ----
  lprior += normal_lpdf(mu0 | 0, 0.4);
  lprior += normal_lpdf(z_mu | 0, 0.3);
  lprior += normal_lpdf(tau_mu | 0, 0.2)
            - normal_lccdf(0 | 0, 0.2);
            
  lprior += normal_lpdf(b | 0, 0.5);

  lprior += normal_lpdf(Intercept_mu2 | 0, 0.01);

  lprior += student_t_lpdf(sigma | 3, 0, 0.4)
            - student_t_lccdf(0 | 3, 0, 0.4);
            
  lprior += gamma_lpdf(nu | 2, 0.1);
  
  lprior += normal_lpdf(gamma10 | 0, 1);
  lprior += normal_lpdf(gamma11 | 0, 0.01);
  lprior += normal_lpdf(gamma30 | 0, 1);
  lprior += normal_lpdf(gamma31 | 0, 0.01);
}

model {
  if (!prior_only) {
    
    for (n in 1:N) {
      
      // gating (mixture-of-experts)
      vector[3] logit_p;
      array[3] real ps;
      vector[3] log_p;

     
      logit_p[1] = gamma10 + gamma11 * x[n];
      logit_p[2] = 0;
      logit_p[3] = gamma30 + gamma31 * x[n];
      
      log_p = log_softmax(logit_p);

      // compute mus on the fly (do not store)
      real min_effect = 0.5;
      real mu1 = - min_effect - log1p_exp(mu0 + tau_mu * z_mu[1] + b[1] * x[n]);
      real mu3 = min_effect + log1p_exp(mu0 + tau_mu * z_mu[2] + b[2] * x[n]);
      real mu2 = Intercept_mu2;
      
      
      // log-likelihood contributions
      ps[1] = log_p[1] + student_t_lpdf(y[n] | nu, mu1, sqrt(se[n]^2 + sigma^2));
      ps[2] = log_p[2] + student_t_lpdf(y[n] | nu, mu2, sqrt(se[n]^2 + sigma^2));
      ps[3] = log_p[3] + student_t_lpdf(y[n] | nu, mu3, sqrt(se[n]^2 + sigma^2));
      
      target += log_sum_exp(ps);
    }
  }

  target += lprior;
}

generated quantities {
  vector[N] r1;
  vector[N] r2;
  vector[N] r3;
  vector[N] y_rep;
  array[N] int z_rep;

  for (n in 1:N) {
    vector[3] logit_p;
    vector[3] log_p;
    vector[3] log_joint;
    vector[3] p;

    real min_effect = 0.5;
    real mu1 = -min_effect - log1p_exp(mu0 + tau_mu * z_mu[1] + b[1] * x[n]);
    real mu2 = Intercept_mu2;
    real mu3 =  min_effect + log1p_exp(mu0 + tau_mu * z_mu[2] + b[2] * x[n]);

    real scale1 = sqrt(se[n]^2 + sigma^2);
    real scale2 = sqrt(se[n]^2 + sigma^2);
    real scale3 = sqrt(se[n]^2 + sigma^2);

    logit_p[1] = gamma10 + gamma11 * x[n];
    logit_p[2] = 0;
    logit_p[3] = gamma30 + gamma31 * x[n];

    log_p = log_softmax(logit_p);
    p = softmax(logit_p);

    // posterior responsibilities for observed y[n]
    log_joint[1] = log_p[1] + student_t_lpdf(y[n] | nu, mu1, scale1);
    log_joint[2] = log_p[2] + student_t_lpdf(y[n] | nu, mu2, scale2);
    log_joint[3] = log_p[3] + student_t_lpdf(y[n] | nu, mu3, scale3);

    {
      vector[3] r = softmax(log_joint);
      r1[n] = r[1];
      r2[n] = r[2];
      r3[n] = r[3];
    }

    // posterior predictive draw
    z_rep[n] = categorical_rng(p);

    if (z_rep[n] == 1) {
      y_rep[n] = student_t_rng(nu, mu1, scale1);
    } else if (z_rep[n] == 2) {
      y_rep[n] = student_t_rng(nu, mu2, scale2);
    } else {
      y_rep[n] = student_t_rng(nu, mu3, scale3);
    }
  }
}