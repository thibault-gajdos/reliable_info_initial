//*  ---------------------------------------------------
//               FUNCTIONS
//*  ---------------------------------------------------

functions {

  // ===== Transform parameters =====
  vector transform_params(vector mu_pr, vector sigma_pr, row_vector param_raw_n) {
    vector[7] p;
    for (j in 1:5)
      p[j] = 2 * Phi_approx(mu_pr[j] + sigma_pr[j] * param_raw_n[j]);
    p[6]  = 10 * Phi_approx(mu_pr[6] + sigma_pr[6] * param_raw_n[6]);
    p[7]  = mu_pr[7] + sigma_pr[7] * param_raw_n[7];
    return p;
  }

  vector transform_group_means(vector mu_pr) {
    vector[7] g;
    for (j in 1:5)
      g[j] = 2 * Phi_approx(mu_pr[j]);
    g[6]  = 10 * Phi_approx(mu_pr[6]);
    g[7]  = mu_pr[7];
    return g;
  }

  // ===== Mapping =====
  real map_linear_logodds(real p, real alpha, real beta, int fix_distortion) {
    if (fix_distortion == 1) return logit(p);
    real pp = fmin(fmax(p, 1e-6), 1 - 1e-6);
    return alpha * logit(pp) + beta;
  }

  // ===== Sequential weight =====
  real get_weight(int s, vector w5, int fix_weights) {
    if (fix_weights == 1) return 1.0;
    if (s <= 5) return w5[s];
    return 1.0;
  }

  // ===== reduce_sum worker =====
  real partial_sum(array[] int slice_indices,
                   int start, int end,
                   vector mu_pr, vector sigma_pr,
                   array[] int Tsubj,
                   array[,] int sample,
                   array[,,] int color,
                   array[,,] real proba,
                   array[,] int choice,
                   matrix param_raw,
                   int fix_weights,
                   int fix_distortion) {

    real lp = 0;

    for (i in 1:size(slice_indices)) {
      int n = slice_indices[i];

      vector[7] params = transform_params(mu_pr, sigma_pr, param_raw[n]);

      vector[5] w5 = params[1:5];
      real alpha   = params[6];
      real beta    = params[7];

      for (t in 1:Tsubj[n]) {
        int S = sample[n, t];
        if (S < 1) continue;

        int ch = choice[n, t];
        if (ch < 1 || ch > 2) continue;

        vector[2] evidence = rep_vector(0.0, 2);

        for (s in 1:S) {
          real p = proba[n, t, s];
          int  c = color[n, t, s];

          if (c < 1 || c > 2) continue;
          if (p <= 0 || p >= 1) continue;

          real m  = map_linear_logodds(p, alpha, beta, fix_distortion);
          real ws = get_weight(s, w5, fix_weights);

          evidence[c] += ws * m;
        }

        real ev_diff = evidence[1] - evidence[2];
        lp += bernoulli_logit_lpmf(ch == 1 | ev_diff);
      }
    }

    return lp;
  }
}

//*  ---------------------------------------------------
//               DATA
//*  ---------------------------------------------------

data {
  int<lower=1> N;
  int<lower=1> T_max;
  int<lower=1> I_max;

  array[N] int<lower=1> Tsubj;
  array[N, T_max] int<lower=-1> sample;
  array[N, T_max, I_max] int<lower=-1, upper=2> color;
  array[N, T_max, I_max] real<lower=-1, upper=1> proba;

  array[N, T_max] int<lower=-1, upper=2> choice;

  int<lower=1> grainsize;

  // Model variant flags
  int<lower=0, upper=1> fix_sequential_weights;     // 1 = w1-w5 all fixed to 1
  int<lower=0, upper=1> fix_probability_distortion; // 1 = alpha=1, beta=0
}

//*  ---------------------------------------------------
//               PARAMETERS
//*  ---------------------------------------------------

parameters {
  vector[7] mu_pr;
  vector<lower=0, upper=10>[7] sigma_pr;
  matrix[N, 7] param_raw;
}

//*  ---------------------------------------------------
//               MODEL
//*  ---------------------------------------------------

model {
  mu_pr               ~ std_normal();
  sigma_pr            ~ normal(0, 0.5);
  to_vector(param_raw) ~ std_normal();

  // Pin unused parameters to avoid wasting sampler effort and
  // inflating effective parameter count (p_loo) for restricted models.
  if (fix_sequential_weights == 1) {
    mu_pr[1:5]                  ~ normal(0, 0.01);
    sigma_pr[1:5]               ~ normal(0, 0.01);
    to_vector(param_raw[, 1:5]) ~ normal(0, 0.01);
  }
  if (fix_probability_distortion == 1) {
    mu_pr[6:7]                  ~ normal(0, 0.01);
    sigma_pr[6:7]               ~ normal(0, 0.01);
    to_vector(param_raw[, 6:7]) ~ normal(0, 0.01);
  }

  array[N] int indices;
  for (n in 1:N) indices[n] = n;

  target += reduce_sum(partial_sum, indices, grainsize,
                       mu_pr, sigma_pr,
                       Tsubj, sample, color, proba, choice,
                       param_raw,
                       fix_sequential_weights,
                       fix_probability_distortion);
}

//*  ---------------------------------------------------
//               GENERATED QUANTITIES
//*  ---------------------------------------------------

generated quantities {

  matrix[N, 7] params;

  array[N, T_max] int  y_pred        = rep_array(-1,   N, T_max);
  array[N, T_max] real pred_proba    = rep_array(-1.0, N, T_max);
  array[N, T_max] real diff_evidence = rep_array(-1.0, N, T_max);

  vector[sum(Tsubj)] log_lik;

  vector[7] g    = transform_group_means(mu_pr);
  real mu_w1     = g[1];
  real mu_w2     = g[2];
  real mu_w3     = g[3];
  real mu_w4     = g[4];
  real mu_w5     = g[5];
  real mu_alpha  = g[6];
  real mu_beta   = g[7];

  int k = 0;
  for (n in 1:N) {

    params[n] = (transform_params(mu_pr, sigma_pr, param_raw[n]))';

    for (t in 1:Tsubj[n]) {
      k += 1;

      int S  = sample[n, t];
      int ch = choice[n, t];

      if (S < 1 || ch < 1 || ch > 2) {
        log_lik[k] = 0;
        continue;
      }

      array[I_max] int  color_trial;
      array[I_max] real proba_trial;
      for (ii in 1:I_max) {
        color_trial[ii] = color[n, t, ii];
        proba_trial[ii] = proba[n, t, ii];
      }

      vector[2] ev = rep_vector(0.0, 2);
      for (s in 1:S) {
        real p = proba_trial[s];
        int  c = color_trial[s];
        if (c < 1 || c > 2) continue;
        if (p <= 0 || p >= 1) continue;
        real m  = map_linear_logodds(p, params[n, 6], params[n, 7],
                                     fix_probability_distortion);
        real ws = get_weight(s, params[n, 1:5]', fix_sequential_weights);
        ev[c] += ws * m;
      }

      real ev_diff = ev[1] - ev[2];
      real p_blue  = inv_logit(ev_diff);

      y_pred[n, t]        = bernoulli_rng(p_blue) ? 1 : 2;
      pred_proba[n, t]    = p_blue;
      diff_evidence[n, t] = ev_diff;

      log_lik[k] = bernoulli_logit_lpmf(ch == 1 | ev_diff);
    }
  }
}
