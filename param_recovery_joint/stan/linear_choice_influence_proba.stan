//*  ---------------------------------------------------
//               FUNCTIONS
//*  ---------------------------------------------------

functions {

  // ===== Transform parameters =====
  // Returns a vector of 10 transformed subject-level parameters:
  //   1-5: w1..w5 in [0, 2]
  //   6:   alpha in [0, 10]
  //   7:   beta (unbounded)
  //   8:   a_infl (unbounded)
  //   9:   b_infl (unbounded)
  //  10:   sigma_infl > 0
  vector transform_params(vector mu_pr, vector sigma_pr, row_vector param_raw_n) {
    vector[10] p;
    for (j in 1:5)
      p[j] = 2 * Phi_approx(mu_pr[j] + sigma_pr[j] * param_raw_n[j]);
    p[6]  = 10 * Phi_approx(mu_pr[6] + sigma_pr[6] * param_raw_n[6]);
    p[7]  = mu_pr[7] + sigma_pr[7] * param_raw_n[7];
    p[8]  = mu_pr[8] + sigma_pr[8] * param_raw_n[8];
    p[9]  = mu_pr[9] + sigma_pr[9] * param_raw_n[9];
    p[10] = exp(mu_pr[10] + sigma_pr[10] * param_raw_n[10]);
    return p;
  }

  vector transform_group_means(vector mu_pr) {
    vector[10] g;
    for (j in 1:5)
      g[j] = 2 * Phi_approx(mu_pr[j]);
    g[6]  = 10 * Phi_approx(mu_pr[6]);
    g[7]  = mu_pr[7];
    g[8]  = mu_pr[8];
    g[9]  = mu_pr[9];
    g[10] = exp(mu_pr[10]);
    return g;
  }

  // ===== Simple mapping: m(p) = alpha * logit(p) + beta =====
  real map_linear_logodds(real p, real alpha, real beta) {
    return alpha * logit(p) + beta;
  }

  // ===== Sequential weight (w6 = 1 fixed) =====
  real get_weight(int s, vector w5) {
    if (s <= 5) return w5[s];
    return 1.0;
  }

  // ===== Trial helper =====
  // Returns: [evidence_blue, evidence_red, sub_evidence]
  vector compute_evidence_and_sub(int sample_size,
                                  array[] int color_data,
                                  array[] real proba_data,
                                  array[] int influence_sample_data,
                                  vector w5,
                                  real alpha, real beta) {
    vector[3] out = rep_vector(0.0, 3);

    for (s in 1:sample_size) {
      real p = proba_data[s];
      int c  = color_data[s];

      if (c < 1 || c > 2) continue;
      if (p <= 0 || p >= 1) continue;

      real m  = map_linear_logodds(p, alpha, beta);
      real ws = get_weight(s, w5);

      out[c] += ws * m;

      if (influence_sample_data[s] == 1) {
        real d_color = (c == 1) ? 1.0 : -1.0;
        out[3] += ws * m * d_color;
      }
    }
    return out;
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
                   array[,,] int influence_sample,
                   array[,] real influence,
                   matrix param_raw) {
    real lp = 0;

    for (i in 1:size(slice_indices)) {
      int n = slice_indices[i];

      vector[10] params = transform_params(mu_pr, sigma_pr, param_raw[n]);

      // locals for speed
      vector[5] w5 = params[1:5];
      real alpha      = params[6];
      real beta       = params[7];
      real a_infl     = params[8];
      real b_infl     = params[9];
      real sigma_infl = params[10];

      for (t in 1:Tsubj[n]) {
        int S = sample[n, t];
        if (S < 1) continue;

        int ch = choice[n, t];
        if (ch < 1 || ch > 2) continue;

        vector[2] evidence = rep_vector(0.0, 2);
        real sub_evidence = 0.0;

        for (s in 1:S) {
          real p = proba[n, t, s];
          int c  = color[n, t, s];

          if (c < 1 || c > 2) continue;
          if (p <= 0 || p >= 1) continue;

          real m  = map_linear_logodds(p, alpha, beta);
          real ws = get_weight(s, w5);

          evidence[c] += ws * m;

          if (influence_sample[n, t, s] == 1) {
            real d_color = (c == 1) ? 1.0 : -1.0;
            sub_evidence += ws * m * d_color;
          }
        }

        // --- Choice likelihood ---
        real ev_diff = evidence[1] - evidence[2];
        lp += bernoulli_logit_lpmf(ch == 1 | ev_diff);

        // --- Influence report likelihood ---
        real total_ev_choice = (ch == 1) ? ev_diff : -ev_diff;
        real sub_ev_choice   = (ch == 1) ? sub_evidence : -sub_evidence;
        real non_sub_ev_choice = total_ev_choice - sub_ev_choice;

        real p_with    = inv_logit(total_ev_choice);
        real p_without = inv_logit(non_sub_ev_choice);
        real delta_p   = p_with - p_without;

        real infl_obs = influence[n, t];
        if (infl_obs > -90) {
          real mu_infl = a_infl * delta_p + b_infl;
          lp += normal_lpdf(infl_obs | mu_infl, sigma_infl);
        }
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

  // 1=Blue, 2=Red, -1 padding/missing
  array[N, T_max] int<lower=-1, upper=2> choice;

  // Influence data
  array[N, T_max, I_max] int<lower=-1, upper=1> influence_sample;
  array[N, T_max] real influence;  // on [-0.5, 0.5], use -99 for missing

  int<lower=5> grainsize;
}

//*  ---------------------------------------------------
//               PARAMETERS
//*  ---------------------------------------------------

parameters {
  vector[10] mu_pr;
  vector<lower=0, upper=10>[10] sigma_pr;
  matrix[N, 10] param_raw;
}

//*  ---------------------------------------------------
//               MODEL
//*  ---------------------------------------------------

model {
  mu_pr ~ std_normal();
  sigma_pr ~ normal(0, 0.5);
  to_vector(param_raw) ~ std_normal();

  array[N] int indices;
  for (n in 1:N) indices[n] = n;

  target += reduce_sum(partial_sum, indices, grainsize,
                       mu_pr, sigma_pr,
                       Tsubj, sample, color, proba, choice,
                       influence_sample, influence,
                       param_raw);
}

//*  ---------------------------------------------------
//               GENERATED QUANTITIES
//*  ---------------------------------------------------

generated quantities {
  matrix[N, 10] params;
  array[N, T_max] int y_pred = rep_array(-1, N, T_max);
  array[N, T_max] real pred_proba = rep_array(-1.0, N, T_max);
  array[N, T_max] real diff_evidence = rep_array(-1.0, N, T_max);
  array[N, T_max] real delta_p_out = rep_array(-99.0, N, T_max);
  array[N, T_max] real sub_ev_out = rep_array(-99.0, N, T_max);
  array[N, T_max] real pred_influence = rep_array(-99.0, N, T_max);

  vector[sum(Tsubj)] log_lik;
  vector[sum(Tsubj)] log_lik_choice;
  vector[sum(Tsubj)] log_lik_influence;

  // --- group-level (interpretable) means ---
  vector[10] g = transform_group_means(mu_pr);
  real mu_w1          = g[1];
  real mu_w2          = g[2];
  real mu_w3          = g[3];
  real mu_w4          = g[4];
  real mu_w5          = g[5];
  real mu_alpha       = g[6];
  real mu_beta        = g[7];
  real mu_a_infl      = g[8];
  real mu_b_infl      = g[9];
  real mu_sigma_infl  = g[10];

  int k = 0;
  for (n in 1:N) {
    // subject-level transformed params
    params[n] = (transform_params(mu_pr, sigma_pr, param_raw[n]))';

    for (t in 1:Tsubj[n]) {
      k += 1;

      int S  = sample[n, t];
      int ch = choice[n, t];

      if (S < 1 || ch < 1 || ch > 2) {
        y_pred[n, t] = -1;
        pred_proba[n, t] = -1.0;
        diff_evidence[n, t] = -1.0;
        log_lik[k] = 0;
        log_lik_choice[k] = 0;
        log_lik_influence[k] = 0;
        continue;
      }

      // pack trial arrays for helper
      array[I_max] int color_trial;
      array[I_max] real proba_trial;
      array[I_max] int infl_sample_trial;
      for (ii in 1:I_max) {
        color_trial[ii] = color[n, t, ii];
        proba_trial[ii] = proba[n, t, ii];
        infl_sample_trial[ii] = influence_sample[n, t, ii];
      }

      vector[3] ev = compute_evidence_and_sub(
        S,
        color_trial, proba_trial,
        infl_sample_trial,
        params[n, 1:5]',
        params[n, 6], params[n, 7]
      );

      // --- Choice ---
      real ev_diff = ev[1] - ev[2];
      real p_blue = inv_logit(ev_diff);

      y_pred[n, t] = bernoulli_rng(p_blue) ? 1 : 2;
      pred_proba[n, t] = p_blue;
      diff_evidence[n, t] = ev_diff;

      real ll_choice = bernoulli_logit_lpmf(ch == 1 | ev_diff);

      // --- Influence report ---
      real total_ev_choice = (ch == 1) ? ev_diff : -ev_diff;
      real sub_ev_choice   = (ch == 1) ? ev[3] : -ev[3];
      real non_sub_ev_choice = total_ev_choice - sub_ev_choice;

      real p_with    = inv_logit(total_ev_choice);
      real p_without = inv_logit(non_sub_ev_choice);
      real delta_p   = p_with - p_without;

      sub_ev_out[n, t] = sub_ev_choice;
      delta_p_out[n, t] = delta_p;

      real ll_infl = 0;
      real infl_obs = influence[n, t];

      if (infl_obs > -90) {
        real mu_infl = params[n, 8] * delta_p + params[n, 9];
        ll_infl = normal_lpdf(infl_obs | mu_infl, params[n, 10]);

        // Posterior predictive
        pred_influence[n, t] = normal_rng(mu_infl, params[n, 10]);
      }

      log_lik_choice[k] = ll_choice;
      log_lik_influence[k] = ll_infl;
      log_lik[k] = ll_choice + ll_infl;
    }
  }
}
