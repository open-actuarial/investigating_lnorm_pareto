functions {
    vector lnormpareto_rng(int N, real alpha, real sigma, real theta) {
        real lnorm_val;
        real pareto_val;

        real mu;
        real mixt_prop;
        real tmp;

        real rng_val;

        vector[N] sample_val;
        vector[3] param_vals;

        param_vals = lnormpareto_vals(alpha, sigma, theta);

        mu        = param_vals[1];
        mixt_prop = param_vals[2];

        lnorm_val = theta + 1;

        for(i in 1:N) {
            rng_val = uniform_rng(0, 1);

            while(lnorm_val > theta) lnorm_val  = lognormal_rng(mu, sigma);

            pareto_val = pareto_type_2_rng(0.0, theta, alpha);

            sample_val[i] = (rng_val < mixt_prop) ? lnorm_val : pareto_val;
        }

        return(sample_val);
    }

    vector lnormpareto_vals(real alpha, real sigma, real theta) {
        vector[3] distrib_vals;

        real tmp;
        real mixt_prop;
        real mu = log(theta) - alpha * sigma * sigma;

        tmp = ( sqrt(2 * pi())
                * alpha
                * sigma
                * Phi(alpha * sigma)
                * exp(0.5 * pow(alpha * sigma, 2)) );

        mixt_prop = tmp / (tmp + 1);

        distrib_vals[1] = mu;
        distrib_vals[2] = mixt_prop;
        distrib_vals[3] = tmp;

        return(distrib_vals);
    }
}



data {
    int<lower=0> N;
    real x[N];
}

parameters {
    real<lower=0> alpha;
    real<lower=0> sigma;
    real<lower=0> theta;
}

transformed parameters {
    real mu = log(theta) - alpha * sigma * sigma;
}

model {
    real lpa[N];
    real tmp;
    real log_tmp;
    real neg_log1p_temp;
    real log_Phi_alpha_times_sigma;

    theta ~ gamma(.001, .001);
    alpha ~ gamma(.001, .001);
    sigma ~ gamma(.001, .001);


    tmp = ( sqrt(2 * pi())
            * alpha
            * sigma
            * Phi(alpha * sigma)
            * exp(0.5 * pow(alpha * sigma, 2)) );

    log_tmp = log(tmp);
    neg_log1p_temp = - log1p(tmp);
    log_Phi_alpha_times_sigma = log(Phi(alpha * sigma));

    for (i in 1:N) {
        if (x[i] < theta)
            lpa[i] = ( log_tmp
                       + neg_log1p_temp
                       - log_Phi_alpha_times_sigma
                       + lognormal_lpdf(x[i] | mu, sigma) );
        else
            lpa[i] = ( neg_log1p_temp
                       + pareto_type_2_lpdf(x[i] | 0.0, theta, alpha) );
    }

    target += sum(lpa);
}

