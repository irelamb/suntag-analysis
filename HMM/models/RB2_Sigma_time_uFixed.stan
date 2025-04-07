// originally model3 in HMM_models.py
// RB2 = Ribosome-Based 2 Conditions (CTRL, Perturbed)




// Fix b0

// Fix u (mature protein intensity)

data{
    int<lower = 0> S; // SunTag length (599aa)
    int<lower = 0> L; // total reporter length
    int ncond; // number of different conditions
    int ntraces; // total number of traces
    int ndata; // total number of data points (traces x length) for all conditions
    
    array[ntraces] int<lower = 1, upper = ndata> starts; // trace starts
    array[ncond] int<lower = 1, upper = ndata> cond_starts; // conditions starts
    
    vector[ndata] I; // measured intensities (traces)
    vector[ndata] t; // elapsed time after treatment (not imaging)
    
    vector[ncond] f_inactive; // fraction of inactive (non-translated) traces (estimate) for both conditions
    
    int N_max; // maximum number of ribosomes
    
    real<lower=0> time_step; // imaging timestep (in sec)
    int<lower=1> l; // particle size in aa
    
    int<lower=0, upper=N_max> max_step_size; // maximum step size (N->N-max_step_size)
    real<lower=0> b0; // intensity measurement offset
    real<lower=0> sigma0;

    real<lower=0> u; // intensity of one mature protein from calibration
    }                                           
transformed data{
    
    // Define ends arrays
    array[ntraces] int<lower = 1, upper = ndata> ends;
    
    for (p in 1 : ntraces){
        ends[p] = p < ntraces ? starts[p + 1] - 1 : ndata;
        }
    
    // Define array with condition per each data point
    array[ndata] int<lower = 1, upper = ncond> cc;
    
    for (c in 1 : ncond){
        int cond_end = c < ncond ? cond_starts[c + 1] - 1 : ndata;
        cc[cond_starts[c] : cond_end] = rep_array(c, cond_end - cond_starts[c] + 1);
        print(cond_starts[c]);
        print(cond_end);
        }
    }
parameters{
    
    // observation parameters:
    //real<lower=0> u; // mature protein intensity (common parameter)
    
    // kinetics parameters:
    vector<lower=5>[ncond] T_i; // initiation waiting times
    vector<lower=0>[ncond] k_e; // elongation rate 
    
    // probability of inactive trace
    vector<lower=0, upper=1>[ncond] p_off;
    }
transformed parameters{

    vector<lower=0>[ncond] k_i = 1 / T_i;

    // transition probabilities
    
    array[ncond] simplex[max_step_size + 1] theta;
    
    for (c in 1 : ncond){
        theta[c][1] = exp(- k_i[c] * time_step);
        
        for (d in 1 : max_step_size - 1){
            theta[c][d + 1] = k_i[c] * time_step * theta[c][d] / d;
            }
        
        real bias = k_e[c] * time_step / L; // correction bias due to run-off
        theta[c][1] = theta[c][1] - bias;
        theta[c][2] = theta[c][2] + bias;
        
        theta[c][max_step_size + 1] = 1 - sum(theta[c][1 : max_step_size]);
        }
    
    //
    
    array[ncond] vector[max_step_size + 1] log_theta;
    for (c in 1 : ncond){
        log_theta[c] = log(theta[c]); // cache log
        }

    // correction factor (continuous approximation)
    
    vector<lower=0>[ndata] f;
    {
        real f0 = 1 - S * 1.0 / (2 * L);// value at t = 0
        // if there is a delay this is not completely exact
        
        for (n in 1 : ndata){
            f[n] = 1 - (1 - f0) * exp(-1.4 * pow((k_e[cc[n]] * t[n] + 206) / S, 3)); 
        }
    }

    // continuous approximation for sigma

    vector<lower=0>[ndata] sigma;
    {   
        real a = 5.4; // average of optimized parameters over different ribosome densities
        real b = 76;

        real var0 = exp(pow(sigma0, 2)) - 1; // variance given sigma0, normalized by mean intensity
        
        real N_avg;
        vector[ncond] s0;

        for (c in 1 : ncond){
            N_avg = L / (T_i[c] * k_e[c]); // average number of ribosomes, given the parameters
            s0[c] = sqrt(log(1 + var0 + 1 / N_avg * S / 3.0 * (4 * L - 3 * S) / pow(2 * L - S, 2))); // value at t=0
        }

        for (n in 1 : ndata){
            sigma[n] = sigma0 + (s0[cc[n]] - sigma0) * exp(-a * pow((k_e[cc[n]] * t[n] + b) / S, 2)); 
        }
    }
    print("sigma[1 : 10] ", sigma[1 : 10]);
}
model{
    
    // priors
    
    //u ~ normal(20, 10);
    
    // priors on kinetic params
    T_i ~ normal(100, 30);
    k_e ~ normal(3, 2);
    
    //
    for (c in 1 : ncond){
        p_off[c] ~ normal(f_inactive[c], 0.1);
        }

    
    vector[max_step_size + 1] acc;
    array[ndata, N_max + 1] real Gamma;
    
    // probability of initial state

    for (s in starts){
        
        int c = cc[s];
        
        real lambda = k_i[c] * (L / k_e[c] - t[s]);
        
        Gamma[s, 1] = log(p_off[c]) + lognormal_lpdf(I[s] | log(b0), sigma0);
        
        for (N in 1 : N_max){
            Gamma[s, N + 1] = log1m(p_off[c]) + poisson_lpmf(N | lambda) - log1m_exp(-lambda)
                              + lognormal_lpdf(I[s] | log(b0 + u * N * f[s]), sigma[s]);
            }
        }
    
    for (p in 1 : ntraces){  
        for (k in (starts[p] + 1) : ends[p]){
            
            int c = cc[k];
            
            ///////// N = 0 /////////
            
            acc = rep_vector(lognormal_lpdf(I[k] | log(b0), sigma0), max_step_size + 1);
            
            real q = 0;
            
            for (d in 0 : max_step_size - 1){
                acc[d + 1] += Gamma[k - 1, d + 1] + log1m(q);   
                q += theta[c][d + 1];
                }
            // d = max_step_size
            acc[max_step_size + 1] += Gamma[k - 1, max_step_size + 1] + log_theta[c][max_step_size + 1];
            
            Gamma[k, 1] = log_sum_exp(acc);
            
            //////// N > 0 /////////
            
            for (N in 1 : N_max){
                
                acc = rep_vector(negative_infinity(), max_step_size + 1);
                
                for (d in 0 : max_step_size){
                        
                    if (N+d < N_max + 1){
                        
                        acc[d + 1] = Gamma[k - 1, N + d + 1] + log_theta[c][d + 1]
                                     + lognormal_lpdf(I[k] | log(b0 + u * N * f[k]), sigma[k]);                 
                        }
                    }
                
                Gamma[k, N + 1] = log_sum_exp(acc);
                }
            }
        target += log_sum_exp(Gamma[ends[p]]);
        }
    }
generated quantities{
    
    // Viterbi algorithm
    
    array[ndata] int<lower = 0, upper = N_max> N_star;
    real logp_N_star;
    
    {
        array[ndata, N_max + 1] int back_ptr;
        array[ndata, N_max + 1] real best_logp;
        real best_total_logp;
        
        
        for (s in starts){
                
            int c = cc[s];
                
            real lambda = k_i[c] * (L / k_e[c] - t[s]);
            
            best_logp[s, 1] = log(p_off[c]) + lognormal_lpdf(I[s] | log(b0), sigma0);
                
            for (N in 1 : N_max){
                best_logp[s, N + 1] = log1m(p_off[c]) + poisson_lpmf(N | lambda) - log1m_exp(-lambda)
                                      + lognormal_lpdf(I[s] | log(b0 + u * N * f[s]), sigma[s]);
                }
            }
        
        for (p in 1 : ntraces){
            for (k in starts[p] + 1 : ends[p]){
                
                int c = cc[k];
                
                ///////// N = 0 /////////
                
                best_logp[k, 1] = negative_infinity();
                real logp = negative_infinity();

                real q = 0;
                
                for (d in 0 : max_step_size - 1){
                    logp = best_logp[k-1, d+1] + log1m(q) + lognormal_lpdf(I[k] | log(b0), sigma0); 
                    q += theta[c][d + 1];
                    
                    if (logp > best_logp[k, 1]) {
                        back_ptr[k, 1] = d;
                        best_logp[k, 1] = logp;
                        }
                    }
                
                // d = max_step_size
                logp = best_logp[k-1, max_step_size+1] + log_theta[c][max_step_size + 1] 
                       + lognormal_lpdf(I[k] | log(b0), sigma0);
                if (logp > best_logp[k, 1]) {
                    back_ptr[k, 1] = max_step_size;
                    best_logp[k, 1] = logp;
                    }
                
                //////// N > 0 /////////
                
                
                for (N in 1 : N_max){
                    
                    best_logp[k, N+1] = negative_infinity();
                    logp = negative_infinity();
                    
                    for (d in 0 : max_step_size){
                            
                        if (N+d < N_max + 1){
                            logp = best_logp[k-1, N+d+1] + log_theta[c][d + 1]
                                   + lognormal_lpdf(I[k] | log(b0 + u * N * f[k]), sigma[k]);                 
                            
                            if (logp > best_logp[k, N+1]) {
                                back_ptr[k, N+1] = N+d;
                                best_logp[k, N+1] = logp;
                                }
                            }
                        }
                    }
                }
            }
        
        for (end in ends){
            logp_N_star = max(best_logp[end]);
            
            for (N in 0 : N_max) {
                if (best_logp[end, N+1] == logp_N_star) {
                    N_star[end] = N;
                    }
                }
            }
        
        for (p in 1 : ntraces){
            for (k in 1 : (ends[p] - starts[p])){
                print("N_star index: ", ends[p] - k + 1);
                print("N_star: ", N_star[ends[p] - k + 1]);
                N_star[ends[p] - k] = back_ptr[ends[p] - k + 1, N_star[ends[p] - k + 1] + 1];
                }
            }
        }
    }
