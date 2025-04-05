#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 15:51:33 2023

@author: Irene

v4) change the standard deviation on i0 prior using input value, cyto_std
v4) testing a model with exponential prior on intensity of the spot (gaussian_plane_exp)
"""

# Models

# ---------------------------------------------------------------------------#

# (2023-02-07)standard deviation for i0 is fixed to 20au, which fits
# both the raw and corrected images (instead of cyto_avg/10).

poisson = """
data {
    int<lower=0> n;
    vector[n] x;
    vector[n] y;
    int intensity[n,n];
    real x_red;
    real y_red;
    real cyto_avg;
    real cyto_std;
}
transformed data {}
parameters {
    real x0;
    real y0;
    real<lower=0> I;
    real<lower=0> sigma;
    
    real i0;
    
    real<lower=0> lambda;
    real<lower=0> tau;
}
transformed parameters {}

model {
    lambda ~ cauchy(0, 1);
    tau ~ cauchy(0, 50);
    
    sigma ~ inv_gamma(80, 16);
    I ~ normal(0, lambda * tau);
    x0 ~ normal(x_red, 0.2);
    y0 ~ normal(y_red, 0.2);
    
    i0 ~ normal(cyto_avg, 20);
    
    for (i in 1:n) {
        for (j in 1:n){
            intensity[i,j] ~ poisson(i0 + I/(2*pi()*pow(sigma,2)) * exp(- pow(x[i]-x0, 2) / (2*pow(sigma, 2)) - pow(y[j]-y0, 2)/(2*pow(sigma, 2))));
        }
    }
}

generated quantities {}
"""

# ---------------------------------------------------------------------------#

# (2023-02-07)standard deviation for i0 is fixed to 20au, which fits
# both the raw and corrected images (instead of cyto_avg/10).

poisson_plane = """
data {
    int<lower=0> n;
    vector[n] x;
    vector[n] y;
    int intensity[n,n];
    real x_red;
    real y_red;
    real cyto_avg;
    real cyto_std;
}
transformed data {}
parameters {
    real x0;
    real y0;
    real<lower=0> I;
    real<lower=0> sigma;
    
    real a;
    real b; 
    real i0;
    
    real<lower=0> lambda;
    real<lower=0> tau;
}
transformed parameters {}

model {
    lambda ~ cauchy(0, 1);
    tau ~ cauchy(0, 50);
    
    sigma ~ inv_gamma(80, 16);
    I ~ normal(0, lambda * tau);
    x0 ~ normal(x_red, 0.2);
    y0 ~ normal(y_red, 0.2);
    
    a ~ normal(0,5);
    b ~ normal(0,5);
    i0 ~ normal(cyto_avg, 20);
    
    for (i in 1:n) {
        for (j in 1:n){
            intensity[i,j] ~ poisson(a * x[i] + b * y[j] + i0 + I/(2*pi()*pow(sigma,2)) * exp(- pow(x[i]-x0, 2) / (2*pow(sigma, 2)) - pow(y[j]-y0, 2)/(2*pow(sigma, 2))));
        }
    }
}

generated quantities {}
"""

# ---------------------------------------------------------------------------#

# (2023-02-07)standard deviation for i0 is fixed to 20au, which fits
# both the raw and corrected images (instead of cyto_avg/10).

gaussian = """
data {
    int<lower=0> n;
    vector[n] x;
    vector[n] y;
    int intensity[n,n];
    real x_red;
    real y_red;
    real cyto_avg;
    real cyto_std;
}
transformed data {}
parameters {
    real x0;
    real y0;
    real<lower=0> I;
    real<lower=0> sigma;
    
    real i0;
    
    real<lower=0> lambda;
    real<lower=0> tau;
}
transformed parameters {}

model {
    lambda ~ cauchy(0, 1);
    tau ~ cauchy(0, 50);
    
    sigma ~ inv_gamma(80, 16);
    I ~ normal(0, lambda * tau);
    x0 ~ normal(x_red, 0.2);
    y0 ~ normal(y_red, 0.2);
    
    i0 ~ normal(cyto_avg, 20);
    
    for (i in 1:n) {
        for (j in 1:n){
            intensity[i,j] ~ normal(i0 + I/(2*pi()*pow(sigma,2)) * exp(- pow(x[i]-x0, 2) / (2*pow(sigma, 2)) - pow(y[j]-y0, 2)/(2*pow(sigma, 2))), cyto_std);
        }
    }
}

generated quantities {}
"""


# ---------------------------------------------------------------------------#

# Gaussian noise model
# Background fitted as a 2D plane (3 parameters)
# Standard deviation of the gaussian is fixed (input parameter)

# (2023-02-07)standard deviation for i0 is fixed to 20au, which fits
# both the raw and corrected images (instead of cyto_avg/10).

# (2023-05-16) standard deviation for i0 is an input parameter


gaussian_plane = """
data {
    int<lower=0> n;
    vector[n] x;
    vector[n] y;
    int intensity[n,n];
    real x_red;
    real y_red;
    real cyto_avg;
    real cyto_std;
}
transformed data {}
parameters {
    real x0;
    real y0;
    real<lower=0> I;
    real<lower=0> sigma;
    
    real a;
    real b; 
    real<lower=0> i0;
    
    real<lower=0> lambda;
    real<lower=0> tau;
}
transformed parameters {}

model {
    lambda ~ cauchy(0, 1);
    tau ~ cauchy(0, 50);
    
    sigma ~ inv_gamma(80, 16);
    I ~ normal(0, lambda * tau);
    x0 ~ normal(x_red, 0.2);
    y0 ~ normal(y_red, 0.2);
    
    a ~ normal(0, 5);
    b ~ normal(0, 5);
    i0 ~ normal(cyto_avg, cyto_std);
    
    for (i in 1:n) {
        for (j in 1:n){
            intensity[i,j] ~ normal(a * x[i] + b * y[j] + i0 + I/(2*pi()*pow(sigma,2)) * exp( - pow(x[i]-x0, 2)/(2*pow(sigma, 2)) - pow(y[j]-y0, 2)/(2*pow(sigma, 2)) ), cyto_std);
        }
    }
}

generated quantities {}
"""

# ---------------------------------------------------------------------------#

# similar to gaussian_plane except that cyto_std represents only the 
# standard deviation on the i0 prior, while the standard deviation
# of the gaussian noise model is calculated independently and it is fixed.

gaussian_plane_v1 = """
data {
    int<lower=0> n;
    vector[n] x;
    vector[n] y;
    int intensity[n,n];
    real x_red;
    real y_red;
    real cyto_avg;
    real cyto_std;
}
transformed data {}
parameters {
    real x0;
    real y0;
    real<lower=0> I;
    real<lower=0> sigma;
    
    real a;
    real b; 
    real<lower=0> i0;
    
    real<lower=0> lambda;
    real<lower=0> tau;
}
transformed parameters {
    real normal_std = 55.0;
}
model {
    lambda ~ cauchy(0, 1);
    tau ~ cauchy(0, 50);
    
    sigma ~ inv_gamma(80, 16);
    I ~ normal(0, lambda * tau);
    x0 ~ normal(x_red, 0.2);
    y0 ~ normal(y_red, 0.2);
    
    a ~ normal(0, 5);
    b ~ normal(0, 5);
    i0 ~ normal(cyto_avg, cyto_std);
    
    for (i in 1:n) {
        for (j in 1:n){
            intensity[i,j] ~ normal(a * x[i] + b * y[j] + i0 + I/(2*pi()*pow(sigma,2)) * exp( - pow(x[i]-x0, 2)/(2*pow(sigma, 2)) - pow(y[j]-y0, 2)/(2*pow(sigma, 2)) ), normal_std);
        }
    }
}
generated quantities {}
"""


# ---------------------------------------------------------------------------#

# Gaussian noise model
# Background fitted as a 2D plane (3 parameters)
# Standard deviation of the gaussian is inferred.

# (2023-02-07)standard deviation for i0 is fixed to 20au, which fits
# both the raw and corrected images (instead of cyto_avg/10).


gaussian_std_plane = """
data {
    int<lower=0> n;
    vector[n] x;
    vector[n] y;
    int intensity[n,n];
    real x_red;
    real y_red;
    real cyto_avg;
    real cyto_std;
}
transformed data {}
parameters {
    real x0;
    real y0;
    real<lower=0> I;
    real<lower=0> sigma;
    
    real a;
    real b; 
    real<lower=0> i0;
    
    real std;
    
    real<lower=0> lambda;
    real<lower=0> tau;
}
transformed parameters {}

model {
    lambda ~ cauchy(0, 1);
    tau ~ cauchy(0, 50);
    
    sigma ~ inv_gamma(80, 16);
    I ~ normal(0, lambda * tau);
    x0 ~ normal(x_red, 0.2);
    y0 ~ normal(y_red, 0.2);
    
    a ~ normal(0,5);
    b ~ normal(0,5);
    i0 ~ normal(cyto_avg, 20);
    
    std ~ normal(cyto_std, 10);
    
    for (i in 1:n) {
        for (j in 1:n){
            intensity[i,j] ~ normal(a * x[i] + b * y[j] + i0 + I/(2*pi()*pow(sigma,2)) * exp(- pow(x[i]-x0, 2) / (2*pow(sigma, 2)) - pow(y[j]-y0, 2)/(2*pow(sigma, 2))), std);
        }
    }
}

generated quantities {}
"""


# ---------------------------------------------------------------------------#


# (2023-02-07)standard deviation for i0 is fixed to 20au, which fits
# both the raw and corrected images (instead of cyto_avg/10).

gaussian_prop_plane = """
data {
    int<lower=0> n;
    vector[n] x;
    vector[n] y;
    int intensity[n,n];
    real x_red;
    real y_red;
    real cyto_avg;
    real cyto_std;
}
transformed data {}
parameters {
    real x0;
    real y0;
    real<lower=0> I;
    real<lower=0> sigma;
    
    real a;
    real b; 
    real<lower=0> i0;
    
    real<lower=0> var_;
    
    real<lower=0> lambda;
    real<lower=0> tau;
}
transformed parameters {}

model {
    lambda ~ cauchy(0, 1);
    tau ~ cauchy(0, 50);
    
    sigma ~ inv_gamma(80, 16);
    I ~ normal(0, lambda * tau);
    x0 ~ normal(x_red, 0.2);
    y0 ~ normal(y_red, 0.2);
    
    a ~ normal(0,5);
    b ~ normal(0,5);
    i0 ~ normal(cyto_avg, 20);
    
    var_ ~ normal(8*i0-350, 250);
    
    for (i in 1:n) {
        for (j in 1:n){
            intensity[i,j] ~ normal(a * x[i] + b * y[j] + i0 + I/(2*pi()*pow(sigma,2)) * exp(- pow(x[i]-x0, 2) / (2*pow(sigma, 2)) - pow(y[j]-y0, 2)/(2*pow(sigma, 2))), sqrt(var_));
        }
    }
}

generated quantities {}
"""

# ---------------------------------------------------------------------------#
#v4

#real<lower=0> I0;

gaussian_plane_exp = """
data {
    int<lower=0> n;
    vector[n] x;
    vector[n] y;
    int intensity[n,n];
    real x_red;
    real y_red;
    real cyto_avg;
    real cyto_std;
    real<lower=0> I0;
}
transformed data {}
parameters {
    real x0;
    real y0;
    real<lower=0> I;
    real<lower=0> sigma;
    
    real a;
    real b; 
    real<lower=0> i0;
    
    real<lower=0> lambda;
    real<lower=0> tau;
}
transformed parameters {}

model { 
    sigma ~ inv_gamma(80, 16);
    I ~ exponential(1/I0);
    x0 ~ normal(x_red, 0.2);
    y0 ~ normal(y_red, 0.2);
    
    a ~ normal(0, 5);
    b ~ normal(0, 5);
    i0 ~ normal(cyto_avg, cyto_std);
    
    for (i in 1:n) {
        for (j in 1:n){
            intensity[i,j] ~ normal(a * x[i] + b * y[j] + i0 + I/(2*pi()*pow(sigma,2)) * exp( - pow(x[i]-x0, 2)/(2*pow(sigma, 2)) - pow(y[j]-y0, 2)/(2*pow(sigma, 2)) ), cyto_std);
        }
    }
}

generated quantities {}
"""


# -------------- Log-normal distribution for pixel intensity ------------------ #

# test a log-normal model for pixel intensity
# the sigma parameter of the log-normal distribution should be computed independently
# from the images

lognorm = """
data {
    int<lower=0> n;
    vector[n] x;
    vector[n] y;
    int intensity[n,n];
    real x_red;
    real y_red;
    real cyto_avg;
    real cyto_std;
}
transformed data {}
parameters {
    real x0;
    real y0;
    real<lower=0> I;
    real<lower=0> sigma;
    
    real a;
    real b; 
    real<lower=0> i0;
    
    real<lower=0> lambda;
    real<lower=0> tau;
}
transformed parameters {
    real sigma_ln = 0.5;
}
model {
    lambda ~ cauchy(0, 1);
    tau ~ cauchy(0, 50);
    
    sigma ~ inv_gamma(80, 16);
    I ~ normal(0, lambda * tau);
    x0 ~ normal(x_red, 0.2);
    y0 ~ normal(y_red, 0.2);
    
    a ~ normal(0, 5);
    b ~ normal(0, 5);
    i0 ~ normal(cyto_avg, cyto_std);
    
    for (i in 1:n) {
        for (j in 1:n){
            intensity[i,j] ~ lognormal(log( a * x[i] + b * y[j] + i0 + I/(2*pi()*pow(sigma,2)) * exp( - pow(x[i]-x0, 2)/(2*pow(sigma, 2)) - pow(y[j]-y0, 2)/(2*pow(sigma, 2)) ) ), sigma_ln);
        }
    }
}
generated quantities {}
"""




# --------------- ONLY FOR CALIBRATION EXPERIMENTS WITH MULTIPLE LASER POWERS USED IN THE SAME ACQUISITION ----------------#



# Trying to adapt the gaussian_plane model, to fit data taken at different laser 
# powers (like calibration), where I cannot set one single value for the background intensity
# or standard deviation, and I need to use very broad priors.
# we prefere to put the standard deviation of the gaussian as a parameter, because it changes a lot
# depending on the laser power


gaussian_plane_cal =  """
data {
    int<lower=0> n;
    vector[n] x;
    vector[n] y;
    int intensity[n,n];
    real x_red;
    real y_red;
    real cyto_avg;
    real cyto_std;
    real std;
}
transformed data {
    // trying to fix the standard deviation
    //real std = 200;
}
parameters {
    real x0;
    real y0;
    real<lower=0> I;
    real<lower=0> sigma;
    
    real a;
    real b; 
    real<lower=0> i0;
    
    //real<lower=0> std;
    
    real<lower=0> lambda;
    real<lower=0> tau;
}
transformed parameters {}

model {
    lambda ~ cauchy(0, 1);
    tau ~ cauchy(0, 50);
    
    sigma ~ inv_gamma(80, 16);
    I ~ normal(0, lambda * tau);
    x0 ~ normal(x_red, 0.2);
    y0 ~ normal(y_red, 0.2);
    
    a ~ normal(0,5);
    b ~ normal(0,5);
    i0 ~ normal(0, 1000); // hopefully uninformative prior
    
    //std ~ normal(0, 100); // also uninformative
    
    for (i in 1:n) {
        for (j in 1:n){
            intensity[i,j] ~ normal(a * x[i] + b * y[j] + i0 + I/(2*pi()*pow(sigma,2)) * exp(- pow(x[i]-x0, 2) / (2*pow(sigma, 2)) - pow(y[j]-y0, 2)/(2*pow(sigma, 2))), std);
        }
    }
}

generated quantities {}
"""

