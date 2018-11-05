library(tidyverse)
library(rstan)
library(tidybayes)
library(bayesplot)


Sys.setenv(USE_CXX14 = 1)

source('fire.data.R')


options(width = 80L
        ,warn  = 1
        ,mc.cores = parallel::detectCores()
)

rstan_options(auto_write = TRUE)


fire_data_lst <- list(x = x
                     ,N = length(x)
                     )


### Run the BUGS stan model from the BUGS Book Vol 3
fire_stanmodel <- stan_model('fire.stan', verbose = TRUE)

fire_stanfit <- sampling(object = fire_stanmodel
                        ,data   = fire_data_lst
                        ,seed   = 42
                        ,show_messages = TRUE
)



### Fit our personal lnorm pareto Stan model
# lnorm_pareto_stanmodel <- stan_model('lnorm_pareto.stan')
#
# lnorm_pareto_stanfit <- sampling(object = lnorm_pareto_stanmodel
#                                 ,data   = fire_data_lst
#                                 ,seed   = 42
#                                 ,show_messages = TRUE
#                                  )

