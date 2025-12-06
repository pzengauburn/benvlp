######################################################################
# Example 1: compare RoLEM with RLMM and LEM  
######################################################################

rm(list = ls()); 
library(benvlp);                  # tested in version 0.1-0

######################################################################
## parameters in simulation 
######################################################################

case = "ar1";                     # correlation structure for generating data 
n = 100;                          # number of subjects, 50, 100, or 200 
J = 5;                            # number of time points, 5 or 10 
iter = 3;                         # number of iterations 

r = 20;                           # number of responses 
p = 30;                           # number of covariates 
u = 3;                            # envelope dimension 
rho = 0.5;                        # temporal correlation 

cor_type = case;                  # correlation structure for fitting data 

size = 1e+3;                      # increase size to 1e+5 
burn_in = 1e+3;                   # increase burn_in to 1e+4
thinning = 20; 
info = 1000;                      # set info = -1 to suppress details 

######################################################################
# hyperparameters  
######################################################################

hyper_pars_rolem = list(eta_mean = matrix(0, nrow = r, ncol = p),    
                        eta_cov_inv = diag(rep(1e-6, p), nrow = p, ncol = p),   
                        Omega_df = u + 1,     
                        Omega_scale = diag(rep(1e-6, u), nrow = u, ncol = u), 
                        Omega0_df = r - u + 1,   
                        Omega0_scale = diag(rep(1e-6, r-u), nrow = r-u, ncol = r-u),
                        P_Mmat = diag(rep(1e-6, r)),        
                        nu_shape = 1.4, 
                        nu_rate = 0.04);

hyper_pars_lem = hyper_pars_rolem;

hyper_pars_rlmm = list(xi = matrix(0, nrow = r, ncol = p),    
                       Hmat = diag(rep(1e-6, p), nrow = p, ncol = p),   
                       df = r + 1,     
                       Psi = diag(rep(1e-6, r), nrow = r, ncol = r), 
                       nu_shape = 1.4,    
                       nu_rate = 0.04);

######################################################################
# tuning parameters 
######################################################################

if((n == 200) && (J == 10))
{
    pro_pars_lem =   list(Pmat_sigma2 = 1e-7,   rho_delta = 0.025); 
    pro_pars_rlmm =  list(                      rho_delta = 0.025, nu_delta = 0.06); 
    pro_pars_rolem = list(Pmat_sigma2 = 7e-8, rho_delta = 0.025, nu_delta = 1.4); 
} else if((n == 200) && (J == 5))
{
    pro_pars_lem =   list(Pmat_sigma2 = 1.4e-7, rho_delta = 0.025); 
    pro_pars_rlmm =  list(                      rho_delta = 0.025, nu_delta = 0.06); 
    pro_pars_rolem = list(Pmat_sigma2 = 1.1e-7, rho_delta = 0.025, nu_delta = 1.5); 
} else if((n == 100) && (J == 10))
{
    pro_pars_lem =   list(Pmat_sigma2 = 1.2e-7, rho_delta = 0.025); 
    pro_pars_rlmm =  list(                      rho_delta = 0.025, nu_delta = 0.09); 
    pro_pars_rolem = list(Pmat_sigma2 = 1.0e-7, rho_delta = 0.025, nu_delta = 0.12); 
} else if((n == 100) && (J == 5))
{
    pro_pars_lem =   list(Pmat_sigma2 = 2.2e-7, rho_delta = 0.03); 
    pro_pars_rlmm =  list(                      rho_delta = 0.03, nu_delta = 0.08); 
    pro_pars_rolem = list(Pmat_sigma2 = 1.4e-7, rho_delta = 0.03, nu_delta = 0.12);  
} else if((n == 50) && (J == 10))
{
    pro_pars_lem =   list(Pmat_sigma2 = 3e-7, rho_delta = 0.03); 
    pro_pars_rlmm =  list(                    rho_delta = 0.03, nu_delta = 0.1); 
    pro_pars_rolem = list(Pmat_sigma2 = 2e-7, rho_delta = 0.03, nu_delta = 0.15); 
} else if((n == 50) && (J == 5))
{
    pro_pars_lem =   list(Pmat_sigma2 = 3.5e-7, rho_delta = 0.045); 
    pro_pars_rlmm =  list(                      rho_delta = 0.045, nu_delta = 0.08); 
    pro_pars_rolem = list(Pmat_sigma2 = 2.4e-7, rho_delta = 0.045, nu_delta = 0.15); 
} else stop("unsupported n or J.\n"); 

######################################################################
# start simulation  
######################################################################

performance = NULL; 
HPDs = NULL; 

time_start = date(); 
cat("start ==>", time_start, "\n");
for(i in 1:iter)
{
    cat(date(), ": i =", i, "\n");

    mydata = synthetic_data(n_subject = n, J = J, r = r, p = p, u = u, 
                    rho = rho, cor_type = case, distr = "t");

    ######################################################################
    # rolem 
    ######################################################################

    ans1 = rolem(mydata$y, mydata$x, mydata$subject, u,  
                     Umat = mydata$Umat, cor_type = cor_type, normal = FALSE,  
                     hyper_pars = hyper_pars_rolem, pro_pars = pro_pars_rolem, 
                     size = size, burn_in = burn_in, thinning = thinning, info = info); 

    perf = get_performance(ans1, mydata);
    perf$cor_type = cor_type; 
    perf$func = "rolem"; 
    perf$iter = i; 
    performance = rbind(performance, perf); 

    tmp = get_more_performance(ans1, mydata, 0.05); 
    tmp$cor_type = cor_type;
    tmp$func = "rolem";
    tmp$iter = i; 
    HPDs = rbind(HPDs, tmp); 

    ######################################################################
    # rlmm 
    ######################################################################

    ans2 = rlmm(mydata$y, mydata$x, mydata$subject, 
                     cor_type = cor_type,
                     normal = FALSE, 
                     hyper_pars = hyper_pars_rlmm, pro_pars = pro_pars_rlmm, 
                     size = size, burn_in = burn_in, thinning = thinning, info = info); 

    perf = get_performance(ans2, mydata);
    perf$cor_type = cor_type; 
    perf$func = "rlmm"; 
    perf$iter = i; 
    performance = rbind(performance, perf); 

    tmp = get_more_performance(ans2, mydata, 0.05); 
    tmp$cor_type = cor_type;
    tmp$func = "rlmm";
    tmp$iter = i; 
    HPDs = rbind(HPDs, tmp); 

    ######################################################################
    # lem 
    ######################################################################

    ans3 = lem(mydata$y, mydata$x, mydata$subject, u,  
                     Umat = mydata$Umat, cor_type = cor_type,
                     hyper_pars = hyper_pars_lem, pro_pars = pro_pars_lem, 
                     size = size, burn_in = burn_in, thinning = thinning, info = info); 

    perf = get_performance(ans3, mydata);
    perf$cor_type = cor_type; 
    perf$func = "lem"; 
    perf$iter = i; 
    performance = rbind(performance, perf); 

    tmp = get_more_performance(ans3, mydata, 0.05); 
    tmp$cor_type = cor_type;
    tmp$func = "lem";
    tmp$iter = i; 
    HPDs = rbind(HPDs, tmp); 

}
time_stop = date(); 
cat("stop ==>", time_stop, "\n");

######################################################################
# save results 
######################################################################

hyper_pars_rolem = ans1$hyper_pars;
pro_pars_rolem = ans1$pro_pars;
     
hyper_pars_rlmm = ans2$hyper_pars; 
pro_pars_rlmm = ans2$pro_pars;

hyper_pars_lem = ans3$hyper_pars; 
pro_pars_lem = ans3$pro_pars;

output_filename = paste0("ex-1-", case, "-", n, "-", J, ".Rdata"); 

save(time_start, time_stop, case, cor_type, performance, HPDs, 
     hyper_pars_rolem, pro_pars_rolem, 
     hyper_pars_rlmm, pro_pars_rlmm, 
     hyper_pars_lem, pro_pars_lem, 
     n, J, p, r, rho, u, iter, size, burn_in, thinning,  
     file = output_filename); 

q('no');

######################################################################
# THE END 
######################################################################
