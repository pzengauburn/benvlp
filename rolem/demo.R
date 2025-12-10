######################################################################
# RoLEM for data analysis 
######################################################################

rm(list = ls()); 
library(benvlp);                  # tested in version 0.1-0

######################################################################
## parameters in simulation 
######################################################################

case = "cs";                      # correlation structure for generating data 
n = 100;                          # number of subjects, 50, 100, or 200 
J = 5;                            # number of time points, 5 or 10 

r = 5;                            # number of responses 
p = 6;                            # number of covariates 
u = 2;                            # envelope dimension 
rho = 0.5;                        # temporal correlation 

cor_type = case;                  # correlation structure for fitting data 

size = 1e+3;                      # increase size to 1e+5 
burn_in = 1e+3;                   # increase burn_in to 1e+4
thinning = 20; 
info = 1000;                      # set info = -1 to suppress details 

######################################################################
# start simulation  
######################################################################

mydata = synthetic_data(n_subject = n, J = J, r = r, p = p, u = u, 
                    rho = rho, cor_type = case, distr = "t");
 
fit = rolem(mydata$y, mydata$x, mydata$subject, u,  
            Umat = mydata$Umat, cor_type = cor_type, normal = FALSE,  
            size = size, burn_in = burn_in, thinning = thinning, info = info); 

print(fit$Amat_accept)            # acceptance rate for A or equivalently P 
print(fit$rho_accept)             # acceptance rate for rho 
print(fit$nu_accept)              # acceptance rate for nu 
print(fit$posterior_mean)         # posterior means 

# par(mfrow = c(1, 2)); 
plot(fit$samples$Amat[1, 1, ], type = "l"); 
acf(fit$samples$Amat[1, 1, ]); 

ans = summary(fit); 

print(ans$BIC);                   # BIC 
print(ans$WAIC);                  # WAIC 
print(ans$beta$mean);             # posterior mean 
print(ans$beta$HPD_low);          # 95% HPD intervals, lower limit 
print(ans$beta$HPD_upp);          # 95% HPD intervals, upper limit  

######################################################################
# THE END 
######################################################################

