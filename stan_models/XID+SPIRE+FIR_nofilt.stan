//Full Bayesian inference fit  XID
data {
  int<lower=0> nsrc;//number of sources
  //----PSW----
  int<lower=0> npix_psw;//number of pixels
  int<lower=0> nnz_psw; //number of non neg entries in A
  vector[npix_psw] db_psw;//flattened map
  vector[npix_psw] sigma_psw;//flattened uncertianty map (assuming no covariance between pixels)
  real bkg_prior_psw;//prior estimate of background
  real bkg_prior_sig_psw;//sigma of prior estimate of background
  vector[nnz_psw] Val_psw;//non neg values in image matrix
  int Row_psw[nnz_psw];//Rows of non neg valies in image matrix
  int Col_psw[nnz_psw];//Cols of non neg values in image matrix
  //----PMW----
  int<lower=0> npix_pmw;//number of pixels
  int<lower=0> nnz_pmw; //number of non neg entries in A
  vector[npix_pmw] db_pmw;//flattened map
  vector[npix_pmw] sigma_pmw;//flattened uncertianty map (assuming no covariance between pixels)
  real bkg_prior_pmw;//prior estimate of background
  real bkg_prior_sig_pmw;//sigma of prior estimate of background
  vector[nnz_pmw] Val_pmw;//non neg values in image matrix
  int Row_pmw[nnz_pmw];//Rows of non neg valies in image matrix
  int Col_pmw[nnz_pmw];//Cols of non neg values in image matrix
  //----PLW----
  int<lower=0> npix_plw;//number of pixels
  int<lower=0> nnz_plw; //number of non neg entries in A
  vector[npix_plw] db_plw;//flattened map
  vector[npix_plw] sigma_plw;//flattened uncertianty map (assuming no covariance between pixels)
  real bkg_prior_plw;//prior estimate of background
  real bkg_prior_sig_plw;//sigma of prior estimate of background
  vector[nnz_plw] Val_plw;//non neg values in image matrix
  int Row_plw[nnz_plw];//Rows of non neg valies in image matrix
  int Col_plw[nnz_plw];//Cols of non neg values in image matrix


}

transformed data {
  real k_B;
  real c;
  real h;
  k_B=1.3806488E-23;
  c=3.0E8;
  h=6.62606957E-34;

}


parameters {
  real bkg_psw;//background
  real bkg_pmw;//background
  real bkg_plw;//background
  real<lower=0.0> sigma_conf_psw;
  real<lower=0.0> sigma_conf_pmw;
  real<lower=0.0> sigma_conf_plw;

  real<lower=0.001,upper=5> z[nsrc];

  real<lower=8, upper=14> Nbb[nsrc];
  real<lower=-1, upper=5> T[nsrc];
  real<lower=1.0,upper=2.5> Beta[nsrc];
  real<lower=1.0,upper=2.6> alpha[nsrc];
  //real<lower=0.001,upper=5> z[nsrc];
}

transformed parameters {
  vector[nsrc] f_vec_psw;//vector of source fluxes
  vector[nsrc] f_vec_pmw;//vector of source fluxes
  vector[nsrc] f_vec_plw;//vector of source fluxes
  vector[nsrc] L_dust;

  for (i in 1:nsrc){


    L_dust[i]=1.0E-23*LIR(1.0,alpha[i],Beta[i],exp(T[i]),c,h,k_B)/3.826E33;
    f_vec_psw[i]=pow(10.0,Nbb[i])*SPIRE_obs_no_filt(1.0,alpha[i],Beta[i],exp(T[i]),z[i],c,h,k_B,250.0)*1E3/L_dust[i];
    f_vec_pmw[i]=pow(10.0,Nbb[i])*SPIRE_obs_no_filt(1.0,alpha[i],Beta[i],exp(T[i]),z[i],c,h,k_B,350.0)*1E3/L_dust[i];
    f_vec_plw[i]=pow(10.0,Nbb[i])*SPIRE_obs_no_filt(1.0,alpha[i],Beta[i],exp(T[i]),z[i],c,h,k_B,500.0)*1E3/L_dust[i];
}
print("Nbb",Nbb)
}



model {
  vector[npix_psw] db_hat_psw;//model of map
  vector[npix_pmw] db_hat_pmw;//model of map
  vector[npix_plw] db_hat_plw;//model of map


  vector[npix_psw] sigma_tot_psw;
  vector[npix_pmw] sigma_tot_pmw;
  vector[npix_plw] sigma_tot_plw;


 //Prior on background 
  bkg_psw ~normal(bkg_prior_psw,bkg_prior_sig_psw);
  bkg_pmw ~normal(bkg_prior_pmw,bkg_prior_sig_pmw);
  bkg_plw ~normal(bkg_prior_plw,bkg_prior_sig_plw); 

 //Prior on conf
  sigma_conf_psw ~normal(0,5);
  sigma_conf_pmw ~normal(0,5);
  sigma_conf_plw ~normal(0,5);

   
  // Create model maps (i.e. db_hat = A*f) using sparse multiplication
  for (k in 1:npix_psw) {
    db_hat_psw[k] <- bkg_psw;
    sigma_tot_psw[k]<-sqrt(square(sigma_psw[k])+square(sigma_conf_psw));
  }
  for (k in 1:nnz_psw) {
    db_hat_psw[Row_psw[k]+1] <- db_hat_psw[Row_psw[k]+1] + Val_psw[k]*f_vec_psw[Col_psw[k]+1];
      }

  for (k in 1:npix_pmw) {
    db_hat_pmw[k] <-  bkg_pmw;
    sigma_tot_pmw[k]<-sqrt(square(sigma_pmw[k])+square(sigma_conf_pmw));
  }
  for (k in 1:nnz_pmw) {
    db_hat_pmw[Row_pmw[k]+1] <- db_hat_pmw[Row_pmw[k]+1] + Val_pmw[k]*f_vec_pmw[Col_pmw[k]+1];
      }

  for (k in 1:npix_plw) {
    db_hat_plw[k] <- bkg_plw;
    sigma_tot_plw[k]<-sqrt(square(sigma_plw[k])+square(sigma_conf_plw));
  }
  for (k in 1:nnz_plw) {
    db_hat_plw[Row_plw[k]+1] <- db_hat_plw[Row_plw[k]+1] + Val_plw[k]*f_vec_plw[Col_plw[k]+1];
      }

  // likelihood of observed map|model map
  db_psw ~ normal(db_hat_psw,sigma_tot_psw);
  db_pmw ~ normal(db_hat_pmw,sigma_tot_pmw);
  db_plw ~ normal(db_hat_plw,sigma_tot_plw);

    }
