
data {
  int Ntrans_250; // number of transmission points in filter
  vector[Ntrans_250] filt_trans_250;
  vector[Ntrans_250] filt_wave_250;
  int Ntrans_350; // number of transmission points in filter
  vector[Ntrans_350] filt_trans_350;
  vector[Ntrans_350] filt_wave_350;
  int Ntrans_500; // number of transmission points in filter
  vector[Ntrans_500] filt_trans_500;
  vector[Ntrans_500] filt_wave_500;
  // real LUV;
  // real Beta_UV;
  // real sigma;
}

transformed data {
  real k_B;
  real c;
  real h;
  real BCratio;
  vector[Ntrans_250] filt_nu_250;
  vector[Ntrans_350] filt_nu_350;
  vector[Ntrans_500] filt_nu_500;


  
  k_B=1.3806488E-23;
  c=3.0E8;
  h=6.62606957E-34;


  for (i in 1:Ntrans_250) {
    filt_nu_250[i]=c/(filt_wave_250[i]*1E-6);
  }
  for (i in 1:Ntrans_350) {
    filt_nu_350[i]=c/(filt_wave_350[i]*1E-6);
  }
  for (i in 1:Ntrans_500) {
    filt_nu_500[i]=c/(filt_wave_500[i]*1E-6);
  }
}


parameters {
  real<lower=0.0,upper=10.0> Nbb;
  real<lower=20,upper=90> T;
  real<lower=1.6,upper=1.7> Beta;
  real<lower=1.5,upper=2.5> alpha;
  real<lower=2.9,upper=3.1> z;
  real<lower=-10,upper=-3> C_0;
  real<lower=0.0> C_1;
  real<lower=9.7,upper=11.6> LFUV;
  real<lower=-3,upper=1> Beta_UV;
}
transformed parameters {
  real Ldust;
  real S_250;
  real S_350;
  real S_500;
  real IRX;
  #real Beta_UV;
  Ldust = log10(1.0E-23*LIR(Nbb,alpha,Beta,T,c,h,k_B)/3.826E33);//convert to solar luminosities
  S_250 = SPIRE_obs(Nbb,alpha,Beta,T,filt_wave_250,filt_nu_250, filt_trans_250,Ntrans_250,z,c,h,k_B,250);
  S_350 = SPIRE_obs(Nbb,alpha,Beta,T,filt_wave_350,filt_nu_350, filt_trans_350,Ntrans_350,z,c,h,k_B,350);
  S_500 = SPIRE_obs(Nbb,alpha,Beta,T,filt_wave_500,filt_nu_500, filt_trans_500,Ntrans_500,z,c,h,k_B,500);
  //Beta_UV = (2.5*log10(0.595*pow(10.0,Ldust-LFUV)+1)-C_0)/C_1;
  IRX=Ldust-LFUV;
}
model {
  C_0~normal(-3.0,0.1);
  C_1~normal(2,0.2);
  Beta_UV~normal(C_0+C_1*IRX,0.1);

  }

