//Casey 2012 MBB+powerlaw model
functions {
  real z_to_DL(real z) {
    // Luminosity distance from redshift, using polynomial fit
    real D_l;
    #D_l = -84.13 + 5208.03*z+ 1898.61*pow(z,2)-247.82*pow(z,3);
    D_l = -0.211 + 4431.95*z + 3329.47*pow(z,2)-1047.52*pow(z,3);

    return D_l*3.086E24;
    
  }

  real trapz(vector y, vector x, int N) {
    real integral;
    integral=0.0;
    for (i in 1:N-1) {
      integral =integral+ (x[i+1]-x[i])*(y[i]+y[i+1])/2.0;
    }
    return integral;
  }
  real Casey(real Nbb, real alpha, real beta, real T, real lam, real c, real h, real k_B) {
    real S_bb;
    real L_a_t;
    real lam_c;
    real A1;
    real A2;
    real B;
    real N_pl;
    real S_pl;
    
    S_bb=Nbb*pow(c/(lam*1.0E-6),3.0+beta)/(exp(h*c/(((1.0E-6)*lam)*k_B*T))-1.0);
    L_a_t=pow((pow(26.68+6.246*alpha,-2.0)+((1.905E-4)+(7.243E-5)*alpha)*T),-1.0);
    lam_c=3.0*L_a_t/4.0;
    A1 = pow((c/(lam_c*1.0E-6)),(beta+3.0));
    A2 = exp((h*c)/((1.0E-6)*lam_c*k_B*T))-1.0;
    B = pow(lam_c,alpha);
    N_pl = Nbb*A1/(A2*B);
    S_pl=N_pl*pow(lam,alpha)*exp(-1.0*pow(lam/lam_c,2.0));
    return (S_bb+S_pl); //erg/s/Hz
  }

  real LIR(real Nbb, real alpha, real beta, real T, real c, real h, real k_B) {
      real a;
      real b;
      real integral;
      integral=0.0;
      for (i in 1:496) {
	a=8+(i-1.0)*2.0;
	b=8+i*2.0;
	integral=integral+1E-6*(b-a)*(Casey(Nbb,alpha,beta,T,b,c,h,k_B)+Casey(Nbb,alpha,beta,T,a,c,h,k_B))/2.0;
      }
      return integral;
	
    }
  
  real filter(vector f_nu, vector filt_nu, vector filt_trans, real nu_0, int npoints) {
    // function to convolve and integrate BB with filter transmission
    real num;
    real denom;
    num=trapz(f_nu .*filt_trans,filt_nu, npoints);
    denom=trapz(nu_0*filt_trans ./ filt_nu, filt_nu,npoints);
    return num/denom;
    
  }


}

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
  vector[3] SPIRE_fluxes;
  vector[3] SPIRE_flux_uncert;
  real z;
  // real LUV;
  // real Beta_UV;
  // real sigma;
}

transformed data {
  real k_B;
  real c;
  real h;
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
  real<lower=0.0> Nbb;
  real<lower=20.0,upper=100.0> T;
  real<lower=1.0,upper=2.0> Beta;
  real<lower=0.5,upper=3.0> alpha;
  //real<lower=0.0> C_0;
  //real<lower=0.0> C_1;
  //real Ldust;
}

transformed parameters {
  real Ldust;
  vector[496] S_nu;

  Ldust = LIR(Nbb,alpha,Beta,T,c,h,k_B)/3.826E33;//convert to solar luminosities
  for ( i in 1:496) {
    S_nu[i] = Casey(Nbb,alpha,Beta,T,6.0+2.0*i,c,h,k_B);//*square(6+2*i*1.0E-6/(1.0+z))/(c*4.0*pi() * square(z_to_DL(z)))*1E23; //1E23convert to Jy
  }


  
}

model{
  vector[Ntrans_250] S_nu_250;
  vector[Ntrans_350] S_nu_350;
  vector[Ntrans_500] S_nu_500;
  vector[3] SPIRE_model_fluxes;
  for (i in 1:Ntrans_250){
    S_nu_250[i]=(1.0+z)*Casey(Nbb,alpha,Beta,T,filt_wave_250[i]/(1.0+z),c,h,k_B)*square(filt_wave_250[i]*1.0E-6/(1.0+z))/(c*4.0*pi() * square(z_to_DL(z)))*1E23; //1E23convert to Jy
  }

  for (i in 1:Ntrans_350){
    S_nu_350[i]=(1.0+z)*Casey(Nbb,alpha,Beta,T,filt_wave_350[i]/(1.0+z),c,h,k_B)*square(filt_wave_350[i]*1.0E-6/(1.0+z))/(c*4.0*pi() * square(z_to_DL(z)))*1E23; //1E23convert to Jy
  }

  for (i in 1:Ntrans_500){
    S_nu_500[i]=(1.0+z)*Casey(Nbb,alpha,Beta,T,filt_wave_500[i]/(1.0+z),c,h,k_B)*square(filt_wave_500[i]*1.0E-6/(1.0+z))/(c*4.0*pi() * square(z_to_DL(z)))*1E23; //1E23convert to Jy
  }

  SPIRE_model_fluxes[1] = filter(S_nu_250,filt_nu_250, filt_trans_250, c/(250.0*1E-6), Ntrans_250);
  SPIRE_model_fluxes[2] = filter(S_nu_350,filt_nu_350, filt_trans_350, c/(350.0*1E-6), Ntrans_350);
  SPIRE_model_fluxes[3] = filter(S_nu_500,filt_nu_500, filt_trans_500, c/(500.0*1E-6), Ntrans_500);
  SPIRE_fluxes ~ normal(SPIRE_model_fluxes,SPIRE_flux_uncert);

}

