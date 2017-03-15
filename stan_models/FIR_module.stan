//Casey 2012 MBB+powerlaw model
functions {
  real z_to_DL(real z) {
    // Luminosity distance from redshift, using polynomial fit
    real D_l;
    #D_l = -84.13 + 5208.03*z+ 1898.61*pow(z,2)-247.82*pow(z,3);
    if (z < 1.0) {
      D_l = -0.211 + 4431.95*z + 3329.47*pow(z,2)-1047.52*pow(z,3);
      return D_l*3.086E22;
    }
    else
      D_l=-1724.8+ 7887.81*z +534.3*pow(z,2)-26.7*pow(z,3);
    return D_l*3.086E22;
  }

  real trapz(vector y, vector x, int N) {
    vector[N-1] integral;;
    for (i in 1:N-1) {
      integral[i]=(x[i+1]-x[i])*(y[i]+y[i+1])/2.0;
    }
    return sum(integral);
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
    B = pow(lam_c*1E-6,alpha);
    N_pl = A1/(A2*B);
    S_pl=N_pl*pow(lam*1E-6,alpha)*exp(-1.0*pow(lam/lam_c,2.0));
    return (S_bb+S_pl); //erg/s/Hz
  
  }

  real LIR(real Nbb, real alpha, real beta, real T, real c, real h, real k_B) {
      real a;
      real b;
      vector[497] L_nu;
      vector[496] integral;
      for (i in 1:497) {
	a=1000.0-(i-1)*2.0;
	L_nu[i]=Casey(Nbb,alpha,beta,T,a,c,h,k_B);
      }
      for (i in 1:496) {
	a=1000.0-(i-1)*2.0;
	b=1000.0-i*2.0;
	//integral[i]=1E6*c*(1.0/b - 1.0/a)*(Casey(Nbb,alpha,beta,T,b,c,h,k_B)+Casey(Nbb,alpha,beta,T,a,c,h,k_B))/2.0;
	integral[i]=1E6*c*(1.0/b - 1.0/a)*(L_nu[i+1]+L_nu[i])/2.0;

      }
      return sum(integral);
	
    }
  
  real filter(vector f_nu, vector filt_nu, vector filt_trans, real nu_0, int npoints) {
    // function to convolve and integrate BB with filter transmission
    real num;
    real denom;
    num=trapz(f_nu .*filt_trans,filt_nu, npoints);
    denom=trapz(nu_0*filt_trans ./ filt_nu, filt_nu,npoints);
    return num/denom;
    
  }
  real SPIRE_obs(real Nbb, real alpha, real Beta, real T, vector filt_wave_250, vector filt_nu_250, vector filt_trans_250, int Ntrans_250, real z,real c, real h, real k_B, real band){
    vector[Ntrans_250] S_nu_250;
    real S_250;
    for (i in 1:Ntrans_250){
    S_nu_250[i]=(1.0+z)*Casey(Nbb,alpha,Beta,T,filt_wave_250[i]/(1.0+z),c,h,k_B)/(4.0*pi() * square(100.0*z_to_DL(z))); //1E23convert to Jy
    }
    S_250=filter(S_nu_250,filt_nu_250, filt_trans_250, c/(band*1.0E-6), Ntrans_250);
    return S_250;
  }
 
  real SPIRE_obs_no_filt(real Nbb, real alpha, real Beta, real T, real z,real c, real h, real k_B, real band){
    real S_250;

    S_250=(1.0+z)*Casey(Nbb,alpha,Beta,T,band/(1.0+z),c,h,k_B)/(4.0*pi() * square(100.0*z_to_DL(z))); //1E23convert to Jy

    return S_250;
}}