#ifndef _opac_util_
#define _opac_util_

#include <fundamental_constants.H>


//const Real o_chi_eV[4][5] = {{13.54,-1.0,-1.0,-1.0,-1.0},{24.48,54.17,-1.0,-1.0,-1.0},{11.26,24.38,47.89,64.49,-1.0},{13.62,35.12,54.93,77.41,113.9}};//ionization potentials

const Real o_chi_eV[2][2] = {{13.54,-1.0},{24.48,54.17}};
const Real o_gfac[2][2] = {{1.0,0.0},{4.0,1.0}};
//const Real o_gfac[4][5] = {{1.0,0.0,0.0,0.0,0.0},{4.0,1.0,0.0,0.0,0.0},{}}; //degeneracy factors
const int o_Z[2] = {1,2};
//const int o_Z[2] = {1,2,6,8}; //atomic number

//constexpr Real eps_root = 1.e-3;//relative convergence threshold

struct o_brent_args
{
  Real f[3]; //H and He only so no more than 2 anyways
  int k; //0=Hydrogen, 1=Helium
  Real n; //Number density of species k
  Real T; // temperature
  Real* ion_frac;
};

//boltzmann function f(T) for jth ion of species k
inline
Real o_fkj(int k, int j, Real T)
{
  Real kT = C::k_B*T;
  Real lf = 2.0_rt*3.14159265359_rt*C::m_e*kT;
  lf /= (C::hplanck*C::hplanck);
  Real lam_fac = pow(lf,1.5_rt);
  Real exp_fac = exp(-o_chi_eV[k-1][j]*C::ev2erg/kT);
  return o_gfac[k-1][j]*lam_fac*exp_fac;
}

// partition ion frac term
inline
Real o_b_ki(int k, int i, Real o_Zbar, Real n, Real* f)
{
  Real prod_sum = 1.0;

  for(int j = 1; j <= i; ++j) prod_sum *= f[j-1];

  return prod_sum /= pow(o_Zbar*n,i);
}


// transcendental equation for o_Zbar for each species
// // see Morozova 2015 (SNEC Code)
inline
Real o_theta(int k, Real n, Real* f, Real o_Zbar, Real& ion_0)
{
  ion_0 = 0.0;
  Real U = 1.0;


  for(int i = 1; i<=o_Z[k-1]; ++i)
  {
    Real beta = o_b_ki(k,i,o_Zbar,n,f);

    U += beta;
    ion_0 += i*beta;
  }

  ion_0 = 1.0_rt/ion_0;
  ion_0 *= o_Zbar;


  return 1.0_rt - ion_0*U;

}

// Brent-Dekker root finding algorithm for average charge o_Zbar=z_s
inline
int o_find_root(Real z_a, Real z_b, o_brent_args* args, Real& z_s)
{

  // fill boltzmann distribution terms
  for(int j = 0; j < args->k; ++j) args->f[j] = o_fkj(args->k,j,args->T);

  Real this_z_a = z_a;
  Real this_z_b = z_b;


  int mflag = 1;

  //
  Real this_ion_frac = args->ion_frac[0];

//  std::cout << "k = " << args->k << ", n = " << args->n << ", f0 = " << args->f[0] << ", ifrac0 = " << this_ion_frac <<"\n";

  Real f_a = o_theta(args->k,args->n,args->f,this_z_a,this_ion_frac);
  Real f_b =  o_theta(args->k,args->n,args->f,this_z_b,this_ion_frac);

 // std::cout << "f_a = " << f_a << "\n";
  //std::cout << "f_b = " << f_b << "\n";

  if (f_a*f_b >= 0 ) return -1;

  if(std::abs(f_a)<std::abs(f_b)){
    this_z_a = z_b;
    this_z_b = z_a;
	Real f_tmp = f_a;
	f_a = f_b;
	f_b = f_tmp;
  }

  Real z_c = this_z_a;
  Real f_c = f_a;
  Real z_d;
  Real f_s;

  Real delta = 1.e-6;//not sure how delta is different from epsilon..

  int MAXITER = 1000;
  int n_iter = 0;

  // main convergence loop until relative o_Zbar error is below threshold
  while((std::abs(this_z_b-this_z_a)/this_z_b > delta) && (n_iter < MAXITER))
  {

    if((f_a != f_c) && (f_b != f_c)){ // inverse quadratic interpolation

      z_s = this_z_a*f_b*f_c/(f_a-f_b)/(f_a-f_c);
      z_s += this_z_b*f_a*f_c/(f_b-f_a)/(f_b-f_c);
      z_s += z_c*f_a*f_b/(f_c-f_a)/(f_c-f_b);
    }
    else z_s = this_z_b-f_b*(this_z_b-this_z_a)/(f_b-f_a); // secant

    if (((z_s < (3.*this_z_a+this_z_b)/4.) || z_s > this_z_b)
       || (mflag > 0  && std::abs(z_s-this_z_b)>=0.5_rt*std::abs(this_z_b-z_c))
       || (mflag <= 0  && std::abs(z_s-this_z_b) >= std::abs(z_c-z_d)/2.0)
       || (mflag > 0 && std::abs(this_z_b-z_c)<delta)
       || (mflag <= 0 && std::abs(z_c-z_d)<delta)){
         z_s = 0.5_rt*(this_z_a+this_z_b);
       }
       else mflag = -1;

       f_s = o_theta(args->k,args->n,args->f,z_s,this_ion_frac);

       z_d = z_c;
       z_c = this_z_b;

       if(f_a*f_s < 0) this_z_b = z_s;
       else this_z_a = z_s;
       if(std::abs(f_a)<std::abs(f_b)){
         Real tmp = this_z_a;
         this_z_a = this_z_b;
         this_z_b = tmp;
       }

       n_iter++;

  }

  if(n_iter >= MAXITER) return 0;

  args->ion_frac[0] = this_ion_frac;

  args->ion_frac[1] = args->ion_frac[0]/(z_s*args->n)*args->f[0];

  if(args->k>1) args->ion_frac[2] = args->ion_frac[1]/(z_s*args->n)*args->f[1];

  

  return 1;

}

#endif
