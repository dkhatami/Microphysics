#ifndef actual_opacity_H
#define actual_opacity_H

#include <AMReX.H>

#include <extern_parameters.H>
#include <prob_parameters.H>
#include <opac_util.H>
#include <fundamental_constants.H>

inline
Real inv_gau(Real x,Real l, Real m)
{
	return 0.4*pow(l,0.5)*pow(x,-1.5)*exp(-l*(x-m)*(x-m)/(2.*m*m*x));
}

inline
Real H_mix_abs(Real rho, Real T)
{
	Real rfac = rho/1.e-10;
	Real dT_a = -3.7-0.04*log10(rfac);
	Real a = 7.0*pow(rfac,0.8)*inv_gau(log10(T)+dT_a,3.2,0.33);
	Real dT_b = -4.2-0.04*log10(rfac);
	Real b = 0.2*pow(rfac,0.8)*inv_gau(log10(T)+dT_b,4.1,0.38);
	return max(1.e-5,max(a,b));
}

inline
Real He_mix_abs(Real rho, Real T)
{
	Real rfac = rho/1.e-10;
	Real dT_a = -3.9-0.04*log10(rfac);
	Real a = 0.16*pow(rfac,0.8)*inv_gau(log10(T)+dT_a,2.6,0.32);
	Real dT_b = -4.2-0.04*log10(rfac);
	Real b = 0.3*pow(rfac,0.8)*inv_gau(log10(T)+dT_b,4.1,0.38);

	return max(1.e-5,max(a,b));


}

inline
Real CO_mix_abs(Real rho, Real T)
{
	Real rfac = rho/1.e-10;
	Real dT = -3.3-0.04*log10(rfac);
	return max(1.e-5,0.5*pow(rfac,0.8)*inv_gau(log10(T)+dT,18.0,1.0));


}


// free-free opacit from Rybicki & Lightman Eq. (5.20)
// there should be a Z^2 and gaunt factor here
// and one of the n_e's should technically
// be replaced by n_ion
// but this should be fine since we just want some
// absorption in case compton scattering is too small
inline
Real ff_abs(Real rhoYe, Real T)
{
	Real n_e = rhoYe/C::m_p;
	return 1.7e-25_rt*pow(T,-3.5)*n_e*n_e;
}




AMREX_INLINE
void actual_opacity_init () {}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_opacity (Real& kp, Real& kr, Real rho, Real temp, Real rhoYe, Real nu,
                     bool get_Planck_mean, bool get_Rosseland_mean)
{
    const Real Ksc = 0.4e0; // Thomson scattering
    const Real fac = 1.e-3; // Planck mean is assumed to be fac * Ksc
    const Real Z_frac = 0.05;
	const Real K_lo = 1.e-4;
	const Real sigma_t = 6.65e-25;

	if(problem::do_grey_es){

		Real kappa_es = problem::kap_grey;

		if(get_Planck_mean) kp = kappa_es*rho*problem::eps;
		if(get_Rosseland_mean) kr = kappa_es*rho;
	}
	else if(temp < problem::T_cold){
		if(get_Planck_mean) kp = K_lo*rho;
		if(get_Rosseland_mean) kr = K_lo*rho;
	}
	else{
    Real kap_abs = 0.0_rt;


	if(problem::opac_use_real_ff && temp < problem::T_ff && rho > problem::rho_ff) kap_abs += ff_abs(rhoYe,temp)/rho;

//     Real t_fac = pow(temp,-3.5);
 //    kap_ff = 1.e23 *rhoYe*t_fac;
  //   kap_bf = 4.34e25*Z_frac*rho*t_fac;
  // }

    if (get_Planck_mean) {
        kp = rhoYe/rho*Ksc*problem::eps;
	//	kp = rho*Ksc*fac;

//        if(temp > 4.5e4) kp += (kap_ff+kap_bf)*rho;
//
		kp += kap_abs;

		kp *= rho;

		kp = max(kp,K_lo*rho);
    }

    if (get_Rosseland_mean) {
	    kr = rhoYe/rho*Ksc;
//		kr = rho*Ksc;

 //       if(temp > 4.5e4) kr += (kap_ff+kap_bf)*rho;
		kr += kap_abs;

		kr *= rho;

		kr = max(kr,K_lo*rho);

    }

	}

}

#endif
