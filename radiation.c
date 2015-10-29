

/* 

model-independent radiation-related utilities.

*/

#include "decs.h"

double Bnu_inv(double nu, double Thetae)
{

	double x;

	x = HPL * nu / (ME * CL * CL * Thetae);

	if (x < 1.e-3)		/* Taylor expand */
		return ((2. * HPL / (CL * CL)) /
			(x / 24. * (24. + x * (12. + x * (4. + x)))));
	else
		return ((2. * HPL / (CL * CL)) / (exp(x) - 1.));
}

double jnu_inv(double nu, double Thetae, double Ne, double B, double theta)
{
	double j;

	j = jnu_synch(nu, Ne, Thetae, B, theta);

	return (j / (nu * nu));
}

/* return Lorentz invariant scattering opacity */
double alpha_inv_scatt(double nu, double Thetae, double Ne)
{
	double kappa;

	kappa = kappa_es(nu, Thetae);

	return (nu * kappa * Ne * MP);
}

#include <gsl/gsl_sf.h>
double kappa_I_abs(double nu, double Thetae, double Ne, double B,
                   double theta
                  )
{
  double electron_charge = EE;
  double mass_electron   = ME;
  double speed_light     = CL;
  double obs_angle       = theta;
  double n_e             = Ne;
  double kappa_width     = Thetae;
  double kappa           = 3.5;


  double nu_c = (electron_charge * B)
    /(2. * M_PI * mass_electron * speed_light);
  double nu_w = pow(kappa_width*kappa, 2.)*nu_c*sin(obs_angle);
  double X_k = nu/nu_w;
  double prefactor = n_e*electron_charge/(B*sin(obs_angle));
  double a = kappa - 1./3.;
  double b = kappa + 1.;
  double c = kappa + 2./3.;
  double z = -kappa*kappa_width;
  double hyp2f1 = pow(1.-z, -a)*tgamma(c)*tgamma(b-a)/(tgamma(b)*tgamma(c-a))
    *gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z))+pow(1.-z, -b)
    *tgamma(c)*tgamma(a-b)/(tgamma(a)*tgamma(c-b))
    *gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z));
  double Nlow = pow(3., 1./6.)*(10./41.)*pow(2.*M_PI, 2.)
    /pow(kappa_width*kappa, 16./3.-kappa)*(kappa-2.)*(kappa-1.)
    *kappa/(3.*kappa-1.)*tgamma(5./3.)*hyp2f1;
  double Nhigh = 2.*pow(M_PI, 5./2.)/3.*(kappa-2.)*(kappa-1.)*kappa
    /pow(kappa_width*kappa, 5.)*(2*tgamma(2.+kappa/2.)
        /(2.+kappa)-1.)*(pow(3./kappa, 19./4.)+3./5.);
  double x = pow(-7./4. + 8.*kappa/5., -43./50.);
  double ans = prefactor*Nlow*pow(X_k, -5./3.)
    *pow(1.+pow(X_k, x*(3.*kappa-1.)/6.)*pow(Nlow/Nhigh, x), -1./x);
  return ans;
}

#define THERMAL (0)
#define KAPPA   (1)
#define ABSORPTIVITY KAPPA
/* return Lorentz invariant absorption opacity */
double alpha_inv_abs(double nu, double Thetae, double Ne, double B,
		     double theta)
{
	double j, bnu;
#if (ABSORPTIVITY == THERMAL)
	j = jnu_inv(nu, Thetae, Ne, B, theta);
	bnu = Bnu_inv(nu, Thetae);

	return (j / (bnu + 1.e-100));
#elif (ABSORPTIVITY == KAPPA)
  return (nu*kappa_I_abs(nu, Thetae, Ne, B, theta));
#endif
}


/* return electron scattering opacity, in cgs */
double kappa_es(double nu, double Thetae)
{
	double Eg;

	/* assume pure hydrogen gas to 
	   convert cross section to opacity */
	Eg = HPL * nu / (ME * CL * CL);
	return (total_compton_cross_lkup(Eg, Thetae) / MP);
}

/* get frequency in fluid frame, in Hz */
double get_fluid_nu(double X[4], double K[4], double Ucov[NDIM])
{
	double ener, nu;

	/* this is the energy in electron rest-mass units */
	ener = -(K[0] * Ucov[0] +
		 K[1] * Ucov[1] + K[2] * Ucov[2] + K[3] * Ucov[3]);

	nu = ener * ME * CL * CL / HPL;

	if (isnan(ener)) {
		fprintf(stderr, "isnan get_fluid_nu, K: %g %g %g %g\n",
			K[0], K[1], K[2], K[3]);
		fprintf(stderr, "isnan get_fluid_nu, X: %g %g %g %g\n",
			X[0], X[1], X[2], X[3]);
		fprintf(stderr, "isnan get_fluid_nu, U: %g %g %g %g\n",
			Ucov[0], Ucov[1], Ucov[2], Ucov[3]);
	}

	return nu;

}

/* return angle between magnetic field and wavevector */
double get_bk_angle(double X[NDIM], double K[NDIM], double Ucov[NDIM],
		    double Bcov[NDIM], double B)
{

	double k, mu;

	if (B == 0.)
		return (M_PI / 2.);

	k = fabs(K[0] * Ucov[0] + K[1] * Ucov[1] + K[2] * Ucov[2] +
		 K[3] * Ucov[3]);

	/* B is in cgs but Bcov is in code units */
	mu = (K[0] * Bcov[0] + K[1] * Bcov[1] + K[2] * Bcov[2] +
	      K[3] * Bcov[3]) / (k * B / B_unit);

	if (fabs(mu) > 1.)
		mu /= fabs(mu);

	return (acos(mu));
}
