#include "jla.h"

#include <class.h>
#include <iomanip>


/*
 * Fills CLASS-Code background structure with the best fit JLA LCDM
 * parameters.
 */
void default_pba(
                 struct background *pba
                 ){
  double sigma_B; /**< Stefan-Boltzmann constant in W/m^2/K^4 = Kg/K^4/s^3 */

  sigma_B = 2. * pow(_PI_,5) * pow(_k_B_,4) / 15. / pow(_h_P_,3) / pow(_c_,2);

  /** - background structure */
      
  pba->h = 0.70;
  pba->H0 = pba->h * 1.e5 / _c_;
  pba->T_cmb = 2.7255;
  pba->Omega0_g = (4.*sigma_B/_c_*pow(pba->T_cmb,4.)) / (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
  pba->Omega0_ur = 3.046*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
  pba->Omega0_b = 0.02253/0.704/0.704;
  pba->Omega0_cdm = 0.295 - pba->Omega0_b;
  pba->N_ncdm = 0;
  pba->Omega0_ncdm_tot = 0.;
  pba->ksi_ncdm_default = 0.;
  pba->ksi_ncdm = NULL;
  pba->T_ncdm_default = pow(4.0/11.0,1.0/3.0);
  pba->T_ncdm = NULL;
  pba->deg_ncdm_default = 1.;
  pba->deg_ncdm = NULL;
  pba->ncdm_psd_parameters = NULL;
  pba->ncdm_psd_files = NULL;

  pba->Omega0_k = 0.;
  pba->K = 0.;
  pba->sgnK = 0;
  pba->Omega0_lambda = 1.-pba->Omega0_k-pba->Omega0_g-pba->Omega0_ur-pba->Omega0_b-pba->Omega0_cdm-pba->Omega0_ncdm_tot;
  pba->Omega0_fld = 0.;     
  pba->a_today = 1.;       
  pba->w0_fld=-1.;
  pba->wa_fld=0.;
  pba->cs2_fld=1.;
}

using namespace std;

int main()
{
  printf("\n*\n* Testing the complete version of the JLALikelihood\n*\n\n");
  // Loading the JLA likelihood
  int verbosity = 3;
  JLALikelihood likelihood(verbosity);
  likelihood.read("data/jla.dataset");

  // Setting up CLASS code
  struct precision pr;
  struct background ba;
  input_default_precision(&pr);
  default_pba(&ba);
  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  // Getting distance modulii from CLASS
  double * z = likelihood.getZ();
  double mu[likelihood.size()];
  for (int i=0; i < likelihood.size(); i++)
    {
      double tau;
      double vec[ba.bg_size];
      int last_index;
      background_tau_of_z(&ba, z[i], &tau);
      background_at_tau(&ba, tau, ba.long_info,
                        ba.inter_normal, &last_index,
                        vec);
      mu[i] = 5*log10(vec[ba.index_bg_lum_distance]) + 25;
    }

  // JLA likelihood evaluation
  double nuisance_parameters[] = {0.141, 3.101, -19.05, -0.070};
  likelihood.computeLikelihood(mu, nuisance_parameters);
  
  // Cleaning
  free(z);

  printf("\n*\n* Testing the fast version of the JLALikelihood\n*\n\n");
  
  // Loading the simplified JLA likelihood

  SimplifiedJLALikelihood fast_likelihood(verbosity);
  fast_likelihood.read("data/jla_simple.dataset");

  // Getting distance modulii from CLASS
  double * z_bins = fast_likelihood.getZ();
  double mu_binned[fast_likelihood.size()];
  for (int i=0; i < fast_likelihood.size(); i++)
    {
      double tau;
      double vec[ba.bg_size];
      int last_index;
      background_tau_of_z(&ba, z_bins[i], &tau);
      background_at_tau(&ba, tau, ba.long_info,
                        ba.inter_normal, &last_index,
                        vec);
      mu_binned[i] = 5*log10(vec[ba.index_bg_lum_distance]) + 25;
    }

  double nuisance_parameters_small[] = {-.00219};
  fast_likelihood.computeLikelihood(mu_binned, nuisance_parameters_small);

  // Cleanning
  free(z_bins);
  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n", ba.error_message);
    return _FAILURE_;
    }

  return _SUCCESS_;
}
