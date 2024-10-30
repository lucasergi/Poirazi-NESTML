#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "CA1.h"
#include "Na_Kdr_spine.h"

CA1_params::CA1_params() :
  //V_peak(0.0),            // mV
  //V_reset(-60.0),         // mV
  //g_L(30.0),              // ns
  //g_L(0.0001),            //
  g_L(0.1),            //
  //C_m(281.0),             // pF
  C_m(1.0),                // pF
  //E_L(-70.6),             // mV
  E_L(-70.0),             // mV
  //Delta_T(2.0),           // mV
  //tau_w(144.0),           // ms
  //a(4.0),
  //b(80.5),                // pA
  //V_th(-50.4),            // mV
  I_e(50.0),               // pA
  tau_nmda_rise(2.0),     // ms    Rise time for NMDA channels
  tau_nmda_fast(10.0),    // ms    Fast decay time for NMDA channels
  tau_nmda_slow(45.0),    // ms    Slow decay time for NMDA channels
  alpha_nmda(0.527),      //       ratio of fast decay NMDA conducance over total NMDA conductance 
  refr_T(2.0),            // ms    Duration of refractory period
  //E_exc(0.0),             // mV    Excitatory synaptic reversal potential
  g_Ca_nmda(25.0),        // mS/cm2
  g_syn_nmda(0.3),        // mS/cm2
  E_Ca_nmda(140.0),       // mV   Calcium current reversal potential
  E_syn_nmda(0.0),        // mV   NMDA channel total current reversal potential
  Mg_mol(2.0)             // mM
{
  h_nmda_fast_const = (tau_nmda_fast - tau_nmda_rise) / (tau_nmda_fast * tau_nmda_rise);
  h_nmda_slow_const = (tau_nmda_slow - tau_nmda_rise) / (tau_nmda_slow * tau_nmda_rise);
  Na_Kdr_spine_obj.p_.gnabar = 7.0;
  Na_Kdr_spine_obj.p_.gkbar = 0.867;
  Na_Kdr_spine_obj.init(E_L);
  
}

double CA1_params::I_Ca_nmda_func(double s_nmda, double v_spine)
{
  double m_Ca_nmda = m_Ca_nmda_func(v_spine);
  double I_Ca_nmda = -g_Ca_nmda*s_nmda*m_Ca_nmda*(v_spine - E_Ca_nmda);
  return I_Ca_nmda;
}
  
double CA1_params::I_syn_nmda_func(double s_nmda, double v_spine)
{
  double m_syn_nmda = m_syn_nmda_func(v_spine);
  double I_syn_nmda = -g_syn_nmda*s_nmda*m_syn_nmda*(v_spine - E_syn_nmda);
  return I_syn_nmda;
}
double CA1_params::m_Ca_nmda_func(double v_spine)
{
  double m_Ca_nmda =  1.0 / (1.0 + 0.3*Mg_mol*exp(-0.124*v_spine));
  return m_Ca_nmda;
}

double CA1_params::m_syn_nmda_func(double v_spine)
{
  double m_syn_nmda = 1.0 / (1.0 + 0.3*Mg_mol*exp(-0.062*v_spine));
  return m_syn_nmda;
}

double CA1_params::g_nmda_func(double g_nmda_fast, double g_nmda_slow)
{
  double g_nmda = alpha_nmda*g_nmda_fast + (1.0 - alpha_nmda)*g_nmda_slow;
  return g_nmda;
}


int CA1_dynamics(double t, const double y[], double f[], void *params)
{
  CA1_params& np = *(reinterpret_cast< CA1_params* >(params));
  const double& V = y[0];
  //const double& w = y[1];
  const double& h_nmda = y[2];
  const double& g_nmda_fast = y[3];
  const double& g_nmda_slow = y[4];
  
  //double I_spike = np.g_L*np.Delta_T*std::exp((V - np.V_th)/np.Delta_T);
  
  //double I_nmda = g_nmda * ( V - np.E_exc );
  double g_nmda = np.g_nmda_func(g_nmda_fast, g_nmda_slow);
  double I_syn_nmda = np.I_syn_nmda_func(g_nmda, V);
  double I_Na = -np.Na_Kdr_spine_obj.ina;
  double I_K = -np.Na_Kdr_spine_obj.ik;
  
  //f[0] = ( -np.g_L*(V - np.E_L) + I_spike + I_syn_nmda - w + np.I_e) / np.C_m;
  //f[0] = 0.5;
  //f[0] = -V-2*V;
  f[0] = (0.05*np.I_e-np.g_L*(V - np.E_L)+I_Na+I_K)/ np.C_m;//( -np.g_L*(V - np.E_L) + I_Na + I_K + I_syn_nmda + np.I_e) / np.C_m;
  f[1] = 0.0; //( np.a*(V - np.E_L) - w ) / np.tau_w;
  f[2] = -h_nmda / np.tau_nmda_rise;
  f[3] = -g_nmda_fast / np.tau_nmda_fast + np.h_nmda_fast_const*h_nmda;
  f[4] = -g_nmda_slow / np.tau_nmda_slow + np.h_nmda_slow_const*h_nmda;

  return GSL_SUCCESS;
}

// costruttore della classe del modello di neurone CA1
CA1::CA1()
{
  //dim_ = 5; // numero di variabili di stato

  // inizializzazione degli oggetti usati dall'ODE
  // (Ordinary Differential Equations) solver
  // Il metodo di soluzione è Runge Kutta al 4o ordine adattivo
  const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
  s_ = gsl_odeiv2_step_alloc(T, dim_);
  c_ = gsl_odeiv2_control_yp_new(1.0e-6, 1.0e-6);
  e_ = gsl_odeiv2_evolve_alloc(dim_);
  sys_ = {CA1_dynamics, NULL, dim_, &params_};
  d_ = gsl_odeiv2_driver_alloc_y_new (&sys_, T, 1.0e-6, 1.0e-6, 1.0e-6);
  gsl_odeiv2_step_set_driver(s_, d_);
  h_ = 1.0e-2; // valore iniziale del passo di integrazione adattivo

  // valori iniziali delle variabili di stato
  state_var_[0] = params_.E_L; // membrane potential at equilibrium
  state_var_[1] = 0.0;
  state_var_[2] = 0.0;
  state_var_[3] = 0.0;
  state_var_[4] = 0.0;
  syn_input_ = 0.0;

  //print_flag_[0] = false;
  //print_flag_[1] = false;
  //print_flag_[2] = false;
  //print_flag_[3] = false;
  //print_flag_[4] = false;
}

// metodo di update della dinamica tra il tempo t e t + t_step
int CA1::Update(double t, double t_step)
{
  double tf = t + t_step; // tempo alla fine del passo temporale
  // La corrente sinaptica viene incrementata di una quantità pari all'input
  // sinaptico
  state_var_[2] += syn_input_;
  syn_input_ = 0.0;
  
  // Il solutore delle equazioni differenziali stima le variabili di stato
  // con un passo di integrazione h_ adattivo, regolato in base alla tolleranza
  // specificata. Se l'errore sulla stima è maggiore della tolleranza il
  // passo viene ridotto altrimenti viene incrementato
  // Ciclo while che si ripete finché non si raggiunge
  // il tempo finale dell'intervallo considerato
  while (t < tf) {
    // funzione di evoluzione temporale delle variabili di stato
    int status = gsl_odeiv2_evolve_apply (e_, c_, s_, &sys_, &t, tf, &h_,
					  state_var_);
     // controlla se viene superata la soglia di emissione dello spike
    //if (state_var_[0] > params_.V_peak_ ) {
    //Spike();
       // in questo caso il potenziale di membrana viene resettato a V_reset_
      //state_var_[0] = params_.V_reset_;
      //state_var_[1] += params_.b_; // spike-driven adaptation
      // anche gli oggetti usati dall'ODE solver vengono resettati
      //gsl_odeiv2_step_reset(s_);
      //gsl_odeiv2_evolve_reset(e_);
    //}	
    //else
    if (status != GSL_SUCCESS) {
      printf ("error, return value=%d\n", status);
      exit(0);
    }
  }
  
  params_.Na_Kdr_spine_obj.Update(t, t_step);
  params_.Na_Kdr_spine_obj.v = state_var_[0];
  
  //printf("%.1f\t%.6e\t%.6e\n", t, state_var_[0], state_var_[1]);
  //for (uint i=0; i<dim_; i++) {
  //  if (print_flag_[i]) {
  //    printf("\t%.6e", state_var_[i]);
  //  }
  //}
  
  return 0;
}







