#ifndef CA1H
#define CA1H

#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "Na_Kdr_spine.h"

class CA1_params
{
public:
  //double V_peak;
  //double V_reset;
  double g_L;
  double C_m;
  double E_L;
  //double Delta_T;
  //double tau_w;
  //double a;
  //double b;
  //double V_th;
  double I_e;
  //double tau_syn;
  double tau_nmda_rise; // 2.0 ms  Rise time for NMDA channels
  double tau_nmda_fast; // 10.0 ms Fast decay time for NMDA channels
  double tau_nmda_slow; // 45.0 ms Slow decay time for NMDA channels
  double alpha_nmda;    // 0.527   ratio of fast decay NMDA conducance over total NMDA conductance 
  double refr_T;        // 2 ms    Duration of refractory period
  //double E_exc;         // 0 mV    Excitatory synaptic reversal potential
  double h_nmda_fast_const;  // 1 / ms
  double h_nmda_slow_const;  // 1 / ms
  double g_Ca_nmda;     // 25.0 mS/cm2
  double g_syn_nmda;    // 0.3 mS/cm2
  double E_Ca_nmda;     // 140.0 mV Calcium current reversal potential
  double E_syn_nmda;    // 0.0 mV   NMDA channel tot current reversal potential
  double Mg_mol;        // 2.0 mM

  // Objects for different currents
  Na_Kdr_spine Na_Kdr_spine_obj;
  
  CA1_params(); // default constructor

  double I_Ca_nmda_func(double s_nmda, double v_spine);
  double I_syn_nmda_func(double s_nmda, double v_spine);
  double m_Ca_nmda_func(double v_spine);
  double m_syn_nmda_func(double v_spine);
  double g_nmda_func(double g_nmda_fast, double g_nmda_slow);

};

class CA1
{
 public:
  static const size_t dim_ = 5; // numero di variabili di stato
  CA1_params params_; // parametri


  // oggetti usati dall'ODE solver
  gsl_odeiv2_step *s_;
  gsl_odeiv2_control *c_;
  gsl_odeiv2_evolve *e_;
  gsl_odeiv2_system sys_;
  gsl_odeiv2_driver * d_;
  double h_; // passo di integrazione adattivo
  double state_var_[dim_]; // variabili di stato
  //bool print_flag_[dim_];
  double syn_input_;
  
  CA1(); // costruttore

  // metodo di update della dinamica tra il tempo t e t + t_step
  int Update(double t, double t_step);

};

int CA1_dynamics(double t, const double y[], double f[], void *params);

#endif
