#ifndef NA_KDR_SPINE_H
#define NA_KDR_SPINE_H

#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

class Na_Kdr_spine_params
{
public:
  double a0r; // 0.0003 (ms)
  double b0r; // 0.0003 (ms)
  double zetar; // 12;
  double zetas; // 12;   
  double gmr; // 0.2;
  double ar2; // 1.0;     // initialized parameter for location-dependent
                          // Na-conductance attenuation, "s", (ar=1 -> zero attenuation)
  double taumin; // 3.0 (ms) min activation time for "s" attenuation system
  double vvs; // 2.0 (mV) slope for "s" attenuation system
  double vhalfr; // -60.0 (mV) half potential for "s" attenuation system
  double W; // 0.016 (/mV) this 1/61.5 mV
  double gnabar; // 0.2 (mho/cm2)  :suggested conductance values
  double gkbar; // 0.12 (mho/cm2)
  double gl; // 0.0001 (mho/cm2)
  // double gnabar; // 0 (mho/cm2)  initialized conductances
  // double gkbar; // 0 (mho/cm2)  :actual values set in cell-setup.hoc
  // double gl; // 0 (mho/cm2)
  double ena; // 60.0 (mV) Na reversal potential (also reset in
  double ek; // -77.0 (mV) K reversal potential  cell-setup.hoc)
  double el; // -70.0 (mV) steady state 
  double celsius; // 34.0 (degC)
  
  Na_Kdr_spine_params(); // default constructor

  double alpv(double v, double vh); // used in "s" activation system infinity calculation
  
  double varss(double v, int i); // steady state values
  
  double alpr(double v);       // used in "s" activation system tau
  
  double betr(double v);       // used in "s" activation system tau
  
  double vartau(double v, int i); // estimate tau values

};

class Na_Kdr_spine
{
public:
  //static const size_t dim_ = 0; // numero di variabili di stato risolte con ODE solver
  Na_Kdr_spine_params p_; // parametri
  // oggetti usati dall'ODE solver
  //gsl_odeiv2_step *s_;
  //gsl_odeiv2_control *c_;
  //gsl_odeiv2_evolve *e_;
  //gsl_odeiv2_system sys_;
  //gsl_odeiv2_driver * d_;
  //double h_; // passo di integrazione adattivo
  //double state_var_[dim_]; // variabili di stato risolte con ODE solver
  //bool print_flag_[dim_];

  double dt;   // (ms) time resolution
  double v;    // (mV) membrane potential, contant for now
  
  double m;
  double h;
  double n;
  double s;
  
  double ina;  // (mA/cm2)
  double ik;   // (mA/cm2)
  double il;   // (mA/cm2)
  double inf[4];
  double fac[4];
  double tau[4];
  
  Na_Kdr_spine(); // costruttore
  void init(double v0); // initialization method
  
  // metodo di update della dinamica tra il tempo t e t + t_step
  int Update(double t, double t_step);

  void calcg();
  void mhn(double v);
  void states();
};

#endif


