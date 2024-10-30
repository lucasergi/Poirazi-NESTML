#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "Na_Kdr_spine.h"

Na_Kdr_spine_params::Na_Kdr_spine_params()
{
  a0r = 0.0003; // (ms)
  b0r = 0.0003; // (ms)
  zetar = 12;
  zetas = 12;   
  gmr = 0.2;
  ar2 = 1.0;       // initialized parameter for location-dependent
                          // Na-conductance attenuation, "s", (ar=1 -> zero attenuation)
  taumin = 3.0;    // (ms) min activation time for "s" attenuation system
  vvs = 2.0;       // (mV) slope for "s" attenuation system
  vhalfr = -60.0;  // (mV) half potential for "s" attenuation system
  W = 0.016;       // (/mV) this 1/61.5 mV
  gnabar = 0.2; // (mho/cm2)  :suggested conductance values
  gkbar = 0.12; // (mho/cm2)
  gl = 0.1;  //  (mho/cm2)
  // gnabar = 0;      // (mho/cm2)  initialized conductances
  // gkbar = 0;       // (mho/cm2)  :actual values set in cell-setup.hoc
  // gl = 0;          // (mho/cm2)
  ena = 60.0;      // (mV) Na reversal potential (also reset in
  ek = -77.0;      // (mV) K reversal potential  cell-setup.hoc)
  el = -70.0;      // (mV) steady state 
  celsius = 34.0;    // (degC)

}



// inizializza l'oggetto della classe
void Na_Kdr_spine::init(double v0)
{
  // valori iniziali delle variabili di stato
  dt = 0.01;
  v = v0;
  mhn(v);
  //m = 0.0;
  //h = 1.0;
  //n = 0.0;
  //s = 1.0;
  m = inf[0];  // Na activation variable
  h = inf[1];  // Na inactivation variable
  n = inf[2];  // K activation variable
  s = inf[3];  // Na attenuation variable

  ina = p_.gnabar*m*m*h*s*(v - p_.ena);
  ik = p_.gkbar*n*n*(v - p_.ek);
  il = p_.gl*(v - p_.el);  
}

// costruttore della classe del modello
Na_Kdr_spine::Na_Kdr_spine()
{
  init(-70.0);
}

// metodo di update della dinamica tra il tempo t e t + t_step
int Na_Kdr_spine::Update(double t, double t_step)
{
  dt = t_step;
  //double tf = t + t_step; // tempo alla fine del passo temporale

  //printf("%.3f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", t, m, h, s, n, ina, ik, il);

  states();
  ina = p_.gnabar*m*m*h*s*(v - p_.ena);
  ik = p_.gkbar*n*n*(v - p_.ek);
  il = p_.gl*(v - p_.el);
  
  return 0;
}

void Na_Kdr_spine::calcg() {
  mhn(v);
  m = m + fac[0]*(inf[0] - m);  // Na activation variable
  h = h + fac[1]*(inf[1] - h);  // Na inactivation variable
  n = n + fac[2]*(inf[2] - n);  // K activation variable
  s = s + fac[3]*(inf[3] - s);  // Na attenuation variable
}	

void Na_Kdr_spine::states() {	// exact when v held constant
  calcg();
}

double Na_Kdr_spine_params::alpv(double v, double vh) { // used in "s" activation system infinity calculation
  double alpv_val = (1.0 + ar2*exp((v-vh)/vvs))/(1.0+exp((v-vh)/vvs));
  return alpv_val;
}


double Na_Kdr_spine_params::varss(double v, int i) { // steady state values
  double varss_val;
  if (i==0) {
    varss_val = 1.0 / (1.0 + exp((v + 40.0)/(-3.0)));    // Na activation
  }
  else if (i==1) {
    varss_val = 1.0 / (1.0 + exp((v + 45.0)/(3)));   // Na inactivation 
  }
  else if (i==2) {	
    varss_val = 1.0 / (1.0 + exp((v + 42.0)/(-2.0))); // K activation
  }
  else { //activation system for spike attenuation - Migliore 96 model
    varss_val = alpv(v,vhalfr);
  }
  return varss_val;
}


double Na_Kdr_spine_params::alpr(double v) {       // used in "s" activation system tau
  double alpr_val = exp(1.e-3*zetar*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius)));
  return alpr_val;
}

double Na_Kdr_spine_params::betr(double v) {       // used in "s" activation system tau
  double betr_val = exp(1.e-3*zetar*gmr*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius)));
  return betr_val;
}

double Na_Kdr_spine_params::vartau(double v, int i) { // estimate tau values
  double tmp;
  double vartau_val;
  
  if (i==0) {
    vartau_val = 0.05;  // Na activation tau
  }
  else if (i==1) {
    vartau_val = 0.5;   // Na inactivation tau
  }
  else if (i==2) {
    vartau_val = 2.2;   // K activation
  } else {
    tmp = betr(v)/(a0r+b0r*alpr(v));
    if (tmp<taumin) {
      tmp=taumin;
    }
    vartau_val = tmp;   // s activation tau
  }
  return vartau_val;
}

void Na_Kdr_spine::mhn(double v) {
  // double a, b; // rest = -70
  //       TABLE infinity, tau, fac DEPEND dt, celsius FROM -100 TO 100 WITH 200
  for(int i=0; i<=3; i++)  {
    tau[i] = p_.vartau(v,i);
    inf[i] = p_.varss(v,i);
    fac[i] = (1.0 - exp(-dt/tau[i]));
  }
}

