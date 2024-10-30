#include <iostream>
#include <cstdio>
#include "CA1.h"

int main()
{
  //size_t dim = 5;

  double t_step = 0.0005; // time step in ms
  double t0 = 0.0;
  double t1 = 300.0; // end of integration interval in ms
  //double t = t0; // time
  //double h = 1.0e-6;
  //double h = 1.0e-3;
     
  CA1 CA1_obj;
  //CA1_obj.print_flag_[0] = true;
  //CA1_obj.print_flag_[2] = true;
  //CA1_obj.print_flag_[3] = true;
  //CA1_obj.print_flag_[4] = true;

  double w = 1.0;
  
  int n = (int)round((t1 - t0) / t_step);

  for (int i=0; i<n; i++) {
    double t = t0 + t_step*i;
    //double E_L = CA1_obj.params_.E_L;
    double V = CA1_obj.state_var_[0];
    double g_nmda_fast = CA1_obj.state_var_[3];
    double g_nmda_slow = CA1_obj.state_var_[4];
    double g_nmda = CA1_obj.params_.g_nmda_func(g_nmda_fast, g_nmda_slow);
    double I_syn_nmda = CA1_obj.params_.I_syn_nmda_func(g_nmda, V);
    double I_Ca_nmda = CA1_obj.params_.I_Ca_nmda_func(g_nmda, V);
    double I_Na_Kdr_spine = CA1_obj.params_.Na_Kdr_spine_obj.ina;
    double ik = CA1_obj.params_.Na_Kdr_spine_obj.ik;
    double il = CA1_obj.params_.Na_Kdr_spine_obj.il;
    double m = CA1_obj.params_.Na_Kdr_spine_obj.m;
    double h = CA1_obj.params_.Na_Kdr_spine_obj.h;
    double n = CA1_obj.params_.Na_Kdr_spine_obj.n;
    double s = CA1_obj.params_.Na_Kdr_spine_obj.s;

    printf("%.3f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", t, V, 
	   g_nmda_fast, g_nmda_slow, g_nmda, I_syn_nmda, I_Ca_nmda,
	   I_Na_Kdr_spine,ik,il,m,h,n,s);

    unsigned int n_spikes; // = gsl_ran_poisson(r, mu);
    if (i==110) {
      n_spikes = 1;
    }
    else {
      n_spikes = 0;
    }
    
    CA1_obj.syn_input_ = w*n_spikes;
    CA1_obj.Update(t, t_step);    
  }

  return 0;
}

