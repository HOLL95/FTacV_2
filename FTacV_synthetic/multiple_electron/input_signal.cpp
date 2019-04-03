#include <math.h>
#include <iostream>
#include <exception>
#include <vector>
#include </auto/users/henney/Documents/Oxford/FTacV/multiple_electron/input_signal.h>
using namespace std;
double classical::et(time_param &time_init, float t){
	double E_dc;
	double E_t;
	if (t<time_init.tr){
		E_dc=time_init.E_start+(time_init.v*t);
	}else {
		E_dc=time_init.E_reverse-(time_init.v*(t-time_init.tr));
	}

	 E_t=E_dc+(time_init.delta_E*(sin((time_init.omega*t)+time_init.phase)));

	return E_t;
}

double classical::dEdt(time_param &time_init, float t){
	double E_dc;
	double dedt;
	if (t < time_init.tr){
		 E_dc=time_init.v;
	}else {
		 E_dc=-time_init.v;
	}
   	dedt=E_dc+(time_init.delta_E*time_init.omega*cos(time_init.omega*t+time_init.phase));

	return dedt;
}
double noramp::et(time_param &time_init, float t){
	double e_time=(time_init.E_start+time_init.amp)+time_init.amp*(std::sin((time_init.omega*t)+time_init.phase));
	return e_time;
}

double noramp::dEdt(time_param &time_init, float t){
   	double dedt=(time_init.amp*time_init.omega*std::cos(time_init.omega*t+time_init.phase));
	return dedt;
}
