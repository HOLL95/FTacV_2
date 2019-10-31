#include <math.h>
#include <iostream>
#include <exception>
#include <vector>
#include </users/henney/Documents/Oxford/FTacV/single_electron/input_signal.h>
using namespace std;
double classical::et(time_param &time_init, double t){
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

double classical::dEdt(time_param &time_init, double t){
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
double noramp::et(time_param &time_init, double t){
	double amp=abs(time_init.E_start-time_init.E_reverse)/2;
	double E_t;


	E_t=(time_init.E_start+amp)+amp*(std::sin((time_init.omega*t)+time_init.phase));

	return E_t;
}

double noramp::dEdt(time_param &time_init, double t){
	double  amp=abs(time_init.E_start-time_init.E_reverse)/2;
	double dedt;

   	dedt=(amp*time_init.omega*std::cos(time_init.omega*t+time_init.phase));

	return dedt;
}
