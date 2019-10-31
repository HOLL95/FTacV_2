#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <math.h>
#include <iostream>
#include <exception>
#include <vector>
#include <time.h>
#include "boost/tuple/tuple.hpp"
#include <boost/math/tools/roots.hpp>
#include </users/henney/Documents/Oxford/C++_libraries/pybind11/include/pybind11/pybind11.h>
#include </users/henney/Documents/Oxford/C++_libraries/pybind11/include/pybind11/stl.h>
#include </users/henney/Documents/Oxford/FTacV/single_electron/input_signal.h>
#include <string>
namespace py = pybind11;
//using namespace boost::python;

struct param1{
	double delta_E;
	double E_start;
	double E_reverse;
	double omega;
	double phase;
	double a;
	double v;
	double tr;
	double E0_sigma;
	double k0_sigma;
	int time_end;
	int duration;
};
struct single_e{
	double E,dE;
	double Cdl;
	double CdlE1;
	double CdlE2;
	double CdlE3;
	double E0_mean;
	double Ru;
	double R;
	double k0_mean;
	double alpha;
	double I_0;
	double theta_0;
	double dt;
	double gamma;
	double exp11, exp12;
	double dexp11, dexp12;
	double dtheta_1;
	double theta_1;
	double Cdlp;

	single_e (
		const double E,
		const double dE,
		const double Cdl,
		const double CdlE1,
		const double CdlE2,
		const double CdlE3,
		const double E0_mean,
		const double Ru,
		const double R,
		const double k0_mean,
		const double alpha,
		const double I_0,
		const double theta_0,
		const double dt,
		const double gamma):

		E(E),dE(dE),Cdl(Cdl),CdlE1(CdlE1),CdlE2(CdlE2),CdlE3(CdlE3),E0_mean(E0_mean),Ru(Ru),R(R),k0_mean(k0_mean),alpha(alpha),I_0(I_0),theta_0(theta_0),dt(dt),gamma(gamma) { }


boost::math::tuple<double,double> operator()(const double I_1) {
	val_iterate(I_1);
        return boost::math::make_tuple(residual(I_1),residual_gradient(I_1));
	}


void val_iterate(double I_1){
		double Er = E - (Ru*I_1);
		double expval1 = Er - E0_mean;
		exp11 = std::exp((1.0-alpha)*expval1);
		exp12 = std::exp(-alpha*expval1);
		dexp11 = -Ru*(1.0-alpha)*exp11;
		dexp12 = Ru*alpha*exp11;
		double theta_1_top = dt*k0_mean*exp11 + theta_0;
		double dtheta_1_top = dt*k0_mean*dexp11;
		double denom = (dt*k0_mean*exp11 + dt*k0_mean*exp12 + 1);
		double ddenom = dt*k0_mean*(dexp11 + dexp12);
		double tmp = 1.0/denom;
		double tmp2 = pow(tmp,2);
		dtheta_1 = -(theta_1_top*ddenom + dtheta_1_top*denom)*tmp2;
		theta_1 = theta_1_top*tmp;
		double Er2=Er*Er;
		double Er3=Er2*Er;
		Cdlp=Cdl*(1+CdlE1*Er+CdlE2*Er2+CdlE3*Er3);
	}
double residual(double I_1) {
        	return Cdlp*(dt*dE-Ru*(I_1-I_0))+dt*R*(E-Ru*I_1) - dt*I_1 +gamma*(theta_1-theta_0)  ;//
	}
double residual_gradient(double I_1) {
	        return -Cdlp*Ru-dt*R*Ru-dt +gamma*dtheta_1 ;//

	}
};




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





py::object I_tot_solver(double Cdl, double CdlE1, double CdlE2, double CdlE3,double omega,double v,double alpha ,double E_start, double E_reverse, \
	double delta_E, double Ru, double dt, std::vector<double> time_vec, double gamma, double E0_mean, double k0_mean, double phase){
	noramp tf;
	time_param times;
	double E0_sigma =1;
	double k0_sigma=1;
	times.v=v;
	times.delta_E=delta_E;
	times.E_start=E_start;
	times.E_reverse=E_reverse;
	times.omega=omega;
	times.phase=phase;
	times.duration=time_vec.size();
	times.tr=((times.E_reverse-times.E_start)/times.v);
	times.time_vec=time_vec;
	double Itot_bound = 10000;
	double R=0;
	float t=0;
	double E=tf.et(times,t);
	double dE=tf.dEdt(times, t+0.5*dt);
	double theta_0=1.0;
	double Cdlp=Cdl*(1.0+(CdlE1*E)+(CdlE2*pow(E,2))+(CdlE3*pow(E,3)));
	double I_0=Cdlp*dE;
	double I_1=I_0;
	std::vector<std::vector<double>> I_matrix(2, std::vector<double>(times.duration, 0));
	const int digits_accuracy = std::numeric_limits<double>::digits;
	I_matrix[0][0]=0;
	for(int j=0; j<times.duration; j++){
		while(t<time_vec[j]){
			I_0=I_1;
			E=tf.et(times,t);
			dE=tf.dEdt(times,t+0.5*dt);//
			single_e init(E,dE,Cdl,CdlE1,CdlE2,CdlE3,E0_mean,Ru,R,k0_mean,alpha,I_0 ,theta_0,dt,gamma);
			I_1=boost::math::tools::newton_raphson_iterate(init, I_0,I_0-Itot_bound,I_0+Itot_bound, digits_accuracy);
			init.val_iterate(I_1);
			theta_0=init.theta_1;
			t+=dt;
		}
		I_matrix[1][j] = (I_1 - I_0) * (time_vec[j] - t + dt) / dt + I_0;
		I_matrix[0][j]=time_vec[j];
	}
	return py::cast(I_matrix);
}

///////////////////////////////////////////////////////////////////////////////DISPERSION SUBROUTINES//////////////////////////////////////////////////////////////////////////////////////
/*/
std::vector<std::vector<double>>weightmatrix(param1& single_e_param,int bins, std::vector<double> E0_disp,std::vector<double> k0_disp ){
	double weight_E;
	double weight_k;
	double x1;
	boost::math::normal_distribution<double> E0(single_e_param.E0_mean, single_e_param.E0_sigma);
	boost::math::lognormal_distribution<double> k0(single_e_param.k0_mean, single_e_param.k0_sigma);
	std::vector< std::vector<double>> weight_matrix((bins), std::vector< double >((bins)));
	for(int i=0; i<bins; i++){
		if(i==0){
			weight_E=cdf(E0,E0_disp[0]);
		}else{
			weight_E=cdf(E0,E0_disp[i])-cdf(E0,E0_disp[i-1]);
		}

		for(int j=0; j<bins; j++){


			if(j==0){
				weight_k=cdf(k0, k0_disp[0]);
			}else{
				weight_k=cdf(k0,k0_disp[j])-cdf(k0, k0_disp[j-1]);
			}
			x1=x1+weight_E*weight_k;
			weight_matrix[i][j]=weight_E*weight_k;
		}

	}

	return weight_matrix;

}

std::vector<double> dispersion_solver(param1& single_e_param, std::vector<std::vector<double>>weight_matrix,std::vector<double> E0_disp,std::vector<double> k0_disp ,int bins){
	std::vector<double> I_matrix(single_e_param.duration);
	std::vector<double> I_disp(single_e_param.duration);
	for(int i=0; i<(bins); i++){
		single_e_param.E0=E0_disp[i];
		for(int j=0; j<bins; j++){
		 	single_e_param.k0=k0_disp[j];
			I_matrix=non_linear_I_solver(single_e_param);
	 		for(int k=0; k<single_e_param.duration; k++){
				I_disp[k]=I_disp[k]+(I_matrix[k]*weight_matrix[i][j]);
			}
		}
	}
	return I_disp;
}





py::object I_tot_solver(double Cdl, double CdlE1, double CdlE2, double CdlE3,double omega,double v,double alpha ,double E_start, double E_reverse, double delta_E, double Ru, double dt, double time_end, int time_length, double E0_mean=0, double k0_mean=0, double E0_sigma=0.1, double k0_sigma=1) {
	int bins=16;
	std::vector<double>I_disp(single_e_param.duration,0);
	std::vector<double>E0_disp(bins,0);
	std::vector<double>k0_disp(bins,0);
	float E0_interval=(abs(-0.4-0.4))/bins;
	float k0_interval=(0+50)/bins;
	E0_disp[0]=-0.3;
	k0_disp[0]=0;
	for(int i=1; i<bins; i++){
		E0_disp[i]=E0_disp[i-1]+E0_interval;
		k0_disp[i]=k0_disp[i-1]+k0_interval;
	}
	std::vector< std::vector< double > > weight_matrix( bins, std::vector< double >( bins ) );
	//weight_matrix=weightmatrix(single_e_param, bins, E0_disp,k0_disp);
	//I_disp=dispersion_solver(single_e_param, weight_matrix, E0_disp, k0_disp, bins);
	I_disp=non_linear_I_solver(single_e_param);
	return py::cast(I_disp);

}
/*/
PYBIND11_MODULE(isolver, m) {
	m.def("I_tot_solver", &I_tot_solver, "solve for I_tot with dispersion");

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(){
	//std::vector<double>I_disp(800,0);
	//I_disp=I_tot_solver(0.000133878548046, 0.000653657774506,0.000245772700637,1.10053945995e-06,boost::math::constants::pi<double>()*2,0,0,0.1,1,40,800);

	//clock_t t;
	//t = clock();
	//I_matrix=non_linear_I_solver(single_e_pa);
	//for(int i=0; i<800; i++){
	//	std::cout<<I_disp[i] << "\n";
	//}
	//t= clock() - t;
	//std::cout<<(t*1.0/CLOCKS_PER_SEC)<<"\n";	//


}
