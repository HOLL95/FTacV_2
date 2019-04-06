#include </users/henney/Documents/Oxford/C++_libraries/pybind11/include/pybind11/pybind11.h>
#include </users/henney/Documents/Oxford/C++_libraries/pybind11/include/pybind11/stl.h>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>
#include <math.h>
#include <iostream>
#include <exception>
namespace py = pybind11;

struct e_surface_fun {
    double E,dE,Edc;
    double Cdl,CdlE,CdlE2,CdlE3;
    double E0;
    double Ru,R;
    double k0;
    double alpha;
    double In0,u1n0;
    double dt;
    double gamma;

    double exp11,exp12;
    double dexp11,dexp12;
    double u1n1;
    double du1n1;
    double Cdlp;

    e_surface_fun (
                    const double E,
                    const double Edc,
                    const double dE,
                    const double Cdl,
                    const double CdlE,
                    const double CdlE2,
                    const double CdlE3,
                    const double E0,
                    const double Ru,
                    const double R,
                    const double k0,
                    const double alpha,
                    const double In0,
                    const double u1n0,
                    const double dt,
                    const double gamma

                    ) :
        E(E),Edc(Edc),dE(dE),Cdl(Cdl),CdlE(CdlE),CdlE2(CdlE2),CdlE3(CdlE3),E0(E0),Ru(Ru),R(R),k0(k0),alpha(alpha),In0(In0),u1n0(u1n0),dt(dt),gamma(gamma) { }

    boost::math::tuple<double,double> operator()(const double In1) {
        update_temporaries(In1);
        return boost::math::make_tuple(residual(In1),residual_gradient(In1));
    }

    double residual(const double In1) const {
        return Cdlp*(dt*dE-Ru*(In1-In0)) + dt*R*(E-Ru*In1) - dt*In1 + gamma*(u1n1-u1n0);
        //return Cdlp*(dt*dE) - dt*In1 + (u1n1-u1n0) + Ru*E*dt;
    }
    double residual_gradient(const double In1) const {
        return -Cdlp*Ru - dt*R*Ru - dt + gamma*du1n1;
        //return -Cdlp*Ru - dt + du1n1;
    }

    void update_temporaries(const double In1) {
        const double Ereduced = E - Ru*In1;
        //const double Ereduced = E;
        const double Ereduced2 = pow(Ereduced,2);
        const double Ereduced3 = Ereduced*Ereduced2;
        const double expval1 = Ereduced - E0;
        exp11 = std::exp((1.0-alpha)*expval1);
        exp12 = std::exp(-alpha*expval1);

        dexp11 = -Ru*(1.0-alpha)*exp11;
        dexp12 = Ru*alpha*exp11;

        const double u1n1_top = dt*k0*exp11 + u1n0;
        const double du1n1_top = dt*k0*dexp11;
        const double denom = (dt*k0*exp11 + dt*k0*exp12 + 1);
        const double ddenom = dt*k0*(dexp11 + dexp12);
        const double tmp = 1.0/denom;
        const double tmp2 = pow(tmp,2);
        u1n1 = u1n1_top*tmp;
        du1n1 = -(u1n1_top*ddenom + du1n1_top*denom)*tmp2;

        Cdlp = Cdl*(1.0 + CdlE*Ereduced + CdlE2*Ereduced2 + CdlE3*Ereduced3);
        //Cdlp = Cdl*(1.0 + CdlE*Edc+ CdlE2*pow(Edc,2)+ CdlE3*pow(Edc,3));
    }
};
double et(double E_start, double omega, double phase, double delta_E,double t){
	double E_t=(E_start+delta_E)+delta_E*(std::sin((omega*t)+phase));

	return E_t;
}

double dEdt(double omega, double phase, double delta_E, double t){
  double dedt=(delta_E*omega*std::cos(omega*t+phase));
	return dedt;
}
py::object martin_surface(double Cdl, double CdlE, double CdlE2, double CdlE3, double omega, double phase, double pi, double alpha, double Estart, double Ereverse, double delta_E, double Ru, double gamma,double E0, double k0, double final_val, std::vector<double> t) {
    const double R = 0;
    const int Ntim = 200.0;

    const int digits_accuracy = std::numeric_limits<double>::digits*0.5;
    const double max_iterations = 100;
    const double dt = (1.0/Ntim)*2*pi/omega;
    const double Tmax = final_val;
    const int Nt = Tmax/dt;
    std::vector<double> Itot;
    if(t.size() == 0){
      Itot.resize(Nt,0);
      t.resize(Nt);
      for (int i=0; i<Nt; i++) {
          t[i] = i*dt;
          }
        } else {
      #ifndef NDEBUG
              //std::cout << "\thave "<<times.size()<<" samples from "<<times[0]<<" to "<<times[times.size()-1]<<std::endl;
      #endif
              Itot.resize(t.size(),0);
          }
    double Itot0,Itot1;
    double u1n0;
    double t1 = 0.0;
    u1n0 = 1.0;

    const double E = et(Estart, omega, phase,delta_E ,t1+dt);
    const double dE = dEdt(omega, phase,delta_E , t1+0.5*dt);
    const double Cdlp = Cdl*(1.0 + CdlE*E + CdlE2*pow(E,2)+ CdlE3*pow(E,2));
    const double Itot_bound =100000;//std::max(10*Cdlp*delta_E*omega/Nt,1.0);
    //std::cout << "Itot_bound = "<<Itot_bound<<std::endl;

    Itot0 =Cdlp*dE;
    Itot1 = Itot0;
    for (int n_out = 0; n_out < t.size(); n_out++) {
        while (t1 < t[n_out]) {
            Itot0 = Itot1;
            const double E = et(Estart, omega, phase,delta_E ,t1+dt);
            const double dE = dEdt(omega, phase,delta_E , t1+0.5*dt);
            const double Edc = 0.0;
            e_surface_fun bc(E,Edc,dE,Cdl,CdlE,CdlE2,CdlE3,E0,Ru,R,k0,alpha,Itot0,u1n0,dt,gamma);
            boost::uintmax_t max_it = max_iterations;
            Itot1 = boost::math::tools::newton_raphson_iterate(bc, Itot0,Itot0-Itot_bound,Itot0+Itot_bound, digits_accuracy, max_it);
            if (max_it == max_iterations) throw std::runtime_error("non-linear solve for Itot[n+1] failed, max number of iterations reached");

            bc.update_temporaries(Itot1);
            u1n0 = bc.u1n1;
            t1 += dt;
        }
        Itot[n_out] = (Itot1-Itot0)*(t[n_out]-t1+dt)/dt + Itot0;
    }
    return py::cast(Itot);
}
PYBIND11_MODULE(isolver_martin, m) {
	m.def("martin_surface", &martin_surface, "solve for I_tot with dispersion");
}
