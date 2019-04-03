#include <iostream>
#include <vector>
#include </auto/users/henney/Oxford/C++_libraries/Eigen/Eigen/Eigen>
#include </auto/users/henney/Oxford/FTacV/multiple_electron/input_signal.h>
#include </auto/users/henney/Oxford/C++_libraries/pybind11/include/pybind11/pybind11.h>
#include </auto/users/henney/Documents/Oxford/C++_libraries/pybind11/include/pybind11/stl.h>
#include <ctime>
//#include </users/henney/Documents/Single_electron/pybind11/include/pybind11/pybind11.h>
//#include </users/henney/Documents/Single_electron/pybind11/include/pybind11/stl.h>
//namespace py = pybind11;
using namespace std;
using namespace Eigen;
namespace py = pybind11;
/////////////////////////////////////////////////////////////////////////////////STRUCTS////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct param_class{
  double Cdl;
  double CdlE1;
  double CdlE2;
  double CdlE3;
  double omega;
  double v;
  double alpha;
  double E_start;
  double E_reverse;
  double Ru;
  double dt;
  double gamma;
  double Er;
  int num_species;
};
/////////////////////////////////////////////////////////////////////////////////FUNCTIONS///////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<vector<double>> calculate_K(const vector<double>& k_t_ox_vector, const vector<double>& k_t_red_vector, int num_species){
  vector<vector<double>> K_matrix(num_species, vector<double>(num_species, 0));
  K_matrix[0][0]=-k_t_red_vector[0];
  for(int i=0; i<num_species; i++){
    for(int j=0; j<num_species; j++){
      if(i==(num_species-1)){
        K_matrix[i][j]=K_matrix[i][j]-k_t_ox_vector[num_species-1];
      }
      if(((i-j)==0) && (i != 0)){
        K_matrix[i][j]=K_matrix[i][j]-k_t_ox_vector[i-1]-k_t_red_vector[i];
      }
      else if((i-j)==1){
        K_matrix[i][j]=K_matrix[i][j]+k_t_red_vector[j];
      }
      else if((i-j)==-1){
        K_matrix[i][j]=K_matrix[i][j]+k_t_ox_vector[i];
      }
    }
  }
  return K_matrix;
}
vector<double> calculate_e(const vector<double>& k_t_ox_vector,const vector<double>& k_t_red_vector, int num_species){
  vector<double> e_vector(num_species, 0);
  e_vector[0]=-k_t_ox_vector[num_species-1]-k_t_red_vector[0];
  for(int i=1; i<num_species; i++){
    e_vector[i]=k_t_ox_vector[i-1]-k_t_red_vector[i]-k_t_ox_vector[num_species-1];
  }
  return e_vector;
}
vector<double> update_theta(param_class &params,const vector<vector<double>>& K_matrix, const vector<double>& prev_theta, int num_species){
  MatrixXd K_mat(num_species, num_species);
  for (int i = 0; i < num_species; i++){
    K_mat.row(i) = VectorXd::Map(&K_matrix[i][0],K_matrix[i].size());
  }
  VectorXd theta_0 = VectorXd::Map(prev_theta.data(), prev_theta.size());
  VectorXd c=VectorXd::Zero(num_species);
  c(num_species-1)=(-K_matrix[num_species-1][0]);
  //cout<<c(1,num_species)<<"\n";
  MatrixXd Ident=MatrixXd::Identity(num_species, num_species);
  MatrixXd A=Ident - params.dt*K_mat;
  VectorXd update_theta=A.colPivHouseholderQr().solve((theta_0+(params.dt*c)));
  vector<double> theta_1(update_theta.data(), update_theta.data() + update_theta.rows() * update_theta.cols());
  return theta_1;
}

double Capacitance(param_class &params, double I_0, double E){
  //double Er =  E - (params.Ru*I_0);
  double Er2=params.Er*params.Er;
  double Er3=Er2*params.Er;
  double Cdlp=params.Cdl*(1+params.CdlE1*params.Er+params.CdlE2*Er2+params.CdlE3*Er3);
  return Cdlp;
}
vector<vector<double>> reaction_rates(param_class &params,vector<double>& k0_vector, vector<double>& E0_vector, int num_species){
  vector<vector<double>> rate_matrix(2, vector<double>(num_species, 0)); //first row is oxidation rates, second row is reduction
  double expval;
  double exp_ox;
  double exp_red;
  for(int i=0; i<num_species; i++){
    expval=params.Er-E0_vector[i];
    exp_red=exp(-params.alpha*expval);
    exp_ox=exp((1-params.alpha)*expval);
    rate_matrix[0][i]=k0_vector[i]*exp_ox;
    rate_matrix[1][i]=k0_vector[i]*exp_red;
  }
  return rate_matrix;
}
double Faradaic_current(vector<double> e_vector, vector<double> next_theta, double last_k_ox){
  VectorXd E_vec = VectorXd::Map(e_vector.data(), e_vector.size());
  VectorXd theta_1 = VectorXd::Map(next_theta.data(), next_theta.size());
  double updated_e=(E_vec.transpose()*theta_1)+last_k_ox;
  return updated_e;
}
double update_total_current(param_class &params,double I_0, double E, double dE, double Cdlp, double faradaic_current){
  double I_1=Cdlp*(dE+params.Ru*I_0/params.dt);
  I_1=I_1+params.gamma*faradaic_current;
  I_1=I_1/(1+(Cdlp*params.Ru/params.dt));
  return I_1;
}
int vector_printer(const vector<double>& vec, int size){
  for (int j=0; j<size; j++){
    cout<<vec[j]<< " ";
  }
  cout<<"\n";
  return 0;
}
int matrix_printer(const vector<vector<double>>& mat, int size){
  for (int j=0; j<size; j++){
    for (int i=0; i<size; i++){
      cout<<mat[j][i]<< " ";
    }
  cout<<"\n";
  }
  return 0;
}
py::object current_solver(double Cdl, double CdlE1, double CdlE2, double CdlE3,double omega,double v,double alpha ,double E_start, double E_reverse, \
  double delta_E, double Ru, double dt, int duration, double gamma, double phase,int num_species, vector<double> k0_vector, vector<double> E0_vector){
  noramp tf;
  param_class params;
  time_param times;
  times.amp=abs(E_start-E_reverse)/2;
  times.v=v;
  times.delta_E=delta_E;
  times.E_start=E_start;
  times.E_reverse=E_reverse;
  times.tr=abs((times.E_reverse-times.E_start)/times.v);
  times.omega=omega;
  times.phase=phase;
  times.duration=duration;
  //times.time_vec=time_vec;
  params.Cdl=Cdl;
  params.CdlE1=CdlE1;
  params.CdlE2=CdlE2;
  params.CdlE3=CdlE3;
  params.omega=omega;
  params.v=v;
  params.alpha=alpha;
  params.E_start=E_start;
  params.E_reverse=E_reverse;
  params.Ru=Ru;
  params.dt=dt;
  params.gamma=gamma;
  double t=0;
  double E=tf.et(times, t);
	double dE=tf.dEdt(times, t+0.5*dt);
	double Cdlp =params.Cdl * (1.0 + params.CdlE1* E + params.CdlE2 * pow(E, 2) + params.CdlE3 * pow(E, 3));
  std::vector<std::vector<double>> I_matrix(2, std::vector<double>(times.duration, 0));
  double I_f;
  double I_0=0//Cdlp*dE;
  double I_1;
  vector<vector<double>> rate_matrix(2, vector<double>(num_species, 0));
  vector<vector<double>> K_matrix(num_species, vector<double>(num_species, 0));
  vector<double> E_vector(num_species);
  vector<double> theta_0(num_species, 0);
  vector<double> theta_1(num_species, 0);
  theta_0[2]=1;
  for (int i=0; i<times.duration; i++){
    t=dt*i;
    E=tf.et(times,t);
    dE=tf.dEdt(times, t+0.5*dt);
    I_matrix[0][i]=I_0;
    I_matrix[1][i]=t;

    params.Er=E-(params.Ru*I_0);
    rate_matrix=reaction_rates(params, k0_vector, E0_vector, num_species);
    K_matrix=calculate_K(rate_matrix[0], rate_matrix[1], num_species);
    E_vector=calculate_e(rate_matrix[0], rate_matrix[1], num_species);
    theta_1=update_theta(params, K_matrix, theta_0,num_species);
    I_f=0;
    I_f=Faradaic_current(E_vector, theta_1, rate_matrix[0][num_species-1]);
    Cdlp=Capacitance(params, I_0,  E);
    I_1=update_total_current(params, I_0, E, dE, Cdlp, I_f);
    I_0=I_1;
    theta_0=theta_1;
    }
  return py::cast(I_matrix);

}

PYBIND11_MODULE(seq_electron, m){
  	m.def("current_solver", &current_solver, "solve for total current for any number of sequential electron steps");

}


int main(){
};
/*double Cdl, double CdlE1, double CdlE2, double CdlE3,double omega,double v,double alpha ,double E_start, double E_reverse, \
	double delta_E, double Ru, double dt, std::vector<double> time_vec, double gamma, double E0_mean, double k0_mean, double phase
for(int i=0; i<num_species; i++){
    std::cout<<e_matrix[i]<<" "<< "\n";
    std::clock_t start;
    std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
    start= std::clock();
    double Cdl, double CdlE1, double CdlE2, double CdlE3,double omega,double v,double alpha ,double E_start, double E_reverse, \
    	double delta_E, double Ru, double dt, vector<double> time_vec, double gamma, double phase,int num_species, vector<double> k0_vector, vector<double> E0_vector
}*/
