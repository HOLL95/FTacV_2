
#ifndef HEADER_H_INPUT
#define HEADER_H_INPUT

struct time_param
{
    double v;
    double amp;
    double tr;
    double E_start;
    double E_reverse;
    double omega;
    double phase;
    double delta_E;
    double duration;
    double time_end;
    std::vector<double> time_vec;

};
class classical{
  public:
    double et(time_param &time_init, float t);
    double dEdt(time_param &time_init, float t);
};
class noramp{
  public:
    double et(time_param &time_init, float t);
    double dEdt(time_param &time_init, float t);
};


#endif
