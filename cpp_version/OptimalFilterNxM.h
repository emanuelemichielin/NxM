//this class is the optimal filter. An object of this class is given templates by the OFManager and it executes the OF on an event.

#include <complex>
#include <vector>



#ifndef OPTIMALFILTERNXM_H  
#define OPTIMALFILTERNXM_H

class OptimalFilterNxM{


    public:

    double dt;
    double t_pre;
    int n_trig;

    int num_templates;     //number of templates SPECIFIC to this optimal filter object. Contrast with num_templates in OFManagerNxM object.
    int num_channels;
    int num_bins_t;

    std::vector<double> N;
    std::vector<std::vector<std::vector<std::complex<double>>>> UU;
    std::vector<std::vector<std::vector<std::complex<double>>>> U_fft;
    std::vector<std::vector<std::vector<std::complex<double>>>> V_inv;
    std::vector<std::vector<std::vector<std::complex<double>>>> F;  //optimal filter (recall for single channel single template, it is conjugate of template / PSD)
    std::vector<std::vector<double>> P;              ///design matrix
    std::vector<std::vector<double>> P_inv;
    std::vector<std::vector<std::complex<double>>> S_fft;

    //A0, chisq0, E0, t0, A, chisq, E
    //0     1      2   3  4    5    6
    std::vector<std::vector<double>> result;

    

    OptimalFilterNxM(double DT, double T_PRE, std::vector<std::vector<std::vector<std::complex<double>>>> U, std::vector<std::vector<std::vector<std::complex<double>>>> V);
    void Execute(std::vector<std::vector<std::complex<double>>> S);
    void reset();
    void Draw(std::string path, int event_count);



    private:
    std::vector<std::vector<std::complex<double>>> V_inv_n;
    std::vector<std::vector<std::complex<double>>> V_n;


};



#endif

