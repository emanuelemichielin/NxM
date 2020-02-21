#include <complex>
#include "header.h"
#include <vector>

#ifndef TEMPLATEGENERATORNXM_H
#define TEMPLATEGENERATORNXM_H


class TemplateGeneratorNxM{
    

    public:

    int num_channels ;
    int num_bins_t ;
    double E_min ;
    double E_max ;

    std::vector<std::vector<std::complex<double>>> J;

    std::vector<std::vector<std::vector<std::vector<double>>>> map_bins_part;
    std::vector<std::vector<std::vector<std::complex<double>>>> cumul_S_fft;
    std::vector<int> event_counts;

    std::vector<vector_distribution> vector_distributions;



    TemplateGeneratorNxM(std::vector<std::vector<std::vector<std::complex<double>>>> V, double E_MIN, double E_MAX, std::vector<double> vec_r_lim, std::vector<std::vector<double>> mat_theta_lim);
    bool IncludeEvent(std::vector<std::vector<std::complex<double>>> S);
    void GetTemplates(std::vector<std::vector<std::vector<std::complex<double>>>> &ans);
    void GetMapBinsPart(std::vector<std::vector<std::vector<std::vector<double>>>> &map_bins_p);
    void Draw(std::string path);


    private:
    int len_t_lim ;
    int len_r_lim ;


};



#endif