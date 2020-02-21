#include <complex>
#include <vector>
#include "OptimalFilterNxM.h"
    

#ifndef OFMANAGERNXM_H  
#define OFMANAGERNXM_H  


class OFManagerNxM{
    
    public:

    int num_templates;
    int num_channels;
    int num_bins_t;
    
    std::vector<std::vector<std::complex<double>>> J;

    double E_min;
    double E_max;

    std::vector<std::vector<std::vector<std::vector<double>>>> map_bins_part;
    std::vector<OptimalFilterNxM> list_of;

    int last_of_index;

    OFManagerNxM(double DT, double T_PRE, std::vector<std::vector<std::vector<std::complex<double>>>> TEMPLATES, std::vector<std::vector<std::vector<std::complex<double>>>> V, double E_MIN, double E_MAX, std::vector<std::vector<std::vector<std::vector<double>>>> MAP_BINS_PART);
    int ProcessEvent(std::vector<std::vector<std::complex<double>>> S, std::vector<std::vector<double>> &result);
    void Draw(std::string path, int event_count);

};



#endif