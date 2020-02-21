//This class generates all templates possibly needed in the optimal filtering process according to inputted (r,theta) bins, given events from a sample called the "calibration sample"
//The templates in a given bin are an average of the events in that bin from the calibration sample



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

    std::vector<std::vector<std::vector<std::vector<double>>>> map_bins_part;  // this multi-dimensional vector stores bin boundary information, along with the corresponding template index.
                                                                               // elements adjacent in the first index to bins defined by the user (central bins) represent bins that are above or below...
                                                                               //...the central bin in the r parameter, which may have different theta boundaries than the user defined.
                                                                               //This means that for first index i, only those satisfying i%3==0 represent bins defined by the user.



    std::vector<std::vector<std::vector<std::complex<double>>>> cumul_S_fft;  // for each bin represented in map_bins_part, this stores the sum of the FFT of the signals of events that lie within the limits of the bin
    std::vector<int> event_counts;

    std::vector<vector_distribution> vector_distributions;  //each vector_distribution object contains a list of (r,theta) coordinates of events that are within the boundaries of a bin.
                                                            //this vector contains such lists for each bin for which a template was generated, whether or not the bin was defined by user.
                                                            //used for plotting purposes (color coding of events).



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