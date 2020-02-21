#include <complex>
#include <vector>
    

#ifndef NOISEPSDGENERATOR_H  
#define NOISEPSDGENERATOR_H


class NoisePSDGenerator{

    public:

    int num_channels;
    int num_bins_t;

    std::vector<std::vector<std::vector<std::complex<double>>>> cumul_V;

    int event_count;

    NoisePSDGenerator(int NUM_CHANNELS, int NUM_BINS_T);
    void IncludeEvent(std::vector<std::vector<std::complex<double>>> S);
    void CalculateNoisePSDs(std::vector<std::vector<std::vector<std::complex<double>>>> &ans);
    void Draw(std::string path);



    private:
    std::vector<std::vector<std::vector<std::complex<double>>>> noise_psds; 


};

#endif
