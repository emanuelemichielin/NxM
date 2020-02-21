#include "header.h"
#include "NoisePSDGenerator.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <sstream>    
#include <cstddef>
#include <string>
#include <complex>
#include <vector>

#include "TH1F.h"



NoisePSDGenerator::NoisePSDGenerator(int NUM_CHANNELS, int NUM_BINS_T) : num_channels(NUM_CHANNELS), num_bins_t(NUM_BINS_T), event_count(0) {

     initialize(cumul_V, num_channels, num_channels, num_bins_t);
     initialize(noise_psds, num_channels, num_channels, num_bins_t);


    for (int a = 0; a < num_channels; a++){
        for (int b = 0; b < num_channels; b++){
            for (int n = 0; n < num_bins_t; n++){
                cumul_V[a][b][n] = {0.0, 0.0};
            }
        }
    }



}


void NoisePSDGenerator::IncludeEvent(std::vector<std::vector<std::complex<double>>> S){

    std::vector<std::vector<std::complex<double>>> S_fft;
    initialize(S_fft, num_channels, num_bins_t);

    for (int a = 0; a < num_channels; a++){
        fft(S_fft[a], S[a]);
    }
    
    for (int a = 0; a < num_channels; a++){
        for (int b = 0; b < num_channels; b++){
            for (int n = 1; n < num_bins_t; n++){
                cumul_V[a][b][n] += std::conj(S_fft[a][n]) * S_fft[b][n];
            }
        }
    }
    event_count += 1; 

}




void NoisePSDGenerator::CalculateNoisePSDs(std::vector<std::vector<std::vector<std::complex<double>>>> &ans){
    if (event_count == 0){
        std::cout << "No Events Loaded in PSD generator" << std::endl;
    }else{
        for (int a = 0; a < num_channels; a++){
            for (int b = 0; b < num_channels; b++){
                for (int n = 0; n < num_bins_t; n++){
                    noise_psds[a][b][n] = cumul_V[a][b][n]/(double)event_count ;
                    ans[a][b][n] = cumul_V[a][b][n]/(double)event_count ;
                    
                }
            }
        }
        
    }
}


void NoisePSDGenerator::Draw(std::string path){
    
    std::string filename;

    std::vector<std::vector<std::vector<std::complex<double>>>> noise_psds;
    initialize(noise_psds, num_channels, num_channels, num_bins_t);
    CalculateNoisePSDs(noise_psds);

    std::vector<TH1F*> hists(num_channels);

    std::vector<int> B(num_channels);

    
    for (int a = 0; a < num_channels; a++){

        hists[a] = new TH1F(create_name_hist(), "", (num_bins_t+2)/2, 0.0, (float)((num_bins_t+2)/2));
        
        B[a] = 2+a;
    }

    std::cout<<  std::abs(noise_psds[0][1][100]) <<std::endl;
    std::cout<<  std::abs(noise_psds[1][0][100]) <<std::endl;

    for (int a = 0; a < num_channels; a++){
        for (int b = 0; b < num_channels; b++){
            for (int n = 1; n < num_bins_t; n++){
                hists[b]->SetBinContent(n+1, std::abs(noise_psds[a][b][n]));
            }
        } 
    
        filename.clear();
	filename = path+"/NoisePSDGenerator"+"_spec_"+std::to_string(a+1)+".png";
        draw_hists(hists, B, "Frequency_(bin)", "", filename, true);

        for (int ii = 0; ii <hists.size(); ii++){
            hists[ii]->Reset();                 
        }
    }

    for (int ii=0; ii<hists.size(); ii++){
        hists[ii]->Delete();
    }

}

