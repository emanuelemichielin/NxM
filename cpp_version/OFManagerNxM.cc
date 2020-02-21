#include "header.h"
#include "OFManagerNxM.h"

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <sstream>    
#include <cstddef>
#include <complex>
#include <vector>


OFManagerNxM::OFManagerNxM(double DT, double T_PRE, std::vector<std::vector<std::vector<std::complex<double>>>> TEMPLATES, std::vector<std::vector<std::vector<std::complex<double>>>> V, double E_MIN, double E_MAX, std::vector<std::vector<std::vector<std::vector<double>>>> MAP_BINS_PART) : num_templates(TEMPLATES.size()), num_channels(TEMPLATES[0].size()), num_bins_t(TEMPLATES[0][0].size()), E_min(E_MIN), E_max(E_MAX), last_of_index(-1) {

    initialize(J, num_channels, num_bins_t);


    for (int a = 0; a < num_channels; a++){
        for (int n = 0; n < num_bins_t; n++){
            J[a][n] = V[a][a][n];
        }
    }


    int j_clockw, j_anticw, i_bot, i_top;
    std::vector<std::vector<std::vector<std::complex<double>>>> UU;

    for (int i = 0; i<MAP_BINS_PART.size(); i++){
        if(i%3 == 0){
            std::vector<std::vector<std::vector<double>>> temp { {{MAP_BINS_PART[i][0][0][0]}}, {{MAP_BINS_PART[i][1][0][0]}}, {} };
            map_bins_part.push_back(temp);

            for(int j = 0; j<MAP_BINS_PART[i][2].size(); j++){
                if(std::round(MAP_BINS_PART[i][2][j][2]) > -1){
                    map_bins_part[map_bins_part.size()-1][2].push_back( std::vector<double> {MAP_BINS_PART[i][2][j][0], MAP_BINS_PART[i][2][j][1], (double)list_of.size()} );

                    std::vector<int> list_adjacent_bins_indices;

                    if(MAP_BINS_PART[i][2].size()>1){
                        j_clockw = j - 1;
                        if (j==0){
                            j_clockw = MAP_BINS_PART[i][2].size() - 1;
                        }

                        if(std::round(MAP_BINS_PART[i][2][j_clockw][2])>-1){
                            list_adjacent_bins_indices.push_back(std::round(MAP_BINS_PART[i][2][j_clockw][2]));
                            
                        }

                        if(MAP_BINS_PART[i][2].size()>2){
                            j_anticw = (j+1)%MAP_BINS_PART[i][2].size();

                            if(std::round(MAP_BINS_PART[i][2][j_anticw][2]) > -1){
                                list_adjacent_bins_indices.push_back(std::round(MAP_BINS_PART[i][2][j_anticw][2]));
                            }
                        }

                    }

                    if(i>0){
                        i_bot = i-1;

                        if(std::round(MAP_BINS_PART[i_bot][2][j][2]) > -1){
                            list_adjacent_bins_indices.push_back(std::round(MAP_BINS_PART[i_bot][2][j][2]));
                        }

                    }

                    if(i<MAP_BINS_PART.size()-1){
                        i_top = i + 1;

                        if (std::round(MAP_BINS_PART[i_top][2][j][2]) > -1){
                            list_adjacent_bins_indices.push_back(std::round(MAP_BINS_PART[i_top][2][j][2]));
                        }
                    }

                    UU.clear();
                    UU.push_back(TEMPLATES[std::round(MAP_BINS_PART[i][2][j][2])]);
                    for (int ii = 0; ii<list_adjacent_bins_indices.size(); ii++){
                        UU.push_back(TEMPLATES[list_adjacent_bins_indices[ii]]);
                    }

                    OptimalFilterNxM optemp(DT, T_PRE, UU, V);
                    list_of.push_back(optemp);


                }else{
                    map_bins_part[map_bins_part.size()-1][2].push_back( std::vector<double> {MAP_BINS_PART[i][2][j][0], MAP_BINS_PART[i][2][j][1], -1} );
                }


            }
        }
    }

}






// 2 is false, 0 is none and 1 is true
int OFManagerNxM::ProcessEvent(std::vector<std::vector<std::complex<double>>> S, std::vector<std::vector<double>> &result){

    std::vector<double> amps(num_channels);

    std::vector<std::complex<double>> coefs_fft(num_bins_t);
    std::vector<std::complex<double>> IFFT(num_bins_t);

    for (int a = 0; a < num_channels; a++){
        fft(coefs_fft, S[a]);
        for (int n = 0; n < num_bins_t; n++){
            coefs_fft[n] *= (1.0-(J[a][n]/(std::conj(coefs_fft[n])*coefs_fft[n])));
        }
        ifft(IFFT, coefs_fft);
        amps[a] = 0.0;
        for (int n = 0; n < num_bins_t; n++){
            amps[a] += std::real(IFFT[n]);
        }
        
    }
    double E = 0.0;
    for (int a = 0; a < num_channels; a++){
        E += amps[a];
    }

    if (E<E_min || E>E_max){
        return 2;    
    }
    
    bool truf = false;
    double r, theta;
    calc_r(r, amps);
    calc_theta(theta, amps);
    get_angle_std(theta,theta);



    for(int i = 0; i < map_bins_part.size(); i++){
        if(r>map_bins_part[i][0][0][0] && r<map_bins_part[i][1][0][0]){
            for(int  j = 0; j<map_bins_part[i][2].size(); j++){
                check_angle(truf, theta, map_bins_part[i][2][j][0], map_bins_part[i][2][j][1]);
                if(truf == true){
                    last_of_index = std::round(map_bins_part[i][2][j][2]);
                    if(last_of_index > -1){
                        list_of[last_of_index].Execute(S);
                        result = list_of[last_of_index].result;
                        return 1;

                    }else{
                        return 0;
                    }
                }
            }
        }
    }

    last_of_index = -1;
    return 0;

}




void OFManagerNxM::Draw(std::string path, int event_count){
  //if(last_of_index>0){
        list_of[last_of_index].Draw(path, event_count);
	//}
}
