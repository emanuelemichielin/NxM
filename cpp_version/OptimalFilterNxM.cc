#include "header.h"
#include "OptimalFilterNxM.h"

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <sstream>    
#include <cstddef>
#include <string>
#include <complex>
#include <vector>

#include "TH1F.h"


OptimalFilterNxM::OptimalFilterNxM(double DT, double T_PRE, std::vector<std::vector<std::vector<std::complex<double>>>> U, std::vector<std::vector<std::vector<std::complex<double>>>> V) : dt(DT), t_pre(T_PRE), n_trig(std::round(T_PRE/DT)), num_templates(U.size()), num_channels(U[0].size()), num_bins_t(U[0][0].size()), N(num_templates, 0.0) {


    result = {{}, {}, {}, {}, {}, {}, {}};

     initialize(UU, num_templates, num_channels, num_bins_t);
     initialize(U_fft, num_templates, num_channels, num_bins_t);
     initialize(V_inv, num_channels, num_channels, num_bins_t) ;
     initialize(F, num_templates, num_channels, num_bins_t) ;
     initialize(P, num_templates, num_templates) ;
     initialize(P_inv, num_templates, num_templates); 
     initialize(S_fft, num_channels,num_bins_t);
     initialize(V_inv_n, num_channels, num_channels);
     initialize(V_n, num_channels, num_channels);

     UU = U;

    for (int i = 0; i < num_templates; i++){
        for (int a = 0; a < num_channels; a++){
                N[i] += std::real(sum(U[i][a]));
            
        }
    }

    for (int i = 0; i < num_templates; i++){
        for (int a = 0; a < num_channels; a++){
            fft(U_fft[i][a], U[i][a]);
        }
    }
    
    for (int n = 0; n < num_bins_t; n++){

        for (int a = 0; a < num_channels; a++){
            for (int b = 0; b < num_channels; b++){
                V_n [a][b] = V[a][b][n];
            }
        }

        inv(V_inv_n, V_n);

        for (int a = 0; a < num_channels; a++){
            for (int b = 0; b < num_channels; b++){
                V_inv[a][b][n] = V_inv_n[a][b];
            }
        }
    }

    for (int i = 0; i < num_templates; i++){
        for (int a = 0; a < num_channels; a++){
            for (int b = 0; b < num_channels; b++){
                for (int n = 0; n < num_bins_t; n++){
            
                    F[i][a][n] += std::conj(U_fft[i][b][n])*V_inv[b][a][n];
                    
                }
            }
        }
    }

    for (int i = 0; i < num_templates; i++){
        for (int j = 0; j < num_templates; j++){

            for (int a = 0; a < num_channels; a++){
                for (int b = 0; b < num_channels; b++){
                    for (int n = 0; n < num_bins_t; n++){

                        P[i][j] += std::real(std::conj(U_fft[i][a][n]) * U_fft[j][b][n] * V_inv[a][b][n]);
                        
                        
                    }
                }
            }

        }
    }


    inv(P_inv, P);



}

void OptimalFilterNxM::reset(){
        result = {{}, {}, {}, {}, {}, {}, {}};
}



void OptimalFilterNxM::Execute(std::vector<std::vector<std::complex<double>>> S){

    reset();        //resetting result
    
    
    std::vector<std::complex<double>> IFFT(num_bins_t);
    

    for (int a = 0; a < num_channels; a++){
        fft(S_fft[a], S[a]);
    }


    std::vector<std::vector<double>> Q;
    initialize(Q, num_templates, num_bins_t);

    for (int i = 0; i < num_templates; i++){
        std::vector<std::complex<double>> coefs_fft(num_bins_t, 0.0);

        for (int a = 0; a < num_channels; a++){
            for (int n = 0; n < num_bins_t; n++){
                coefs_fft[n] += F[i][a][n] * S_fft[a][n];
                
            }
        }

        ifft(IFFT, coefs_fft);
        for (int n = 0; n < num_bins_t; n++){
            Q[i][n] = std::real(IFFT[n]);
        }

    }

    std::vector<std::vector<double>> funcs_a;
    initialize(funcs_a, num_templates, num_bins_t);

    for (int i = 0; i < num_templates; i++){

        for (int j = 0; j < num_templates; j++){
            for (int n = 0; n < num_bins_t; n++){
                
                funcs_a[i][n] += P_inv[i][j]*Q[j][n];
                
            }
        }
        for (int n = 0; n < num_bins_t; n++){
            funcs_a[i][n] *= (double) num_bins_t;
        }

    }


    for (int i = 0; i < num_templates; i++){
        result[0].push_back(funcs_a[i][0]);
    }

    double chisq_base = 0.0;


    for (int a = 0; a < num_channels; a++){
        for (int b = 0; b < num_channels; b++){
            for (int n = 0; n < num_bins_t; n++){
                chisq_base += std::real(std::conj(S_fft[a][n]) * S_fft[b][n] * V_inv[a][b][n]);
            }
        }
    }


    std::vector<double> func_chisq(num_bins_t);

    for (int n = 0; n < num_bins_t; n++){
        func_chisq[n] = chisq_base;
    }

    for (int i = 0; i < num_templates; i++){
        for (int j = 0; j < num_templates; j++){
            for (int n = 0; n < num_bins_t; n++){
                func_chisq[n] -= funcs_a[i][n] * funcs_a[j][n] * P[i][j];
            }
        }
    }

    result[1].push_back(func_chisq[0]);


    result[2].push_back(0.0);
    for (int i = 0; i < num_templates; i++){
        result[2][0] += (double)(result[0][i] * N[i]);
    }


    int n_min = mindex(func_chisq);

    if (n_min >= num_bins_t - n_trig){
        n_min -= num_bins_t;
    }
    
    if (n_min > -n_trig && n_min < ((num_bins_t-n_trig)-1)){
        int n_prev = n_min-1;
        int n_post = (n_min+1)%num_bins_t;

        std::vector<double> parab_chisq(3);
        interpolate_parab(parab_chisq, func_chisq[cind(n_prev,num_bins_t)], func_chisq[cind(n_min,num_bins_t)], func_chisq[cind(n_post,num_bins_t)]);
        double dn0 = -parab_chisq[1]/(2.0*parab_chisq[0]);

        result[3].push_back((double)((0.5+(double)n_min+dn0)*dt));                             //A0, chisq0, E0, t0, A, chisq, E
                                                                                              //0     1      2   3  4    5    6

        std::vector<double> parab_a(3);
        for (int i = 0; i < num_templates; i++){
            interpolate_parab(parab_a, funcs_a[i][cind(n_prev,num_bins_t)], funcs_a[i][cind(n_min,num_bins_t)], funcs_a[i][cind(n_post,num_bins_t)]);
            result[4].push_back((double)((parab_a[0]*(dn0*dn0))+(parab_a[1]*dn0)+parab_a[2]));
        }

        result[5].push_back((double)(parab_chisq[2]-((parab_chisq[1]*parab_chisq[1])/(4.0*parab_chisq[0]))));
        result[6].push_back(0.0);

        for (int i = 0; i < num_templates; i++){
            result[6][0] +=  (double)(result[4][i]*N[i]);
        }
            

    }else{
        result[3].push_back((double)((0.5+(double)(n_min))*dt));
        
        for (int i = 0; i < num_templates; i++){
            result[4].push_back((double)(funcs_a[i][cind(n_min,num_bins_t)]));
        }

        result[5].push_back( (double)(func_chisq[cind(n_min,num_bins_t)]) );
        result[6].push_back(0.0);

        for (int i = 0; i < num_templates; i++){
            result[6][0] +=  (double)(result[4][i]*N[i]);
        }
    }
    
}



void OptimalFilterNxM::Draw(std::string path, int event_count){
    std::vector<TH1F*> hists(2);


    hists[0] = new TH1F( create_name_hist(), "", num_bins_t, 0.0, (float)( num_bins_t ) );
    hists[1] = new TH1F( create_name_hist(), "", num_bins_t, 0.0, (float)( num_bins_t ) );

    std::vector<int> colors(2);
    colors[0] = 3;
    colors[1] = 2;

    std::vector<std::complex<double>> IFFT(num_bins_t);
    std::vector<std::vector<double>> S;
    initialize(S, num_channels, num_bins_t);

    
    for (int a = 0; a < num_channels; a++){
        ifft(IFFT, S_fft[a]);
        for (int n = 0; n < num_bins_t; n++){
            S[a][n] = std::real(IFFT[n]);
        }
    }

    std::vector<double> AU(num_templates);
    
    for (int a = 0; a < num_channels; a++){
        for (int n = 0; n < num_bins_t; n++){

            for (int i = 0; i < num_templates; i++){
                AU[i] = result[4][i]*std::real(UU[i][a][n]);
            }

            hists[0]->SetBinContent(n+1, (float)(S[a][n]));
            hists[1]->SetBinContent(n+1, sum(AU));

        }

        std::string filename = path+"/OptimalFilterNxM"+"_pulse_"+std::to_string( event_count )+"_"+std::to_string( a+1 )+".png";
        draw_hists( hists, colors, "Time_(bin)", "", filename );

        for (int ii = 0; ii<hists.size(); ii++){
            hists[ii]->Reset();
        }
    }

    for (int ii = 0; ii<hists.size(); ii++){
        hists[ii]->Delete();
    }
}
