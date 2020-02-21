#include "header.h"
#include "TemplateGeneratorNxM.h"

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <sstream>    
#include <cstddef>
#include <string>
#include <complex>
#include <vector>

#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"





TemplateGeneratorNxM::TemplateGeneratorNxM(std::vector<std::vector<std::vector<std::complex<double>>>> V, double E_MIN, double E_MAX, std::vector<double> vec_r_lim, std::vector<std::vector<double>> mat_theta_lim) : len_t_lim(mat_theta_lim.size()), len_r_lim(vec_r_lim.size()), num_channels(V.size()), num_bins_t(V[0][0].size()), E_min(E_MIN), E_max(E_MAX) {

    initialize(J, num_channels, num_bins_t);


    for (int a = 0; a < num_channels; a++){
        for (int n = 0; n < num_bins_t; n++){
            J[a][n] = V[a][a][n];               //auto-power  spectral densities (diagonals of cross-power spectral density matrix)
        }
    }


    std::vector<std::vector<std::complex<double>>> cumul_S_fft_0;
    initialize(cumul_S_fft_0, num_channels, num_bins_t);

    for (int a = 0; a < num_channels; a++){
        for (int n = 0; n < num_bins_t; n++){
            cumul_S_fft_0[a][n] = 0.0;        //a 2D vector of zeros that will be appended for each template, placeholder for all events included in a template's production
        }
    }


    //initializing cumul_S_fft and map_bins_part, described in the header file of this class
    int len_t_lim_i;
    for (int i = 0; i < len_t_lim; i++){
        len_t_lim_i = mat_theta_lim[i].size();
                
        if (i>0){

            std::vector<std::vector<std::vector<double>>> temp { {{vec_r_lim[i-1]}}, {{vec_r_lim[i]}}, {} };

            map_bins_part.push_back(temp);


            for (int j = 0; j< len_t_lim_i; j++){
                if(j==0){
                    map_bins_part[map_bins_part.size()-1][2].push_back( std::vector<double> {mat_theta_lim[i][len_t_lim_i-1], mat_theta_lim[i][j], (double)cumul_S_fft.size()} );
                }else{
                    map_bins_part[map_bins_part.size()-1][2].push_back( std::vector<double> {mat_theta_lim[i][j-1], mat_theta_lim[i][j], (double)cumul_S_fft.size()} );
                }
                cumul_S_fft.push_back(cumul_S_fft_0);
                event_counts.push_back(0);
            }

        }

        std::vector<std::vector<std::vector<double>>> temp { {{vec_r_lim[i]}}, {{vec_r_lim[i+1]}}, {} };
        map_bins_part.push_back(temp);

        for (int j = 0; j< len_t_lim_i; j++){
            if(j==0){
                map_bins_part[map_bins_part.size()-1][2].push_back( std::vector<double> {mat_theta_lim[i][len_t_lim_i-1], mat_theta_lim[i][j], (double)cumul_S_fft.size()} );
            }else{
                map_bins_part[map_bins_part.size()-1][2].push_back( std::vector<double> {mat_theta_lim[i][j-1], mat_theta_lim[i][j], (double)cumul_S_fft.size()} );
            }
            cumul_S_fft.push_back(cumul_S_fft_0);
            event_counts.push_back(0);
        }

        if (i<len_t_lim-1){

            std::vector<std::vector<std::vector<double>>> temp { {{vec_r_lim[i+1]}}, {{vec_r_lim[i+2]}}, {} };

            map_bins_part.push_back(temp);

            for (int j = 0; j< len_t_lim_i; j++){
                if(j==0){
                    map_bins_part[map_bins_part.size()-1][2].push_back( std::vector<double> {mat_theta_lim[i][len_t_lim_i-1], mat_theta_lim[i][j], (double)cumul_S_fft.size()} );
                }else{
                    map_bins_part[map_bins_part.size()-1][2].push_back( std::vector<double> {mat_theta_lim[i][j-1], mat_theta_lim[i][j], (double)cumul_S_fft.size()} );
                }
                cumul_S_fft.push_back(cumul_S_fft_0);
                event_counts.push_back(0);
            }

        }
    }

    vector_distributions = std::vector<vector_distribution>(cumul_S_fft.size()+1);      

    for (int i = 0; i<cumul_S_fft.size()+1; i++){
        vector_distribution Vect;
        vector_distributions[i] = Vect;
    }
}


//function to include an event in the creation of templates (we include all events in the "calibration sample")
bool TemplateGeneratorNxM::IncludeEvent(std::vector<std::vector<std::complex<double>>> S){
    std::vector<std::vector<std::complex<double>>> S_fft;
    initialize(S_fft, num_channels, num_bins_t);
    std::vector<double> amps(num_channels);

    std::vector<std::complex<double>> coefs_fft(num_bins_t);
    std::vector<std::complex<double>> IFFT(num_bins_t);

    for(int a = 0; a < num_channels; a++){

        amps[a] = 0.0;

        fft(S_fft[a], S[a]);

        for (int n = 0; n < num_bins_t; n++){  
            coefs_fft[n] = ( 1.0-( J[ a ][ n ]/( std::conj(S_fft[a][n])*(S_fft[a][n]) ) ) )*S_fft[ a ][ n ];
        } 
        ifft(IFFT, coefs_fft);

        for (int n = 0; n < num_bins_t; n++){
            amps[a] += std::real(IFFT[n]);
        }
    }


    //pre-construction of the quantities we want to reconstruct
    double E = 0.0;

    for (int a = 0; a < num_channels; a++){
        E += amps[a];
    }
    
    if (E<E_min || E>E_max){
        return false;
    }else{

        double r;
        double theta;
        calc_r(r, amps);
        calc_theta(theta, amps);
        get_angle_std(theta, theta);

        bool found = false;

        //searching for all bins (user-defined or otherwise) that the included event lies within
        for(int i = 0; i < map_bins_part.size(); i++){
            if (r > map_bins_part[i][0][0][0] && r < map_bins_part[i][1][0][0]){
                for(int j = 0; j < map_bins_part[i][2].size(); j++){
                    bool truf = false;
                    check_angle(truf, theta, map_bins_part[i][2][j][0], map_bins_part[i][2][j][1]);
                    if (truf == true){
                        int k = std::round(map_bins_part[i][2][j][2]);
                        for (int a = 0; a < num_channels; a++){
                            for (int n = 0; n < num_bins_t; n++){
                                cumul_S_fft[k][a][n] += S_fft[a][n];   //including event in average of traces (templates)
                            }
                        }
                        event_counts[k] += 1;
                        found = true;

                        vector_distributions[k].add(theta,r);           //discussed in header file of this class
                    }
                }
            }
        }

        vector_distributions[vector_distributions.size()-1].add(theta,r);
        return found;


    }

}


//get generated templates
void TemplateGeneratorNxM::GetTemplates(std::vector<std::vector<std::vector<std::complex<double>>>> &ans){

    std::vector<std::complex<double>> IFFT(num_bins_t);
    std::vector<std::vector<std::complex<double>>> TEMPLATE;
    initialize(TEMPLATE, num_channels, num_bins_t);

    for (int i = 0; i < cumul_S_fft.size(); i++){
        if (event_counts[i] > 0){
            for(int a = 0; a < num_channels; a++){

                ifft(IFFT, cumul_S_fft[i][a]); 

                for(int n = 0; n < num_bins_t; n++){
                    TEMPLATE[a][n] = {std::real(IFFT[n])/((double)event_counts[i]), 0.0};
                }
	        }
            ans.push_back(TEMPLATE);   //only save templates which have events in them   
            
        }
    }
}




//get map_bins_part that corresponds to the selected templates in the function above.
void TemplateGeneratorNxM::GetMapBinsPart(std::vector<std::vector<std::vector<std::vector<double>>>> &map_bins_p){

    int template_count = 0;

    for (int i = 0; i<map_bins_part.size(); i++){


        std::vector<std::vector<std::vector<double>>> temp { {{map_bins_part[i][0][0][0]}}, {{map_bins_part[i][1][0][0]}}, {} };
        map_bins_p.push_back(temp);

        for (int j = 0; j<map_bins_part[i][2].size(); j++){
            int k = std::round(map_bins_part[i][2][j][2]);
            if (event_counts[k] > 0){
                map_bins_p[map_bins_p.size()-1][2].push_back( std::vector<double> {map_bins_part[i][2][j][0], map_bins_part[i][2][j][1], (double)template_count} );
                template_count += 1;
            }else{
                map_bins_p[map_bins_p.size()-1][2].push_back( std::vector<double> {map_bins_part[i][2][j][0], map_bins_part[i][2][j][1], -1} );
            }
        }
    }
}



//Draw calibration set in (r,theta) space, and bins for all possible optimal filters. Also draw some templates.
void TemplateGeneratorNxM::Draw(std::string path){

    std::string filename;

    const double pi = 3.14159265358979323846;

    int vds = vector_distributions.size();
    int maps_size = map_bins_part.size();

    std::vector<TGraph*> graphs(vds);

    for (int i = 0; i<vds; i++){
        if (vector_distributions[i].get_size() > 0){
            int n = vector_distributions[i].get_size();
            double * array_x = new double[n];
            double * array_y = new double[n];
            vector_distributions[i].get_array_x(array_x, n);
            vector_distributions[i].get_array_y(array_y, n);
            graphs[i] = new TGraph( n, array_x, array_y );
            delete[] array_x;
            delete[] array_y;
        }else{
            graphs[i] = new TGraph();
        }
    }

    std::vector<double> limits_x {-pi, pi};
    std::vector<double> limits_y {map_bins_part[0][0][0][0], map_bins_part[map_bins_part.size()-1][1][0][0]};
	
    std::vector<TLine*> lines;

    double r_bot, r_top, theta_clockw, theta_anticw;
    int maps_size_i2;

    for (int i = 0; i<maps_size; i++){
        if (i%3 == 0){
            maps_size_i2 = map_bins_part[i][2].size();

            r_bot = map_bins_part[i][0][0][0];
            r_top = map_bins_part[i][1][0][0];

            for (int j = 0; j<maps_size_i2; j++){
                theta_clockw = map_bins_part[i][2][j][0];
                theta_anticw = map_bins_part[i][2][j][1];

                lines.push_back( new TLine(theta_clockw, r_bot, theta_clockw, r_top));
                lines.push_back( new TLine(theta_anticw, r_bot, theta_anticw, r_top));

                if(theta_anticw > theta_clockw){
                    lines.push_back( new TLine( theta_clockw, r_bot, theta_anticw, r_bot ));
                    lines.push_back( new TLine( theta_clockw, r_top, theta_anticw, r_top ));
                }else{
                    lines.push_back( new TLine( theta_clockw, r_bot, pi     , r_bot ));
                    lines.push_back( new TLine( -pi    , r_bot, theta_anticw, r_bot ));
                    lines.push_back( new TLine( theta_clockw, r_top, pi     , r_top ));
                    lines.push_back( new TLine( -pi    , r_top, theta_anticw, r_top ));
                }
            }
            
        }
    }

    for (int i = 0; i<lines.size(); i++){
        lines[i]->SetLineColor(15);
    }

    //num graphs = vds
    for (int i = 0; i<maps_size; i++){
        if(i%3 == 0){
            maps_size_i2 = map_bins_part[i][2].size();

            for (int j = 0; j<maps_size_i2; j++){
                std::vector<TGraph*> graphs_subset {graphs[vds-1], graphs[std::round(map_bins_part[i][2][j][2])]};
                std::vector<int> colors {3,2};

                if(maps_size_i2>1){
                    if(j==0){
                        graphs_subset.push_back(graphs[std::round(map_bins_part[i][2][maps_size_i2-1][2])]);
                        colors.push_back(4);
                    }else{
                        graphs_subset.push_back(graphs[std::round(map_bins_part[i][2][j-1][2])]);
                        colors.push_back(4);
                    }

                    if(maps_size_i2>2){
                        graphs_subset.push_back(graphs[std::round(map_bins_part[i][2][(j+1)%maps_size_i2][2])]);
                        colors.push_back(4);
                    }
                }

                if(i>0){
                    graphs_subset.push_back(graphs[std::round(map_bins_part[i-1][2][j][2])]);
                    colors.push_back(4);
                }

                if(i<(maps_size-1)){
                    graphs_subset.push_back(graphs[std::round(map_bins_part[i+1][2][j][2])]);
                    colors.push_back(4);
                }

		filename.clear();
                filename = path+"/TemplateGeneratorNxM"+"_map_"+std::to_string( ( i/3 )+1 )+'_'+std::to_string( j+1 )+".png" ;
                draw_graphs(graphs_subset, colors, 0.5, "#theta", "r", limits_x, limits_y, filename, lines); 
            }



        }
    }

    for (int ii = 0; ii<graphs.size(); ii++){
        graphs[ii]->Delete();
    }

    for (int ii = 0; ii<lines.size(); ii++){
        delete lines[ii];
    }

    for (int i = 0; i<vds; i++){
        vector_distributions[i].reset();
    }


    //histograms:

    std::vector<std::vector<std::vector<std::complex<double>>>> templates;
    GetTemplates(templates);

    std::vector<TH1F*> hists(num_channels);

    std::vector<int> A(num_channels);

    for(int a = 0; a < num_channels; a++){
        hists[a] = new TH1F(create_name_hist(), "", num_bins_t, 0.0, (float)num_bins_t);
        A[a] = 2+a;
    }

    std::vector<std::vector<std::vector<std::vector<double>>>> mbpg; //map_bins_part_get;
    GetMapBinsPart(mbpg);
    int mbpgs = mbpg.size();

    std::string filename1;
    std::string filename2;
    std::string filename3;

    int k;
    for (int i = 0; i<mbpgs; i++){
        if(i%3 == 0){
            maps_size_i2 = mbpg[i][2].size();
            for(int j = 0; j<maps_size_i2; j++){
                k = std::round(mbpg[i][2][j][2]);

                if (k > -1){
                    for (int a = 0; a < num_channels; a++){
                        for(int n = 0; n < num_bins_t; n++){
                            hists[a]->SetBinContent(n+1, std::real(templates[k][a][n]));
                        }
                    }

		            filename1.clear();
                    filename1 = path+"/TemplateGeneratorNxM"+"_pulse_"+std::to_string( ( i/3 )+1 )+"_"+std::to_string( j+1 )+".png";
                    draw_hists(hists, A, "Time_(bin)", "", filename1);

                    for (int ii = 0; ii<hists.size(); ii++){
                        hists[ii]->Reset();
                    }



                }

                if(i > 0){
                    k = std::round(mbpg[i-1][2][j][2]);
                    if(k > -1){
                        for (int a = 0; a < num_channels; a++){
                            for(int n = 0; n < num_bins_t; n++){
                                hists[a]->SetBinContent(n+1, std::real(templates[k][a][n]));
                            }
                        }

			            filename2.clear();
                        filename2 = path+"/TemplateGeneratorNxM"+"_pulse_"+std::to_string( ( i/3 )+1 )+"_"+std::to_string( j+1 )+"_bot.png";
                        draw_hists(hists, A, "Time_(bin)", "", filename2);

                        for (int ii = 0; ii<hists.size(); ii++){
                            hists[ii]->Reset();
                        }

                    }
                }


                if(i < (mbpgs-1)){
                    k = std::round(mbpg[i+1][2][j][2]);
                    if(k > -1){
                        for (int a = 0; a < num_channels; a++){
                            for(int n = 0; n < num_bins_t; n++){
                                hists[a]->SetBinContent(n+1, std::real(templates[k][a][n]));
                            }
                        }
			            filename3.clear();
                        filename3 = path+"/TemplateGeneratorNxM"+"_pulse_"+std::to_string( ( i/3 )+1 )+"_"+std::to_string( j+1 )+"_top.png";
                        draw_hists(hists, A, "Time_(bin)", "", filename3);

                        for (int ii = 0; ii<hists.size(); ii++){
                            hists[ii]->Reset();
                        }

                    }


                }
            }
        }
    }

    for (int ii = 0; ii < hists.size(); ii++){
        hists[ii]->Delete();
    }    

}


