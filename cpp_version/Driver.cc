#include "header.h"
#include "OFManagerNxM.h"
#include "TemplateGeneratorNxM.h"
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

#include "TFile.h"



void Driver(){

    configure_draw();

    std::string path = "/home/d/dcurtin/ialkh/scdmsPyTools/reco_nxm/reco_cpp_new";
    std::string data_noise = "noise.txt";
    std::string data_calib = "calibration.txt";
    std::string data_test = "test.txt";
    std::string save_psds = "psds.txt";
    std::string save_templates = "templates.txt";

    std::string filename_root_s = path +"/root/reco.root";
    const char * filename_root = filename_root_s.c_str();

    int NUM_EVENTS = 200;
    int NUM_CHANNELS = 6;
    int NUM_BINS_T = 1024;
    double DT = 0.0016;
    double T_PRE = 0.2655;
    double E_MIN = 0.8;
    double E_MAX = 1.2;

    int STEP_MONITOR = 10;
    
    std::vector<double> VEC_R_LIM;
    std::vector<std::vector<double>> MAT_THETA_LIM;
    define_lims(VEC_R_LIM, MAT_THETA_LIM);

    std::vector<std::vector<std::vector<std::complex<double>>>> NOISE_PSDS;
    std::vector<std::vector<std::vector<std::complex<double>>>> TEMPLATES;
    std::vector<std::vector<std::vector<std::vector<double>>>> MAP_BINS_PART;
    std::vector<std::vector<double>> RESULT;
    std::vector<std::vector<std::complex<double>>> temp_S;
    initialize(NOISE_PSDS, NUM_CHANNELS, NUM_CHANNELS, NUM_BINS_T);
    initialize(temp_S, NUM_CHANNELS, NUM_BINS_T);

    double sreal;


    auto gen = new NoisePSDGenerator(NUM_CHANNELS, NUM_BINS_T);
   
    std::ifstream inFile; 
    inFile.open(data_noise);
    for (int i = 0; i < NUM_EVENTS; i++){
        for (int a = 0; a < NUM_CHANNELS; a++){
            for (int n = 0; n < NUM_BINS_T; n++){
                inFile >> sreal;
                temp_S[a][n] = {sreal, 0.0};
            }
        }
        gen->IncludeEvent(temp_S);
    }
    inFile.close();
    
    gen->CalculateNoisePSDs(NOISE_PSDS);
    gen->Draw(path+"/png");
    save_noise_psds(NOISE_PSDS, save_psds);
    delete gen;

    auto gent = new TemplateGeneratorNxM(NOISE_PSDS, E_MIN, E_MAX, VEC_R_LIM, MAT_THETA_LIM);


    inFile.open(data_calib);
    for (int i = 0; i < NUM_EVENTS; i++){
        for (int a = 0; a < NUM_CHANNELS; a++){
            for (int n = 0; n < NUM_BINS_T; n++){
                inFile >> sreal;
                temp_S[a][n] = {sreal, 0.0};
            }
        }
        gent->IncludeEvent(temp_S);
    }
    inFile.close();

    gent->GetTemplates(TEMPLATES);
    gent->Draw(path+"/png");
    gent->GetMapBinsPart(MAP_BINS_PART);


    save_templates_nxm(TEMPLATES, E_MIN, E_MAX, MAP_BINS_PART, save_templates);
    delete gent;


    auto man = new OFManagerNxM(DT, T_PRE, TEMPLATES, NOISE_PSDS, E_MIN, E_MAX, MAP_BINS_PART);
    auto tm_nxm = new tree_manager("NxM");
    int stat;

    inFile.open(data_test);
    for (int i = 0; i < NUM_EVENTS; i++){
        for (int a = 0; a < NUM_CHANNELS; a++){
            for (int n = 0; n < NUM_BINS_T; n++){
                inFile >> sreal;
                temp_S[a][n] = {sreal, 0.0};
            }
        }
        stat = man->ProcessEvent(temp_S, RESULT);

        

        if (stat == 1){

        if((i+1)%STEP_MONITOR == 1){
            man->Draw(path+"/png", i+1);
        }
        
	    tm_nxm->set(RESULT[3][0], RESULT[5][0], RESULT[6][0]);
        }else{
	    tm_nxm->set(-999999.0, -999999.0, -999999.0);
        }

        tm_nxm->Fill();

        RESULT.clear();
    }
    inFile.close();

    auto filepointer = new TFile(filename_root, "recreate");
    tm_nxm->Write();
    filepointer->Close();
    delete tm_nxm;
    delete man;

    std::cout << std::endl << "|------------------------- END PROGRAM -------------------------|" << std::endl;
    std::exit(0);
}
