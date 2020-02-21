#include "header.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <stdio.h>
#include <sstream>    
#include <cstddef>
#include <cblas.h>
#include <string>
#include <iomanip>
#include <complex>
#include<fftw3.h>
#include <lapacke.h>
#include <vector>

#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THStack.h"
#include "TTree.h"
#include "TDirectory.h"


//gPad->GetUymin() disabled in draw_hists to fix ylim problem (too big)


int cind(int n, int size){        
    if (n<0){                 
        return size + n;      
     }else{       
        return n;    
     }      
}  


void define_lims(std::vector<double> &vec_r_lim, std::vector<std::vector<double>> &mat_theta_lim){
    const double pi = 3.14159265358979323846;

    vec_r_lim = std::vector<double> {0 , 10};
    mat_theta_lim = std::vector<std::vector<double>> {{}};

    for (int i = 0; i<mat_theta_lim.size(); i++){
        std::vector<double> sym_theta_lims;
	
        for (int j = 0; j<mat_theta_lim[i].size(); j++){
            sym_theta_lims.push_back( -( mat_theta_lim[ i ][ j ]+( 2.0*pi/3.0 ) ) );
            sym_theta_lims.push_back( mat_theta_lim[ i ][ j ]-( 2.0*pi/3.0 )      );
            sym_theta_lims.push_back( -mat_theta_lim[ i ][ j ]                         );
            sym_theta_lims.push_back( mat_theta_lim[ i ][ j ]+( 2.0*pi/3.0 )      );
            sym_theta_lims.push_back( -( mat_theta_lim[ i ][ j ]-( 2.0*pi/3.0 ) ) );
        }
	
        for (int ii = 0; ii<sym_theta_lims.size(); ii++){ 
            mat_theta_lim[i].push_back(sym_theta_lims[ii]);
        }
	
        mat_theta_lim[ i ].push_back( pi  );
	
        std::sort(mat_theta_lim[i].begin(), mat_theta_lim[i].end());

    }

}





void get_angle_std(double &ans, double angle){

    const double pi = 3.14159265358979323846;

    double angle_std = angle;
    while (angle_std>pi){
        angle_std -= 2.0 * pi;
    }
    while (angle_std<-pi){
        angle_std += 2.0 * pi;
    }

    ans = angle_std;
}


void check_angle(bool &ans, double theta, double limit_theta_clockw, double limit_theta_anticw){
    if (limit_theta_anticw > limit_theta_clockw){
        if (theta > limit_theta_clockw && theta <= limit_theta_anticw){
            ans = true;
        }else{
            ans = false;
        }
    }else{
        if (theta > limit_theta_clockw  ||	theta <= limit_theta_anticw){
            ans = true;
        }else{
            ans = false;
        }
    }

}



void calculate_delta_angle(double &ans, double theta_clockw, double theta_anticw){

    const double pi = 3.14159265358979323846;

    if (theta_anticw > theta_clockw){
        ans = theta_anticw-theta_clockw;
    }else{
        ans = (2.0*pi) - (theta_clockw-theta_anticw);
    }
}


void interpolate_parab(std::vector<double> &ans, double y1, double y2, double y3){

    ans[0] = (0.5*(y1+y3))-y2;
    ans[1] = 0.5*(y3-y1);
    ans[2] = y2;

}



void vector_distribution::add(double x, double y){
    list_x.push_back(x);
    list_y.push_back(y);
}
int vector_distribution::get_size(){
    return list_x.size();
}
void vector_distribution::get_array_x(double* ptr, int n){
    for (int i = 0; i<n; i++){
        ptr[i] = list_x[i];
    }
}
void vector_distribution::get_array_y(double* ptr, int n){
    for (int i = 0; i<n; i++){
        ptr[i] = list_y[i];
    }
}
void vector_distribution::reset(){
    list_x.clear();
    list_y.clear();
}



void fft(std::vector<std::complex<double>> &ans, std::vector<std::complex<double>> x){
    int n = x.size();
    fftw_plan p = fftw_plan_dft_1d(n, (fftw_complex*)x.data(), (fftw_complex*)ans.data(), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
}


void ifft(std::vector<std::complex<double>> &ans, std::vector<std::complex<double>> x){
    int n = x.size();
    fftw_plan p = fftw_plan_dft_1d(n, (fftw_complex*)x.data(), (fftw_complex*)ans.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    for (int i = 0; i<n; i++){
	ans[i] /= n;
    }
}


void inv(std::vector<std::vector<std::complex<double>>> &ans, std::vector<std::vector<std::complex<double>>> MAT){   //intel math kernel library

    int N = MAT.size();

    int *IPIV = new int[N];

    __complex__ double * arr = new __complex__ double [N*N];
    for (int i = 0; i<N; i++){
        for (int j = 0; j<N; j++){
	    int idx = i*N + j;
	    arr[idx] = std::real(MAT[i][j]) + _Complex_I * std::imag(MAT[i][j]);
        }
    }


    LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, arr, N, IPIV);
    LAPACKE_zgetri(LAPACK_ROW_MAJOR, N, arr, N, IPIV);

    for (int i = 0; i<N; i++){
        for (int j = 0; j<N; j++){
	    int idx = i*N + j;
            ans[i][j] = {creal(arr[idx]), cimag(arr[idx])};
        }
    }
    
    delete[] IPIV;
    delete[] arr;
}




void inv(std::vector<std::vector<double>> &ans, std::vector<std::vector<double>> MAT){

    int N = MAT.size();

    int *IPIV = new int[N];

    double * arr = new double[N*N];
    for (int i = 0; i<N; i++){
        for (int j = 0; j<N; j++){
            int idx = i*N + j;
            arr[idx] = MAT[i][j];
        }
    }

    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, arr, N, IPIV);
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, N, arr, N, IPIV);

     for (int i = 0; i<N; i++){
        for (int j = 0; j<N; j++){
	        int idx = i*N + j;
            ans[i][j] = arr[idx];
        }
    }
    
    delete[] IPIV;
    delete[] arr;
}






int mindex(std::vector<double> A){
    int ans = A[0];
    int ind = 0;

    for (int i = 1; i<A.size(); i++){
        if (A[i] < ans){
            ans = A[i];
            ind = i;
        }
    }
    return ind;
}



std::complex<double> sum(std::vector<std::complex<double>> A){
    std::complex<double> ans = {0.0, 0.0};

    for(int i = 0; i<A.size(); i++){
        ans += A[i];
    }

    return ans;
}


double sum(std::vector<double> A){
    double ans = 0.0;

    for(int i = 0; i<A.size(); i++){
        ans += A[i];
    }

    return ans;
}







// #########################################################       CALC       #########################################################


void calc_r(double &ans, std::vector<double> A){
  //double sumw = 0.5 * (A[0]);
    ans = (4.37*A[0])/(A[0]);
}

void calc_theta(double &ans, std::vector<double> A){

    double sumAr = A[1] + A[2] + A[3];
    double X = (1.0 * A[1] - 0.5 * A[2] - 0.5 * A[3])/sumAr;
    double Y = (0.5 * std::sqrt(3.0) * A[2] - 0.5 * std::sqrt(3.0) * A[3])/sumAr;

    ans = std::atan2(Y,X);
    ans = 0;

}



// #########################################################       WIENER       #########################################################



void filter_weiner(std::vector<double> &ans, std::vector<std::complex<double>> S, std::vector<std::complex<double>> J){

    int S_size = S.size();
    std::vector<std::complex<double>> S_fft(S_size);
    std::vector<std::complex<double>> S_ifft(S_size);
    fft(S_fft, S);

    for (int n = 0; n<S_size; n++){
        S_fft[n] *= (1.0 - (J[n]/(S_fft[n]*std::conj(S_fft[n]))));
    }

    ifft(S_ifft, S_fft);

    for (int i = 0; i<S_size; i++){
        ans[i] = std::real(S_ifft[i]);
    }

}



// #########################################################       SAVING AND RETRIEVING       #########################################################

bool save_noise_psds(std::vector<std::vector<std::vector<std::complex<double>>>> noise_psds, std::string filename){
    int num_channels = noise_psds.size();
    int num_bins_t = noise_psds[0][0].size();

    std::ofstream outFile;

    outFile.open(filename);
    outFile << num_channels << " " << num_bins_t << std::endl;

    for (int a = 0; a < num_channels; a++){
        for (int b = 0; b < num_channels; b++){
            for (int n = 0; n < num_bins_t; n++){
                outFile << std::scientific << std::setprecision(14) << std::real(noise_psds[a][b][n]) << " " << std::imag(noise_psds[a][b][n]) << " " << a << " " << b << " " << n << std::endl;
            }
        }
    }

    outFile.close();
    return true;

}



void get_noise_psds(std::vector<std::vector<std::vector<std::complex<double>>>> &ans, std::string filename){
    std::ifstream inFile;
    inFile.open(filename);

    int num_channels;
    int num_bins_t;

    inFile >> num_channels >> num_bins_t;

    double real_part;
    double imag_part;

    int a;
    int b;
    int n;

    while (inFile >> real_part >> imag_part >> a >> b >> n) {
        ans [a][b][n] = {real_part, imag_part};
    }

    inFile.close();

}



bool save_templates_nxm(std::vector<std::vector<std::vector<std::complex<double>>>> templates, double E_min, double E_max, std::vector<std::vector<std::vector<std::vector<double>>>> map_bins_part, std::string filename){

    int num_templates = templates.size();
    int num_channels = templates[0].size();
    int num_bins_t = templates[0][0].size();
    int num_bins_r = map_bins_part.size();
    int num_bins_theta;
    

    std::ofstream outFile;

    outFile.open(filename);
    outFile << std::scientific << std::setprecision(14) << num_templates << " " << num_channels << " " << num_bins_t << " " << E_min << " " << E_max << " " << num_bins_r << std::endl;

    for (int i = 0; i<num_templates; i++){
        for (int a = 0; a<num_channels; a++){
            for (int n = 0; n<num_bins_t; n++){
                outFile << std::scientific << std::setprecision(14) << std::real(templates[i][a][n]) << std::endl;
            }
        }
    }

    for (int i = 0; i<num_bins_r; i++){
        num_bins_theta = map_bins_part[i][2].size();
        outFile << std::scientific << std::setprecision(14) << map_bins_part[i][0][0][0] << " " << map_bins_part[i][1][0][0] << " " << num_bins_theta << std::endl;

        for (int j = 0; j<num_bins_theta; j++){
            outFile << std::scientific << std::setprecision(14) << map_bins_part[i][2][j][0] << " " << map_bins_part[i][2][j][1] << " " << std::round(map_bins_part[i][2][j][2]) << std::endl;
        }
    }

    outFile.close();
    return true;
}



void get_templates_nxm(std::vector<std::vector<std::vector<std::complex<double>>>> &templates, double &E_min, double &E_max, std::vector<std::vector<std::vector<std::vector<double>>>> &map_bins_part, std::string filename){
    std::ifstream inFile;
    inFile.open(filename);

    int num_channels;
    int num_bins_t;
    int num_templates;
    int num_bins_r;
    int num_bins_theta;

    inFile >> num_templates >> num_channels >> num_bins_t >> E_min >> E_max >> num_bins_r;

    std::vector<std::vector<std::complex<double>>> TEMPLATE;
    initialize(TEMPLATE, num_channels, num_bins_t);


    double real_part;
    for (int i = 0; i<num_templates; i++){
        for (int a = 0; a<num_channels; a++){
            for (int n = 0; n<num_bins_t; n++){
                inFile >> real_part;
                TEMPLATE[a][n] = {real_part, 0.0};
            }
        }
        templates.push_back(TEMPLATE);
    }

    double r_bot; 
    double r_top; 
    double theta_clockw; 
    double theta_anticw;
    int template_index;


    int map_bins_part_size;
    for (int i = 0; i<num_bins_r; i++){
        inFile >> r_bot >> r_top >> num_bins_theta;

        std::vector<std::vector<std::vector<double>>> temp { {{r_bot}}, {{r_top}}, {} };
        map_bins_part.push_back(temp);

        map_bins_part_size = map_bins_part.size();  // can discard, because it's i+1, also line 219 in TemplateGeneratorNxM.cc

        for (int j = 0; j<num_bins_theta; j++){
            inFile >> theta_clockw >> theta_anticw >> template_index;
            map_bins_part[map_bins_part_size-1][2].push_back( std::vector<double> {theta_clockw, theta_anticw, (double) template_index} );
        }
    }

    inFile.close();
}



// #########################################################    DRAW UTILS     #########################################################

int global_counter = 1;
const char* create_name_hist(){
    extern int global_counter;
    std::string name_s = "h"+std::to_string(global_counter);
    global_counter += 1;
    const char *name_c = name_s.c_str();
    return name_c;
}


void configure_draw(){
    gROOT->Reset();
    gROOT->SetStyle("Plain");

    gStyle->SetEndErrorSize( 0.);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFont(22, "");
    gStyle->SetTitleSize(0.06, "");
    gStyle->SetTitleX(0.1);
    gStyle->SetTitleW(0.8);
}


void draw_hists(std::vector<TH1F*> hists, std::vector<int> colors,  const char* title_x,  const char* title_y,  std::string filename_s, bool log_y){

    const char *filename = filename_s.c_str();

    if(hists.size() == 0 || colors.size() != hists.size()){
        std::cout << "ERROR in draw_hists" << std::endl;
    }else{
        for (int i = 0; i<hists.size(); i++){
            hists[i]->SetLineColor(colors[i]);
            hists[i]->SetMarkerColor(colors[i]);
        }

        auto c = new TCanvas("c","");

        if(log_y==true){
            gPad->SetLogy(1);
	    gPad->SetLogx(1);
        }

        if (hists.size() == 1){
            hists[ 0 ]->GetXaxis()->SetLabelFont( 22 );
            hists[ 0 ]->GetXaxis()->SetLabelOffset( 0.012 );
            hists[ 0 ]->GetXaxis()->SetLabelSize( 0.035 );
            hists[ 0 ]->GetXaxis()->SetTitle( title_x );
            hists[ 0 ]->GetXaxis()->SetTitleFont( 22 );
            hists[ 0 ]->GetXaxis()->SetTitleOffset( 1.34 );
            hists[ 0 ]->GetXaxis()->SetTitleSize( 0.04 );
    
            hists[ 0 ]->GetYaxis()->SetLabelFont( 22 );
            hists[ 0 ]->GetYaxis()->SetLabelOffset( 0.012 );
            hists[ 0 ]->GetYaxis()->SetLabelSize( 0.035 );
            hists[ 0 ]->GetYaxis()->SetTitle( title_y );
            hists[ 0 ]->GetYaxis()->SetTitleFont( 22 );
            hists[ 0 ]->GetYaxis()->SetTitleOffset( 1.4 );
            hists[ 0 ]->GetYaxis()->SetTitleSize( 0.04 );
    
            hists[ 0 ]->Draw();
            c->SaveAs( filename );

        }else{
            auto hstack = new THStack("hstack", "");

            for (int i = 0; i<hists.size(); i++){
                hstack->Add(hists[i]);
            }

            hstack->Draw("nostack");
            c->Update();

            hstack->GetXaxis()->SetLabelFont( 22 );
            hstack->GetXaxis()->SetLabelOffset( 0.012 );
            hstack->GetXaxis()->SetLabelSize( 0.035 );
            hstack->GetXaxis()->SetTitle( title_x );
            hstack->GetXaxis()->SetTitleFont( 22 );
            hstack->GetXaxis()->SetTitleOffset( 1.34 );
            hstack->GetXaxis()->SetTitleSize( 0.04 );
    
            hstack->GetYaxis()->SetLabelFont( 22 );
            hstack->GetYaxis()->SetLabelOffset( 0.012 );
            hstack->GetYaxis()->SetLabelSize( 0.035 );
            hstack->GetYaxis()->SetTitle( title_y );
            hstack->GetYaxis()->SetTitleFont( 22 );
            hstack->GetYaxis()->SetTitleOffset( 1.4 );
            hstack->GetYaxis()->SetTitleSize( 0.04 );
            
            //hstack->SetMinimum( gPad->GetUymin() );
            //hstack->SetMaximum( gPad->GetUymax() );
    
            hstack->Draw( "nostack" );
            c->SaveAs( filename );
            hstack->Delete();
            
        }

        if(log_y == true){
            gPad->SetLogy(0);
        }

        c->Close();


    }
}




void draw_hist_2d(TH1F* hist,  const char* title_x, const char* title_y,  std::string filename_s){

    const char *filename = filename_s.c_str();

    auto c = new TCanvas( "c", "" );
 
    hist->GetXaxis()->SetLabelFont( 22 );
    hist->GetXaxis()->SetLabelOffset( 0.012 );
    hist->GetXaxis()->SetLabelSize( 0.035 );
    hist->GetXaxis()->SetTitle( title_x );
    hist->GetXaxis()->SetTitleFont( 22 );
    hist->GetXaxis()->SetTitleOffset( 1.34 );
    hist->GetXaxis()->SetTitleSize( 0.04 );
 
    hist->GetYaxis()->SetLabelFont( 22 );
    hist->GetYaxis()->SetLabelOffset( 0.012 );
    hist->GetYaxis()->SetLabelSize( 0.035 );
    hist->GetYaxis()->SetTitle( title_y );
    hist->GetYaxis()->SetTitleFont( 22 );
    hist->GetYaxis()->SetTitleOffset( 1.4 );
    hist->GetYaxis()->SetTitleSize( 0.04 );
 
    hist->Draw( "colz" );
    c->SaveAs( filename );
 
    c->Close();
}



void draw_graphs(std::vector<TGraph*> graphs, std::vector<int> colors, float dotsize,  const char* title_x,  const char* title_y, std::vector<double> limits_x, std::vector<double> limits_y,  std::string filename_s , std::vector<TLine*> lines){

    const char *filename = filename_s.c_str();


    if ( graphs.size() == 0 or colors.size()==0 or limits_x.size() != 2 or limits_y.size() != 2){
        std::cout << "ERROR in draw_graphs" << std::endl;
    }else{
        for (int i = 0; i<graphs.size(); i++){
            graphs[i]->SetMarkerColor(colors[i]);
            graphs[i]->SetMarkerSize(dotsize);
            graphs[i]->SetMarkerStyle(8);
        }

        auto c = new TCanvas("c","");

        int ind = 0;

        while(true==true){
            if (graphs[ind]->GetN() > 0){
                graphs[ ind ]->GetXaxis()->SetLabelFont( 22 );
                graphs[ ind ]->GetXaxis()->SetLabelOffset( 0.012 );
                graphs[ ind ]->GetXaxis()->SetLabelSize( 0.035 );
                graphs[ ind ]->GetXaxis()->SetTitle( title_x );
                graphs[ ind ]->GetXaxis()->SetTitleFont( 22 );
                graphs[ ind ]->GetXaxis()->SetTitleOffset( 1.34 );
                graphs[ ind ]->GetXaxis()->SetTitleSize( 0.04 );
    
                graphs[ ind ]->GetYaxis()->SetLabelFont( 22 );
                graphs[ ind ]->GetYaxis()->SetLabelOffset( 0.012 );
                graphs[ ind ]->GetYaxis()->SetLabelSize( 0.035 );
                graphs[ ind ]->GetYaxis()->SetTitle( title_y );
                graphs[ ind ]->GetYaxis()->SetTitleFont( 22 );
                graphs[ ind ]->GetYaxis()->SetTitleOffset( 1.4 );
                graphs[ ind ]->GetYaxis()->SetTitleSize( 0.04 );
    
                graphs[ ind ]->SetTitle( "" );
    
                graphs[ ind ]->GetXaxis()->SetLimits( limits_x[ 0 ], limits_x[ 1 ] );
    
                graphs[ ind ]->SetMinimum( limits_y[ 0 ] );
                graphs[ ind ]->SetMaximum( limits_y[ 1 ] );
    
                graphs[ ind ]->Draw( "ap" );
                break;
            }
            ind += 1;
        }

        if(ind == graphs.size()){
            std::cout << "ERROR in draw_graphs: NO POINTS" << std::endl;
        }else{
            for (int ii = 0; ii<graphs.size(); ii++){
                if(graphs[ii]->GetN()>0){
                    graphs[ii]->Draw("p");      //python code separates to 2 loops, [:ind] and [ind+1:]. Why?
                }
            }

            for (int ii = 0; ii<lines.size(); ii++){
                lines[ii]->Draw();
            }

            c->SaveAs(filename);
            c->Close();
        }

    }


}


tree_manager::tree_manager(const char* NAME) : array(4){
    tree = new TTree(NAME, NAME);
    tree->Branch("t0", &array[0], "t0/D");
    tree->Branch("A", &array[1], "A/D");
    tree->Branch("chisq", &array[2], "chisq/D");
    tree->Branch("E", &array[3], "E/D");
}
tree_manager::~tree_manager(){
    delete tree;
}
void tree_manager::set(double t0, double A, double chisq, double E){
    array[0] = t0;
    array[1] = A;
    array[2] = chisq;
    array[3] = E;
}
void tree_manager::Fill(){
    tree->Fill();
}
void tree_manager::Write(){
    tree->Write();
}




void initialize(std::vector<std::vector<std::vector<std::complex<double>>>> &vect, int a, int b, int c){
    std::vector<std::vector<std::vector<std::complex<double>>>> vect_t(a);
    for (int i = 0; i<a; i++){
        vect_t[i] = std::vector<std::vector<std::complex<double>>> (b);
    }

    for (int i = 0; i<a; i++){
        for (int j = 0; j<b; j++){
            vect_t[i][j] = std::vector<std::complex<double>> (c, {0.0,0.0});
        }
    }
    vect = vect_t;

}

void initialize(std::vector<std::vector<std::vector<double>>> &vect, int a, int b, int c){
    std::vector<std::vector<std::vector<double>>> vect_t(a);
    for (int i = 0; i<a; i++){
        vect_t[i] = std::vector<std::vector<double>> (b);
    }

    for (int i = 0; i<a; i++){
        for (int j = 0; j<b; j++){
            vect_t[i][j] = std::vector<double> (c,0.0);
        }
    }
    vect = vect_t;

}


void initialize(std::vector<std::vector<std::complex<double>>> &vect, int a, int b){
    std::vector<std::vector<std::complex<double>>> vect_t(a, std::vector<std::complex<double>> (b, {0.0, 0.0}));
    vect = vect_t;

}


void initialize(std::vector<std::vector<double>> &vect, int a, int b){
    std::vector<std::vector<double>> vect_t(a, std::vector<double>(b, 0.0));
    vect = vect_t;

}
