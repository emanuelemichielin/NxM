
#include <complex>
#include <vector>

#include "TGraph.h"
#include "TH1F.h"
#include "THStack.h"
#include "TTree.h"
#include "TLine.h"

#ifndef HEADER_H  
#define HEADER_H


class vector_distribution;
class tree_manager;


class vector_distribution{
    public:
    std::vector<double> list_x;
    std::vector<double> list_y;

    void add(double x, double y);
    int get_size();
    void get_array_x(double* ptr, int n);
    void get_array_y(double* ptr, int n);
    void reset();

};


class tree_manager{
    private:
    TTree * tree;       //t0, chisq, E
    std::vector<double> array;

    public:
    tree_manager(const char* NAME);
    ~tree_manager();
    void set(double t0, double chisq, double E);
    void Fill();
    void Write();
};



int cind(int n, int size);
void define_lims(std::vector<double> &vec_r_lim, std::vector<std::vector<double>> &mat_theta_lim);
const char* create_name_hist();
void get_angle_std(double &ans, double angle);
void check_angle(bool &ans, double theta, double limit_theta_clockw, double limit_theta_anticw);
void calculate_delta_angle(double &ans, double theta_clockw, double theta_anticw);
void interpolate_parab(std::vector<double> &ans, double y1, double y2, double y3);
void fft(std::vector<std::complex<double>> &ans, std::vector<std::complex<double>> x);
void ifft(std::vector<std::complex<double>> &ans, std::vector<std::complex<double>> x);
void inv(std::vector<std::vector<std::complex<double>>> &ans, std::vector<std::vector<std::complex<double>>> MAT);
void inv(std::vector<std::vector<double>> &ans, std::vector<std::vector<double>> MAT);
int mindex(std::vector<double> A);
std::complex<double> sum(std::vector<std::complex<double>> A);
double sum(std::vector<double> A);
void calc_r(double &ans, std::vector<double> A);
void calc_theta(double &ans, std::vector<double> A);
void filter_weiner(std::vector<double> &ans, std::vector<std::complex<double>> S, std::vector<std::complex<double>> J);
bool save_noise_psds(std::vector<std::vector<std::vector<std::complex<double>>>> noise_psds, std::string filename);
void get_noise_psds(std::vector<std::vector<std::vector<std::complex<double>>>> &ans, std::string filename);
bool save_templates_nxm(std::vector<std::vector<std::vector<std::complex<double>>>> templates, double E_min, double E_max, std::vector<std::vector<std::vector<std::vector<double>>>> map_bins_part, std::string filename);
void get_templates_nxm(std::vector<std::vector<std::vector<std::complex<double>>>> &templates, double &E_min, double &E_max, std::vector<std::vector<std::vector<std::vector<double>>>> &map_bins_part, std::string filename);
void configure_draw();
void draw_hists(std::vector<TH1F*> hists, std::vector<int> colors, const char* title_x, const char* title_y,  std::string filename, bool log_y = false);
void draw_hist_2d(TH1F* hist,  const char* title_x,  const char* title_y,  std::string filename);
void draw_graphs(std::vector<TGraph*> graphs, std::vector<int> colors, float dotsize,  const char* title_x,  const char* title_y, std::vector<double> limits_x, std::vector<double> limits_y,  std::string filename , std::vector<TLine*> lines = {});
void Driver();
void initialize(std::vector<std::vector<std::vector<std::complex<double>>>> &vect, int a, int b, int c);
void initialize(std::vector<std::vector<std::vector<double>>> &vect, int a, int b, int c);
void initialize(std::vector<std::vector<std::complex<double>>> &vect, int a, int b);
void initialize(std::vector<std::vector<double>> &vect, int a, int b);

#endif
