
#ifndef INVERSEPROBLEMSOLVER_H
#define INVERSEPROBLEMSOLVER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <ctime>
#include <stdlib.h>
#include <complex>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iomanip>
#include <stdio.h>
#include <omp.h>

#include <nlopt.hpp>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_statistics.h>

#define NUMBER_OF_PARAMETERS 5
#define DELTAT 8
#define TRUNC 5
#define POLY_DEGREE 6
#define NCOEFFS 8
#define NBREAK (NCOEFFS-2)

using namespace std;

typedef std::vector<double> vec;
typedef std::vector<complex<double> > vecc;
typedef std::vector<vector<double> > vvec;

struct Ne_pars
{
	double Nmf2;
	double Hmf2;
	double H_b;
	double H_t;
	double Phi;
};

struct Plasma_pars
{
	double lambda;
	double Te;
	double Tr;
};

struct Data_T
{
    std::vector<vecc> CorMatrix;
    vec far;
	vec height;
    vec knots;
	std::vector<bool> isempty;
    std::map<Plasma_pars, vec> Ssec;
    size_t Clutter_nlags;
    double Pulse_Length;
    unsigned nE;
    unsigned nR;
	bool Global_solver_on;
};

typedef struct {
    size_t i;
	size_t j;
    size_t nE;
    size_t nR;
} constraint_data;

struct Optim_option
{
	double minf;
	unsigned neval;
	double time;
	bool opt_failed;
	unsigned maxeval;
	double stopval;
	double ftol_rel;
	double ftol_abs;
	double xtol_rel;
};

struct OF_data
{
	ofstream * out;
	vec data;
	vec height;
	std::vector<bool> isempty;
	bool Global_solver_on;
	unsigned channel;
	unsigned frequency;
	unsigned pulselength;
};

struct polynom_params
{
	double *A;
	int T;
	int n;
	int trunc;
};

struct data
{
	size_t n;
	size_t trunc;
	double *y;
};

struct DUH_data
{
	unsigned DoY;
	unsigned Hour;
	unsigned Minute;
	double Height;
	double Ne;
	double Te;
	double Tr;
	double P0;
	double P_exp;
	double StoN;
	vec ReCM;
	vec dReCM;

	bool operator==(const DUH_data& a) const
    {
        return (DoY==a.DoY and Hour==a.Hour and Minute==a.Minute and Height==a.Height);
    }
};

double Pulse(int t, int T);
double Ne_profile(double x, void * params);
vec Faraday(vec height, std::vector<bool> isempty, void * params);
vec Power_v(vec Ne);
vec Radeq(vec Far, vec height, std::vector<bool> isempty, double Length);
double Misfit(vec f, vec Y);
double ObjFunc(const vec &x, vec &grad, void *f_data);
double GenObjFunc(const vec &x, vec &grad, void *f_data);
void Global_solver(OF_data Data, double UT, \
	struct Optim_option &option, vec &x, std::string &mesout);
void Local_solver(OF_data Data, double UT, \
	struct Optim_option &option, vec &x, std::string &mesout);

std::map<Plasma_pars, vec> Scatter_section();
vec T_profile(vec coefs, vec knots, vec height);
std::vector<vec> Radeq_ACF(std::vector<vec> B0, vec Far, vec height, size_t Length);
void Get_Tprofiles(vec x, Data_T D, vec &Te_p, vec &Tr_p,\
     std::vector<vec> &B, vec &A);
double Ti_constraint(const vec &x, vec &grad, void *data);
double Tr_constraint(const vec &x, vec &grad, void *data);
double Te_constraint(const vec &x, vec &grad, void *data);
double ObjFuncT(const vec &x, vec &grad, void *f_data);
void Global_solver_T(Data_T data, double UT, \
	Optim_option &option, vec &x, std::string &mesout);
void Local_solver_T(Data_T data, double UT, \
	Optim_option &option, vec &x, std::string &mesout);

void ACF_pardet(double &t_0, double &t_min, double &a_min, double T, int trunc, vec p);
void ACFEST_pardet(double &t_0, double &t_min, double &a_min, double T, int trunc, vec p);
double polynom(double x, void *p);
double polyval(double x, vec A);
int repoly_f(const gsl_vector *x, void *data, gsl_vector *f);
int repoly_df(const gsl_vector *x, void *data, gsl_matrix *J);
// void ReACF_interpolator(double *y, size_t n, int trunc, vec &Ar, vec &Er, double &chi);
void ReACF_interpolator(vec y, const unsigned trunc, vec &Ar, vec &Er, double &chi);
void Get_Temperatures(double &Te, double &Tr, double t_0, double t_min, double A_min);

std::vector<string> LineSplit(std::string line, std::string sep);
double F_spline_fit(vec y, double x);
void B_spline_fit(vec &Cv, vec &Yerr, vec Xv, vec Yv, vec Wv);

void SkyNoise(double Freq, int DoY, map<int, double> &N_sky);


#endif
