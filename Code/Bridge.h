// This .h file links the .cpp files with the .h files.

#pragma once
#include "Constants.h"
#include "Libraries.h"
using namespace std;

// -------------------------------------------------- \\

// Forward Declarations:
// Objects:
class particle;
class NHthermostat;
// Functions:

// -------------------------------------------------- \\

// Class Defintions:

class particle
{
public:
	// Constructor:
	particle(void);
	// Destructor:
	~particle(void);
	// Functions:
	void initialize(double[],double[],double);
    inline double measure_ppm(void) { return (p[0]*p[0]+p[1]*p[1]+p[2]*p[2])/mass; }
	// Data elements:
	double r[3];
	double p[3];
	double mass;
private:
};

class NHthermostat
{
public:
	// Constructor:
	NHthermostat(void);
	// Destructor:
	~NHthermostat(void);
	// Functions:
	void initialize(int,double,double,double,double,double,double,double,double,double,double);
	void set_rgrid_prand(double,double);
	void build_region(void);
	void measure_temp(double&,double&,double&,double&);
	void set_temp(double);
	void build_force_arr(void);
	void tstep(double);
	// Data elements:
    particle *particle_arr;
    double **force_arr;
    int *region;
    int nregion[3];
    double h1;
    double h2;
    double Q1;
    double Q2;
    double Teq1;
    double Teq2;
    double rc2;
    double zeta1;
    double zeta2;
    int n;
    int g1;
    int gm;
    int g2;
    int g;
    double gmag;
private:
};
