#include "Bridge.h"

// -------------------------------------------------- \\

// External global variables:

// -------------------------------------------------- \\

// Constructor:
particle::particle(void)
{
	int i;
	for(i=0;i<3;++i)
	{
		r[i] = 0.0;
		p[i] = 0.0;
	}
	mass = 0.0;
}

// Destructor:
particle::~particle(void)
{

}

// Functions:

void particle::initialize(double rv[3], double pv[3], double massv)
{
    int i;
    for(i=0;i<3;++i)
	{
		r[i] = rv[i];
		p[i] = pv[i];
	}
	mass = massv;
    return;
}
