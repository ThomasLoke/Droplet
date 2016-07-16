#include "Bridge.h"

// -------------------------------------------------- \\

// External global variables:

// -------------------------------------------------- \\

// Constructor:
NHthermostat::NHthermostat(void)
{
    particle_arr = NULL;
    force_arr = NULL;
    region = NULL;
    h1 = 0.0;
    h2 = 0.0;
    Q1 = 0.0;
    Q2 = 0.0;
    Teq1 = 0.0;
    Teq2 = 0.0;
    rc2 = 10000000000;
    zeta1 = 0.0;
    zeta2 = 0.0;
    n = 0;
    g1 = 0;
    gm = 0;
    g2 = 0;
    g = 0;
    gmag = 0;
}

// Destructor:
NHthermostat::~NHthermostat(void)
{
    int i;
    delete[] particle_arr;
    for(i=0;i<n;++i) delete[] force_arr[i];
    delete[] force_arr;
    delete[] region;
}

// Functions:

void NHthermostat::initialize(int nv, double h1v, double h2v, double Q1v, double Q2v, double Teq1v, double Teq2v, double rcv, double zeta1v, double zeta2v, double gmagv)
{
    int i, j;
    particle_arr = new particle[nv];
    force_arr = new double*[nv];
    for(i=0;i<nv;++i) force_arr[i] = new double[3];
    region = new int[nv];
    for(i=0;i<nv;++i)
    {
        for(j=0;j<3;++j) force_arr[i][j] = 0.0;
        region[i] = 0;
    }
    h1 = h1v;
    h2 = h2v;
    Q1 = Q1v;
    Q2 = Q2v;
    Teq1 = Teq1v;
    Teq2 = Teq2v;
    rc2 = rcv*rcv;
    zeta1 = zeta1v;
    zeta2 = zeta2v;
    n = nv;
    g1 = 3*nv;
    gm = 0;
    g2 = 0;
    g = 3*nv;
    gmag = gmagv;
    return;
}

void NHthermostat::set_rgrid_prand(double space, double sfac)
{
    int i, j, k, idx, gridsize, gridsize_xy, gridsize_z;
    double r_start_xy, r_start_z;
    double r_temp[3] = {0.0,0.0,0.0};
    double p_temp[3] = {0.0,0.0,0.0};
    // Initialize random number stuff
    typedef uniform_int<> NumberDistribution;
    typedef mt19937 RandomNumberGenerator;
    typedef variate_generator<RandomNumberGenerator&,NumberDistribution> Generator;

    NumberDistribution distribution(0.0,1.0);
    RandomNumberGenerator generator;
    Generator numberGenerator(generator, distribution);
    generator.seed(time(0));
    // Initialize grid of particles with random momenta
    gridsize = ceil(pow(n,1.0/3.0));
    gridsize_xy = ceil(gridsize * sfac);
    gridsize_z = ceil(gridsize / (sfac*sfac));
    r_start_xy = -space*gridsize_xy/2.0;
    r_start_z = space*gridsize_z;
    idx = 0;
    for(i=0;i<=gridsize_xy-1;++i)
    {
        r_temp[0] = r_start_xy + (i*space);
        for(j=0;j<=gridsize_xy-1;++j)
        {
            r_temp[1] = r_start_xy + (j*space);
            for(k=0;k<=gridsize_z-1;++k)
            {
                r_temp[2] = r_start_z - (k*space);
                p_temp[0] = numberGenerator();
                p_temp[1] = numberGenerator();
                p_temp[2] = numberGenerator();
                /*
                    WARNING: Currently treats all particles as m = 1.0
                */
                particle_arr[idx].initialize(r_temp,p_temp,1.0);
                ++idx;
                if(idx >= n) return;
            }
        }
    }
    return;
}

void NHthermostat::build_region(void)
{
    int i;
    double temp;
    nregion[0] = 0;
    nregion[1] = 0;
    nregion[2] = 0;
    for(i=0;i<n;++i)
    {
        temp = particle_arr[i].r[2];
        if(temp < h1)
        {
            region[i] = 1;
            ++nregion[0];
        } else {
            if(temp < h2)
            {
                region[i] = 2;
                ++nregion[1];
            } else {
                region[i] = 3;
                ++nregion[2];
            }
        }
    }
    g1 = nregion[0] * 3;
    gm = nregion[1] * 3;
    g2 = nregion[2] * 3;
    return;
}

void NHthermostat::measure_temp(double& temp_hot, double& temp_mid, double& temp_cold, double& temp_sys)
{
    int i;
    double temp;
    // Assumes that region has been built prior to calling this
    temp_hot = 0.0;
    temp_mid = 0.0;
    temp_cold = 0.0;
    temp_sys = 0.0;
    for(i=0;i<n;++i)
    {
        temp = particle_arr[i].measure_ppm();
        switch (region[i])
        {
            case 1:
                temp_hot += temp;
                break;
            case 2:
                temp_mid += temp;
                break;
            case 3:
                temp_cold += temp;
                break;
            default:
                cout << "Classification error! Regions improperly built. Terminating program." << endl;
                exit(EXIT_FAILURE);
        }
    }
    temp_sys = (temp_hot + temp_mid + temp_cold) / (k_B * g);
    if(nregion[0] > 0) temp_hot /= k_B * g1;
    if(nregion[1] > 0) temp_mid /= k_B * gm;
    if(nregion[2] > 0) temp_hot /= k_B * g2;
    return;
}

void NHthermostat::set_temp(double T0)
{
    double sval, temp;
    int i;
    temp = 0.0;
    for(i=0;i<n;++i) temp += particle_arr[i].measure_ppm();
    temp /= k_B * g;
    sval = sqrt(T0/temp);
    for(i=0;i<n;++i)
    {
        particle_arr[i].p[0] *= sval;
        particle_arr[i].p[1] *= sval;
        particle_arr[i].p[2] *= sval;
    }
    return;
}

void NHthermostat::build_force_arr(void)
{
    double temp[3];
    double rd2, rd6, fac;
    int i, j, k;
    for(i=0;i<n;++i) for(j=0;j<3;++j) force_arr[i][j] = 0.0;
	#pragma omp parallel for private(i,j,k,rd2,rd6,temp,fac)
    for(i=0;i<n;++i)
    {
        for(j=0;j<i;++j)
        {
			//if(i==j) continue;
            for(k=0;k<3;++k) temp[k] = particle_arr[i].r[k] - particle_arr[j].r[k];
            rd2 = temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2];
            // Compute force if rc*rc > (ri-rj)*(ri-rj), i.e. cutoff distance is greater than distance between particles i and j
            if(rc2 > rd2)
            {
                rd6 = rd2 * rd2 * rd2;
                fac = 24.0 * ((2.0/rd6) - 1.0) / (rd6*rd2);
                for(k=0;k<3;++k) temp[k] *= fac;
				#pragma omp critical
				for(k=0;k<3;++k)
                {
                    force_arr[i][k] = force_arr[i][k] + temp[k];
                    force_arr[j][k] = force_arr[j][k] - temp[k];
                }
            }
        }
    }
    // Apply gravitational effects to z-component of force
    for(i=0;i<n;++i) force_arr[i][2] -= gmag;
    return;
}

void NHthermostat::tstep(double dt)
{
    double dth, expval1, expval2;
    double temp_hot, temp_mid, temp_cold, temp_sys;
    int i, j;
    dth = dt * 0.5;
    expval1 = exp(-zeta1*dth);
    expval2 = exp(-zeta2*dth);
    for(i=0;i<n;++i)
    {
        // Apply exp(D3*dt/2) to p for hot and cold regions
        if(region[i] == 1)
        {
            for(j=0;j<3;++j) particle_arr[i].p[j] *= expval1;
        } else if(region[i] == 3) {
            for(j=0;j<3;++j) particle_arr[i].p[j] *= expval2;
        }
        // Apply exp(D2*dt/2) to r for all regions
        for(j=0;j<3;++j) particle_arr[i].r[j] += particle_arr[i].p[j] * dth / particle_arr[i].mass;
    }
    // Apply exp(Dt*dt/2) to zeta1 and zeta2
    measure_temp(temp_hot,temp_mid,temp_cold,temp_sys);
    if(nregion[0] > 0) zeta1 += (dth*g1*k_B/Q1)*(temp_hot - Teq1);
    if(nregion[2] > 0) zeta2 += (dth*g2*k_B/Q2)*(temp_cold - Teq2);
    // Recalculate forces
    build_force_arr();
    for(i=0;i<n;++i)
    {
        // Apply exp(D1*dt) for all regions
        for(j=0;j<3;++j) particle_arr[i].p[j] += force_arr[i][j] * dt;
        // Apply exp(D2*dt/2) for all regions
        for(j=0;j<3;++j) particle_arr[i].r[j] += particle_arr[i].p[j] * dth / particle_arr[i].mass;
    }
    // Apply exp(Dt*dt/2) to zeta1 and zeta2
    measure_temp(temp_hot,temp_mid,temp_cold,temp_sys);
    if(nregion[0] > 0) zeta1 += (dth*g1*k_B/Q1)*(temp_hot - Teq1);
    if(nregion[2] > 0) zeta2 += (dth*g2*k_B/Q2)*(temp_cold - Teq2);
    expval1 = exp(-zeta1*dth);
    expval2 = exp(-zeta2*dth);
    for(i=0;i<n;++i)
    {
        // Apply exp(D3*dt/2) to p for hot and cold regions
        if(region[i] == 1)
        {
            for(j=0;j<3;++j) particle_arr[i].p[j] *= expval1;
        } else if(region[i] == 3) {
            for(j=0;j<3;++j) particle_arr[i].p[j] *= expval2;
        }
    }
    // Enforce barrier at z = 0 after applying time-propagation operator
    for(i=0;i<n;++i)
    {
        if(particle_arr[i].r[2] < 0.0)
        {
            // Reflect z-components of r and z
            particle_arr[i].r[2] *= -1.0;
            particle_arr[i].p[2] *= -1.0;
        }
    }
    return;
}
