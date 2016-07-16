#include "Bridge.h"

int main(int argc, char** argv)
{
    // Objects
    NHthermostat NHinst;
    // Given parameters
    int N;                              // Number of particles
    double rho;                         // Density, rho = N /
    double sfac;                        // Stretching factor for xy-plane initialization
    double h1;                          // Height (z-axis) dividing hot/normal region, h1 > 0
    double h2;                          // Height (z-axis) dividing normal/cold region, h2 > h1
    double Q1;                          // Mass of thermostat in hot region
    double Q2;                          // Mass of thermostat in cold region
    double Teq1;                        // Equilibrium temperature in hot region, Teq1 > Teq2
    double Teq2;                        // Equilibrium temperature in cold region, Teq2 < Teq1
    double T0;                          // Initial temperature
    double rc;                          // Cutoff distance
    double zeta0_1;                     // Initial value for zeta in hot region
    double zeta0_2;                     // Initial value for zeta in cold region
    double gmag;                        // Magnitude of gravitational pull
    double dt;                          // Time step size
    double tprog;                       // Length of time to propagate
    double tproj;                       // Length of time between each projection
    // Derived parameters
    double psize;                       // Particle box size, psize**3 = V
    double space;                       // Space between initial particles, space = psize / ceiling(L), where L**3 = N
    int Nt;                             // Number of time steps to run
    int Ntp;                            // Ntp = Nt + 1
    int Nproj;                          // Number of time steps between each projection
    // Output variables
    double *cur_t;                      // Series of time values used
    double *temp_hot_t;                 // Time series for instantaneous temperature for hot region
    double *temp_mid_t;                 // Time series for instantaneous temperature for intermediate region
    double *temp_cold_t;                // Time series for instantaneous temperature for cold region
    double *temp_sys_t;                 // Time series for instantaneous temperature for all regions
    double *zeta1_t;                    // Time series for zeta in hot region
    double *zeta2_t;                    // Time series for zeta in cold region
    int **nregion_t;                    // Time series for number of particles in each region
    string ftemp_out, fvars_out, fnreg_out, ftpos_out;
    ofstream fout_file;
    // Timer variables
    double c_start, c_end;
    double e_time;
    // I/O variables
    string finput;
    ifstream finput_file;
    // Misc variables:
    int i, j, projidx, tempvar;
    string sbuf;

    // Initialize given parameters to default values
    N = 500;
    rho = 0.80;
    sfac = 2.0;
    h1 = 0.0;
    h2 = 0.0;
    Q1 = 1.0;
    Q2 = 1.0;
    Teq1 = 4.0;
    Teq2 = 1.0;
    T0 = 0.50;
    rc = 4.0;
    zeta0_1 = 0.0;
    zeta0_2 = 0.0;
    gmag = 0.0;
    dt = 0.01;
    tprog = 5.0;
    tproj = 1.0;

    // Extract filename from command line
    if(argc == 1) cout << "No file specified. Using default parameter values." << endl;
    else
    {
        // Read in parameter values from given input file
        finput = argv[1];
        cout << "Reading parameter values from file " << finput << endl;
        finput_file.open(finput);
        if(finput_file.is_open())
        {
            finput_file >> N >> sbuf;
            finput_file >> rho >> sbuf;
            finput_file >> sfac >> sbuf;
            finput_file >> h1 >> sbuf;
            finput_file >> h2 >> sbuf;
            finput_file >> Q1 >> sbuf;
            finput_file >> Q2 >> sbuf;
            finput_file >> Teq1 >> sbuf;
            finput_file >> Teq2 >> sbuf;
            finput_file >> T0 >> sbuf;
            finput_file >> rc >> sbuf;
            finput_file >> zeta0_1 >> sbuf;
            finput_file >> zeta0_2 >> sbuf;
            finput_file >> gmag >> sbuf;
            finput_file >> dt >> sbuf;
            finput_file >> tprog >> sbuf;
            finput_file >> tproj >> sbuf;
            finput_file.close();
            cout << "Parameter values initialized." << endl;
            //cout << N << " " << rho << " " << sfac << " "  << h1 << " " << h2 << " " << Q1 << " " << Q2 << " " << Teq1 << " " << Teq2 << " " << T0 << " " << rc << " " << zeta0_1 << " " << zeta0_2 << " " << gmag << " " << dt << " " << tprog << " " << tproj << endl;
        } else cout << "Unable to open file " << finput << ". Using default parameter values instead." << endl;
    }
	omp_set_num_threads(omp_get_num_procs());

    // Initialize derived parameters
    psize = pow((N/rho),1.0/3.0);
    space = psize / ceil(pow(N,1.0/3.0));
    Nt = ceil(tprog/dt);
    Ntp = Nt + 1;
    Nproj = ceil(tproj/dt);

    // Initialize objects
    NHinst.initialize(N,h1,h2,Q1,Q2,Teq1,Teq2,rc,zeta0_1,zeta0_2,gmag);
    NHinst.set_rgrid_prand(space,sfac);
    NHinst.build_region();
    // Allocate output variables
    cur_t = new double[Ntp];
    temp_hot_t = new double[Ntp];
    temp_mid_t = new double[Ntp];
    temp_cold_t = new double[Ntp];
    temp_sys_t = new double[Ntp];
    zeta1_t = new double[Ntp];
    zeta2_t = new double[Ntp];
    nregion_t = new int*[Ntp];
    for(i=0;i<Ntp;++i) nregion_t[i] = new int[3];
    // Initialize output variables
    for(i=0;i<Ntp;++i)
    {
        cur_t[i] = i * dt;
        temp_hot_t[i] = 0.0;
        temp_mid_t[i] = 0.0;
        temp_cold_t[i] = 0.0;
        temp_sys_t[i] = 0.0;
        zeta1_t[i] = 0.0;
        zeta2_t[i] = 0.0;
        for(j=0;j<3;++j) nregion_t[i][j] = 0;
    }
    ftemp_out = "temp_t.dat";
    fvars_out = "vars_t.dat";
    fnreg_out = "nreg_t.dat";
    ftpos_out = "tpos_t.dat";
    projidx = 0;
    // Open ftpos_out as a binary file
    fout_file.open(ftpos_out,ios::out|ios::binary);
    // First line of binary file: Number of particles, number of time steps written to file (including t = 0), time between each time step
    fout_file.write((char*)&N,sizeof(int));
    tempvar = 1+floor(double(Nt)/Nproj);
    fout_file.write((char*)&tempvar,sizeof(int));
    fout_file.write((char*)&tproj,sizeof(int));
    fout_file.write((char*)&h1,sizeof(double));
    fout_file.write((char*)&h2,sizeof(double));
    for(i=0;i<N;++i) fout_file.write((char*)NHinst.particle_arr[i].r,3*sizeof(double));

    cout << "Running time-propagation of system..." << endl;
    c_start = omp_get_wtime();
    // Fix initial temperature
    NHinst.set_temp(T0);
    // Run time propagation
    NHinst.measure_temp(temp_hot_t[0],temp_mid_t[0],temp_cold_t[0],temp_sys_t[i]);
    zeta1_t[0] = zeta0_1;
    zeta2_t[0] = zeta0_2;
    for(i=1;i<Ntp;++i)
    {
        NHinst.tstep(dt);
        NHinst.build_region();
        NHinst.measure_temp(temp_hot_t[i],temp_mid_t[i],temp_cold_t[i],temp_sys_t[i]);
        zeta1_t[i] = NHinst.zeta1;
        zeta2_t[i] = NHinst.zeta2;
        for(j=0;j<3;++j) nregion_t[i][j] = NHinst.nregion[j];
        ++projidx;
        if(projidx == Nproj)
        {
            for(j=0;j<N;++j) fout_file.write((char*)NHinst.particle_arr[j].r,3*sizeof(double));
            projidx = 0;
        }
    }
    c_end = omp_get_wtime();
    e_time = c_end - c_start;
    cout << "Time-propagation of system completed." << endl;
    cout << "Time taken = " << e_time << " seconds." << endl;
    // Close ftpos_out
    fout_file.close();

    // Write output to file
    fout_file.open(ftemp_out);
    for(i=0;i<Ntp;++i) fout_file << cur_t[i] << " " << temp_hot_t[i] << " " << temp_mid_t[i]  << " " << temp_cold_t[i]  << " " << temp_sys_t[i] << endl;
    fout_file.close();
    fout_file.open(fvars_out);
    for(i=0;i<Ntp;++i) fout_file << " " << cur_t[i] << " " << zeta1_t[i] << " " << zeta2_t[i] << endl;
    fout_file.close();
    fout_file.open(fnreg_out);
    for(i=0;i<Ntp;++i) fout_file << cur_t[i] << " " << nregion_t[i][0] << " " << nregion_t[i][1] << " " << nregion_t[i][2] << endl;
    fout_file.close();
    cout << "Output of results to files completed." << endl;

    // Deallocate output variables
    delete[] cur_t;
    delete[] temp_hot_t;
    delete[] temp_mid_t;
    delete[] temp_cold_t;
    delete[] temp_sys_t;
    delete[] zeta1_t;
    delete[] zeta2_t;
    for(i=0;i<Ntp;++i) delete[] nregion_t[i];
    delete[] nregion_t;

	return 0;
}
