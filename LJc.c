#include <stdio.h>
#include <stdlib.h>
#include <math.h>
  
#define LOG_FILENAME "LJ.dat"  // output file for log data
#define XYZ_FILENAME "LJ.xyz"  // output file for structure data
#define LOG_STEP 10            // output step interval of log data
#define XYZ_STEP 10            // output step interval of xyz data
#define N 256                  // total number of atoms
#define RAW 4                  // 4*RAW^3 = N
#define BOXLX 100.0            // x-side length of unit cell [angstrom]
#define BOXLY 100.0            // y-side length of unit cell [angstrom]
#define BOXLZ 100.0            // z-side length of unit cell [angstrom]
#define SIGMA 3.405            // Lennard-Jones parameter [angstrom]
#define EPSILON 0.01032        // Lennard-Jones parameter [eV]
#define CUTOFF 3.0             // cutoff length [SIGMA]
#define MASS 39.95             // mass of Ar [amu]
#define T_INIT 300.0           // initial temperaturee [K]
#define TIME_END 50.0          // time length of simulation [ps]
#define TIME_STEP 0.01         // time step [ps]
#define RAND_SEED 10           // seed of random number generator
#define ps 1.0e-12             // 1 pico second [s]
#define angstrom 1.0e-10       // 1 angstrom [m]
#define eV 1.602176e-19        // 1 electron volt [J]
#define amu 1.660539e-27       // 1 atomic mass unit [kg]
#define kb 8.161734e-5         // Boltzmann constant [eV/K]
  
void calc_LJ_potential_force(double *fx, double *fy, double *fz, 
                             double *rx, double *ry, double *rz,
                             double *PE);
  
int main(int argc, char* argv[])
{
  double rx[N],  ry[N],  rz[N];  // coordinates 
  double vx[N],  vy[N],  vz[N];  // velocity
  double fx[N],  fy[N],  fz[N];  // force
  double f0x[N], f0y[N], f0z[N]; // force at the previous step
  double m;                      // mass of atoms
  double rij;                    // inter-atomic distance
  double dxij, dyij, dzij;       // relative position vector
  double t, dt;                  // time, time step
  double KE;                     // kinetic energy [eV]
  double PE;                     // potential energy [eV]
  double latticeconst;           // lattice constant at initial state
  double massunitconv;           // unit conversion factor for mass
  double T;                      // temperature [K]
  int counter;                   // MD simulation step counter
  int i,j,k,l;
  FILE *fp1, *fp2;
  
  fp1 = fopen(LOG_FILENAME,"w");
  if(fp1 == NULL){
    printf("ERROR: File open failed.\n");  exit(EXIT_FAILURE);
  }
  fp2 = fopen(XYZ_FILENAME,"w");                              
  if(fp2 == NULL){
    printf("ERROR: File open failed.\n");  exit(EXIT_FAILURE);
  }
  fprintf(fp1,"t KE PE KE+PE T\n");
  if(CUTOFF*SIGMA > 0.5*BOXLX || CUTOFF*SIGMA > 0.5*BOXLY || CUTOFF*SIGMA > 0.5*BOXLZ){
    printf("ERROR: Unitcell size is too small.\n");
    exit(EXIT_FAILURE);
  }
  
  // Initialize
  dt = TIME_STEP;
  massunitconv = amu * angstrom * angstrom / ps / ps / eV;
  m = MASS * massunitconv;
  srand(RAND_SEED);                                                     
  // Initial configuration (face-centered cubic)
  latticeconst = pow(2.0, 1.0/6.0) * SIGMA * sqrt(2.0);
  i = 0;                                                                       
  for(j=0; j<RAW; ++j){                                              
    for(k=0; k<RAW; ++k){                                        
      for(l=0; l<RAW; ++l){                                          
        rx[i]  = latticeconst * (j - 0.5 * RAW) + 0.5 * BOXLX;
        ry[i]  = latticeconst * (k - 0.5 * RAW) + 0.5 * BOXLY;
        rz[i]  = latticeconst * (l - 0.5 * RAW) + 0.5 * BOXLZ;
        ++i;
        rx[i]  = latticeconst * (0.5 + j - 0.5 * RAW) + 0.5 * BOXLX;
        ry[i]  = latticeconst * (0.5 + k - 0.5 * RAW) + 0.5 * BOXLY;
        rz[i]  = latticeconst * (l - 0.5 * RAW) + 0.5 * BOXLZ;
        ++i;
        rx[i]  = latticeconst * (0.5 + j - 0.5 * RAW) + 0.5 * BOXLX;
        ry[i]  = latticeconst * (k - 0.5 * RAW) + 0.5 * BOXLY;
        rz[i]  = latticeconst * (0.5 + l - 0.5 * RAW) + 0.5 * BOXLZ;
        ++i;
        rx[i]  = latticeconst * (j - 0.5 * RAW) + 0.5 * BOXLX;
        ry[i]  = latticeconst * (0.5 + k - 0.5 * RAW) + 0.5 * BOXLY;
        rz[i]  = latticeconst * (0.5 + l - 0.5 * RAW) + 0.5 * BOXLZ;
        ++i;
      }
    }
  }
  // set initial velocity by the Box-Muller method
  for(i=0; i<N; ++i){
    vx[i] = sqrt(-2.0*kb*T_INIT/m*log((double)rand()/RAND_MAX))
            * cos(2.0*M_PI*(double)rand()/RAND_MAX);
    vy[i] = sqrt(-2.0*kb*T_INIT/m*log((double)rand()/RAND_MAX))
            * cos(2.0*M_PI*(double)rand()/RAND_MAX);
    vz[i] = sqrt(-2.0*kb*T_INIT/m*log((double)rand()/RAND_MAX))
            * cos(2.0*M_PI*(double)rand()/RAND_MAX);
  }
  // initial force calculation
  calc_LJ_potential_force(fx, fy, fz, rx, ry, rz, &PE);
 
  counter = 0;
  for(t = 0.0; t < TIME_END; t += dt){
    // Kinetic energy calculation
    KE = 0.0;
    for(i=0; i<N; ++i){
      KE += 0.5 * m * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    }
    // Temperature calculation
    T = 2.0 / (3.0 * kb) * KE / (double)N;
    if(counter % LOG_STEP == 0){
      // print kinetic, potential and total energies
      fprintf(fp1,"%f %f %f %f %f\n", t, KE, PE, KE + PE, T);
    }
    if(counter % XYZ_STEP == 0){
      // print structure
      fprintf(fp2, "%d\n", N);
      fprintf(fp2, "Lennard-Jones potential\n");
      for(i=0; i<N; ++i){
        fprintf(fp2, "Ar, %f, %f, %f\n", rx[i], ry[i], rz[i]);
      }
    }
    //  velocity-Verlet algorithm
    for(i=0; i<N; ++i){
      rx[i] += vx[i] * dt + 0.5 * fx[i] / m * dt * dt;
      ry[i] += vy[i] * dt + 0.5 * fy[i] / m * dt * dt;
      rz[i] += vz[i] * dt + 0.5 * fz[i] / m * dt * dt;
      f0x[i] = fx[i];
      f0y[i] = fy[i];
      f0z[i] = fz[i];
    }
    calc_LJ_potential_force(fx, fy, fz, rx, ry, rz, &PE);
    for(i=0; i<N; ++i){
      vx[i] += 0.5 * (fx[i] + f0x[i]) / m * dt;
      vy[i] += 0.5 * (fy[i] + f0y[i]) / m * dt;
      vz[i] += 0.5 * (fz[i] + f0z[i]) / m * dt;
    }
    // Periodic boundary condition
    for(i=0; i<N; ++i){
      if(rx[i] > BOXLX) rx[i] -= BOXLX;
      if(rx[i] < 0.0)   rx[i] += BOXLX;
      if(ry[i] > BOXLY) ry[i] -= BOXLY;
      if(ry[i] < 0.0)   ry[i] += BOXLY;
      if(rz[i] > BOXLZ) rz[i] -= BOXLZ;
      if(rz[i] < 0.0)   rz[i] += BOXLZ;
    }
    counter++;
  }
  fclose(fp1);
  fclose(fp2);
  return EXIT_SUCCESS;
}

void calc_LJ_potential_force(double *fx, double *fy, double *fz, 
                             double *rx, double *ry, double *rz,
                             double *PE)
{
  int i,j;
  double rij, fij;         // distance, force 
  double dxij, dyij, dzij; // relative position vector
  
  (*PE) = 0.0;
  for(i=0; i<N; ++i){
    fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0;                                 
  }
  for(i=0; i<N; ++i){
    for(j = i+1; j<N; ++j){
      dxij = rx[i] - rx[j];
      dyij = ry[i] - ry[j];
      dzij = rz[i] - rz[j];
  
      // Minimum image convention
      if(dxij >  0.5 * BOXLX) dxij -= BOXLX;
      if(dxij < -0.5 * BOXLX) dxij += BOXLX;
      if(dyij >  0.5 * BOXLY) dyij -= BOXLY;
      if(dyij < -0.5 * BOXLY) dyij += BOXLY;
      if(dzij >  0.5 * BOXLZ) dzij -= BOXLZ;
      if(dzij < -0.5 * BOXLZ) dzij += BOXLZ;

      rij  = sqrt(dxij * dxij + dyij * dyij + dzij * dzij);
      if(rij < CUTOFF *SIGMA){
        (*PE) += 4.0 * EPSILON * ( pow(SIGMA/rij,12.0)-pow(SIGMA/rij,6.0));
        fij  = - 4.0 * EPSILON * ( -12.0 * pow(SIGMA/rij, 12.0) / rij
                                   + 6.0 * pow(SIGMA/rij, 6.0) / rij);
        fx[i] += fij * dxij / rij;
        fy[i] += fij * dyij / rij;
        fz[i] += fij * dzij / rij;
        fx[j] -= fij * dxij / rij;
        fy[j] -= fij * dyij / rij;
        fz[j] -= fij * dzij / rij;
      }
    }
  }
  return;
} 
