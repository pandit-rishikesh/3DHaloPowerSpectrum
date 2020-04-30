//------------ Parameters and Constants ------------//

//const int N = 160;
//const int N = 40;

//namespace constants{
#define PI 3.141592653589793

const int L = 1500;		   	// Box side in "h^{-1} Mpc"
const int H_0 = 100; 			// Hubble constant at z=0 in "h km/s/Mpc"
const int N_grid = 256; 		// Number of grids
const float H = (float) L/N_grid;		// Grid size

double k_N = PI / H;
double k_F = 2.*PI/L;                             //Fundamental frequency
const int kbins = 64;
double k_min =  0;
double k_max = 0.53197635600790005;                               // Nyquist frequency
double binsize_k = 2*k_F;

//double mu_min = -1;
//double mu_max = 1.;
//const int mubins = 100 ;
//double binsize_mu = (mu_max-mu_min)/mubins;


const int Ndim = pow(N_grid, 3);
const int Nbins = kbins;
double Norm = pow(N_grid, 6)/pow (L,3);            // Normalise Volume-Density

const float omega_m=0.25; // matter density
const float omega_lambda = 1 - omega_m;
const float Z = 0.;
const float Hofz = H_0 * sqrt(omega_m*pow((1+Z), 3.)+omega_lambda);
const float gadget_normalize = sqrt(1/(1+Z));

