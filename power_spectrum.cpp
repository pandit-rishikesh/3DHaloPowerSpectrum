// This code takes ASCII input file with position coordinates, respective velocities and masses and gives out the power spectrum as its output.

#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<valarray>
#include<fstream>
#include<vector>
#include<cmath>
#include<cstring>
#include<string>
#include<iomanip>
#include<cassert>
#include<sstream>
#include<fftw3.h>
#include<complex>
#include"initialize_power.h"

using namespace std;

//===================================================================================================================
//                  Declaring universal variables and constants
//===================================================================================================================

struct Particle{double pos[3]; double mass;};
int N_part;
vector <Particle> particles;
double aux;
double mesh[N_grid][N_grid][N_grid];
double delta[N_grid][N_grid][N_grid];


enum AssignmentType {NGP, CIC, TSC};

//===================================================================================================================
//                  MAS: particle distribution function to form a density grid
//===================================================================================================================


void assignMass(AssignmentType type, int sim){


	string number;
	stringstream convert;
	convert << sim;
	number = convert.str();
	string spacetype= "RedshiftSpace";

	cout << "\33[32;40m" << "Measuring "+spacetype+" power spectrum P(k) for sim "+number << "\33[0m" << endl;

	int numLines = 0;

        ifstream in1(string("/home/rpandit/Dropbox/workspace/data/sim/z_0/1-40/run"+number+".dat").c_str());	
	assert(in1);
	string unused;
	while(getline(in1, unused))
	++numLines;
	in1.close();

	cout << "\33[1;31m" << numLines << "\33[36m" << " particles found in sim "+number << endl;

	N_part = numLines;
	ifstream in2(string("/home/rpandit/Dropbox/workspace/data/sim/z_0/1-40/run"+number+".dat").c_str());
	assert(in2);

	double vz[N_part];

	cout << "\33[1;30m" << "Acquiring positions from the box" << endl;

	for(int iPart=0; iPart < N_part; iPart++){
		Particle readParticles;
		in2 >> readParticles.pos[0] >> readParticles.pos[1] >> readParticles.pos[2] >> aux >> aux >> vz[iPart] >> aux;
		particles.push_back(readParticles);
		particles[iPart].mass = 1.0;
	}

	in2.close();

	//comment following loop for calculating realspace power spectrum

	for (int iPart = 0; iPart < N_part; iPart++){
			vz[iPart]*=(gadget_normalize/Hofz);
			particles[iPart].pos[2] = particles[iPart].pos[2] + vz[iPart];
	}


	cout << "Distributing particles on the density grid using "<< "\33[1;36m" <<"MAS"<< "\33[0;2m" << endl;

// window function for NGP CIC and TSC to make a regular cubic grid of size N_Grid x N_Grid x N_Grid

  for(int iPart =0; iPart<N_part; iPart++){
		int iiMesh = (int) floor(particles[iPart].pos[0]*(1./H));
		int jjMesh = (int) floor(particles[iPart].pos[1]*(1./H));
		int kkMesh = (int) floor(particles[iPart].pos[2]*(1./H));


		switch(type){
			case NGP:
				mesh[iiMesh][jjMesh][kkMesh]+=particles[iPart].mass;
				break;

			case CIC:
				for(int ii=-1; ii<=1; ii++){
					for(int jj=-1; jj<=1; jj++){
						for(int kk=-1; kk<=1; kk++){
								double delta, fraction = 1.0f;
								delta = (1.0 - abs(particles[iPart].pos[0]-(iiMesh+ii+0.5)*H) * (1./H));
								fraction *= delta > 0.0f ? delta : 0.0f;
								delta = (1.0 - abs(particles[iPart].pos[1]-(jjMesh+jj+0.5)*H) * (1./H));
								fraction *= delta > 0.0f ? delta : 0.0f;
								delta = (1.0 - abs(particles[iPart].pos[2]-(kkMesh+kk+0.5)*H) * (1./H));
								fraction *= delta > 0.0f ? delta : 0.0f;
								int iiGrid = iiMesh + ii;
								int jjGrid = jjMesh + jj;
								int kkGrid = kkMesh + kk;

								//PBC
								if (iiGrid < 0) iiGrid += N_grid;
								if (iiGrid >= N_grid) iiGrid -= N_grid;
								if(jjGrid < 0) jjGrid += N_grid;
								if(jjGrid >= N_grid) jjGrid -= N_grid;
								if(kkGrid < 0) kkGrid += N_grid;
								if(kkGrid >= N_grid) kkGrid -= N_grid;
								mesh[iiGrid][jjGrid][kkGrid] += particles[iPart].mass * fraction;
							}
						}
					}

				break;


			//case TSC:


				default: assignMass(CIC, sim);
				break;
			}

		}

	}

int main(int argc, char *argv[]){

	if(argc != 3){
	cout << "\33[2m" << "Please provide the interval of simulations e.g " << "\33[1;32m" << "./redpower " << "\33[3;31m" << " <from>  <to> " << "\33[0m"<< endl;
	return 0;
	}

	int from = atoi(argv[1]);
	int to = atoi(argv[2]);

	
	for(int sim=from; sim<=to; sim++){

		assignMass(CIC, sim);
		//assignMass(NGP, sim);
		
		string number;
		stringstream convert;
	       	convert << sim;
		number = convert.str();

//===================================================================================================================
//                  Density grid in real space ====== >>> Density grid in Fourier space
//===================================================================================================================

	double number_density=0;
	double mean_density = N_part/pow(L, 3);
	double poisson_shot_noise = 1./mean_density;
	//cout << "poisson shot noise is " << poisson_shot_noise << endl;

	cout << "computing the density contrast "<< "\u03B4" << "(r)" << endl;

	for (int i=0; i<N_grid; i++)
		for (int j=0; j<N_grid; j++)
			for (int k=0; k<N_grid; k++){
				mesh[i][j][k]*=pow((1./H), 3);
				//number_density+=mesh[i][j][k];
				delta[i][j][k] = (mesh[i][j][k]/mean_density)-1.0;
				//cout << delta[i][j][k] << endl;
			}


	cout << "\33[1;34m" << "FFTW :"<< "\33[36m" <<"\nplan    " << "\33[33m" << "##\n";


// Loading the default FFTW_FORWARD configuration and making the plan for obtaining delta_k

	fftw_complex *delta_in, *delta_out;

	delta_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ndim);
	delta_out =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ndim);

	fftw_plan plan_forward = fftw_plan_dft_3d(N_grid, N_grid, N_grid, delta_in, delta_out, FFTW_FORWARD, FFTW_ESTIMATE);


// Now loading the arrays of delta_in with the 1 dimensional grid density array 'grid_1d'

	int entry = 0;

	for (int i=0; i<N_grid; i++)
		for (int j=0; j<N_grid; j++)
			for (int k=0; k<N_grid; k++){
				delta_in[entry][0] = delta[i][j][k];
				delta_in[entry][1] = 0;
				entry++;
			}


	cout << "\33[36m" << "execute " << "\33[33m" << "##\n";

	fftw_execute(plan_forward);

	cout << "\33[36m" << "destroy " << "\33[33m" << "##" << "\33[0m" << endl;

	fftw_destroy_plan(plan_forward);


//===================================================================================================================
//                      Computation of Power Spectra in (k), (k_perp k_par), (k mu)
//===================================================================================================================


	cout << "\33[1;32m" << "Computing power spectrum " << "\33[31m;" << " P(k)" << "\33[0m" << endl;

	int p = 2;

	int k1, k2, k3, K, K1, K2;
	double k_mag, k_perp, k_par, mu_k, Legendre_2, Legendre_4, Sampling_Window;
	double P_0[kbins]={}, P_2[kbins]={}, P_4[kbins]={}, Power_K1_K2[kbins][kbins]={};
	double counts_K[kbins]={}, counts_K1_K2[kbins][kbins]={};


//-----------------------------------K-------------------------------------------
	int	element = 0;

	for (int i=0; i<N_grid; i++){
		for (int j=0; j<N_grid; j++){
			for (int k=0; k<N_grid; k++){


				if(i!=0 || j!=0 || k!=0){

					// retrieving the frequencies from the 3D grid
					k1 = i<N_grid/2 ? i : i-N_grid;
					k2 = j<N_grid/2 ? j : j-N_grid;
					k3 = k<N_grid/2 ? k : k-N_grid;

					k_mag = k_F*sqrt(k1*k1 + k2*k2 + k3*k3);
					k_par = k_F*k3;
					mu_k  = k_par/k_mag;
					
					Legendre_2=(3*mu_k*mu_k-1)/2.;
					Legendre_4=(35*pow(mu_k, 4)-30*pow(mu_k, 2)+3)/8.;


					if(k_mag >= k_min && k_mag <= k_max){
						//correcting the fourier modes for the sampling effects from MAS
						Sampling_Window = pow(sin(k_mag*PI/(2*k_N))/((k_mag*PI)/(2*k_N)), p);
						delta_out[element][0] /= Sampling_Window;
						delta_out[element][1] /= Sampling_Window;

						//linear binning for power_spectra
						K = floor((k_mag-k_min)/binsize_k);
						P_0[K] += delta_out[element][0]*delta_out[element][0] + delta_out[element][1]*delta_out[element][1];
						P_2[K] += (delta_out[element][0]*delta_out[element][0] + delta_out[element][1]*delta_out[element][1])*Legendre_2;
						P_4[K] += (delta_out[element][0]*delta_out[element][0] + delta_out[element][1]*delta_out[element][1])*Legendre_4;
						
						counts_K[K]++;
						

					}
				}
				element++;
			}
		}
	}


	for(int i=0; i<kbins; i++){
		if(counts_K[i]!=0){
			P_0[i]/=(Norm*counts_K[i]);
			P_2[i]/=(Norm*counts_K[i]);
			P_4[i]/=(Norm*counts_K[i]);
			P_0[i]-=poisson_shot_noise; //*(1-2./3.*sin(PI*(k_min+(i+0.5)*binsize_k)/(2*k_N))*sin(PI*(k_min+(i+0.5)*binsize_k)/(2*k_N)));
		}

		//else Power_K[i] = -1000. ;

	}



//-----------------------------------K_par_K_perp---------------------------------------

	cout << "\33[1;32m" << "Computing power spectrum " << "\33[31m;" << " P(k||,k_|_)" << "\33[0m" << endl;

	element = 0;

	for (int i=0; i<N_grid; i++){
		for (int j=0; j<N_grid; j++){
			for (int k=0; k<N_grid; k++){

				if(i!=0 || j!=0 || k!=0){
					// retrieving the frequencies from the 3D grid
					k1 = i<N_grid/2 ? i : i-N_grid;
					k2 = j<N_grid/2 ? j : j-N_grid;
					k3 = k<N_grid/2 ? k : k-N_grid;

					k_mag = k_F*sqrt(k1*k1 + k2*k2 + k3*k3);

					k_par = k_F*sqrt(k3*k3);
					k_perp = k_F*sqrt(k1*k1+k2*k2);

					if((k_par >= k_min && k_par <= k_max)&&(k_perp >= k_min && k_perp <= k_max)){
						//correcting the fourier modes for the sampling effects from MAS
						Sampling_Window = pow(sin(k_mag*PI/(2*k_N))/((k_mag*PI)/(2*k_N)), p);
						delta_out[element][0] /= Sampling_Window;
						delta_out[element][1] /= Sampling_Window;

						//linear binning for power_spectra
						K1 = floor((k_par-k_min)/binsize_k);
						K2 = floor((k_perp-k_min)/binsize_k);
						Power_K1_K2[K1][K2] += delta_out[element][0]*delta_out[element][0] + delta_out[element][1]*delta_out[element][1];
						counts_K1_K2[K1][K2]++;
					}
				}
				element++;
			}
		}
	}


	for(int i=0; i<kbins; i++)
		for(int j=0; j<kbins; j++){
			if(counts_K1_K2[i][j]!=0){
				Power_K1_K2[i][j]/=(Norm*counts_K1_K2[i][j]);
				Power_K1_K2[i][j]-=poisson_shot_noise; //*(1-2./3.*sin(PI*(k_min+(i+0.5)*binsize_k)/(2*k_N))*sin(PI*(k_min+(i+0.5)*binsize_k)/(2*k_N)));
			}

			else Power_K1_K2[i][j] = -1000. ;
		}

//===================================================================================================================
//                      Writing output files with k, Pk, counts for all 3 spectra
//===================================================================================================================


	cout << "\33[33m" << "Writing output files "  << "\33[36m"<< "####\n";

	ofstream outfile(string("Pk_zHorizon_direct_sampling_s_Ngrid256_CIC_run"+number+".dat").c_str());
	for(int i=0; i<kbins; i++){

		outfile << std::scientific << setiosflags(ios::uppercase) << k_min+(i+0.5)*binsize_k << "\t" << P_0[i] << "\t" << 5*P_2[i] << "\t" << 9*P_4[i] << "\t" << counts_K[i] << endl;

	}

	outfile.close();



	ofstream outfile2(string("Pk_kpar_kperp_zHorizon_direct_sampling_s_Ngrid256_CIC_run"+number+".dat").c_str());
	for(int i=0; i<kbins; i++)
		for(int j=0; j<kbins; j++){

		outfile2 << k_min+(i+0.5)*binsize_k << "\t" << k_min+(j+0.5)*binsize_k << "\t" << Power_K1_K2[i][j] << "\t" << counts_K1_K2[i][j] << endl;

	}

	outfile2.close();


 	cout << "\33[33m" << "Free-ing arrays and memory "  << "\33[36m"<< "####\n";

	fftw_free(delta_in);

	fftw_free(delta_out);

	particles.clear();

	cout << "\33[1;32m" << "DONE!! "  << "\33[0;36m"<< "####" << "\33[0m"<< endl;

	
	}
	return 0;

}
