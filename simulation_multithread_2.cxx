/*
 * QWI_FEB-26.cxx
 * PY4115 MAJOR RESEARCH PROJECT
 * Copyright 2024 jmccw <jmcc0@DESKTOP-65IEBIH>
 */

#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <cmath>     // for std::abs
using namespace std;
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <complex.h>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Eigenvalues> 
#include <algorithm>
#include <fstream>
#include "gnuplot.cxx"
#include <fftw3.h>
#include <sstream>
#include <map>
#include <functional>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#define complex _Complex
const double pi = M_PI;

//Global Variables
bool over_ride_offsets = false;
vector<double> CBO_override, VBO_override;
double electron_charge = 1.0; //we use eV, so just let this = 1.
//double electric_field = 0.0;
double mass_electron = 1.0; //9.109534e-31; [kg]
double h_c = 1240.00; // eV nm
double hbar_sqaured_2m = 3.81; // [eV Amstrong squared]
double hbar = 1.05457182e-34; // [eV s]
double electron_mass = 9.1093837e-31; //kg
double e = 1.60217663e-19;

const int number_steps = 500;
const int max_matrix = 500;
const int ND = 1024;
const double reduce_matrix = 1.0; //cut down relatively the amount of computations - 4.0 considers only the first octant on matrices - data outside not required - should be checked for wells/materials


void print_matrix(double M[][number_steps]){
	cout << "[";						
	for(int i = 0; i < number_steps; i++){			//loop rows
		for(int j = 0; j < number_steps; j++){ 	//loop row elements
			cout << M[i][j];			//print element
			if(j != number_steps-1) cout << ", ";	// (decoration*)
			else if (j == number_steps-1 && i == number_steps-1) cout << "]";
		}
		cout << "\n ";
	}
}

// Input (destroyed during calculation of eigenvalues and eigenvectors)
double trimatrix_diag [max_matrix];
double trimatrix_subdiag [max_matrix];
// Output
double trimatrix_eigenvalue [max_matrix];
double trimatrix_eigenvector [max_matrix][max_matrix];
// only used during calculation of eigenvalues and eigenvectors
double trimatrix_result[max_matrix][max_matrix];
double* pointer_matrix [max_matrix];

#include "num_recipes_tridiagonal.cxx"

// as defined in example program
//~ void eigenvalues_and_eigenvectors_tridiagonal_matrix (int nn)
//~ {
	//~ int n, m;
	//~ for (n=0; n<nn; n++)
	//~ {
		//~ for (m=0; m<nn; m++)
		//~ {
			//~ if (n==m)
				//~ trimatrix_result [n][m] = 1.0;
			//~ else
				//~ trimatrix_result [n][m] = 0.0;
		//~ }
	//~ }
	//~ for (n=0; n<nn; n++)
		//~ pointer_matrix [n] = &(trimatrix_result[n][0]);
	//~ tqli (trimatrix_diag, trimatrix_subdiag, nn, pointer_matrix);
	//~ for (n=0; n<nn; n++)
	//~ {
		//~ trimatrix_eigenvalue [n] = trimatrix_diag[n];
		//~ for (m=0; m<nn; m++)
			//~ trimatrix_eigenvector [n][m] = trimatrix_result[m][n];
	//~ }
//~ }

//~ void sortTrimatrixEigenvaluesAndEigenvectors(double trimatrix_eigenvalue[], double trimatrix_eigenvector[][max_matrix], int nn) {
    //~ // Create a vector of pairs to store eigenvalue-eigenvector pairs
    //~ std::vector<std::pair<double, std::vector<double>>> eigenPairs(nn);
    
    //~ // Populate the vector of pairs
    //~ for (int i = 0; i < nn; ++i) {
        //~ eigenPairs[i].first = trimatrix_eigenvalue[i];
        //~ eigenPairs[i].second = std::vector<double>(trimatrix_eigenvector[i], trimatrix_eigenvector[i] + nn);
    //~ }
    
    //~ // Sort the vector of pairs based on eigenvalues
    //~ std::sort(eigenPairs.begin(), eigenPairs.end(), [](const std::pair<double, std::vector<double>>& a, const std::pair<double, std::vector<double>>& b) {
        //~ return a.first < b.first;
    //~ });
    
    //~ // Update the trimatrix_eigenvalue and trimatrix_eigenvector arrays with sorted values
    //~ for (int i = 0; i < nn; ++i) {
        //~ trimatrix_eigenvalue[i] = eigenPairs[i].first;
        //~ std::copy(eigenPairs[i].second.begin(), eigenPairs[i].second.end(), trimatrix_eigenvector[i]);
    //~ }
//~ }

std::vector<double> shiftEdgeToCenter(const std::vector<double>& input) {
    std::vector<double> shiftedVector = input;
    int halfSize = input.size() / 2;
    
    // Calculate the number of rotations needed to center the edge values
    int rotations = halfSize % input.size();

    // Perform cyclic rotations
    std::rotate(shiftedVector.begin(), shiftedVector.begin() + rotations, shiftedVector.end());

    return shiftedVector;
}

vector<double> convolution(vector<double> function1, vector<double> function2){
	//cout << "in convolution\n";
	if((int)function1.size() != (int)function2.size()) throw std::runtime_error("ERROR @ convolution() :: sizes not equal, "+to_string((int)function1.size())+" and "+to_string((int)function2.size())+"\n");
	const int N = (int)function1.size(); // Size of input signal -  assumes both vectors are same size
    double* input1 = (double*) fftw_malloc(sizeof(double) * N);
    double* input2 = (double*) fftw_malloc(sizeof(double) * N);
    fftw_complex* output1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* output2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Initialize input signal (e.g., with data from your computational physics project)
    for (int i = 0; i < N; ++i) {
        input1[i] = function1[i]; // Sample data
        input2[i] = function2[i]; // Sample data
    }

    // Create FFTW plan for forward transform
    fftw_plan plan1 = fftw_plan_dft_r2c_1d(N, input1, output1, FFTW_ESTIMATE);
    fftw_plan plan2 = fftw_plan_dft_r2c_1d(N, input2, output2, FFTW_ESTIMATE);

    // Execute forward transform
    fftw_execute(plan1);
    fftw_execute(plan2);
    
    fftw_complex* conv_mult = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    vector<double> result_vec(N);
    for (int i = 0; i < N; ++i) { //BEWARE OF DIVISION BY N HERE !!! I dont fully understand why this has to happen, but it conserves "area"
        conv_mult[i][0] = (output1[i][0] * output2[i][0] - output1[i][1] * output2[i][1])/N; // Real part
        conv_mult[i][1] = (output1[i][0] * output2[i][1] + output1[i][1] * output2[i][0])/N; // Imaginary part 
    }
    
    double* result = (double*) fftw_malloc(sizeof(double) * N);
    fftw_plan plan3 = fftw_plan_dft_c2r_1d(N, conv_mult, result, FFTW_ESTIMATE);
    fftw_execute(plan3);
    
    for (int i = 0; i < N; ++i) {
        result_vec[i] = result[i]; 
    }

    // Cleanup resources
    fftw_destroy_plan(plan1);
    fftw_free(input1);
    fftw_free(output1);
    fftw_destroy_plan(plan2);
    fftw_free(input2);
    fftw_free(output2);
    
    fftw_free(result);
    fftw_free(conv_mult);

    return result_vec;
}

vector<double> convolution_abso(vector<double> function1, vector<double> function2){
	//cout << "in convolution\n";
	if((int)function1.size() != (int)function2.size()) throw std::runtime_error("ERROR @ convolution() :: sizes not equal, "+to_string((int)function1.size())+" and "+to_string((int)function2.size())+"\n");
	const int N = (int)function1.size(); // Size of input signal -  assumes both vectors are same size
    double* input1 = (double*) fftw_malloc(sizeof(double) * N);
    double* input2 = (double*) fftw_malloc(sizeof(double) * N);
    fftw_complex* output1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* output2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Initialize input signal (e.g., with data from your computational physics project)
    for (int i = 0; i < N; ++i) {
        input1[i] = function1[i]; // Sample data
        input2[i] = function2[i]; // Sample data
    }

    // Create FFTW plan for forward transform
    fftw_plan plan1 = fftw_plan_dft_r2c_1d(N, input1, output1, FFTW_ESTIMATE);
    fftw_plan plan2 = fftw_plan_dft_r2c_1d(N, input2, output2, FFTW_ESTIMATE);

    // Execute forward transform
    fftw_execute(plan1);
    fftw_execute(plan2);
    
    fftw_complex* conv_mult = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    vector<double> result_vec(N);
    for (int i = 0; i < N; ++i) { //BEWARE OF DIVISION BY N HERE !!! I dont fully understand why this has to happen, but it conserves "area"
        conv_mult[i][0] = (output1[i][0] * output2[i][0] - output1[i][1] * output2[i][1]); // Real part
        conv_mult[i][1] = (output1[i][0] * output2[i][1] + output1[i][1] * output2[i][0]); // Imaginary part 
    }
    
    double* result = (double*) fftw_malloc(sizeof(double) * N);
    fftw_plan plan3 = fftw_plan_dft_c2r_1d(N, conv_mult, result, FFTW_ESTIMATE);
    fftw_execute(plan3);
    
    for (int i = 0; i < N; ++i) {
        result_vec[i] = result[i]; 
    }

    // Cleanup resources
    fftw_destroy_plan(plan1);
    fftw_free(input1);
    fftw_free(output1);
    fftw_destroy_plan(plan2);
    fftw_free(input2);
    fftw_free(output2);
    
    fftw_free(result);
    fftw_free(conv_mult);

    return result_vec;
}

class Material {
	private:
		std::string name;
		double electronAffinity; // eV
		double bandGap; // eV
		double eEffMass; // kg
		double lhEffMass; // kg
		double hhEffMass; // kg
		double refractive; // unitless

	public:
		// Constructor
		Material(std::string name_, double electronAffinity_, double bandGap_, double eEffMass_, double lhEffMass_, double hhEffMass_, double refractive_)
			: name(name_), electronAffinity(electronAffinity_), bandGap(bandGap_), eEffMass(eEffMass_), lhEffMass(lhEffMass_), hhEffMass(hhEffMass_), refractive(refractive_) {}
		// Default constructor
		Material() : name(""), electronAffinity(0.0), bandGap(0.0), eEffMass(0.0), lhEffMass(0.0), hhEffMass(0.0), refractive(0.0) {}

		// Getter methods
		std::string getName() const { return name; }
		double getBG() const { return bandGap; } //get band gap // eV
		double getEF() const { return electronAffinity; } //get electron affinity // eV
		double getRefractive() const { return refractive; } //get electron affinity // eV
		double getEffectiveMass(int p) const {
			if (p == 0) { // electron
				return eEffMass;
			} else if (p == 1) { // light hole
				return lhEffMass;
			} else if (p == 2) { // heavy hole
				return hhEffMass;
			} else {
				return 0.0;
			}
		}

		void display() const {
			std::cout << "Material: " << name << std::endl;
			std::cout << "Electron Affinity: " << electronAffinity << " eV" << std::endl;
			std::cout << "Band Gap: " << bandGap << " eV" << std::endl;
			std::cout << "n: " << refractive << "" << std::endl;
			std::cout << "Effective Masses" << std::endl;
			for (int i = 0; i <= 2; i++) {
				std::cout << i << " : " << this->getEffectiveMass(i) << " kg" << std::endl;
			}
		}
};

struct Layer {
	private:
		Material material;
		double thickness; // in Amstrong
		
	public:
		// Constructor 
		Layer(Material material_, double thickness_) : material(material_), thickness(thickness_) {}
		
		double getThickness() const { return thickness; }
		Material getMaterial() const { return material; }
		
		void display() {
			cout << "Thickness: " << this->getThickness() << " A" << endl;
			cout << "Material: " << this->getMaterial().getName() << endl;
		}
};

//Default band offset calculations | These functions utilise ^ Materials ^ and are utilised by Heterostructure (next class)
double VBO(Material material1, Material material2){ //conduction band offset default
	return material1.getEF()-material2.getEF();
}
double CBO(Material material1, Material material2){ //valence band offset default
	double calc1 = material1.getEF() + material1.getBG();
    double calc2 = material2.getEF() + material2.getBG();
    return calc2 - calc1;
}

class Heterostructure {
	private:
		//vector<Layer> layers;
		double heterostructure_thickness;
		double avg_refractive;
		double electric_field;
		vector<vector<double>> potential_;
		std::map<std::string, Material> materials; // Map to store materials
		std::vector<Layer> layers;
		
	public:
		//Constructor, reads in material parameters and layer file
		Heterostructure(const std::string& materialFileName, const std::string& layerFileName, double electric_field_) : potential_(3, vector<double>(number_steps)) {
			// Read material information from the input file
			std::ifstream materialFile(materialFileName);
			std::string line;

			// Read materials line by line
			while (std::getline(materialFile, line)) {
				std::istringstream iss(line);
				std::string materialName;
				double bandGap, electronAffinity, eEffMass, lhEffMass, hhEffMass, refractive;
				if (!(iss >> materialName >> bandGap >> electronAffinity >> eEffMass >> lhEffMass >> hhEffMass >> refractive)) {
					throw std::runtime_error("Error: Invalid material file format");
				}
				// Create Material object and store in the map
				Material material(materialName, bandGap, electronAffinity, eEffMass, lhEffMass, hhEffMass, refractive);
				materials[materialName] = material;
			}

			// Read layer information from the input file
			std::ifstream layerFile(layerFileName);
			int numLayers;
			if (!(layerFile >> numLayers)) {
				throw std::runtime_error("Error: Invalid layer file format");
			}

			// Read layers line by line
			for (int i = 0; i < numLayers; ++i) {
				std::string materialName;
				double thickness;
				if (!(layerFile >> materialName >> thickness)) {
					throw std::runtime_error("Error: Invalid layer file format");
				}
				// Retrieve Material object from the map using the material name
				auto it = materials.find(materialName);
				if (it == materials.end()) {
					throw std::runtime_error("Error: Unknown material name in layer file");
				}
				Material material = it->second;
				// Create Layer object and add to layers vector
				layers.push_back(Layer(material, thickness));
			}
			double total_thickness = 0.0;
			for(const Layer& layer : layers) {
				total_thickness += layer.getThickness();
			}
			this->heterostructure_thickness = total_thickness;
			resetPotential();
			double temp = 0.0; //avg refr
			int count = 0;
			for(const Layer& layer : layers) {
				temp += layer.getMaterial().getRefractive();
				count++;
			}
			this->avg_refractive = temp / count;
			this->electric_field=electric_field_;
		}
		
		double getAvgRefractive() {
			return this->avg_refractive;
		}
		
		double getElectricField() { return this->electric_field; }
		void setElectricField(double e_new) { this->electric_field = e_new; }
		
		double getPotential(int particle, double x, int x_int) {
			if (particle == 1 || particle == 2){ // if particle is hole
				return this->potential_[particle][x_int] - electron_charge*(this->electric_field)*(x-this->getThickness()/2.0); // not as efficient as defining a strict function globally but gives ALOT of freedom (only adds +(number of layers)*3 computations at most  
			}
			else if (particle == 0){             // else it is electron
				return this->potential_[particle][x_int] + electron_charge*(this->electric_field)*(x-this->getThickness()/2.0);
			}
			return potential_[particle][x_int];
		}
		
		vector<Layer> getLayers() const {
			return this->layers;
		}
		
		// Getters / Setters
		double getThickness() const {
			return this->heterostructure_thickness;
		}
		
		void resetPotential(){
			double x_set[number_steps];
			double delta_x = this->getThickness()/(number_steps+1);
			for(int p = 0; p < 3; p++){
				for(int i = 0; i < number_steps; i++){
					x_set[i] = delta_x + i*delta_x;
					this->potential_[p][i] = potential(p, x_set[i]);
				}
			}
		}
		
		double eff_mass(int particle, double x) const { //
			double material_threshold = 0.0;
			for(const Layer& layer : layers) {
				material_threshold += layer.getThickness();		//Material threshold determines whether or not we need to add a band offset to the potential at x
				if(x <= material_threshold) {			
					return layer.getMaterial().getEffectiveMass(particle);
				}
			}
			throw std::runtime_error("Error: Unable to find effective mass for the given position."+to_string(x));
			return 0.0; //this will never happen under normal circumstances
		}
		
		double potential(int particle, double x) { // particle = 0,1,2 -> elecron, lh, hh
			double material_threshold = 0.0;
			double U = 0.0, V = 0.0;
			int i = 0;

			for(Layer layer : this->layers) {
				material_threshold += layer.getThickness();		//Material threshold determines whether or not we need to add a band offset to the potential at x
				if(x >= material_threshold) {
					if(over_ride_offsets == false) {				// so if x is not within material of layer[i] add/subtract relevant band offset
						U += CBO(this->layers[i].getMaterial(), this->layers[i+1].getMaterial());
						V += VBO(this->layers[i].getMaterial(), this->layers[i+1].getMaterial());						
					} else if (over_ride_offsets == true) {
						U += CBO_override[i]; //TEMPORARY SOLUTION ??
						V += VBO_override[i]; //TEMPORARY SOLUTION ??
					}
					i++;
				}
				else { // x within new material threshold, return relevant potential - if used correctly, there should never be an error here (one might occur when an x input is out of analysis bounds)
					if (particle == 1 || particle == 2){ // if particle is hole
						return V; // not as efficient as defining a strict function globally but gives ALOT of freedom (only adds +(number of layers)*3 computations at most  
					}
					else if (particle == 0){             // else it is electron
						return U;
					}
				}
			}

			throw std::runtime_error("Error: Unable to find valid potential for the given position."+to_string(x));
			return 0.0; //this will never happen under normal circumstances but leaving this here to make compiler happy.
		}
		
		void intermixPotential(double strength) {
			vector<double> Gauss_y(number_steps), x_het(number_steps);
			vector<vector<double>> Potential_y(3, vector<double>(number_steps));
			double sum = 0, x_0 = this->getThickness()/2.0, sum_pot_init = 0;
			double sigma = strength;//*this->heterostructure_thickness;
			//cout << "Creating Gaussian.\n";
			for(int p = 0; p < 3; p++){
				for(int i = 0; i < number_steps; i++){
					double delta_x = (this->getThickness()/(number_steps+1));
					x_het[i] = delta_x + i*delta_x;
					Potential_y[p][i] = this->potential(p, x_het[i]);
					Gauss_y[i] = 1.0/(sigma*sqrt(2*pi))*exp(-(x_het[i]-x_0)*(x_het[i]-x_0)/(2.0*sigma*sigma));
				}
			}
			//cout << "Done.\n";
			
			// Normalise Gaussian - so that potential is "CTS"
			for(int i = 0; i < number_steps; i++){
				sum += Gauss_y[i];
				sum_pot_init += Potential_y[0][i];
			}

			for(int i = 0; i < number_steps; i++){
				Gauss_y[i] = Gauss_y[i]/abs(sum);///sum;
			}

			//cout << "Conv.\n";
			vector<vector<double>> Potential_QWI(3, vector<double>(number_steps));
			
			Potential_QWI[0] = convolution(Potential_y[0], Gauss_y);
			Potential_QWI[1] = convolution(Potential_y[1], Gauss_y);
			Potential_QWI[2] = convolution(Potential_y[2], Gauss_y);
			
			//cout << "Conv returned.\n";
			Potential_QWI[0] = shiftEdgeToCenter(Potential_QWI[0]);
			Potential_QWI[1] = shiftEdgeToCenter(Potential_QWI[1]);
			Potential_QWI[2] = shiftEdgeToCenter(Potential_QWI[2]);
			
			// Print convolution result
			double sum_pot_QWI = 0;
			//cout << "Conv returned.\n";
			for(int i = 0; i<number_steps; i++){
				sum_pot_QWI += Potential_QWI[0][i];
			}
			//cout << "Area before Intermixing > " << sum_pot_init << endl;
			//cout << "Area after Intermixing > " << sum_pot_QWI << endl;
			for(int p = 0; p < 3; p++){
				this->potential_[p] = Potential_QWI[p];
			}		
		}
		
		void display() { // print heterostructure details
			cout << "Layers: " << endl;
			int i = 1;
			for(Layer layer : layers){
				cout << i << " : " << layer.getMaterial().getName() << " : " << layer.getThickness() << " A"<< endl;
				i++;
			}	
			cout << "Total Thickness : " << this->getThickness() << " A" << endl;
			cout << "\nMaterial properties by layer:"<<endl;	
			for(Layer layer : layers) {
				layer.getMaterial().display();
			}	
		}
};

double relative_energy(double energy, int p, Heterostructure& QW) { // accept a 2 paramters [energy relative to well of..][particle type]
    // corrects an (inital) calculated energy from QW solution for respective band of type particle type p
    //~ double EF_offset = electron_charge*electric_field*QW.getThickness();
    double E_REL;
    
    // this should be okay with and without override on
    double BG = abs(QW.getLayers()[0].getMaterial().getBG());
	for (int i = 1; i < (int)QW.getLayers().size() - 1; ++i) {
		double band_gap = abs(QW.getLayers()[i].getMaterial().getBG());
		BG = min(BG, band_gap);
	}
	
    if(p==0){
		double max_CBO = CBO(QW.getLayers()[0].getMaterial(), QW.getLayers()[1].getMaterial());
		double pos = CBO(QW.getLayers()[0].getMaterial(), QW.getLayers()[1].getMaterial());
		for (int i = 1; i < (int)QW.getLayers().size()-1; ++i) {
			pos += CBO(QW.getLayers()[i].getMaterial(), QW.getLayers()[i+1].getMaterial());
			max_CBO = std::max(abs(max_CBO), abs(pos));
		}
		//cout << "max_CBO" << max_CBO << endl;
		E_REL = energy + BG + max_CBO;
    } else { // i.e. p==1 or p==2
		double max_VBO = VBO(QW.getLayers()[0].getMaterial(), QW.getLayers()[1].getMaterial());
		double pos = VBO(QW.getLayers()[0].getMaterial(), QW.getLayers()[1].getMaterial());
		for (int i = 1; i < (int)QW.getLayers().size()-1; ++i) {
			pos += VBO(QW.getLayers()[i].getMaterial(), QW.getLayers()[i+1].getMaterial());
			max_VBO = std::max(abs(max_VBO), abs(pos));
		}
		//cout << "max_VBO" << max_VBO << endl;
		E_REL = -(max_VBO) - energy;	
	}
    return E_REL;
}

// Define a function to solve the heterostructure
void solve(Heterostructure& heterostructure, std::vector<double>& x_out,
           std::vector<std::vector<double>>& energies,
           std::vector<std::vector<std::vector<double>>>& eigenVectors) {
	    
	double delta_x = heterostructure.getThickness() / (number_steps + 1.0); 
	double x[number_steps];

	x[0] = delta_x; 
	x_out[0] = x[0]; 

	for(int i = 1; i<number_steps; i++){ //initialise x (once)
		x[i] = x[0] + i*delta_x;
		x_out[i] = x[i];
	}

	//implementation of Kowano & Kito : 'Optical Waveguide Analysis' : solution of SE with effective mass approximation for all bands/particles
    for(int p = 0; p<=2; p++){ //for all particle types
		//initialise solution matrix M
		double M[number_steps][number_steps];
		for(int i = 0; i<number_steps; i++){
			for(int j = 0; j<number_steps; j++){
				M[i][j] = 0;
			}
		}

		double alpha_w[number_steps], alpha_e[number_steps], alpha_x[number_steps];
		alpha_w[0] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[0])+heterostructure.eff_mass(p,x[0]));
        alpha_e[0] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[0])+heterostructure.eff_mass(p,x[1]));
        alpha_x[0] = -alpha_w[0]-alpha_e[0];
        
        M[0][0] = alpha_x[0] + heterostructure.getPotential(p,x[0],0);
        M[0][1] = alpha_e[0];

        for(int nr = 1; nr < number_steps-1; nr++){
            alpha_w[nr] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[nr])+heterostructure.eff_mass(p,x[nr-1]));
            alpha_e[nr] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[nr])+heterostructure.eff_mass(p,x[nr+1]));
            alpha_x[nr] = -alpha_w[nr]-alpha_e[nr];

            M[nr][nr-1] = alpha_w[nr];    //sub-diagonal
            M[nr][nr] = alpha_x[nr] + heterostructure.getPotential(p,x[nr],nr); //diagonal
            M[nr][nr+1] = alpha_e[nr];   //upper diagonal   
		}

        alpha_w[number_steps-1] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[number_steps-1])+heterostructure.eff_mass(p,x[number_steps-1-1]));
        alpha_e[number_steps-1] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[number_steps-1])+heterostructure.eff_mass(p,x[number_steps-1])); // assuming m(x_edge-dx) = m(x_edge) as boundary condition
        alpha_x[number_steps-1] = -alpha_w[number_steps-1]-alpha_e[number_steps-1];
        M[number_steps-1][number_steps-2] = alpha_w[number_steps-1];
        M[number_steps-1][number_steps-1] = alpha_x[number_steps-1] + heterostructure.getPotential(p,x[number_steps-1],number_steps-1);
		
		//print_matrix(M);
		
		//~ trimatrix_diag[0] = M[0][0];
		//~ trimatrix_subdiag[0] = M[0][1];
		//~ for (int i = 1; i < number_steps; ++i) {
			//~ // Assign eigenvalue to energies
			//~ trimatrix_subdiag[i] = M[i][i-1];
			//~ trimatrix_diag[i] = M[i][i];
		//~ }
		
		//~ eigenvalues_and_eigenvectors_tridiagonal_matrix (number_steps);
		//~ sortTrimatrixEigenvaluesAndEigenvectors(trimatrix_eigenvalue, trimatrix_eigenvector, number_steps);
		//~ for (int i = 0; i < number_steps; ++i) {
			//~ // Assign eigenvalue to energies
			//~ energies[p][i] = trimatrix_eigenvalue[i];
			//~ // Assign eigenvector to eigenVectors
			//~ for (int j = 0; j < number_steps; ++j) {
				//~ eigenVectors[p][i][j] = trimatrix_eigenvector[i][j];
				//~ //test[j] = eigenvectors(0, j);
			//~ }
		//~ }
		
		
		
		
		
		//cout << "test0 ";
		//~ //solve Matrix (using Eigen)
		Eigen::MatrixXd M_eigen(number_steps, number_steps);
		//cout << "test1 ";
		for (int i = 0; i < number_steps; ++i) { // really not a speed issue here
			for (int j = 0; j < number_steps; ++j) {
				M_eigen(i, j) = M[i][j];
			}
		}
		//cout << "test2 ";
		
		// Solve the matrix using Eigen's EigenSolver
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M_eigen);
		Eigen::VectorXd eigenvalues = solver.eigenvalues();
		//cout << "TEST4 ";
		Eigen::MatrixXd eigenvectors = solver.eigenvectors();
		//cout << "TEST5 ";

		// Print the eigenvalues
		//cout << "test3 ";
		for (int i = 0; i < number_steps; ++i) {
			energies[p][i] = eigenvalues(i);			// Assign eigenvalue to energies
			for (int j = 0; j < number_steps; ++j) {	// Assign eigenvector to eigenVectors
				eigenVectors[p][i][j] = eigenvectors(j, i);
			}
		}
		//cout << "test4 " << endl;
	}
}


vector<vector<vector<double>>> findTransitions(vector<vector<double>> energies){
	vector<vector<vector<double>>> E_gap(number_steps, vector<vector<double>>(number_steps, vector<double>(2))); //Gap between [electron state] and [hole state] of [hole type]
	for (int k=0; k<2; k++){ // for 2 hole types
        for (int i=0; i<(int)energies[0].size(); i++){ // for all electrons
            for (int j=0; j<(int)energies[0].size(); j++){ // and all hole states
                E_gap[i][j][k] = abs(energies[0][i]-energies[k+1][j]);
			}
		}
	}
    return E_gap; // a matrix/array of energy transtions indexed as [electron state][hole state][hole type] 
}



vector<vector<double>> findEnergiesRelative(vector<vector<double>> energies, Heterostructure& heterostructure) {
	vector<vector<double>> energies_relative(3, vector<double>(number_steps));
	for (int k=0; k<3; k++){ // for 3 particle types
        for (int i=0; i<(int)energies[0].size(); i++){ //for all states
            energies_relative[k][i] = relative_energy(energies[k][i], k, heterostructure);
		}
	}
	return energies_relative; //returns [3][number_solutions] matrix of energies relative to their ectual value in well structure
}



vector<double> overlapIntegral(vector<double> vector1, vector<double> vector2){
	if((int)vector1.size() != (int)vector2.size()) throw std::runtime_error("vector sizes not equal, "+to_string((int)vector1.size())+" and "+to_string((int)vector2.size()));
	vector<double> overlap((int)vector1.size()), vector1_dummy((int)vector1.size()), vector2_dummy((int)vector2.size());
	double N1 = 0.0, N2 = 0.0;
	// possibly might need to declare some dummy vectors ?? - this was a python issue
	for(int i = 0; i < (int)overlap.size(); i++){
		N1 += fabs(vector1[i])*fabs(vector1[i]);
		N2 += fabs(vector2[i])*fabs(vector2[i]);
	}
	for(int i = 0; i < (int)overlap.size(); i++){
		vector1_dummy[i] = vector1[i]/N1;
		vector2_dummy[i] = vector2[i]/N2;
		overlap[i] = vector1_dummy[i]*vector2_dummy[i];
	}
	return overlap;
}



double I_squared(vector<double> vector1, vector<double> vector2){
	vector<double> overlap = overlapIntegral(vector1, vector2);
	double I_squared = 0;
	for(int i = 0; i < (int)overlap.size(); i++) I_squared += fabs(overlap[i]);
	I_squared *= I_squared; // square result
	return I_squared;
}



vector<vector<vector<double>>> findOverlapsAll(vector<vector<vector<double>>> wavefunctions) {
	vector<vector<vector<double>>> I_squared_matrix(number_steps, vector<vector<double>>(number_steps, vector<double>(2)));
	// [electron state][hole state][hole type]
	for (int k=0; k<2; k++) { // for 2 hole types
        for (int i=0; i<(int)wavefunctions[0].size(); i++) { // for all electrons
			vector<double> state1 = wavefunctions[0][i];
            for (int j=0; j<(int)wavefunctions[0].size(); j++) { // and all hole states
				vector<double> state2 = wavefunctions[1+k][j];
                I_squared_matrix[i][j][k] = I_squared(state1, state2);
			}
		}
	}
	return I_squared_matrix;
}



vector<double> pad_func_zeros(vector<double> func){
	vector<double> func_padded(2*func.size());
	for(int i = 0; i<(int)(0.25*func_padded.size()); i++) func_padded[i] = 0.0;
	for(int i = (int)(0.75*func_padded.size()); i<(int)(func_padded.size()); i++) func_padded[i] = 0.0;
	int j = 0;
	for(int i = (int)(0.25*func_padded.size()); i<(int)(0.75*func_padded.size()); i++){
		func_padded[i] = func[j];
		j++;
	}
	return func_padded;
}



vector<vector<vector<double>>> wavelengthTransformation(vector<vector<vector<double>>> data_in) {
	//The intention here is that the user passes in E_GAP to find an adjacent matrix in wavelength terms.
	vector<vector<vector<double>>> data_out(number_steps, vector<vector<double>>(number_steps, vector<double>(2)));
	for (std::size_t i = 0; i < data_out.size(); ++i) {
		for (std::size_t j = 0; j < data_out[i].size(); ++j) {
			for (std::size_t k = 0; k < data_out[i][j].size(); ++k) {
				data_out[i][j][k] = h_c / data_in[i][j][k]; //h_c = 1240 eV nm :: so [eV nm / eV] = [nm]
			}
		}
	}
	return data_out;
}



vector<vector<vector<double>>> findOverlapsOnly(vector<vector<vector<double>>> wavefunctions) {
	vector<vector<vector<double>>> I_squared_matrix(number_steps, vector<vector<double>>(number_steps, vector<double>(2)));
	// [electron state][hole state][hole type]
	for (int k=0; k<2; k++) { // for 2 hole types
        for (int i=0; i<(int)(wavefunctions[0].size()/reduce_matrix); i++) { // for all electrons
			vector<double> state1 = wavefunctions[0][i];
            for (int j=0; j<(int)(wavefunctions[0].size()/reduce_matrix); j++) { // and all hole states
				vector<double> state2 = wavefunctions[1+k][j];
                I_squared_matrix[i][j][k] = I_squared(state1, state2);
			}
		}
	}
	return I_squared_matrix;
}



double Rydberg(double m_e, double mass_lh, double mass_hh, int particle, double n_avg) {
       double h=6.62377E-34;
       double q0=1.60219E-19;
       double mass,mu,R;
       double m0 = 9.109534E-31; //[kg]
       //double J_eV = 6.242E+18;
       m_e*=m0;
       mass_lh*=m0;
       mass_hh*=m0;
       double e0=8.8542E-12;
       double eps=n_avg*n_avg*e0;
       //cout << "R1: " << eps <<endl;
       //cout << "n: " << n_avg <<endl;
       if (particle==0) { //lh
		   mass = (4.0 * mass_hh * mass_lh) / (mass_hh + 3.0 * mass_lh);
       } else { //hh
           mass = (4.0 * mass_hh * mass_lh) / (3.0 * mass_hh + mass_lh);
       }
       mu = (m_e * mass) / (m_e + mass);
       //cout << "mass_hh: " << mass_lh <<endl;
       //cout << "mass_hh: " << mass_hh <<endl;
       //cout << "mu: " << mu <<endl;
       R = q0*q0*q0 * mu / (8.0 * eps*eps * h*h); // in eV
       //cout << "R: " << R <<endl;
       return R;
}




int delta(double limit, double in, double precision){
	if(in >= limit-precision && in < limit+precision) {
		//cout << "\n" << limit << "=" << in;
		return 1;
	}
	else return 0;
}



double findYValue(const vector<vector<vector<double>>>& E_GAP, const vector<vector<vector<double>>>& I_squared_Matrix, double x, double tolerance) {
    for (size_t i = 0; i < E_GAP.size(); ++i) {
        for (size_t j = 0; j < E_GAP[i].size(); ++j) {
            for (size_t k = 0; k < E_GAP[i][j].size(); ++k) {
                if (fabs(E_GAP[i][j][k] - x) < tolerance) {
                    return I_squared_Matrix[i][j][k];
                }
            }
        }
    }
    return 0;  // If x value is not found within tolerance, return 0
}

void plot_potential(Heterostructure& QW){
	double potential1[number_steps], potential2[number_steps];
	double x_test[number_steps];
	for(int i = 0; i<number_steps; i++){
		x_test[i] = (QW.getThickness()/(number_steps+1)) + i*(QW.getThickness()/(number_steps+1));
		potential1[i] = relative_energy(QW.getPotential(2, x_test[i],i), 2, QW);
		potential2[i] = relative_energy(QW.getPotential(0, x_test[i],i), 0, QW);
	}
    gnuplot_two_functions ("Potentials", "linespoints", "x [A]", "Effective Potentials",
			x_test, potential1, number_steps, "Conduction", x_test, potential2, number_steps, "Valence");
}




double findGroundStateQWI(Heterostructure& heterostructure, int p) {  
	double delta_x = heterostructure.getThickness() / (number_steps + 1.0); 
	double x[number_steps];
	x[0] = delta_x; 
	//x_out[0] = x[0]; 

	for(int i = 1; i<number_steps; i++){ //initialise x (once)
		x[i] = x[0] + i*delta_x;
		//x_out[i] = x[i];
	}

	//implementation of Kowano & Kito : 'Optical Waveguide Analysis' : solution of SE with effective mass approximation for all bands/particles
		//initialise solution matrix M
		double M[number_steps][number_steps];
		for(int i = 0; i<number_steps; i++){
			for(int j = 0; j<number_steps; j++){
				M[i][j] = 0;
			}
		}

		double alpha_w[number_steps], alpha_e[number_steps], alpha_x[number_steps];
		alpha_w[0] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[0])+heterostructure.eff_mass(p,x[0]));
        alpha_e[0] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[0])+heterostructure.eff_mass(p,x[1]));
        alpha_x[0] = -alpha_w[0]-alpha_e[0];
        
        M[0][0] = alpha_x[0] + heterostructure.getPotential(p,x[0],0);
        M[0][1] = alpha_e[0];

        for(int nr = 1; nr < number_steps-1; nr++){
            alpha_w[nr] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[nr])+heterostructure.eff_mass(p,x[nr-1]));
            alpha_e[nr] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[nr])+heterostructure.eff_mass(p,x[nr+1]));
            alpha_x[nr] = -alpha_w[nr]-alpha_e[nr];

            M[nr][nr-1] = alpha_w[nr];    //sub-diagonal
            M[nr][nr] = alpha_x[nr] + heterostructure.getPotential(p,x[nr],nr); //diagonal
            M[nr][nr+1] = alpha_e[nr];   //upper diagonal   
		}

        alpha_w[number_steps-1] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[number_steps-1])+heterostructure.eff_mass(p,x[number_steps-1-1]));
        alpha_e[number_steps-1] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[number_steps-1])+heterostructure.eff_mass(p,x[number_steps-1])); // assuming m(x_edge-dx) = m(x_edge) as boundary condition
        alpha_x[number_steps-1] = -alpha_w[number_steps-1]-alpha_e[number_steps-1];
        M[number_steps-1][number_steps-2] = alpha_w[number_steps-1];
        M[number_steps-1][number_steps-1] = alpha_x[number_steps-1] + heterostructure.getPotential(p,x[number_steps-1],number_steps-1);
		
		//solve Matrix (using Eigen)
		Eigen::MatrixXd M_eigen(number_steps, number_steps);
		//cout << "test1 ";
		for (int i = 0; i < number_steps; ++i) { // really not a speed issue here
			for (int j = 0; j < number_steps; ++j) {
				M_eigen(i, j) = M[i][j];
			}
		}
		
		// Solve the matrix using Eigen's EigenSolver
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M_eigen);
		Eigen::VectorXd eigenvalues = solver.eigenvalues();
		//cout << "TEST4 ";
		//Eigen::MatrixXd eigenvectors = solver.eigenvectors();
		//cout << "TEST5 ";
		vector<double> energies(number_steps);
		for (int i = 0; i < number_steps; ++i) { // really not a speed issue here
			energies[i] = eigenvalues(i);
		}
		double w = *min_element(energies.begin(), energies.end());
		return relative_energy(w, p, heterostructure);
}



std::mutex mtx;
void findGroundStateQWIAndReturnResult(Heterostructure& heterostructure, int p, double& result) {
    double val = findGroundStateQWI(heterostructure, p);
    // Lock the mutex before accessing shared resources
    mtx.lock();
    result = val;
    // Unlock the mutex after updating the shared resources
    mtx.unlock();
}



void ShiftBandgapQWI(Heterostructure& QW, double shift_amount){ 					// shift amount in nm ?
    std::vector<std::thread> threads;
	QW.resetPotential();
	if(shift_amount < 0.001) return;
	// solve QW with no Electric field
	QW.setElectricField(0.0);
		
	// ------- write code here ---------
	vector<double> energies(3);
    
    // Launch a thread for each value of p
    for (int i = 0; i < 3; ++i) {
        threads.emplace_back(findGroundStateQWIAndReturnResult, std::ref(QW), i, std::ref(energies[i]));
    }
    
    // Join all threads with the main thread
    for (auto& thread : threads) {
        thread.join();
    }
	// Lambda function to be executed by each thread
	cout << "here1 1st";
	// Wait for all tasks to complete
	// ------- write code here ---------

	cout << "done" << endl;
	//cout << energies[0] << endl;
	//cout << energies[1] << endl;
	//cout << energies[2] << endl;
	//cout << "done" << endl;
	
	double initBG_valence;
	if (energies[1] > energies[2]) initBG_valence = energies[1];
	else initBG_valence = energies[2];
	cout << energies[2] << " " << energies[1]<< " " << energies[0] << endl;
	double initBG = h_c / abs(initBG_valence-energies[0]); // [nm] relative_energy(energies[0],0,QW)
	
	
	// Bisection Method inspired routine
	double prox = 0.5; 						// nm
	double sigma = QW.getThickness()/2.0; 	// nm - max possible Intermixing
	double sigma0 = 0;
	double sigma_prev = sigma, temp;
	double BG_valence;
	int max_steps = 100;
	int count = 0;
	
	cout << "Starting QWI routine."<<endl;
	cout << "Initial Bandgap " << initBG << " [nm]" << endl;
	for(int i = 0; i < max_steps; i++) {
		std::vector<std::thread> threads_;
		QW.resetPotential();
		QW.intermixPotential(sigma);

		// Enqueue tasks for each value of p
		for (int p = 0; p < 3; ++p) {
			threads_.emplace_back(findGroundStateQWIAndReturnResult, std::ref(QW), p, std::ref(energies[p]));
		}
		
		// Join all threads with the main thread
		for (auto& thread : threads_) {
			thread.join();
		}
		cout << energies[2] << " " << energies[1]<< " " << energies[0] << endl;
		
		if (energies[1] < energies[2]) BG_valence = energies[2];
		else BG_valence = energies[1];
		double BG_diff = abs(initBG - h_c / abs(BG_valence-energies[0])); // [nm]
		
		temp = sigma;
		if (BG_diff > shift_amount) {
			sigma = (sigma0+sigma)/2.0;
			sigma_prev = temp;
			// don't change sigma0
		} else if (BG_diff < shift_amount) { 
			sigma0 = sigma;
			sigma = (sigma_prev+sigma)/2.0;
		}
		cout << "Progress towards Bandgap " << BG_diff << " / " << shift_amount << endl;
		if ((BG_diff < shift_amount+prox) && (BG_diff > shift_amount-prox)) {
			cout << "New Bandgap " << h_c / abs(BG_valence-energies[0]) << " [nm]"<<endl;
			break;
		}
		
		count++;
	}
	cout << "Found solution in " << count << " iterations."<<endl;
}



struct SimulationParameters {
    bool intermixingEnabled;
    double targetBandgapShift;
    int numElectricFields;
    double maxElectricField;
};



SimulationParameters readSimulationParameters(const std::string& filename) {
    SimulationParameters params;

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Failed to open file.. returning default parameters " << filename << std::endl;
        // might want to handle this error condition appropriately
        // For now, just return default-initialized parameters
        return params;
    }

    // Read line 1: intermixing parameters
    std::string intermixing;
    double targetBandgapShift;
    if (!(file >> intermixing >> targetBandgapShift)) {
        std::cerr << "Error: Failed to read intermixing parameters" << std::endl;
        return params;
    }
    params.intermixingEnabled = (intermixing == "TRUE");
    params.targetBandgapShift = targetBandgapShift;

    // Read line 2: number of electric fields
    if (!(file >> params.numElectricFields)) {
        std::cerr << "Error: Failed to read number of electric fields" << std::endl;
        return params;
    }

    // Read line 3: max electric field
    if (!(file >> params.maxElectricField)) {
        std::cerr << "Error: Failed to read max electric field" << std::endl;
        return params;
    }

    return params;
}

vector<vector<double>> interpolateY(const std::vector<double>& x_values, const std::vector<double>& y_values,
                                 double x_min, double x_max, size_t num_points, double precision) {
    std::vector<double> x_vector;
    std::vector<double> y_vector;
	
    // Discretize the x-axis
    for (size_t i = 0; i < num_points; ++i) {
        double x = x_min + (x_max - x_min) * i / (num_points - 1);
        x_vector.push_back(x);

        // Find the closest x-value in the given x_values
        size_t closest_index = 0;
        double min_distance = std::abs(x_values[0] - x);
        for (size_t j = 1; j < x_values.size(); ++j) {
            double distance = std::abs(x_values[j] - x);
            if (distance < min_distance) {
                min_distance = distance;
                closest_index = j;
            }
        }

        // Interpolate y-value if x is close to a given x-value, otherwise set y to 0
        if (min_distance < precision) {  // adjust threshold as needed
            y_vector.push_back(y_values[closest_index]);
            cout << "check"<<endl;
        } else {
            y_vector.push_back(0.0);
        }
    }
    vector<vector<double>> result(ND, vector<double>(2));
	for(int i = 0; i < ND; i++){
		result[i][0] = x_vector[i];
		result[i][1] = y_vector[i];
	}
	cout << "test";
	return result;
}

struct Point {
    double x;
    double y;
};

void generateXYVector(std::vector<Point>& xyPairs, double xRangeStart, double xRangeEnd, int numPoints, std::vector<double>& xVector, std::vector<double>& yVector) {
    // Step 1: Define the range of the x-axis
    double step = (xRangeEnd - xRangeStart) / (numPoints - 1);

    // Step 2: Initialize y-axis vector with zeros
    xVector.clear();
    yVector.clear();
    for (int i = 0; i < numPoints; ++i) {
        xVector.push_back(xRangeStart + i * step);
        yVector.push_back(0.0);
    }

    // Step 3: For each x-y pair, find the closest point in the x-axis vector
    for (const auto& point : xyPairs) {
        int closestIndex = 0;
        double minDistance = std::abs(xVector[0] - point.x);
        for (int i = 1; i < numPoints; ++i) {
            double distance = std::abs(xVector[i] - point.x);
            if (distance < minDistance) {
                minDistance = distance;
                closestIndex = i;
            }
        }
        yVector[closestIndex] = point.y;
    }
}

void script(int j, Heterostructure& QW_init, double electric_field){
	
	//electric_field = 0.0;
	//double delta_field = params.maxElectricField/(double)params.numElectricFields; //in put in
	//delta_field *= 0.0001; //conversion to V/Amstrong
	
	Heterostructure QW = QW_init;
	QW.setElectricField(electric_field);
	
	vector<double> x(number_steps); //length element
	vector<vector<double>> energies(3, vector<double>(number_steps)); //particle, energy_level
	vector<vector<vector<double>>> eigenVectors(3, vector<vector<double>>(number_steps, vector<double>(number_steps)));

	vector<double> x_limits;
	//electric_field = j*delta_field;
	cout << "E-Field = " << electric_field << endl;
	cout << "Run " << j << endl;
	
	vector<vector<vector<double>>> RESULT;
	
	//if(j==0) electric_field = 0.0001;
	solve(QW, x, energies, eigenVectors);
	cout << "energies[0][0][0] " << electric_field << " " << energies[0][0] <<endl;

	// sort results
	vector<vector<double>> energies_relative(3, vector<double>(number_steps));
	energies_relative = findEnergiesRelative(energies, QW);
	//cout << "energies_relative:" << energies_relative[0][0] << " " << energies_relative[2][0] << endl;
	//cout << "energies_relative1-1:" << energies_relative[0][1] << " " << energies_relative[2][1] << endl;
	//cout << "energies_relative2-2:" << energies_relative[0][2] << " " << energies_relative[2][2] << endl;
	//cout << "energies_relative0-2:" << energies_relative[0][0] << " " << energies_relative[2][2] << endl;
	// Absorption Routine (write to files for python)
	cout << "Sorting transitions.\n";
	vector<vector<vector<double>>> E_GAP = findTransitions(energies_relative); //Gap between [electron state] and [hole state] of [hole type]
	cout << "E_GAP[0][0][1] " << electric_field << " " << E_GAP[0][0][1] <<endl;
	//cout << "Converting to wavelength.\n";
	vector<vector<vector<double>>> E_GAP_WL = wavelengthTransformation(E_GAP);
	cout << "E_GAP[0][0][1] " << electric_field << " " << E_GAP_WL[0][0][1] <<endl;
	//cout << "transitions:" << E_GAP_WL[0][0][0] << " " << E_GAP_WL[0][0][1] << endl;
	//~ for(int j = 0; j<(int)E_GAP[0].size(); j++){
		//~ for(int i = 0; i<(int)E_GAP[0].size(); i++){
			//~ cout << i << " " << j << " : EGAP " << E_GAP[i][j][0]<<endl;;
		//~ }
	//~ }
	
	//~ double x_array_2[number_steps], y_array_2[number_steps];
	//~ for (int i = 0; i < number_steps; i++) {
		//~ y_array_2[i] = eigenVectors[0][0][i];
		//~ x_array_2[i] = i;
	//~ }
	//~ gnuplot_two_functions ("Numerical Eigenfunctions", "linespoints", "x", "Eigenfunction |E_1> (x)",
		//~ x_array_2, y_array_2, number_steps, "Ground", x_array_2, y_array_2, number_steps, "5th excited state");
	
	vector<vector<vector<double>>> I_squared_matrix = findOverlapsOnly(eigenVectors);
	// [electron state][hole state][hole type]
	
	vector<double> delta_elements_I_lh;
	vector<double> delta_elements_E_GAP_lh;
	vector<double> delta_elements_I_hh;
	vector<double> delta_elements_E_GAP_hh;
	
		double Nq_con = 4800; // cm-1 eV-1 // from Whitehead paper - determined experimentally and applicable to this material system
	double x_center = QW.getThickness()/2.0; // add refractive indeces
	double n_avg = QW.getAvgRefractive();
	
	//for(int i = 0; i < (int)delta_E_GAP_lh.size(); i++){
		//cout << "E_GAP_lh" << i << " " << delta_E_GAP_lh[i];
		//cout << " I_lh" << i << " " << delta_I_lh[i];
		//cout << "E_GAP_hh" << i << " " << delta_E_GAP_hh[i];
		//cout << " I_hh" << i << " " << delta_I_hh[i] << endl;
	//}

	//~ double delta_x = abs(x_range[1]-x_range[0])/(double)(ND-1);
	//~ double precision = delta_x*0.5;	
	double R_0 = Rydberg(QW.eff_mass(0,x_center), QW.eff_mass(1,x_center), QW.eff_mass(2,x_center), 0, n_avg); // Rydberg(double m_e, double mass_lh, double mass_hh, int particle, double n_avg)
	double R_1 = Rydberg(QW.eff_mass(0,x_center), QW.eff_mass(1,x_center), QW.eff_mass(2,x_center), 1, n_avg); // Rydberg(double m_e, double mass_lh, double mass_hh, int particle, double n_avg)
	//~ //cout << R_0<<endl;
	//~ //cout << R_1<<endl;
	
	double q_ex_0 = R_0*12*Nq_con;
	double q_ex_1 = R_1*12*Nq_con;
	
	std::vector<Point> xyPairs;
	
	int k = 0;
	for(int j = 0; j<30; j++){
		for(int i = 0; i<30; i++){
			if(i==0 && j ==0) xyPairs.push_back({E_GAP_WL[i][j][k],I_squared_matrix[i][j][k]*0.5});
			else xyPairs.push_back({E_GAP_WL[i][j][k],I_squared_matrix[i][j][k]});
		}
	}
	k = 1;
	for(int j = 0; j<30; j++){
		for(int i = 0; i<30; i++){
			if(i==0 && j ==0) xyPairs.push_back({E_GAP_WL[i][j][k],I_squared_matrix[i][j][k]*1.5});
			else xyPairs.push_back({E_GAP_WL[i][j][k],I_squared_matrix[i][j][k]});
		}
	}
	

	
	x_limits = {1200, 1700};
    // Specified range for x-axis
    double xRangeStart = 1300;
    double xRangeEnd = 1700;
    // Number of points in x-axis vector
    int numPoints = ND;

    std::vector<double> xVector, yVector;
    generateXYVector(xyPairs, xRangeStart, xRangeEnd, numPoints, xVector, yVector);

    vector<vector<double>> out(ND, vector<double>(2));
	for(int i = 0; i < ND; i++){
		out[i][0] = xVector[i];
		out[i][1] = yVector[i];
	}
    
    
	//cout << E_GAP[0][0].size()<<"<size"<<endl;
	cout << "Calculating overlaps.\n";
	
	//vector<double> delta_E_GAP_lh, vector<double> delta_I_lh, vector<double> delta_E_GAP_hh, vector<double> delta_I_hh
	
	// vector<vector<double>> out = createDeltaSpace_two(delta_elements_E_GAP_lh, delta_elements_I_lh, delta_elements_E_GAP_hh, delta_elements_I_hh,  QW, x_limits);

		//~ double x_array_2[ND], y_array_2[ND];
		//~ for (int i = 0; i < ND; i++) {
			//~ y_array_2[i] = out[i][1];
			//~ x_array_2[i] = out[i][0];
		//~ }
		//~ gnuplot_two_functions ("Numerical Eigenfunctions", "linespoints", "x", "Eigenfunction |E_1> (x)",
			//~ x_array_2, y_array_2, ND, "Ground", x_array_2, y_array_2, ND, "5th excited state");
			
			
	//~ vector<vector<double>> absorption_data = absorption(E_GAP_WL, I_squared_matrix, x_limits); // , {900, 1300}
	
	// Write I_squared_matrix to a file
	//~ ofstream I_squared_file("I_squared_matrix_"+to_string(j)+".txt");
	//~ if(I_squared_file.is_open()){
		//~ for(const auto& row : I_squared_matrix){
			//~ for(const auto& col : row){
				//~ for(double val : col){
					//~ I_squared_file << val << " ";
				//~ }
				//~ I_squared_file << "\n";
			//~ }
			//~ I_squared_file << "\n";
		//~ }
		//~ I_squared_file.close();
	//~ } else {
		//~ cerr << "Unable to open I_squared_matrix file for writing.\n";
	//~ }

	//~ // Write E_GAP_WL to a file
	//~ ofstream E_GAP_WL_file("E_GAP_WL_"+to_string(j)+".txt");
	//~ if(E_GAP_WL_file.is_open()){
		//~ for(const auto& row : E_GAP_WL){
			//~ for(const auto& col : row){
				//~ for(double val : col){
					//~ E_GAP_WL_file << val << " ";
				//~ }
				//~ E_GAP_WL_file << "\n";
			//~ }
			//~ E_GAP_WL_file << "\n";
		//~ }
		//~ E_GAP_WL_file.close();
	//~ } else {
		//~ cerr << "Unable to open E_GAP_WL file for writing.\n";
	//~ }
	
	//~ // Write absorption to a file
	ofstream absorption_data_file("absorption_"+to_string(j)+".txt");
	if (absorption_data_file.is_open()) {
		for (const auto& row : out) {
			for (double val : row) {
				absorption_data_file << val << " ";
			}
			absorption_data_file << "\n";
		}
		absorption_data_file.close();
		//std::cout << "Matrix written to file: " << absorption_data_file << std::endl;
	} else {
		std::cerr << "Unable to open file: " << "absorption_"<<to_string(j)<<".txt" << std::endl;
	}
}

class ThreadPool {
public:
    ThreadPool(size_t numThreads) : stop(false) {
        for (size_t i = 0; i < numThreads; ++i) {
            workers.emplace_back([this] {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(queueMutex);
                        condition.wait(lock, [this] { return stop || !tasks.empty(); });
                        if (stop && tasks.empty())
                            return;
                        task = std::move(tasks.front());
                        tasks.pop();
                    }
                    task();
                }
            });
        }
    }

    template<class F>
    void enqueue(F&& f) {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            tasks.emplace(std::forward<F>(f));
        }
        condition.notify_one();
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread& worker : workers)
            worker.join();
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queueMutex;
    std::condition_variable condition;
    bool stop;
};

int main(int argc, char **argv)
{
	SimulationParameters params = readSimulationParameters("simulation_parameters.txt");
    std::cout << "Intermixing Enabled: " << (params.intermixingEnabled ? "TRUE" : "FALSE") << std::endl;
    std::cout << "Target Bandgap Shift: " << params.targetBandgapShift << std::endl;
    std::cout << "Number of Electric Fields: " << params.numElectricFields << std::endl;
    std::cout << "Max Electric Field: " << params.maxElectricField << std::endl;

	double delta_field = params.maxElectricField/(double)params.numElectricFields; //in put in
	delta_field *= 0.0001; //conversion to V/Amstrong

	Heterostructure QW_init("materials.txt","input.txt", 0.0);
	ShiftBandgapQWI(QW_init, (double)params.targetBandgapShift);
	QW_init.display();
	//plot_potential(QW);

    const int numSimulations = params.numElectricFields;
    const int numThreads = 5; // Adjust as needed
    ThreadPool pool(numThreads);
    std::vector<std::thread> threads;

    // Create threads for each simulation
    for (int i = 0; i < numSimulations; ++i) { //(int)numSimulations/
		double electric_field = i * delta_field; // Assuming 'delta_field' is defined
		if(i==0) electric_field = 0.00005;
		//cout << electric_field << endl;
        pool.enqueue(std::bind(script, i, QW_init, electric_field));
    }

    return 0;
}
