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
double h_c = 1240; // eV nm
double hbar_sqaured_2m = 3.81; // [eV Amstrong squared]

const int number_steps = 500;
const int max_matrix = 500;
const int ND = 1024;
const double reduce_matrix = 4.0; //cut down relatively the amount of computations - 4.0 considers only the first octant on matrices - data outside not required - should be checked for wells/materials


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
void eigenvalues_and_eigenvectors_tridiagonal_matrix (int nn)
{
	int n, m;
	for (n=0; n<nn; n++)
	{
		for (m=0; m<nn; m++)
		{
			if (n==m)
				trimatrix_result [n][m] = 1.0;
			else
				trimatrix_result [n][m] = 0.0;
		}
	}
	for (n=0; n<nn; n++)
		pointer_matrix [n] = &(trimatrix_result[n][0]);
	tqli (trimatrix_diag, trimatrix_subdiag, nn, pointer_matrix);
	for (n=0; n<nn; n++)
	{
		trimatrix_eigenvalue [n] = trimatrix_diag[n];
		for (m=0; m<nn; m++)
			trimatrix_eigenvector [n][m] = trimatrix_result[m][n];
	}
}

void sortTrimatrixEigenvaluesAndEigenvectors(double trimatrix_eigenvalue[], double trimatrix_eigenvector[][max_matrix], int nn) {
    // Create a vector of pairs to store eigenvalue-eigenvector pairs
    std::vector<std::pair<double, std::vector<double>>> eigenPairs(nn);
    
    // Populate the vector of pairs
    for (int i = 0; i < nn; ++i) {
        eigenPairs[i].first = trimatrix_eigenvalue[i];
        eigenPairs[i].second = std::vector<double>(trimatrix_eigenvector[i], trimatrix_eigenvector[i] + nn);
    }
    
    // Sort the vector of pairs based on eigenvalues
    std::sort(eigenPairs.begin(), eigenPairs.end(), [](const std::pair<double, std::vector<double>>& a, const std::pair<double, std::vector<double>>& b) {
        return a.first < b.first;
    });
    
    // Update the trimatrix_eigenvalue and trimatrix_eigenvector arrays with sorted values
    for (int i = 0; i < nn; ++i) {
        trimatrix_eigenvalue[i] = eigenPairs[i].first;
        std::copy(eigenPairs[i].second.begin(), eigenPairs[i].second.end(), trimatrix_eigenvector[i]);
    }
}

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
				return this->potential_[particle][x_int] - electron_charge*(this->electric_field)*(x-this->getThickness()); // not as efficient as defining a strict function globally but gives ALOT of freedom (only adds +(number of layers)*3 computations at most  
			}
			else if (particle == 0){             // else it is electron
				return this->potential_[particle][x_int] + electron_charge*(this->electric_field)*x;
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
    double EF_offset = electron_charge*QW.getElectricField()*QW.getThickness();
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
		E_REL = energy + BG + max_CBO;
    } else { // i.e. p==1 or p==2
		double max_VBO = VBO(QW.getLayers()[0].getMaterial(), QW.getLayers()[1].getMaterial());
		double pos = VBO(QW.getLayers()[0].getMaterial(), QW.getLayers()[1].getMaterial());
		for (int i = 1; i < (int)QW.getLayers().size()-1; ++i) {
			pos += VBO(QW.getLayers()[i].getMaterial(), QW.getLayers()[i+1].getMaterial());
			max_VBO = std::max(abs(max_VBO), abs(pos));
		}
		E_REL = -(max_VBO) + EF_offset - energy;	
	}
    return E_REL;
}

void findGroundState(Heterostructure& heterostructure,  int p) {  
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
		
		//~ //cout << "test0 ";
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
		std::string str = "";
		str += std::to_string(p);
		str += ":";
		str += std::to_string(relative_energy(w, p, heterostructure));
		str += "\n";
		cout << str;
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

void ShiftBandgapQWI(Heterostructure& QW, double shift_amount){ 					// shift amount in nm ?
	if(abs(shift_amount) < 0.01) {
		cout << "no QWI specified";
		return;
	}
	QW.resetPotential();
	// solve QW with no Electric field
	//QW.getElectricField()
	vector<double> x(number_steps); 												//length element
	vector<vector<double>> energies(3, vector<double>(number_steps)); 				//particle, energy_level
	vector<vector<vector<double>>> eigenVectors(3, vector<vector<double>>(number_steps, vector<double>(number_steps)));
	//solve(QW, x, energies, eigenVectors);
	
	double initBG_valence;
	if (relative_energy(energies[1][0],1,QW) < relative_energy(energies[2][0],2,QW)) initBG_valence = relative_energy(energies[2][0],2,QW);
	else initBG_valence = relative_energy(energies[1][0],1,QW);
	double initBG = h_c / abs(initBG_valence-relative_energy(energies[0][0],0,QW)); // [nm]
	
	// Bisection Method inspired routine
	double prox = 0.1; 						// nm
	double sigma = QW.getThickness(); 	// nm - max possible Intermixing
	double sigma0 = 0;
	double sigma_prev = sigma, temp;
	double BG_valence;
	int max_steps = 100;
	int count = 0;
	
	cout << "Starting QWI routine."<<endl;
	cout << "Initial Bandgap " << initBG << " [nm]" << endl;
	for(int i = 0; i < max_steps; i++) {
		QW.resetPotential();
		QW.intermixPotential(sigma);
		//solve(QW, x, energies, eigenVectors);
		
		if (relative_energy(energies[1][0],1,QW) < relative_energy(energies[2][0],2,QW)) BG_valence = relative_energy(energies[2][0],2,QW);
		else BG_valence = h_c / relative_energy(energies[1][0],1,QW);
		double BG_diff = initBG - h_c / abs(BG_valence-relative_energy(energies[0][0],0,QW)); // [nm]
		
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
			cout << "New Bandgap " << h_c / abs(BG_valence-relative_energy(energies[0][0],0,QW)) << " [nm]"<<endl;
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
	Heterostructure QW_init("materials.txt","input.txt", 0.0000);
	//ShiftBandgapQWI(QW_init, (double)params.targetBandgapShift);
	//plot_potential(QW);
	
    //const int numSimulations = params.numElectricFields;
    const int numThreads = 4; // Adjust as needed
    ThreadPool pool(numThreads);
    std::vector<std::thread> threads;

    // Create threads for each simulation
    for (int i = 0; i < 3; ++i) { //(int)numSimulations/
        pool.enqueue(std::bind(findGroundState, QW_init, i));
    }

    return 0;
}
