#include <iostream>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TChain.h"
#include <string>
using namespace std;
fstream fin;
void calulate_cross_sec(){
//initialize vectors for counting events
vector<double> vect_target_inc_carbon{0,0,0,0};//indices: 400/550/650/800 AMeV
vector<double> vect_target_twim_carbon{0,0,0,0};
vector<double> vect_target_particles_mw3{0,0,0,0};
vector<double> vect_target_11c10c{0,0,0,0};
vector<double> vect_target_survived_twim{0,0,0,0};

vector<double> vect_empty_inc_carbon{0,0,0,0};
vector<double> vect_empty_twim_carbon{0,0,0,0};
vector<double> vect_empty_particles_mw3{0,0,0,0};
vector<double> vect_empty_11c10c{0,0,0,0};
vector<double> vect_empty_survived_twim{0,0,0,0};
string line, word;

//here I check the runnumber to map it to the right target/energy
//emtpy target,400amev
string empty_400amev[2] = {"0069_0001","0074_0001"};
//carbon target 2.198cm,400amev
string carbon_2198_400amev[1] = {"0183_0001"};
for (int i = 0; i < (sizeof(carbon_2198_400amev)/sizeof(carbon_2198_400amev[0])); i++){
	fin.open("/scratch5/ge37liw/macros_s444/analysis_p2p/400amev/output_carbon_2198_400amev_run_"+carbon_2198_400amev[i]+".csv",ios::in);
	getline(fin, line);
	getline(fin, line);
	stringstream s(line);
	getline(s, word, ',');
	vect_target_inc_carbon[0] += stod(word);
	getline(s, word, ',');
	vect_target_twim_carbon[0] += stod(word);
	getline(s, word, ',');
	vect_target_particles_mw3[0] += stod(word);
	getline(s, word, ',');
	vect_target_11c10c[0] += stod(word);
	getline(s, word, ',');
	vect_target_survived_twim[0] += stod(word);
	fin.close();
	
	}
for (int i = 0; i < (sizeof(empty_400amev)/sizeof(empty_400amev[0])); i++){
	fin.open("/scratch5/ge37liw/macros_s444/analysis_p2p/400amev/output_empty_400amev_run_"+empty_400amev[i]+".csv",ios::in);
	getline(fin, line);
	getline(fin, line);
	stringstream s(line);
	getline(s, word, ',');
	vect_empty_inc_carbon[0] += stod(word);
	getline(s, word, ',');
	vect_empty_twim_carbon[0] += stod(word);
	getline(s, word, ',');
	vect_empty_particles_mw3[0] += stod(word);
	getline(s, word, ',');
	vect_empty_11c10c[0] += stod(word);
	getline(s, word, ',');
	vect_empty_survived_twim[0] += stod(word);
	fin.close();
	}

//empty target,800 amev
string empty_800amev[2] = {"0164_0001","0173_0001"};
//carbon target 2.198cm,800amev
string carbon_2198_800amev[1] = {"0170_0001"};
for (int i = 0; i < (sizeof(carbon_2198_800amev)/sizeof(carbon_2198_800amev[0])); i++){

	fin.open("/scratch5/ge37liw/macros_s444/analysis_p2p/800amev/output_carbon_2198_800amev_run_"+carbon_2198_800amev[i]+".csv",ios::in);
	getline(fin, line);
	getline(fin, line);
	stringstream s(line);
	getline(s, word, ',');
	vect_target_inc_carbon[3] += stod(word);
	getline(s, word, ',');
	vect_target_twim_carbon[3] += stod(word);
	getline(s, word, ',');
	vect_target_particles_mw3[3] += stod(word);
	getline(s, word, ',');
	vect_target_11c10c[3] += stod(word);
	getline(s, word, ',');
	vect_target_survived_twim[3] += stod(word);
	fin.close();
	}
for (int i = 0; i < (sizeof(empty_800amev)/sizeof(empty_800amev[0])); i++){

	fin.open("/scratch5/ge37liw/macros_s444/analysis_p2p/800amev/output_empty_800amev_run_"+empty_800amev[i]+".csv",ios::in);
	getline(fin, line);
	getline(fin, line);
	stringstream s(line);
	getline(s, word, ',');
	vect_empty_inc_carbon[3] += stod(word);
	getline(s, word, ',');
	vect_empty_twim_carbon[3] += stod(word);
	getline(s, word, ',');
	vect_empty_particles_mw3[3] += stod(word);
	getline(s, word, ',');
	vect_empty_11c10c[3] += stod(word);
	getline(s, word, ',');
	vect_empty_survived_twim[3] += stod(word);
	fin.close();
	}


//empty target,650amev
string empty_650amev[2] = {"0124_0001","0122_0001"};
//carbon target 2.198cm,650amev
string carbon_2198_650amev[2] = {"0130_0001","0134_0001"};
for (int i = 0; i < (sizeof(carbon_2198_650amev)/sizeof(carbon_2198_650amev[0])); i++){

	fin.open("/scratch5/ge37liw/macros_s444/analysis_p2p/650amev/output_carbon_2198_650amev_run_"+carbon_2198_650amev[i]+".csv",ios::in);
	getline(fin, line);
	getline(fin, line);
	stringstream s(line);
	getline(s, word, ',');
	vect_target_inc_carbon[2] += stod(word);
	getline(s, word, ',');
	vect_target_twim_carbon[2] += stod(word);
	getline(s, word, ',');
	vect_target_particles_mw3[2] += stod(word);
	getline(s, word, ',');
	vect_target_11c10c[2] += stod(word);
	getline(s, word, ',');
	vect_target_survived_twim[2] += stod(word);
	fin.close();
	}
for (int i = 0; i < (sizeof(empty_650amev)/sizeof(empty_650amev[0])); i++){

	fin.open("/scratch5/ge37liw/macros_s444/analysis_p2p/650amev/output_empty_650amev_run_"+empty_650amev[i]+".csv",ios::in);
	getline(fin, line);
	getline(fin, line);
	stringstream s(line);
	getline(s, word, ',');
	vect_empty_inc_carbon[2] += stod(word);
	getline(s, word, ',');
	vect_empty_twim_carbon[2] += stod(word);
	getline(s, word, ',');
	vect_empty_particles_mw3[2] += stod(word);
	getline(s, word, ',');
	vect_empty_11c10c[2] += stod(word);
	getline(s, word, ',');
	vect_empty_survived_twim[2] += stod(word);
	fin.close();
		}

//empty target,550 amev
//string empty_550amev[3] = {"0121_0001","0121_0002","0117_0001"};
string empty_550amev[3] = {"0094_0001","0095_0001","0096_0001"};
//carbon target 2.198cm,550amev
string carbon_2198_550amev[1] = {"0103_0001"};
for (int i = 0; i < (sizeof(carbon_2198_550amev)/sizeof(carbon_2198_550amev[0])); i++){

	fin.open("/scratch5/ge37liw/macros_s444/analysis_p2p/550amev/output_carbon_2198_550amev_run_"+carbon_2198_550amev[i]+".csv",ios::in);
	getline(fin, line);
	getline(fin, line);
	stringstream s(line);
	getline(s, word, ',');
	vect_target_inc_carbon[1] += stod(word);
	getline(s, word, ',');
	vect_target_twim_carbon[1] += stod(word);
	getline(s, word, ',');
	vect_target_particles_mw3[1] += stod(word);
	getline(s, word, ',');
	vect_target_11c10c[1] += stod(word);
	getline(s, word, ',');
	vect_target_survived_twim[1] += stod(word);
	fin.close();
	}
for (int i = 0; i < (sizeof(empty_550amev)/sizeof(empty_550amev[0])); i++){

	fin.open("/scratch5/ge37liw/macros_s444/analysis_p2p/550amev/output_empty_550amev_run_"+empty_550amev[i]+".csv",ios::in);
	getline(fin, line);
	getline(fin, line);
	stringstream s(line);
	getline(s, word, ',');
	vect_empty_inc_carbon[1] += stod(word);
	getline(s, word, ',');
	vect_empty_twim_carbon[1] += stod(word);
	getline(s, word, ',');
	vect_empty_particles_mw3[1] += stod(word);
	getline(s, word, ',');
	vect_empty_11c10c[1] += stod(word);
	getline(s, word, ',');
	vect_empty_survived_twim[1] += stod(word);
	fin.close();
	}
//now final calculation of the cross sections.... yeii
const double prefactor = 1./(4.034752*(1./12.0107)*6.02214*pow(10,23));
//400AMeV
Double_t a_400 = vect_empty_twim_carbon[0]/vect_empty_inc_carbon[0];
Double_t bc_400 = vect_empty_particles_mw3[0]/vect_empty_twim_carbon[0];
Double_t clean_incoming_12c_400 = vect_target_inc_carbon[0];
Double_t n_reacted_400 = clean_incoming_12c_400 - vect_target_twim_carbon[0] - (1-a_400)*clean_incoming_12c_400 ;//+ 1.02*vect_target_11c10c[0]*(1./bc_400); //here I use a simple sumup of the reactions in mw23
Double_t cross_sec_400 = (n_reacted_400/clean_incoming_12c_400)*prefactor*pow(10,27);
cout << "cross section 400 AMeV:\t" << cross_sec_400 << "mbarn" <<  "   and a:\t" << a_400 << endl;
cout << "N_survived/N_in :\t" << vect_target_twim_carbon[0]/clean_incoming_12c_400 << endl;

//550 AMeV

Double_t a_550 = vect_empty_twim_carbon[1]/vect_empty_inc_carbon[1];
Double_t bc_550 = vect_empty_particles_mw3[1]/vect_empty_twim_carbon[1];
Double_t clean_incoming_12c_550 = vect_target_inc_carbon[1];
Double_t n_reacted_550 = clean_incoming_12c_550 - vect_target_twim_carbon[1] - (1-a_550)*clean_incoming_12c_550 ;//+ 1.02*vect_target_11c10c[1]*(1./bc_550);
Double_t cross_sec_550 = (n_reacted_550/clean_incoming_12c_550)*prefactor*pow(10,27);
cout << "cross section 550 AMeV:\t" << cross_sec_550 << "mbarn" <<  "   and a:\t" << a_550 << endl;
cout << "N_survived/N_in :\t" << vect_target_twim_carbon[1]/clean_incoming_12c_550 << endl;

//650 AMeV

Double_t a_650 = vect_empty_twim_carbon[2]/vect_empty_inc_carbon[2];
Double_t bc_650 = vect_empty_particles_mw3[2]/vect_empty_twim_carbon[2];
Double_t clean_incoming_12c_650 = vect_target_inc_carbon[2];
Double_t n_reacted_650 = clean_incoming_12c_650 - vect_target_twim_carbon[2] - (1-a_650)*clean_incoming_12c_650 ;//+ 1.02*vect_target_11c10c[2]*(1./bc_650);
Double_t cross_sec_650 = (n_reacted_650/clean_incoming_12c_650)*prefactor*pow(10,27);
cout << "cross section 650 AMeV:\t" << cross_sec_650 << "mbarn" <<  "   and a:\t" << a_650 << endl;
cout << "N_survived/N_in :\t" << vect_target_twim_carbon[2]/clean_incoming_12c_650 << endl;
//cout << "now just by hand:\t" << (vect_target_twim_carbon[2]-1273)/clean_incoming_12c_650 << endl;

//800 AMeV

Double_t a_800 = vect_empty_twim_carbon[3]/vect_empty_inc_carbon[3];
Double_t bc_800 = vect_empty_particles_mw3[3]/vect_empty_twim_carbon[3];
Double_t clean_incoming_12c_800 = vect_target_inc_carbon[3];
Double_t n_reacted_800 = clean_incoming_12c_800 - vect_target_twim_carbon[3] - (1-a_800)*clean_incoming_12c_800 ;//+ 1.02*vect_target_11c10c[3]*(1./bc_800);
Double_t cross_sec_800 = (n_reacted_650/clean_incoming_12c_800)*prefactor*pow(10,27);
cout << "cross section 800 AMeV:\t" << cross_sec_800 << "mbarn"  <<  "   and a:\t" << a_800 << endl;
cout << "N_survived/N_in :\t" << vect_target_twim_carbon[3]/clean_incoming_12c_800 << endl;


//Further calculation with transmission method:

Double_t n_reacted_transmission_400 =(1 - (vect_target_twim_carbon[0]/vect_target_inc_carbon[0])/(vect_empty_twim_carbon[0]/vect_empty_inc_carbon[0]))*vect_target_inc_carbon[0];
Double_t n_reacted_transmission_550 =(1 - (vect_target_twim_carbon[1]/vect_target_inc_carbon[1])/(vect_empty_twim_carbon[1]/vect_empty_inc_carbon[1]))*vect_target_inc_carbon[1];
Double_t n_reacted_transmission_650 =(1 - (vect_target_twim_carbon[2]/vect_target_inc_carbon[2])/(vect_empty_twim_carbon[2]/vect_empty_inc_carbon[2]))*vect_target_inc_carbon[2];
Double_t n_reacted_transmission_800 =(1 - (vect_target_twim_carbon[3]/vect_target_inc_carbon[3])/(vect_empty_twim_carbon[3]/vect_empty_inc_carbon[3]))*vect_target_inc_carbon[3];

Double_t cross_sec_400_transmission = (n_reacted_transmission_400/clean_incoming_12c_400)*prefactor*pow(10,27);
Double_t cross_sec_550_transmission = (n_reacted_transmission_550/clean_incoming_12c_550)*prefactor*pow(10,27);
Double_t cross_sec_650_transmission = (n_reacted_transmission_650/clean_incoming_12c_650)*prefactor*pow(10,27);
Double_t cross_sec_800_transmission = (n_reacted_transmission_800/clean_incoming_12c_800)*prefactor*pow(10,27);

cout << "Cross section measured with transmission method:" << endl;
cout << "Cross section at 400 AMeV:\t" << cross_sec_400_transmission << endl;
cout << "Cross section at 550 AMeV:\t" << cross_sec_550_transmission << endl;
cout << "Cross section at 650 AMeV:\t" << cross_sec_650_transmission << endl;
cout << "Cross section at 800 AMeV:\t" << cross_sec_800_transmission << endl;

cout << "Cross section with transmission and also considering 10C/11C in some way... " << endl;
Double_t cross_sec_400_transmission_inclusive = ((n_reacted_transmission_400 + 1.02*vect_target_11c10c[0]*(1./bc_400))/clean_incoming_12c_400)*prefactor*pow(10,27);
Double_t cross_sec_550_transmission_inclusive = ((n_reacted_transmission_550 + 1.02*vect_target_11c10c[1]*(1./bc_550))/clean_incoming_12c_550)*prefactor*pow(10,27);
Double_t cross_sec_650_transmission_inclusive = ((n_reacted_transmission_650 + 1.02*vect_target_11c10c[2]*(1./bc_650))/clean_incoming_12c_650)*prefactor*pow(10,27);
Double_t cross_sec_800_transmission_inclusive = ((n_reacted_transmission_800 + + 1.02*vect_target_11c10c[3]*(1./bc_800))/clean_incoming_12c_800)*prefactor*pow(10,27);
cout << "INCLUSIVE Cross section at 400 AMeV:\t" << cross_sec_400_transmission_inclusive << endl;
cout << "INCLUSIVE Cross section at 550 AMeV:\t" << cross_sec_550_transmission_inclusive << endl;
cout << "INCLUSIVE Cross section at 650 AMeV:\t" << cross_sec_650_transmission_inclusive << endl;
cout << "INCLUSIVE Cross section at 800 AMeV:\t" << cross_sec_800_transmission_inclusive << endl;


cout << "<<<<<<<<<<<<<<<now correct transition method!! <<<<<<<<" << endl;
Double_t corr_400_exp = (vect_target_twim_carbon[0]/vect_target_inc_carbon[0])/(vect_empty_twim_carbon[0]/vect_empty_inc_carbon[0]);
Double_t corr_550_exp = (vect_target_twim_carbon[1]/vect_target_inc_carbon[1])/(vect_empty_twim_carbon[1]/vect_empty_inc_carbon[1]);
Double_t corr_650_exp = (vect_target_twim_carbon[2]/vect_target_inc_carbon[2])/(vect_empty_twim_carbon[2]/vect_empty_inc_carbon[2]);
Double_t corr_800_exp = (vect_target_twim_carbon[3]/vect_target_inc_carbon[3])/(vect_empty_twim_carbon[3]/vect_empty_inc_carbon[3]);

Double_t corr_cross_trans_400 = (-prefactor)*log(corr_400_exp)*pow(10,27);
Double_t corr_cross_trans_550 = (-prefactor)*log(corr_550_exp)*pow(10,27);
Double_t corr_cross_trans_650 = (-prefactor)*log(corr_650_exp)*pow(10,27);
Double_t corr_cross_trans_800 = (-prefactor)*log(corr_800_exp)*pow(10,27);
cout << "corrected cross section transition method, 400amev:\t" << corr_cross_trans_400 << endl;
cout << "corrected cross section transition method, 550amev:\t" << corr_cross_trans_550 << endl;
cout << "corrected cross section transition method, 650amev:\t" << corr_cross_trans_650 << endl;
cout << "corrected cross section transition method, 800amev:\t" << corr_cross_trans_800 << endl;

cout << "now even with some kind of correction:" << endl;
Double_t corr_factor_400 = ((vect_target_twim_carbon[0]-1.02*vect_target_11c10c[0]*(1./bc_400))/vect_target_inc_carbon[0])/(vect_empty_twim_carbon[0]/vect_empty_inc_carbon[0]);
Double_t corr_factor_550 = ((vect_target_twim_carbon[1]-1.02*vect_target_11c10c[1]*(1./bc_550))/vect_target_inc_carbon[1])/(vect_empty_twim_carbon[1]/vect_empty_inc_carbon[1]);
Double_t corr_factor_650 = ((vect_target_twim_carbon[2]-1.02*vect_target_11c10c[2]*(1./bc_650))/vect_target_inc_carbon[2])/(vect_empty_twim_carbon[2]/vect_empty_inc_carbon[2]);
Double_t corr_factor_800 = ((vect_target_twim_carbon[3]-1.02*vect_target_11c10c[3]*(1./bc_800))/vect_target_inc_carbon[3])/(vect_empty_twim_carbon[3]/vect_empty_inc_carbon[3]);

Double_t cross_sec_incl_trans_400 = (-prefactor)*log(corr_factor_400)*pow(10,27);
Double_t cross_sec_incl_trans_550 = (-prefactor)*log(corr_factor_550)*pow(10,27);
Double_t cross_sec_incl_trans_650 = (-prefactor)*log(corr_factor_650)*pow(10,27);
Double_t cross_sec_incl_trans_800 = (-prefactor)*log(corr_factor_800)*pow(10,27);
cout << "cross sec inclusive transision method, 400 amev:\t" << cross_sec_incl_trans_400 << endl;
cout << "cross sec inclusive transision method, 550 amev:\t" << cross_sec_incl_trans_550 << endl;
cout << "cross sec inclusive transision method, 650 amev:\t" << cross_sec_incl_trans_650 << endl;
cout << "cross sec inclusive transision method, 800 amev:\t" << cross_sec_incl_trans_800 << endl;


//surviving ranges---------------
cout << "-----------surviving rates--------------" << endl;
cout << "Efficiency,400 AMeV,empty:\t" << vect_empty_survived_twim[0]/vect_empty_inc_carbon[0] << endl;
cout << "Efficiency,550 AMeV,empty:\t" << vect_empty_survived_twim[1]/vect_empty_inc_carbon[1] << endl;
cout << "Efficiency,650 AMeV,empty:\t" << vect_empty_survived_twim[2]/vect_empty_inc_carbon[2] << endl;
cout << "Efficiency,800 AMeV,empty:\t" << vect_empty_survived_twim[3]/vect_empty_inc_carbon[3] << endl;

cout << "--------....and with target...--------------------" << endl;
cout << "Efficiency,400 AMeV,target:\t" << vect_target_survived_twim[0]/vect_target_inc_carbon[0] << endl;
cout << "Efficiency,550 AMeV,target:\t" << vect_target_survived_twim[1]/vect_target_inc_carbon[1] << endl;
cout << "Efficiency,650 AMeV,target:\t" << vect_target_survived_twim[2]/vect_target_inc_carbon[2] << endl;
cout << "Efficiency,800 AMeV,target:\t" << vect_target_survived_twim[3]/vect_target_inc_carbon[3] << endl;

}
