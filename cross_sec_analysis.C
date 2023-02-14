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

char fname[500];


void cross_sec_analysis(char const count_i[50]){
cout << "this is count_i " << count_i << endl;
cout << "and this is count_i[50] " << count_i[50] << endl;
stringstream ss;
ss << (count_i);
cout << "stringstream is:\t" << ss.str() << endl;
string input_run = string(count_i);
cout << "conversion to string:\t" << input_run << endl;

//here I check the runnumber to map it to the right target/energy
//emtpy target,400amev
string empty_400amev[2] = {"0069_0001","0074_0001"};
//carbon target 2.198cm,400amev
string carbon_2198_400amev[1] = {"0183_0001"};
//empty target,800 amev
string empty_800amev[2] = {"0164_0001","0173_0001"};
//carbon target 2.198cm,800amev
string carbon_2198_800amev[1] = {"0170_0001"};
//empty target,650amev
string empty_650amev[2] = {"0124_0001","0122_0001"};
//carbon target 2.198cm,650amev
string carbon_2198_650amev[2] = {"0130_0001","0134_0001"};
//empty target,550 amev
//string empty_550amev[3] = {"0121_0001","0121_0002","0117_0001"};
string empty_550amev[3] = {"0094_0001","0095_0001","0096_0001"};
//carbon target 2.198cm,550amev
string carbon_2198_550amev[1] = {"0103_0001"};

Double_t max_mw0x;
Double_t min_mw0x;
Double_t max_mw0y;
Double_t min_mw0y;
Double_t max_music_charge;
Double_t min_music_charge;
Double_t max_twim_charge;
Double_t min_twim_charge;
char f_out_name[500];
char f_par_file[500];
char f_cut_file[500];
char f_cut_file_all[500];
bool enable_cut = false;
bool enable_cut_all = false;
ofstream par_file;
ofstream csv_file;
char f_csv_file[500];
char canv_name[500];
char canv_beam[500];
char canv_mult_twim[500];


TFile* file_tcut_11c;
TCutG* cut_11c;
//for inclusive cut (11c +10c)
TFile* file_tcut_all;
TCutG* cut_all;

//400 amev runs--------------------------
for (int i = 0; i < (sizeof(carbon_2198_400amev)/sizeof(carbon_2198_400amev[0])); i++){
	if (input_run.compare(carbon_2198_400amev[i]) == 0){
		//insert here the cuts and names for the files
		sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_2198_400amev/stitched_and_unpacked_main%s.root",count_i);
		sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/carbon_target_2198_400amev_run_%s.root", count_i);
		sprintf(f_cut_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/cut_carbon_2198_400amev_%s.root",count_i);
		sprintf(f_cut_file_all,"/scratch5/ge37liw/macros_s444/analysis_p2p/cut_all_2198_400amev_%s.root",count_i);
		sprintf(f_csv_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/400amev/output_carbon_2198_400amev_run_%s.csv",count_i);
		csv_file.open(f_csv_file);
		sprintf(canv_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/400amev/carbon_target_charge_distr_%s.png",count_i);	
		sprintf(canv_beam,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/carbon_target_400_charge_distr_%s.png",count_i);	
		sprintf(canv_mult_twim,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/carbon_target_400_twim_mult_%s.png",count_i);	

		ifstream cut_file;
		cut_file.open(f_cut_file);
		if (cut_file){
			cout << "file to cut on 11C exists...." << endl;
			file_tcut_11c = TFile::Open(f_cut_file,"READ");
			cut_11c = (TCutG*)file_tcut_11c->Get("cut_11c");
			enable_cut = true;
			}
		ifstream cut_file_all;
		cut_file_all.open(f_cut_file_all);
		if (cut_file_all){
			cout << "file for inclusive 10C & 11C  exists..." << endl;
			file_tcut_11c = TFile::Open(f_cut_file_all);
			cut_all = (TCutG*)file_tcut_11c->Get("cut_all");
			enable_cut_all = true;
			}
		max_mw0x = -3.54069;
		min_mw0x = -5.54069;
		max_mw0y = -10.5;
		min_mw0y = -13;
		max_music_charge = 6.19862;
		min_music_charge = 5.79862;
		max_twim_charge = 7.7;
		min_twim_charge = 6.07;
		sprintf(f_par_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/400amev/output_carbon_2198_400amev_run_%s.txt",count_i);
		par_file.open(f_par_file);
		par_file << "****************CARBON RUN 21.98 mm, 400 AMeV, Runnr:  " << count_i << "   ******************" << endl;
		par_file << "CUTS:" << endl;
		par_file << "-----------Before target region:---------" << endl;
		par_file << "max_mw0x:\t" << max_mw0x << "\t" << "min_mw0x:\t" << min_mw0x << endl;
		par_file << "max_mw0y:\t" << max_mw0y << "\t" << "min_mw0y:\t" << min_mw0y << endl;
		par_file << "max_music_charge:\t" << max_music_charge << "\t" << "min_music_charge:\t" << min_music_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "-----------After target region:---------" << endl;
		par_file << "max_twim_charge:\t" << max_twim_charge << "\t" << "min_twim_charge:\t" << min_twim_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "------------OUTPUT------------------" << endl;
		}
	}

for (int i = 0; i < (sizeof(empty_400amev)/sizeof(empty_400amev[0])); i++){
	if (input_run.compare(empty_400amev[i]) == 0){
		//insert here the cuts and names for the files
		sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/empty_400amev/stitched_and_unpacked_main%s.root",count_i);
		sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/empty_target_2198_400amev_run_%s.root", count_i);
		sprintf(f_csv_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/400amev/output_empty_400amev_run_%s.csv",count_i);
        csv_file.open(f_csv_file);
		sprintf(canv_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/400amev/empty_target_charge_distr_%s.png",count_i);
		sprintf(canv_beam,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/empty_target_400_charge_distr_%s.png",count_i);	
		sprintf(canv_mult_twim,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/empty_target_400_twim_mult_%s.png",count_i);	
		max_mw0x = -2.15499; 
		min_mw0x = -4.15499;
		max_mw0y = -10.5;
		min_mw0y = -14.1386;
		max_music_charge = 6.15281;
		min_music_charge = 5.83281;
		max_twim_charge = 7.;
		min_twim_charge = 5.35;
		sprintf(f_par_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/400amev/output_empty_400amev_run_%s.txt",count_i);
		par_file.open(f_par_file);
		par_file << "****************EMPTY RUN, 400 AMeV, Runnr:  " << count_i << "   ******************" << endl;
		par_file << "CUTS:" << endl;
		par_file << "-----------Before target region:---------" << endl;
		par_file << "max_mw0x:\t" << max_mw0x << "\t" << "min_mw0x:\t" << min_mw0x << endl;
		par_file << "max_mw0y:\t" << max_mw0y << "\t" << "min_mw0y:\t" << min_mw0y << endl;
		par_file << "max_music_charge:\t" << max_music_charge << "\t" << "min_music_charge:\t" << min_music_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "-----------After target region:---------" << endl;
		par_file << "max_twim_charge:\t" << max_twim_charge << "\t" << "min_twim_charge:\t" << min_twim_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "------------OUTPUT------------------" << endl;
		}
	}
//end of 400 amev runs -------------------

//800 amev runs -------------------------

for (int i = 0; i < (sizeof(carbon_2198_800amev)/sizeof(carbon_2198_800amev[0])); i++){
	if (input_run.compare(carbon_2198_800amev[i]) == 0){
		//insert here the cuts and names for the files
		sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_2198_800amev/stitched_and_unpacked_main%s.root",count_i);
		sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/carbon_target_2198_800amev_run_%s.root", count_i);
		sprintf(f_cut_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/cut_carbon_2198_800amev_%s.root",count_i);
		sprintf(f_cut_file_all,"/scratch5/ge37liw/macros_s444/analysis_p2p/cut_all_2198_800amev_%s.root",count_i);
		sprintf(f_csv_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/800amev/output_carbon_2198_800amev_run_%s.csv",count_i);
        csv_file.open(f_csv_file);
		sprintf(canv_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/800amev/carbon_target_charge_distr_%s.png",count_i);
		sprintf(canv_beam,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/carbon_target_800_charge_distr_%s.png",count_i);	
		sprintf(canv_mult_twim,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/carbon_target_800_twim_mult_%s.png",count_i);	
		ifstream cut_file;
		cut_file.open(f_cut_file);
		if (cut_file){
			cout << "file to cut on 11C exists...." << endl;
			file_tcut_11c = TFile::Open(f_cut_file,"READ");
			cut_11c = (TCutG*)file_tcut_11c->Get("cut_11c");
			enable_cut = true;
			}
		ifstream cut_file_all;
		cut_file_all.open(f_cut_file_all);
		if (cut_file_all){
			cout << "file for inclusive 10C & 11C  exists..." << endl;
			file_tcut_11c = TFile::Open(f_cut_file_all);
			cut_all = (TCutG*)file_tcut_11c->Get("cut_all");
			enable_cut_all = true;
			}
		max_mw0x = -1.;
		min_mw0x = -3.5;
		max_mw0y = -9.6;
		min_mw0y = -15.;
		max_music_charge = 5.96448;
		min_music_charge = 5.60448;
		max_twim_charge = 7;
		min_twim_charge = 5.09;
		sprintf(f_par_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/800amev/output_carbon_2198_800amev_run_%s.txt",count_i);
		par_file.open(f_par_file);
		par_file << "****************CARBON RUN 21.98 mm, 800 AMeV, Runnr:  " << count_i << "   ******************" << endl;
		par_file << "CUTS:" << endl;
		par_file << "-----------Before target region:---------" << endl;
		par_file << "max_mw0x:\t" << max_mw0x << "\t" << "min_mw0x:\t" << min_mw0x << endl;
		par_file << "max_mw0y:\t" << max_mw0y << "\t" << "min_mw0y:\t" << min_mw0y << endl;
		par_file << "max_music_charge:\t" << max_music_charge << "\t" << "min_music_charge:\t" << min_music_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "-----------After target region:---------" << endl;
		par_file << "max_twim_charge:\t" << max_twim_charge << "\t" << "min_twim_charge:\t" << min_twim_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "------------OUTPUT------------------" << endl;
		}
	}

for (int i = 0; i < (sizeof(empty_800amev)/sizeof(empty_800amev[0])); i++){
	if (input_run.compare(empty_800amev[i]) == 0){
		//insert here the cuts and names for the files
		sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/empty_800amev/stitched_and_unpacked_main%s.root",count_i);
		sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/empty_target_2198_800amev_run_%s.root", count_i);
		sprintf(f_csv_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/800amev/output_empty_800amev_run_%s.csv",count_i);
        csv_file.open(f_csv_file);
		sprintf(canv_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/800amev/empty_target_charge_distr_%s.png",count_i);
		sprintf(canv_beam,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/empty_target_800_charge_distr_%s.png",count_i);	
		sprintf(canv_mult_twim,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/empty_target_800_twim_mult_%s.png",count_i);	
		max_mw0x = -1.; 
		min_mw0x = -3.3;
		max_mw0y = -9.6;
		min_mw0y = -15.;
		max_music_charge = 5.96817;
		min_music_charge = 5.60817;
		max_twim_charge = 6.5;
		min_twim_charge = 5.02;
		sprintf(f_par_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/800amev/output_empty_800amev_run_%s.txt",count_i);
		par_file.open(f_par_file);
		par_file << "****************EMPTY RUN, 800 AMeV, Runnr:  " << count_i << "   ******************" << endl;
		par_file << "CUTS:" << endl;
		par_file << "-----------Before target region:---------" << endl;
		par_file << "max_mw0x:\t" << max_mw0x << "\t" << "min_mw0x:\t" << min_mw0x << endl;
		par_file << "max_mw0y:\t" << max_mw0y << "\t" << "min_mw0y:\t" << min_mw0y << endl;
		par_file << "max_music_charge:\t" << max_music_charge << "\t" << "min_music_charge:\t" << min_music_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "-----------After target region:---------" << endl;
		par_file << "max_twim_charge:\t" << max_twim_charge << "\t" << "min_twim_charge:\t" << min_twim_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "------------OUTPUT------------------" << endl;
		}
	}
//end of 800 amev runs ------------------

//650 amev runs -------------------------

for (int i = 0; i < (sizeof(carbon_2198_650amev)/sizeof(carbon_2198_650amev[0])); i++){
	if (input_run.compare(carbon_2198_650amev[i]) == 0){
		//insert here the cuts and names for the files
		sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_2198_650amev/stitched_and_unpacked_main%s.root",count_i);
		sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/carbon_target_2198_650amev_run_%s.root", count_i);
		sprintf(f_cut_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/cut_carbon_2198_650amev_%s.root",count_i);
		sprintf(f_cut_file_all,"/scratch5/ge37liw/macros_s444/analysis_p2p/cut_all_2198_650amev_%s.root",count_i);
		sprintf(f_csv_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/650amev/output_carbon_2198_650amev_run_%s.csv",count_i);
        csv_file.open(f_csv_file);
		sprintf(canv_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/650amev/carbon_target_charge_distr_%s.png",count_i);
		sprintf(canv_beam,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/carbon_target_650_charge_distr_%s.png",count_i);	
		sprintf(canv_mult_twim,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/carbon_target_650_twim_mult_%s.png",count_i);	
		ifstream cut_file;
		cut_file.open(f_cut_file);
		if (cut_file){
			cout << "file to cut on 11C exists...." << endl;
			file_tcut_11c = TFile::Open(f_cut_file,"READ");
			cut_11c = (TCutG*)file_tcut_11c->Get("cut_11c");
			enable_cut = true;
			}
		ifstream cut_file_all;
		cut_file_all.open(f_cut_file_all);
		if (cut_file_all){
			cout << "file for inclusive 10C & 11C  exists..." << endl;
			file_tcut_11c = TFile::Open(f_cut_file_all);
			cut_all = (TCutG*)file_tcut_11c->Get("cut_all");
			enable_cut_all = true;
			}
		max_mw0x = -2.3;
		min_mw0x = -5;
		max_mw0y = -9.7;
		min_mw0y = -15.;
		max_music_charge = 5.73287;
		min_music_charge = 5.43286;
		max_twim_charge = 6.6;
		min_twim_charge = 4.82;
		sprintf(f_par_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/650amev/output_carbon_2198_650amev_run_%s.txt",count_i);
		par_file.open(f_par_file);
		par_file << "****************CARBON RUN 21.98 mm, 650 AMeV, Runnr:  " << count_i << "   ******************" << endl;
		par_file << "CUTS:" << endl;
		par_file << "-----------Before target region:---------" << endl;
		par_file << "max_mw0x:\t" << max_mw0x << "\t" << "min_mw0x:\t" << min_mw0x << endl;
		par_file << "max_mw0y:\t" << max_mw0y << "\t" << "min_mw0y:\t" << min_mw0y << endl;
		par_file << "max_music_charge:\t" << max_music_charge << "\t" << "min_music_charge:\t" << min_music_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "-----------After target region:---------" << endl;
		par_file << "max_twim_charge:\t" << max_twim_charge << "\t" << "min_twim_charge:\t" << min_twim_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "------------OUTPUT------------------" << endl;
		}
	}

for (int i = 0; i < (sizeof(empty_650amev)/sizeof(empty_650amev[0])); i++){
	if (input_run.compare(empty_650amev[i]) == 0){
		//insert here the cuts and names for the files
		sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/empty_650amev/stitched_and_unpacked_main%s.root",count_i);
		sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/empty_target_2198_650amev_run_%s.root", count_i);
		sprintf(f_csv_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/650amev/output_empty_650amev_run_%s.csv",count_i);
        csv_file.open(f_csv_file);
		sprintf(canv_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/650amev/empty_target_charge_distr_%s.png",count_i);
		sprintf(canv_beam,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/empty_target_650_charge_distr_%s.png",count_i);	
		sprintf(canv_mult_twim,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/empty_target_650_twim_mult_%s.png",count_i);	
		max_mw0x = -2.3; 
		min_mw0x = -5;
		max_mw0y = -9.7;
		min_mw0y = -15.;
		max_music_charge = 5.73341;
		min_music_charge = 5.43340;
		max_twim_charge = 6.6;
		min_twim_charge = 4.74;
		sprintf(f_par_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/650amev/output_empty_650amev_run_%s.txt",count_i);
		par_file.open(f_par_file);
		par_file << "****************EMPTY RUN, 650 AMeV, Runnr:  " << count_i << "   ******************" << endl;
		par_file << "CUTS:" << endl;
		par_file << "-----------Before target region:---------" << endl;
		par_file << "max_mw0x:\t" << max_mw0x << "\t" << "min_mw0x:\t" << min_mw0x << endl;
		par_file << "max_mw0y:\t" << max_mw0y << "\t" << "min_mw0y:\t" << min_mw0y << endl;
		par_file << "max_music_charge:\t" << max_music_charge << "\t" << "min_music_charge:\t" << min_music_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "-----------After target region:---------" << endl;
		par_file << "max_twim_charge:\t" << max_twim_charge << "\t" << "min_twim_charge:\t" << min_twim_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "------------OUTPUT------------------" << endl;
		}
	}
//end of 650 amev runs ------------------

//550 amev runs -------------------------

for (int i = 0; i < (sizeof(carbon_2198_550amev)/sizeof(carbon_2198_550amev[0])); i++){
	if (input_run.compare(carbon_2198_550amev[i]) == 0){
		//insert here the cuts and names for the files
		sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_2198_550amev/stitched_and_unpacked_main%s.root",count_i);
		sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/carbon_target_2198_550amev_run_%s.root", count_i);
		sprintf(f_cut_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/cut_carbon_2198_550amev_%s.root",count_i);
		sprintf(f_cut_file_all,"/scratch5/ge37liw/macros_s444/analysis_p2p/cut_all_2198_550amev_%s.root",count_i);
		sprintf(f_csv_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/550amev/output_carbon_2198_550amev_run_%s.csv",count_i);
        csv_file.open(f_csv_file);
		sprintf(canv_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/550amev/carbon_target_charge_distr_%s.png",count_i);
		sprintf(canv_beam,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/carbon_target_550_charge_distr_%s.png",count_i);	
		sprintf(canv_mult_twim,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/carbon_target_550_twim_mult_%s.png",count_i);	
		ifstream cut_file;
		cut_file.open(f_cut_file);
		if (cut_file){
			cout << "file to cut on 11C exists...." << endl;
			file_tcut_11c = TFile::Open(f_cut_file,"READ");
			cut_11c = (TCutG*)file_tcut_11c->Get("cut_11c");
			enable_cut = true;
			}
		ifstream cut_file_all;
		cut_file_all.open(f_cut_file_all);
		if (cut_file_all){
			cout << "file for inclusive 10C & 11C  exists..." << endl;
			file_tcut_11c = TFile::Open(f_cut_file_all);
			cut_all = (TCutG*)file_tcut_11c->Get("cut_all");
			enable_cut_all = true;
			}
		max_mw0x = -2.5;
		min_mw0x = -2.98;
		max_mw0y = -9.59;
		min_mw0y = -14.8;
		max_music_charge = 5.86289;
		min_music_charge = 5.54289;
		max_twim_charge = 6.7;
		min_twim_charge = 5.02;
		sprintf(f_par_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/550amev/output_carbon_2198_550amev_run_%s.txt",count_i);
		par_file.open(f_par_file);
		par_file << "****************CARBON RUN 21.98 mm, 550 AMeV, Runnr:  " << count_i << "   ******************" << endl;
		par_file << "CUTS:" << endl;
		par_file << "-----------Before target region:---------" << endl;
		par_file << "max_mw0x:\t" << max_mw0x << "\t" << "min_mw0x:\t" << min_mw0x << endl;
		par_file << "max_mw0y:\t" << max_mw0y << "\t" << "min_mw0y:\t" << min_mw0y << endl;
		par_file << "max_music_charge:\t" << max_music_charge << "\t" << "min_music_charge:\t" << min_music_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "-----------After target region:---------" << endl;
		par_file << "max_twim_charge:\t" << max_twim_charge << "\t" << "min_twim_charge:\t" << min_twim_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "------------OUTPUT------------------" << endl;
		}
	}

for (int i = 0; i < (sizeof(empty_550amev)/sizeof(empty_550amev[0])); i++){
	if (input_run.compare(empty_550amev[i]) == 0){
		//insert here the cuts and names for the files
		sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/empty_550amev/stitched_and_unpacked_main%s.root",count_i);
		sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/empty_target_2198_550amev_run_%s.root", count_i);
		sprintf(f_csv_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/550amev/output_empty_550amev_run_%s.csv",count_i);
        csv_file.open(f_csv_file);
		sprintf(canv_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/550amev/empty_target_charge_distr_%s.png",count_i);
		sprintf(canv_beam,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/empty_target_550_charge_distr_%s.png",count_i);	
		sprintf(canv_mult_twim,"/scratch5/ge37liw/macros_s444/analysis_p2p/intensities/empty_target_550_twim_mult_%s.png",count_i);	
		max_mw0x = -0.5; 
		min_mw0x = -2.98;
		max_mw0y = -9.59;
		min_mw0y = -14.8;
		max_music_charge = 5.85558;
		min_music_charge = 5.53557;
		max_twim_charge = 6.6;
		min_twim_charge = 4.92;
		sprintf(f_par_file,"/scratch5/ge37liw/macros_s444/analysis_p2p/550amev/output_empty_550amev_run_%s.txt",count_i);
		par_file.open(f_par_file);
		par_file << "****************EMPTY RUN, 550 AMeV, Runnr:  " << count_i << "   ******************" << endl;
		par_file << "CUTS:" << endl;
		par_file << "-----------Before target region:---------" << endl;
		par_file << "max_mw0x:\t" << max_mw0x << "\t" << "min_mw0x:\t" << min_mw0x << endl;
		par_file << "max_mw0y:\t" << max_mw0y << "\t" << "min_mw0y:\t" << min_mw0y << endl;
		par_file << "max_music_charge:\t" << max_music_charge << "\t" << "min_music_charge:\t" << min_music_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "-----------After target region:---------" << endl;
		par_file << "max_twim_charge:\t" << max_twim_charge << "\t" << "min_twim_charge:\t" << min_twim_charge << endl;
		par_file << "-----------------------------------------" <<endl;
		par_file << endl;
		par_file << "------------OUTPUT------------------" << endl;
		}
	}
//end of 550 amev runs ------------------

//MY HISTOS
char hist_name[500];

TH2F* h2_xmw2_vs_xmw3_survived12c;
sprintf(hist_name, "X MW2 vs MW3 in mm survived  12c");
h2_xmw2_vs_xmw3_survived12c = new TH2F(hist_name,hist_name,800,-400,400,1000,-500,500);
h2_xmw2_vs_xmw3_survived12c->GetXaxis()->SetTitle("x position MW2");
h2_xmw2_vs_xmw3_survived12c->GetYaxis()->SetTitle("x position MW3");
h2_xmw2_vs_xmw3_survived12c->GetXaxis()->CenterTitle(true);
h2_xmw2_vs_xmw3_survived12c->GetYaxis()->CenterTitle(true);
h2_xmw2_vs_xmw3_survived12c->GetYaxis()->SetLabelSize(0.045);
h2_xmw2_vs_xmw3_survived12c->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_xmw2_vs_xmw3_mw0cut;
sprintf(hist_name, "MW2 vs MW3 in mm with cut on mw0_x and cut around mean r3bmusic +-0.2");
h2_xmw2_vs_xmw3_mw0cut = new TH2F(hist_name,hist_name,800,-400,400,1000,-500,500);
h2_xmw2_vs_xmw3_mw0cut->GetXaxis()->SetTitle("x position MW2");
h2_xmw2_vs_xmw3_mw0cut->GetYaxis()->SetTitle("x position MW3");
h2_xmw2_vs_xmw3_mw0cut->GetXaxis()->CenterTitle(true);
h2_xmw2_vs_xmw3_mw0cut->GetYaxis()->CenterTitle(true);
h2_xmw2_vs_xmw3_mw0cut->GetYaxis()->SetLabelSize(0.045);
h2_xmw2_vs_xmw3_mw0cut->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_xmw2_vs_xmw3_11c;
sprintf(hist_name, "MW2 vs MW3 in mm for all survived 11C's ");
h2_xmw2_vs_xmw3_11c = new TH2F(hist_name,hist_name,800,-400,400,1000,-500,500);
h2_xmw2_vs_xmw3_11c->GetXaxis()->SetTitle("x position MW2");
h2_xmw2_vs_xmw3_11c->GetYaxis()->SetTitle("x position MW3");
h2_xmw2_vs_xmw3_11c->GetXaxis()->CenterTitle(true);
h2_xmw2_vs_xmw3_11c->GetYaxis()->CenterTitle(true);
h2_xmw2_vs_xmw3_11c->GetYaxis()->SetLabelSize(0.045);
h2_xmw2_vs_xmw3_11c->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_xmw2_vs_xmw3_inclusive;
sprintf(hist_name, "MW2 vs MW3 in mm for all survived 11C's and 10C's (inclusive) ");
h2_xmw2_vs_xmw3_inclusive = new TH2F(hist_name,hist_name,800,-400,400,1000,-500,500);
h2_xmw2_vs_xmw3_inclusive->GetXaxis()->SetTitle("x position MW2");
h2_xmw2_vs_xmw3_inclusive->GetYaxis()->SetTitle("x position MW3");
h2_xmw2_vs_xmw3_inclusive->GetXaxis()->CenterTitle(true);
h2_xmw2_vs_xmw3_inclusive->GetYaxis()->CenterTitle(true);
h2_xmw2_vs_xmw3_inclusive->GetYaxis()->SetLabelSize(0.045);
h2_xmw2_vs_xmw3_inclusive->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_charge_val;
sprintf(hist_name, "Charge Value Twim(for condition entries_start ==2 && entries_mw0 == 1 && softwimhitdata[0] && musichitdata[0])");
h1_charge_val = new TH1F(hist_name,hist_name,600,2,8);
h1_charge_val->GetXaxis()->SetTitle("Charge Value (uncalibrated)");
h1_charge_val->GetYaxis()->SetTitle("Counts");
h1_charge_val->GetXaxis()->CenterTitle(true);
h1_charge_val->GetYaxis()->CenterTitle(true);
h1_charge_val->GetYaxis()->SetLabelSize(0.045);
h1_charge_val->GetYaxis()->SetTitleSize(0.045);


TH1F* h1_charge_val_music;
sprintf(hist_name, "Charge Value r3bmusic (for condition entries_start ==2 && entries_mw0 == 1 && softwimhitdata[0] && musichitdata[0])");
h1_charge_val_music = new TH1F(hist_name,hist_name,600,2,8);
h1_charge_val_music->GetXaxis()->SetTitle("Charge Value (uncalibrated)");
h1_charge_val_music->GetYaxis()->SetTitle("Counts");
h1_charge_val_music->GetXaxis()->CenterTitle(true);
h1_charge_val_music->GetYaxis()->CenterTitle(true);
h1_charge_val_music->GetYaxis()->SetLabelSize(0.045);
h1_charge_val_music->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw0_x_pos;
sprintf(hist_name, "x position on mw0 ");
h1_mw0_x_pos = new TH1F(hist_name,hist_name,600,-100,100);
h1_mw0_x_pos->GetXaxis()->SetTitle("x value [mm]");
h1_mw0_x_pos->GetYaxis()->SetTitle("Counts");
h1_mw0_x_pos->GetXaxis()->CenterTitle(true);
h1_mw0_x_pos->GetYaxis()->CenterTitle(true);
h1_mw0_x_pos->GetYaxis()->SetLabelSize(0.045);
h1_mw0_x_pos->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw0_y_pos;
sprintf(hist_name, "y position on mw0 ");
h1_mw0_y_pos = new TH1F(hist_name,hist_name,600,-100,100);
h1_mw0_y_pos->GetXaxis()->SetTitle("x value [mm]");
h1_mw0_y_pos->GetYaxis()->SetTitle("Counts");
h1_mw0_y_pos->GetXaxis()->CenterTitle(true);
h1_mw0_y_pos->GetYaxis()->CenterTitle(true);
h1_mw0_y_pos->GetYaxis()->SetLabelSize(0.045);
h1_mw0_y_pos->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_clean_incoming_12c;
sprintf(hist_name, "clean incoming 12C ");
h1_clean_incoming_12c = new TH1F(hist_name,hist_name,600,-100,100);
h1_clean_incoming_12c->GetXaxis()->SetTitle("clean incoming 12c");
h1_clean_incoming_12c->GetYaxis()->SetTitle("Counts");
h1_clean_incoming_12c->GetXaxis()->CenterTitle(true);
h1_clean_incoming_12c->GetYaxis()->CenterTitle(true);
h1_clean_incoming_12c->GetYaxis()->SetLabelSize(0.045);
h1_clean_incoming_12c->GetYaxis()->SetTitleSize(0.045);


TH2F* h2_xmw0_vs_xmw2;
sprintf(hist_name, "MW0 vs MW2 in mm for all incoming particles, to select good ones...X ");
h2_xmw0_vs_xmw2 = new TH2F(hist_name,hist_name,600,-100,100,800,-400,400);
h2_xmw0_vs_xmw2->GetXaxis()->SetTitle("x position MW0");
h2_xmw0_vs_xmw2->GetYaxis()->SetTitle("x position MW2");
h2_xmw0_vs_xmw2->GetXaxis()->CenterTitle(true);
h2_xmw0_vs_xmw2->GetYaxis()->CenterTitle(true);
h2_xmw0_vs_xmw2->GetYaxis()->SetLabelSize(0.045);
h2_xmw0_vs_xmw2->GetYaxis()->SetTitleSize(0.045);


TH2F* h2_xmw0_vs_ymw0;
sprintf(hist_name, "x vs y MW0 in mm for all incoming particles...");
h2_xmw0_vs_ymw0 = new TH2F(hist_name,hist_name,800,-100,100,800,-100,100);
h2_xmw0_vs_ymw0->GetXaxis()->SetTitle("x position MW0");
h2_xmw0_vs_ymw0->GetYaxis()->SetTitle("y position MW0");
h2_xmw0_vs_ymw0->GetXaxis()->CenterTitle(true);
h2_xmw0_vs_ymw0->GetYaxis()->CenterTitle(true);
h2_xmw0_vs_ymw0->GetYaxis()->SetLabelSize(0.045);
h2_xmw0_vs_ymw0->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_xmw0_vs_ymw0_no_twim;
sprintf(hist_name, "x vs y MW0 in mm for all incoming particles...,no twim");
h2_xmw0_vs_ymw0_no_twim = new TH2F(hist_name,hist_name,800,-100,100,800,-100,100);
h2_xmw0_vs_ymw0_no_twim->GetXaxis()->SetTitle("x position MW0");
h2_xmw0_vs_ymw0_no_twim->GetYaxis()->SetTitle("y position MW0");
h2_xmw0_vs_ymw0_no_twim->GetXaxis()->CenterTitle(true);
h2_xmw0_vs_ymw0_no_twim->GetYaxis()->CenterTitle(true);
h2_xmw0_vs_ymw0_no_twim->GetYaxis()->SetLabelSize(0.045);
h2_xmw0_vs_ymw0_no_twim->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_xmw0_vs_ymw0_no_twim_no_mw_cut;
sprintf(hist_name, "x vs y MW0 in mm for all incoming particles...,no twim, no mw cut at all...");
h2_xmw0_vs_ymw0_no_twim_no_mw_cut = new TH2F(hist_name,hist_name,800,-100,100,800,-100,100);
h2_xmw0_vs_ymw0_no_twim_no_mw_cut->GetXaxis()->SetTitle("x position MW0");
h2_xmw0_vs_ymw0_no_twim_no_mw_cut->GetYaxis()->SetTitle("y position MW0");
h2_xmw0_vs_ymw0_no_twim_no_mw_cut->GetXaxis()->CenterTitle(true);
h2_xmw0_vs_ymw0_no_twim_no_mw_cut->GetYaxis()->CenterTitle(true);
h2_xmw0_vs_ymw0_no_twim_no_mw_cut->GetYaxis()->SetLabelSize(0.045);
h2_xmw0_vs_ymw0_no_twim_no_mw_cut->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_xmw2_vs_ymw2;
sprintf(hist_name, "x vs y MW2 in mm for all incoming particles...");
h2_xmw2_vs_ymw2 = new TH2F(hist_name,hist_name,800,-100,100,800,-100,100);
h2_xmw2_vs_ymw2->GetXaxis()->SetTitle("x position MW2");
h2_xmw2_vs_ymw2->GetYaxis()->SetTitle("y position MW2");
h2_xmw2_vs_ymw2->GetXaxis()->CenterTitle(true);
h2_xmw2_vs_ymw2->GetYaxis()->CenterTitle(true);
h2_xmw2_vs_ymw2->GetYaxis()->SetLabelSize(0.045);
h2_xmw2_vs_ymw2->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_ymw0_vs_ymw2;
sprintf(hist_name, "MW0 vs MW2 in mm for all incoming particles, to select good ones...Y ");
h2_ymw0_vs_ymw2 = new TH2F(hist_name,hist_name,600,-100,100,800,-400,400);
h2_ymw0_vs_ymw2->GetXaxis()->SetTitle("y position MW0");
h2_ymw0_vs_ymw2->GetYaxis()->SetTitle("y position MW2");
h2_ymw0_vs_ymw2->GetXaxis()->CenterTitle(true);
h2_ymw0_vs_ymw2->GetYaxis()->CenterTitle(true);
h2_ymw0_vs_ymw2->GetYaxis()->SetLabelSize(0.045);
h2_ymw0_vs_ymw2->GetYaxis()->SetTitleSize(0.045);


TH1F* h1_mw0_radius;
sprintf(hist_name, "Radius of clean incoming particles in MW0 ");
h1_mw0_radius = new TH1F(hist_name,hist_name,1000,0,200);
h1_mw0_radius->GetXaxis()->SetTitle("radius mw0 in mm");
h1_mw0_radius->GetYaxis()->SetTitle("Counts");
h1_mw0_radius->GetXaxis()->CenterTitle(true);
h1_mw0_radius->GetYaxis()->CenterTitle(true);
h1_mw0_radius->GetYaxis()->SetLabelSize(0.045);
h1_mw0_radius->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw0_radius_no_twim;
sprintf(hist_name, "Radius of clean incoming particles in MW0 no hit in twim ");
h1_mw0_radius_no_twim = new TH1F(hist_name,hist_name,1000,0,200);
h1_mw0_radius_no_twim->GetXaxis()->SetTitle("radius mw0 in mm");
h1_mw0_radius_no_twim->GetYaxis()->SetTitle("Counts");
h1_mw0_radius_no_twim->GetXaxis()->CenterTitle(true);
h1_mw0_radius_no_twim->GetYaxis()->CenterTitle(true);
h1_mw0_radius_no_twim->GetYaxis()->SetLabelSize(0.045);
h1_mw0_radius_no_twim->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw2_radius;
sprintf(hist_name, "Radius of clean incoming particles in MW2 ");
h1_mw2_radius = new TH1F(hist_name,hist_name,1200,0,300);
h1_mw2_radius->GetXaxis()->SetTitle("radius mw2 in mm");
h1_mw2_radius->GetYaxis()->SetTitle("Counts");
h1_mw2_radius->GetXaxis()->CenterTitle(true);
h1_mw2_radius->GetYaxis()->CenterTitle(true);
h1_mw2_radius->GetYaxis()->SetLabelSize(0.045);
h1_mw2_radius->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_radius_mw0_mw2;
sprintf(hist_name, "Radius MW0 vs Radius MW2");
h2_radius_mw0_mw2 = new TH2F(hist_name,hist_name,500,0,200,600,0,300);
h2_radius_mw0_mw2->GetXaxis()->SetTitle("Radius MW0");
h2_radius_mw0_mw2->GetYaxis()->SetTitle("Radius MW2");
h2_radius_mw0_mw2->GetXaxis()->CenterTitle(true);
h2_radius_mw0_mw2->GetYaxis()->CenterTitle(true);
h2_radius_mw0_mw2->GetYaxis()->SetLabelSize(0.045);
h2_radius_mw0_mw2->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_wr_td_previous_current;
sprintf(hist_name, "WR t-diff between current time and previous one");
h1_wr_td_previous_current = new TH1F(hist_name,hist_name,1000,0,1000);
h1_wr_td_previous_current->GetXaxis()->SetTitle("wr time diff in us");
h1_wr_td_previous_current->GetYaxis()->SetTitle("Counts");
h1_wr_td_previous_current->GetXaxis()->CenterTitle(true);
h1_wr_td_previous_current->GetYaxis()->CenterTitle(true);
h1_wr_td_previous_current->GetYaxis()->SetLabelSize(0.045);
h1_wr_td_previous_current->GetYaxis()->SetTitleSize(0.045);


TH1F* h1_wr_td_previous_current_only_start;
sprintf(hist_name, "WR t-diff between current time and previous one only start detector");
h1_wr_td_previous_current_only_start = new TH1F(hist_name,hist_name,1000,0,1000);
h1_wr_td_previous_current_only_start->GetXaxis()->SetTitle("wr time diff in us");
h1_wr_td_previous_current_only_start->GetYaxis()->SetTitle("Counts");
h1_wr_td_previous_current_only_start->GetXaxis()->CenterTitle(true);
h1_wr_td_previous_current_only_start->GetYaxis()->CenterTitle(true);
h1_wr_td_previous_current_only_start->GetYaxis()->SetLabelSize(0.045);
h1_wr_td_previous_current_only_start->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_wr_td_previous_current_music_start;
sprintf(hist_name, "WR t-diff between current time and previous one only start and music detector");
h1_wr_td_previous_current_music_start = new TH1F(hist_name,hist_name,1000,0,1000);
h1_wr_td_previous_current_music_start->GetXaxis()->SetTitle("wr time diff in us");
h1_wr_td_previous_current_music_start->GetYaxis()->SetTitle("Counts");
h1_wr_td_previous_current_music_start->GetXaxis()->CenterTitle(true);
h1_wr_td_previous_current_music_start->GetYaxis()->CenterTitle(true);
h1_wr_td_previous_current_music_start->GetYaxis()->SetLabelSize(0.045);
h1_wr_td_previous_current_music_start->GetYaxis()->SetTitleSize(0.045);


TH1F* h1_mw0_x_mult;
sprintf(hist_name, "Multiplicity of hits in mapped level MW0,X");
h1_mw0_x_mult = new TH1F(hist_name,hist_name,100,-0.5,99.5);
h1_mw0_x_mult->GetXaxis()->SetTitle("Multiplicity");
h1_mw0_x_mult->GetYaxis()->SetTitle("Counts");
h1_mw0_x_mult->GetXaxis()->CenterTitle(true);
h1_mw0_x_mult->GetYaxis()->CenterTitle(true);
h1_mw0_x_mult->GetYaxis()->SetLabelSize(0.045);
h1_mw0_x_mult->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw0_x_mult_selected;
sprintf(hist_name, "Multiplicity of hits in mapped level MW0,X selected events");
h1_mw0_x_mult_selected = new TH1F(hist_name,hist_name,100,-0.5,99.5);
h1_mw0_x_mult_selected->GetXaxis()->SetTitle("Multiplicity");
h1_mw0_x_mult_selected->GetYaxis()->SetTitle("Counts");
h1_mw0_x_mult_selected->GetXaxis()->CenterTitle(true);
h1_mw0_x_mult_selected->GetYaxis()->CenterTitle(true);
h1_mw0_x_mult_selected->GetYaxis()->SetLabelSize(0.045);
h1_mw0_x_mult_selected->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw0_y_mult;
sprintf(hist_name, "Multiplicity of hits in mapped level MW0,Y");
h1_mw0_y_mult = new TH1F(hist_name,hist_name,100,-0.5,99.5);
h1_mw0_y_mult->GetXaxis()->SetTitle("Multiplicity");
h1_mw0_y_mult->GetYaxis()->SetTitle("Counts");
h1_mw0_y_mult->GetXaxis()->CenterTitle(true);
h1_mw0_y_mult->GetYaxis()->CenterTitle(true);
h1_mw0_y_mult->GetYaxis()->SetLabelSize(0.045);
h1_mw0_y_mult->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_mw0_y_mult_selected;
sprintf(hist_name, "Multiplicity of hits in mapped level MW0,Y selected events");
h1_mw0_y_mult_selected = new TH1F(hist_name,hist_name,100,-0.5,99.5);
h1_mw0_y_mult_selected->GetXaxis()->SetTitle("Multiplicity");
h1_mw0_y_mult_selected->GetYaxis()->SetTitle("Counts");
h1_mw0_y_mult_selected->GetXaxis()->CenterTitle(true);
h1_mw0_y_mult_selected->GetYaxis()->CenterTitle(true);
h1_mw0_y_mult_selected->GetYaxis()->SetLabelSize(0.045);
h1_mw0_y_mult_selected->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_r3bmusic_mult;
sprintf(hist_name, "R3B Music Multiplicity of hits in mapped level");
h1_r3bmusic_mult = new TH1F(hist_name,hist_name,100,-0.5,99.5);
h1_r3bmusic_mult->GetXaxis()->SetTitle("Multiplicity");
h1_r3bmusic_mult->GetYaxis()->SetTitle("Counts");
h1_r3bmusic_mult->GetXaxis()->CenterTitle(true);
h1_r3bmusic_mult->GetYaxis()->CenterTitle(true);
h1_r3bmusic_mult->GetYaxis()->SetLabelSize(0.045);
h1_r3bmusic_mult->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_r3bmusic_mult_selected;
sprintf(hist_name, "R3B Music Multiplicity of hits in mapped level,selected events");
h1_r3bmusic_mult_selected = new TH1F(hist_name,hist_name,100,-0.5,99.5);
h1_r3bmusic_mult_selected->GetXaxis()->SetTitle("Multiplicity");
h1_r3bmusic_mult_selected->GetYaxis()->SetTitle("Counts");
h1_r3bmusic_mult_selected->GetXaxis()->CenterTitle(true);
h1_r3bmusic_mult_selected->GetYaxis()->CenterTitle(true);
h1_r3bmusic_mult_selected->GetYaxis()->SetLabelSize(0.045);
h1_r3bmusic_mult_selected->GetYaxis()->SetTitleSize(0.045);
Double_t x_low = -0.5;
Double_t x_high = 19.5;
Int_t bin_number_size = 20;
TH1F* h1_twim_music_mult_selected;
sprintf(hist_name, "TWIM Music Multiplicity of hits in Cal level,selected events");
h1_twim_music_mult_selected = new TH1F(hist_name,hist_name,bin_number_size,x_low,x_high);
h1_twim_music_mult_selected->GetXaxis()->SetTitle("Multiplicity");
h1_twim_music_mult_selected->GetYaxis()->SetTitle("Counts");
h1_twim_music_mult_selected->GetXaxis()->CenterTitle(true);
h1_twim_music_mult_selected->GetYaxis()->CenterTitle(true);
h1_twim_music_mult_selected->GetYaxis()->SetLabelSize(0.045);
h1_twim_music_mult_selected->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_twim_music_mult_rel;
sprintf(hist_name, "Relative TWIM Music Multiplicity of hits in Cal level,selected events");
h1_twim_music_mult_rel = new TH1F(hist_name,hist_name,bin_number_size,x_low,x_high);
h1_twim_music_mult_rel->GetXaxis()->SetTitle("Multiplicity");
h1_twim_music_mult_rel->GetYaxis()->SetTitle("Counts");
h1_twim_music_mult_rel->GetXaxis()->CenterTitle(true);
h1_twim_music_mult_rel->GetYaxis()->CenterTitle(true);
h1_twim_music_mult_rel->GetYaxis()->SetLabelSize(0.045);
h1_twim_music_mult_rel->GetYaxis()->SetTitleSize(0.045);

//

TChain* chain = new TChain("evt");
chain->Reset();
chain->Add(fname);

//TClonesArrays containing the TObjects R3BSofToFWTcalData,... which allow access to data over function calls
//
TClonesArray* Mwpc3HitData = new TClonesArray("R3BMwpcHitData",5);
R3BMwpcHitData** sofmwpc3hitdata;
TBranch *branchMwpc3HitData = chain->GetBranch("Mwpc3HitData");
branchMwpc3HitData->SetAddress(&Mwpc3HitData);
//
TClonesArray* SofSciTcalData = new TClonesArray("R3BSofSciTcalData",2);
R3BSofSciTcalData** sofscitcaldata;
TBranch *branchSofSciTcalData = chain->GetBranch("SofSciTcalData");
branchSofSciTcalData->SetAddress(&SofSciTcalData);
//
TClonesArray* SofTofWTcalData = new TClonesArray("R3BSofTofWTcalData",2);
R3BSofTofWTcalData** softofwtcaldata;
TBranch *branchSofTofWTcalData = chain->GetBranch("SofTofWTcalData");
branchSofTofWTcalData->SetAddress(&SofTofWTcalData);
//
TClonesArray* Mwpc0HitData = new TClonesArray("R3BMwpcHitData",5);
R3BMwpcHitData** sofmwpc0hitdata;
TBranch *branchMwpc0HitData = chain->GetBranch("Mwpc0HitData");
branchMwpc0HitData->SetAddress(&Mwpc0HitData);
//
TClonesArray* Mwpc0CalData = new TClonesArray("R3BMwpcCalData",5);
R3BMwpcCalData** sofmwpc0caldata;
TBranch *branchMwpc0CalData = chain->GetBranch("Mwpc0CalData");
branchMwpc0CalData->SetAddress(&Mwpc0CalData);
//
TClonesArray* Mwpc1HitData = new TClonesArray("R3BMwpcHitData",5);
R3BMwpcHitData** sofmwpc1hitdata;
TBranch *branchMwpc1HitData = chain->GetBranch("Mwpc1HitData");
branchMwpc1HitData->SetAddress(&Mwpc1HitData);
//
TClonesArray* Mwpc2HitData = new TClonesArray("R3BMwpcHitData",5);
R3BMwpcHitData** sofmwpc2hitdata;
TBranch *branchMwpc2HitData = chain->GetBranch("Mwpc2HitData");
branchMwpc2HitData->SetAddress(&Mwpc2HitData);
//
TClonesArray* Mwpc2CalData = new TClonesArray("R3BMwpcCalData",5);
R3BMwpcCalData** sofmwpc2caldata;
TBranch *branchMwpc2CalData = chain->GetBranch("Mwpc2CalData");
branchMwpc2CalData->SetAddress(&Mwpc2CalData);
//
TClonesArray* Mwpc3CalData = new TClonesArray("R3BMwpcCalData",5);
R3BMwpcCalData** sofmwpc3caldata;
TBranch *branchMwpc3CalData = chain->GetBranch("Mwpc3CalData");
branchMwpc3CalData->SetAddress(&Mwpc3CalData);
//
TClonesArray* SofTofWMappedData = new TClonesArray("R3BSofTofWMappedData",2);
R3BSofTofWMappedData** softofwmappeddata;
TBranch *branchSofTofWMappedData = chain->GetBranch("SofTofWMappedData");
branchSofTofWMappedData->SetAddress(&SofTofWMappedData);
//
TClonesArray* TwimHitData = new TClonesArray("R3BTwimHitData",2);
R3BTwimHitData** softwimhitdata;
TBranch *branchTwimHitData = chain->GetBranch("TwimHitData");
branchTwimHitData->SetAddress(&TwimHitData);
//

TClonesArray* TwimCalData = new TClonesArray("R3BTwimCalData",2);
R3BTwimHitData** softwimcaldata;
TBranch *branchTwimCalData = chain->GetBranch("TwimCalData");
branchTwimCalData->SetAddress(&TwimCalData);
//
TClonesArray* CalifaHitData = new TClonesArray("R3BCalifaHitData",3);
R3BCalifaHitData** califahitdata;
TBranch *branchCalifaHitData = chain->GetBranch("CalifaHitData");
branchCalifaHitData->SetAddress(&CalifaHitData);
//
TClonesArray* MusicHitData = new TClonesArray("R3BMusicHitData",1);
R3BMusicHitData** musichitdata;
TBranch *branchMusicHitData = chain->GetBranch("MusicHitData");
branchMusicHitData->SetAddress(&MusicHitData);
//
TClonesArray* MusicMappedData = new TClonesArray("R3BMusicMappedData",1);
R3BMusicHitData** musicmappeddata;
TBranch *branchMusicMappedData = chain->GetBranch("MusicMappedData");
branchMusicMappedData->SetAddress(&MusicMappedData);
//
R3BEventHeader* DataCA = new R3BEventHeader();
TBranch* branchData = chain->GetBranch("EventHeader.");
branchData->SetAddress(&DataCA);


//Counters
ULong64_t count_12c_incoming = 0;
ULong64_t count_twim_survived_all = 0;
ULong64_t count_twim_survived_carbon = 0;
ULong64_t count_mw2_survived_all = 0;
ULong64_t count_mw3_survived_all = 0;
ULong64_t count_mw3_survived_carbon = 0;
ULong64_t count_mw3_detected_11c = 0;
ULong64_t count_mw3_detected_inclusive = 0;
uint64_t current_timestamp = 0;
uint64_t previous_timestamp = 0;
uint64_t time_diff_current_previous = 0;

uint64_t current_timestamp_only_start = 0;
uint64_t previous_timestamp_only_start = 0;
uint64_t time_diff_current_previous_only_start = 0;

uint64_t current_timestamp_music_start = 0;
uint64_t previous_timestamp_music_start = 0;
uint64_t time_diff_current_previous_music_start = 0;
Long64_t nevents = chain->GetEntries(); //number of events in file with name "fname"
chain->GetEvent(1);
uint64_t time_first_event = DataCA->GetTimeStamp();
cout << "time of first event:\t" << time_first_event << endl;
chain->GetEvent(nevents-1);
uint64_t time_last_event = DataCA->GetTimeStamp();
cout << "Time last event:\t" << time_last_event << endl;
uint64_t time_diff_first_last = (time_last_event - time_first_event)/(pow(10,9));
cout << "Time diff in seconds:\t" << time_diff_first_last << endl;
Int_t bin_nr = int(time_diff_first_last) * 2;

TH1F* h1_beam_freq;
sprintf(hist_name, " Beam Intensity ");
h1_beam_freq = new TH1F(hist_name,hist_name,bin_nr,0,int(time_diff_first_last));
h1_beam_freq->GetXaxis()->SetTitle("Time (0.5s bins)");
h1_beam_freq->GetYaxis()->SetTitle("Counts");
h1_beam_freq->GetXaxis()->CenterTitle(true);
h1_beam_freq->GetYaxis()->CenterTitle(true);
h1_beam_freq->GetYaxis()->SetLabelSize(0.045);
h1_beam_freq->GetYaxis()->SetTitleSize(0.045);

for(Long64_t i=0; i< nevents;i++){
    Long64_t evtnr = i;

    if (i%100000==0)
        cout<<"Processing event "<<i<<endl;

	chain->GetEvent(i); //event_i_call
	Int_t entries_mw0 = Mwpc0HitData->GetEntries();
    Int_t entries_mw1 = Mwpc1HitData->GetEntries();
    Int_t entries_mw2 = Mwpc2HitData->GetEntries();
    Int_t entries_mw3 = Mwpc3HitData->GetEntries();
    Int_t entries_start = SofSciTcalData->GetEntriesFast();
    Int_t entries_tofw = SofTofWTcalData->GetEntriesFast();
    Int_t entries_califa = CalifaHitData->GetEntries();
	Int_t entries_music = MusicHitData->GetEntries();
	Int_t entries_twim = TwimHitData->GetEntries();
	Int_t entries_mw2_cal = Mwpc2CalData->GetEntries(); 
	Int_t entries_mw3_cal = Mwpc3CalData->GetEntries(); 
	Int_t entries_mw0_cal = Mwpc0CalData->GetEntries();
	Int_t entries_music_cal = MusicMappedData->GetEntries();
	Int_t entries_twim_cal = TwimCalData->GetEntries();

	//this is a pre step to check the start detector
	if (entries_start == 2 && entries_music == 1){
		current_timestamp_music_start = DataCA->GetTimeStamp();
		time_diff_current_previous_music_start = current_timestamp_music_start-previous_timestamp_music_start;
		//fill histo
		h1_wr_td_previous_current_music_start->Fill(double(time_diff_current_previous_music_start)/1000.);
		previous_timestamp_music_start = current_timestamp_music_start;
		//cout << "time diff to first:\t" << int64_t(current_timestamp_music_start-time_first_event) << endl;
		//cout << "intensity:\t" << (double(current_timestamp_music_start-time_first_event)/pow(10,9)) << endl;
		h1_beam_freq->Fill((double(current_timestamp_music_start-time_first_event)/pow(10,9)));
		}
	if (entries_start){
		current_timestamp_only_start = DataCA->GetTimeStamp();
		time_diff_current_previous_only_start = current_timestamp_only_start-previous_timestamp_only_start;
		//fill histo
		h1_wr_td_previous_current_only_start->Fill(double(time_diff_current_previous_only_start)/1000.);
		previous_timestamp_only_start = current_timestamp_only_start;
		}

	//this is the new part
	if (entries_start ==2 && entries_music == 1 && entries_mw0 == 1 ){
		//check wr timestamps
		current_timestamp = DataCA->GetTimeStamp();
		time_diff_current_previous = current_timestamp-previous_timestamp;
		h1_wr_td_previous_current->Fill(double(time_diff_current_previous)/1000.);
		previous_timestamp = current_timestamp;
		//fill music, no cuts 
		musichitdata = new R3BMusicHitData*[1];
		musichitdata[0] = (R3BMusicHitData*)MusicHitData->At(0);
		Double_t charge_music = musichitdata[0]->GetZcharge();
		h1_charge_val_music->Fill(charge_music);
		//fill mw0, no cuts
		sofmwpc0hitdata = new R3BMwpcHitData*[1];
		sofmwpc0hitdata[0] = (R3BMwpcHitData*)Mwpc0HitData->At(0);
		Double_t xMW0 = sofmwpc0hitdata[0]->GetX();
		Double_t yMW0 = sofmwpc0hitdata[0]->GetY();
		h1_mw0_x_pos->Fill(xMW0);
		h1_mw0_y_pos->Fill(yMW0);
		
		//in the following check on multiplicities for MW0(X &Y) and R3BMusic
		//MW0
		sofmwpc0caldata = new R3BMwpcCalData*[entries_mw0_cal];		
		Int_t count_mw0_x = 0;
		Int_t count_mw0_y = 0;
		for (Int_t l = 0; l < entries_mw0_cal;l++){
			sofmwpc0caldata[l] = (R3BMwpcCalData*)Mwpc0CalData->At(l);
			UInt_t plane = sofmwpc0caldata[l]->GetPlane();
			if (plane == 1) count_mw0_x++;
			if (plane == 3) count_mw0_y++;
			}
		delete [] sofmwpc0caldata;
		h1_mw0_x_mult->Fill(count_mw0_x);
		h1_mw0_y_mult->Fill(count_mw0_y);
		h1_r3bmusic_mult->Fill(entries_music_cal);
		//if (entries_music_cal > 20) cout << "This is eventnr with large r3bmusic multiplicity:\t" << evtnr << endl;
		

		
		//following doing check on the incoming particles. I just want to use the ones which are in a good focus
		if(entries_start ==2 && entries_music == 1 && entries_mw0 == 1 && entries_mw2 == 1){
			sofmwpc2hitdata = new R3BMwpcHitData*[1];
			sofmwpc2hitdata[0] = (R3BMwpcHitData*)Mwpc2HitData->At(0);
			Double_t xMW2 = sofmwpc2hitdata[0]->GetX();
			Double_t yMW2 = sofmwpc2hitdata[0]->GetY();
			h2_xmw0_vs_xmw2->Fill(xMW0,xMW2);
			h2_ymw0_vs_ymw2->Fill(yMW0,yMW2);
			h2_xmw0_vs_ymw0->Fill(xMW0,yMW0);
			h2_xmw2_vs_ymw2->Fill(xMW2,yMW2);
			h1_mw0_radius->Fill(sqrt(xMW0*xMW0+yMW0*yMW0));
			h1_mw2_radius->Fill(sqrt(xMW2*xMW2+yMW2*yMW2));
			h2_radius_mw0_mw2->Fill(sqrt(xMW0*xMW0+yMW0*yMW0),sqrt(xMW2*xMW2+yMW2*yMW2));
			delete [] sofmwpc2hitdata;
			}
			//check what's happening with events that do not hit TWIM
			if(charge_music < max_music_charge && charge_music > min_music_charge && xMW0 < max_mw0x && xMW0 > min_mw0x && yMW0 < max_mw0y && yMW0 > min_mw0y && entries_twim == 0){
				h1_mw0_radius_no_twim->Fill(sqrt(xMW0*xMW0+yMW0*yMW0));
				h2_xmw0_vs_ymw0_no_twim->Fill(xMW0,yMW0);		
				}
			if(charge_music < max_music_charge && charge_music > min_music_charge && entries_twim == 0){
				h2_xmw0_vs_ymw0_no_twim_no_mw_cut->Fill(xMW0,yMW0);	
				}

		//now cutting on mw0 and music to get clean incoming 12C beam, also cut on the multiplicity of MW0 and R3B Music
		if(charge_music < max_music_charge && charge_music > min_music_charge && xMW0 < max_mw0x && xMW0 > min_mw0x && yMW0 < max_mw0y && yMW0 > min_mw0y && count_mw0_x < 6 && count_mw0_y < 6 && entries_music_cal == 10){
			count_12c_incoming += 1;
			h1_mw0_x_mult_selected->Fill(count_mw0_x);
			h1_mw0_y_mult_selected->Fill(count_mw0_y);
			h1_r3bmusic_mult_selected->Fill(entries_music_cal);
			//check how many events do we see in twim
			if (entries_twim == 1 && entries_twim_cal < 17){
				//if (entries_twim_cal ==17) cout << "eventnr with mult_twim == 17:\t" << evtnr << endl;
				h1_twim_music_mult_selected->Fill(entries_twim_cal);
				softwimhitdata = new R3BTwimHitData*[1];
				softwimhitdata[0] = (R3BTwimHitData*)TwimHitData->At(0);
				count_twim_survived_all +=1;
				Double_t charge_twim = softwimhitdata[0]->GetZcharge();
				h1_charge_val->Fill(charge_twim);
				//cut on charge == 6 on twim
				if (charge_twim < max_twim_charge && charge_twim > min_twim_charge){
					count_twim_survived_carbon +=1;
					//look at mw2 and mw3
					if (entries_mw2 == 1){
						//check multiplicity of mw2x
						sofmwpc2caldata = new R3BMwpcCalData*[entries_mw2_cal];		
						vector <int> v_pads;
						for (Int_t l = 0; l < entries_mw2_cal;l++){
							sofmwpc2caldata[l] = (R3BMwpcCalData*)Mwpc2CalData->At(l);
							UInt_t plane = sofmwpc2caldata[l]->GetPlane();
							UInt_t pad = sofmwpc2caldata[l]->GetPad();
							if( plane ==1 || plane ==2){
								v_pads.push_back(pad);	
								}
							}
						delete [] sofmwpc2caldata;
						sort( v_pads.begin(), v_pads.end() );
						v_pads.erase( unique( v_pads.begin(), v_pads.end() ), v_pads.end() );
						Int_t size_hits_cal = v_pads.size();
						
						if (size_hits_cal < 5){
						sofmwpc2hitdata = new R3BMwpcHitData*[1];
						sofmwpc2hitdata[0] = (R3BMwpcHitData*)Mwpc2HitData->At(0);
						Double_t xMW2 = sofmwpc2hitdata[0]->GetX();
						Double_t yMW2 = sofmwpc2hitdata[0]->GetY();
						count_mw2_survived_all +=1;		
						if (entries_mw3 == 1){
							//check multiplicity of mw3x
							
							sofmwpc3caldata = new R3BMwpcCalData*[entries_mw3_cal];		
							vector <int> v_pads3;
							for (Int_t l = 0; l < entries_mw3_cal;l++){
								sofmwpc3caldata[l] = (R3BMwpcCalData*)Mwpc3CalData->At(l);
								UInt_t plane = sofmwpc3caldata[l]->GetPlane();
								UInt_t pad = sofmwpc3caldata[l]->GetPad();
								if( plane ==1 || plane ==2){
									v_pads3.push_back(pad);	
									}
								}
							delete [] sofmwpc3caldata;
							sort( v_pads3.begin(), v_pads3.end() );
							v_pads3.erase( unique( v_pads3.begin(), v_pads3.end() ), v_pads3.end() );
							Int_t size_hits_cal3 = v_pads3.size();


							if (size_hits_cal3 < 6){
							sofmwpc3hitdata = new R3BMwpcHitData*[1];
        					sofmwpc3hitdata[0] = (R3BMwpcHitData*)Mwpc3HitData->At(0);
        					Double_t xMW3 = sofmwpc3hitdata[0]->GetX();
        					Double_t yMW3 = sofmwpc3hitdata[0]->GetY();
							count_mw3_survived_all += 1;
							h2_xmw2_vs_xmw3_mw0cut->Fill(xMW2,xMW3);
							//if (!(xMW2< 400 && xMW2> -400)) cout << "this is event with xmw2 larger 400:\t" << evtnr <<  "\t and correpsonding x value:\t" << xMW2 << endl;
							if (xMW3<-300 ) cout << "this is event with xmw3 smaller -300:\t" << evtnr << endl;
							if (enable_cut){
								if(cut_11c->IsInside(xMW2,xMW3)){
									count_mw3_detected_11c += 1;	
									h2_xmw2_vs_xmw3_11c->Fill(xMW2,xMW3);
									}	
								}	
							if (enable_cut_all){
								if(cut_all->IsInside(xMW2,xMW3)){
									count_mw3_detected_inclusive += 1;	
									h2_xmw2_vs_xmw3_inclusive->Fill(xMW2,xMW3);
									}	
								}	
						

							delete [] sofmwpc3hitdata;
							}
							}
							
							delete [] sofmwpc2hitdata;
							}
						}
					}
				delete [] softwimhitdata;	
				}
			}
	    delete [] musichitdata;
		delete [] sofmwpc0hitdata;
		}
	//-------------------- end of new part
}


//insert counting numbers to par file
par_file << "Clean incoming 12C:\t" << count_12c_incoming << endl;
par_file << "Number of survived hits in TWIM:\t" << count_twim_survived_all << endl;
par_file << "Number of carbon in TWIM:\t" << count_twim_survived_carbon << endl; 
par_file << "Number of survived hits in MW2:\t" << count_mw2_survived_all << endl;
par_file << "Number of survived hits in MW3:\t" << count_mw3_survived_all << endl;
if (enable_cut) par_file << " Number of survived 11C's in MW3:\t" << count_mw3_detected_11c << endl;
if (enable_cut_all) par_file << "Number of survived 11C's & 10C's (inclusive) in MW3:\t" << count_mw3_detected_inclusive << endl;
par_file << "*******END OF FILE **************" << endl;
par_file.close();

//write important number to csv file
csv_file << "Clean incoming 12C (R3BMusic),carbon_twim,hits_mw3,particles_cut10_11C_mw2_mw3,count_twim_survived_all" << endl;
csv_file << count_12c_incoming << "," << count_twim_survived_carbon << "," << count_mw3_survived_all << "," << count_mw3_detected_inclusive << "," << count_twim_survived_all << endl;
csv_file.close();

//create canvas with charge cuts
TCanvas charge_dist_canv("charge_dist_canv", "charge_dist_canv",2000,1600);
charge_dist_canv.Divide(2,1);
charge_dist_canv.cd(1);
gPad-> SetLogy();
h1_charge_val_music->Draw();
Float_t ymax = h1_charge_val_music->GetMaximum();
TLine* line_r3b_music_low = new TLine(min_music_charge,0,min_music_charge,ymax);
TLine* line_r3b_music_high = new TLine(max_music_charge,0,max_music_charge,ymax);
line_r3b_music_low->SetLineColor(kRed);
line_r3b_music_high->SetLineColor(kRed);
line_r3b_music_low->Draw();
line_r3b_music_high->Draw();
charge_dist_canv.cd(2);
gPad-> SetLogy();
h1_charge_val->Draw();
Float_t ymax_twim = h1_charge_val->GetMaximum();
TLine* line_twim_music_low = new TLine(min_twim_charge,0,min_twim_charge,ymax_twim);
TLine* line_twim_music_high = new TLine(max_twim_charge,0,max_twim_charge,ymax_twim);
line_twim_music_low->SetLineColor(kRed);
line_twim_music_high->SetLineColor(kRed);
line_twim_music_low->Draw();
line_twim_music_high->Draw();
charge_dist_canv.Modified();
charge_dist_canv.Update();
charge_dist_canv.Print(canv_name);

TCanvas intensity_canv("intensity_canv", "intensity_canv",2000,1600);
h1_beam_freq->Draw();
intensity_canv.Modified();
intensity_canv.Update();
intensity_canv.Print(canv_beam);

TFile * f = new TFile(f_out_name,"RECREATE");
TList *l = new TList();
l->Add(h1_mw0_x_pos);
l->Add(h1_mw0_y_pos);
l->Add(h2_xmw0_vs_xmw2);
l->Add(h2_ymw0_vs_ymw2);
l->Add(h2_xmw0_vs_ymw0);
l->Add(h2_xmw2_vs_ymw2);
l->Add(h1_mw0_radius);
l->Add(h1_mw2_radius);
l->Add(h2_radius_mw0_mw2);
l->Add(h2_xmw2_vs_xmw3_mw0cut);
l->Add(h2_xmw2_vs_xmw3_11c);
l->Add(h2_xmw2_vs_xmw3_inclusive);
l->Add(h1_charge_val_music);
l->Add(h1_charge_val);
l->Add(h1_clean_incoming_12c);
l->Add(h2_xmw0_vs_ymw0_no_twim);
l->Add(h1_mw0_radius_no_twim);
l->Add(h2_xmw0_vs_ymw0_no_twim_no_mw_cut);
l->Add(h1_wr_td_previous_current);
l->Add(h1_wr_td_previous_current_only_start);
l->Add(h1_wr_td_previous_current_music_start);
l->Add(h1_mw0_x_mult);
l->Add(h1_mw0_x_mult_selected);
l->Add(h1_mw0_y_mult);
l->Add(h1_mw0_y_mult_selected);
l->Add(h1_r3bmusic_mult);
l->Add(h1_r3bmusic_mult_selected);
l->Add(h1_twim_music_mult_selected);
l->Add(h1_beam_freq);
Double_t tot_entries_twim = h1_twim_music_mult_selected->GetEntries();
for (Int_t i = 1; i<= bin_number_size; i++){
	Double_t bin_content = Double_t((h1_twim_music_mult_selected->GetBinContent(i)))/tot_entries_twim;
	h1_twim_music_mult_rel->SetBinContent(i,bin_content);
	} 
TCanvas twim_mult_canv("twim_mult_canv", "twim_mult_canv",2000,1600);
h1_twim_music_mult_rel->Draw();
twim_mult_canv.Modified();
twim_mult_canv.Update();
twim_mult_canv.Print(canv_mult_twim);
l->Add(h1_twim_music_mult_rel);

l->Write("histlist", TObject::kSingleKey);

delete Mwpc3HitData;
delete SofSciTcalData;
delete SofTofWTcalData;
delete Mwpc0HitData;
delete Mwpc1HitData;
delete Mwpc2HitData;
delete SofTofWMappedData;
delete TwimHitData;
delete CalifaHitData;
delete MusicHitData;
cout << "end of process, successful" <<endl;
}



