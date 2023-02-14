#include <iostream>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "find_z.C"
#include "find_isotopes.C"
#include "all_params.C"
#include "find_psi_par.C"
#include <string>
#include<ctime>
#include <cmath>


#define MASS_PROTON 938.270
#define PI 3.14159265359
using namespace std;
//sort tlorentzvector according to energy
bool sortcol_lorentz( const TLorentzVector& v1,const TLorentzVector& v2 ) {
                                    return v1[3] > v2[3];
                                    }

void measure_cross_sec_two_protons(char const count_i[50]){
srand(time(0));
char hist_name[500];
char f_out_name[500];
char f_events_name[500];
char fname[500];
char output_values_stream[500];
sprintf(output_values_stream,"/home/ge37liw/new_vec_for_cross_sec/values_run_%s.txt", count_i);
ofstream mystream;
mystream.open(output_values_stream);
const Double_t beta_given = 0.714549;//with ingoing 12C 400 MeV
bool simulation = false;
//if ((string(count_i)).compare("sim_400") == 0 || (string(count_i)).compare("sim_550") == 0 || (string(count_i)).compare("sim_650") == 0 || (string(count_i)).compare("sim_800") == 0){
//	simulation = true;
//	sprintf(fname,"/scratch5/ge37liw/updated_r3b/R3BRoot/macros/r3b/califa/%s.root",count_i);
//	}
//newer version of above code ...
if ((string(count_i)).find("sim") != string::npos){
	simulation = true;
	sprintf(fname,"/home/ge37liw/%s.root",count_i);
        if(string(inputfile).find("400") != string::npos){
            if(string(inputfile).find("CH2Target_1229") != string::npos){
                cout << "sim" << "\t" << "CH2Target_1229" << "\t" << "400" << "\t" << f1->GetParameter(1) << "\t" << f1->GetParameter(2)  << endl;
				f_out_name = "/home/ge37liw/califa_cross_sec_400_CH2Target_1229_sim.root";
                }
            if(string(inputfile).find("CH2Target_2453") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_400_CH2Target_2453_sim.root";
                }
            if(string(inputfile).find("CTarget_54") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_400_CTarget_54_sim.root";
                }
            if(string(inputfile).find("CTarget_1086") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_400_CTarget_1086_sim.root";
                }
            if(string(inputfile).find("CTarget_2198") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_400_CTarget_2198_sim.root";
                }
            }
        if(string(inputfile).find("550") != string::npos){
            if(string(inputfile).find("CH2Target_1229") != string::npos){
                cout << "sim" << "\t" << "CH2Target_1229" << "\t" << "550" << "\t" << f1->GetParameter(1) << "\t" << f1->GetParameter(2)  << endl;
				f_out_name = "/home/ge37liw/califa_cross_sec_550_CH2Target_1229_sim.root";
                }
            if(string(inputfile).find("CH2Target_2453") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_550_CH2Target_2453_sim.root";
                }
            if(string(inputfile).find("CTarget_54") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_550_CTarget_54_sim.root";
                }
            if(string(inputfile).find("CTarget_1086") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_550_CTarget_1086_sim.root";
                }
            if(string(inputfile).find("CTarget_2198") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_550_CTarget_2198_sim.root";
                }
            }
        if(string(inputfile).find("650") != string::npos){
            if(string(inputfile).find("CH2Target_1229") != string::npos){
                cout << "sim" << "\t" << "CH2Target_1229" << "\t" << "650" << "\t" << f1->GetParameter(1) << "\t" << f1->GetParameter(2)  << endl;
				f_out_name = "/home/ge37liw/califa_cross_sec_650_CH2Target_1229_sim.root";
                }
            if(string(inputfile).find("CH2Target_2453") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_650_CH2Target_2453_sim.root";
                }
            if(string(inputfile).find("CTarget_54") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_650_CTarget_54_sim.root";
                }
            if(string(inputfile).find("CTarget_1086") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_650_CTarget_1086_sim.root";
                }
            if(string(inputfile).find("CTarget_2198") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_650_CTarget_2198_sim.root";
                }
            }
        if(string(inputfile).find("800") != string::npos){
            if(string(inputfile).find("CH2Target_1229") != string::npos){
                cout << "sim" << "\t" << "CH2Target_1229" << "\t" << "800" << "\t" << f1->GetParameter(1) << "\t" << f1->GetParameter(2)  << endl;
				f_out_name = "/home/ge37liw/califa_cross_sec_800_CH2Target_1229_sim.root";
                }
            if(string(inputfile).find("CH2Target_2453") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_800_CH2Target_2453_sim.root";
                }
            if(string(inputfile).find("CTarget_54") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_800_CTarget_54_sim.root";
                }
            if(string(inputfile).find("CTarget_1086") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_800_CTarget_1086_sim.root";
                }
            if(string(inputfile).find("CTarget_2198") != string::npos){
				f_out_name = "/home/ge37liw/califa_cross_sec_800_CTarget_2198_sim.root";
                }
            }
		
	}

//write important numbers to .txt file
ofstream myfile ("data_cross_section_s444.txt", myfile.out | myfile.app);


//sprintf(f_events_name,"/home/ge37liw/vectors_for_cross_sec/cross_sec_p2p_events_%s.root", count_i);
sprintf(f_events_name,"/home/ge37liw/new_vec_for_cross_sec/cross_sec_p2p_events_%s.root", count_i);
TFile* eventFile = TFile::Open(f_events_name,"READ");
vector<Long64_t> v_events;
vector<Long64_t> vec_interating_loop;
TFile * f_events;
TTree *t_events;
if (!eventFile){
	cout << "hello, I have to create new file..." << endl;
    f_events = TFile::Open(f_events_name,"RECREATE");
    t_events = new TTree("tvec","Tree with vectors");
    t_events->Branch("v_events",&v_events);
    }
if (eventFile){
	cout << "FILE EXISTS; NOW IT GOES FAST!!!!!!!!!!! " << endl;
    TTreeReader myReader("tvec", eventFile);
    TTreeReaderValue<vector<Long64_t> >vec_good_events(myReader,"v_events");
    myReader.Next();
    vec_interating_loop = *vec_good_events;
    eventFile->Close();
    }

//SIMULATION-----------
//sprintf(f_out_name,"/home/ge37liw/califa_opening_angle_simulation.root"); //simulation run -<<
//sprintf(fname,"/scratch5/ge37liw/updated_r3b/R3BRoot/macros/r3b/califa/clean_event_sim_out.root"); //simulation run-<<
//sprintf(fname,"/scratch5/ge37liw/updated_r3b/R3BRoot/macros/r3b/califa/sim_out.root"); //simulation run-<<
//
int64_t beam_energy;
string target_type;

//write important numbers to file -->empty runs!
if((string(count_i)).compare("0069_0001") == 0 || (string(count_i)).compare("0074_0001") == 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/empty_400amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 400;
	target_type = "empty";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_400_empty_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0096_0001") == 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/empty_550amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 550;
	target_type = "empty";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_550_empty_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0124_0001") == 0 || (string(count_i)).compare("0122_0001") == 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/empty_650amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 650;
	target_type = "empty";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_650_empty_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0173_0001") == 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/empty_800amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 800;
	target_type = "empty";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_800_empty_exp_%s.root",count_i);
	}
//write important numbers to file -->target runs!
//PLASTIC 12.29mm
if((string(count_i)).compare("0180_0001") == 0 || (string(count_i)).compare("0188_0001") == 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/plastic_target_1229_400amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 400;
	target_type = "ch2_1229";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_400_CH2Target_1229_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0100_0001") == 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/plastic_target_1229_550amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 550;
	target_type = "ch2_1229";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_550_CH2Target_1229_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0125_0001")== 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/plastic_target_1229_650amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 650;
	target_type = "ch2_1229";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_650_CH2Target_1229_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0161_0001")== 0 ||(string(count_i)).compare("0167_0001")== 0 ){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/plastic_target_1229_800amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 800;
	target_type = "ch2_1229";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_800_CH2Target_1229_exp_%s.root",count_i);
	}
//PLASTIC 24.53mm
if((string(count_i)).compare("0018_0001")== 0 || (string(count_i)).compare("0182_0001")== 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/plastic_target_2453_400amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 400;
	target_type = "ch2_2453";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_400_CH2Target_2453_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0102_0001")== 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/plastic_target_2453_550amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 550;
	target_type = "ch2_2453";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_550_CH2Target_2453_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0128_0001")== 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/plastic_target_2453_650amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 650;
	target_type = "ch2_2453";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_650_CH2Target_2453_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0169_0001")== 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/plastic_target_2453_800amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 800;
	target_type = "ch2_2453";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_800_CH2Target_2453_exp_%s.root",count_i);
	}
//CARBON 5.4mm
if((string(count_i)).compare("0179_0001")== 0 || (string(count_i)).compare("0192_0001")== 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_5_400amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 400;
	target_type = "c_54";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_400_CTarget_54_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0099_0001")== 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_5_550amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 550;
	target_type = "c_54";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_550_CTarget_54_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0126_0001")== 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_5_650amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 650;
	target_type = "c_54";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_650_CTarget_54_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0166_0001")== 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_5_800amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 800;
	target_type = "c_54";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_800_CTarget_54_exp_%s.root",count_i);
	}
//CARBON 10.86mm
if((string(count_i)).compare("0101_0001")== 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_1086_550amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 550;
	target_type = "c_1086";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_550_CTarget_1086_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0127_0001")== 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_1086_650amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 650;
	target_type = "c_1086";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_650_CTarget_1086_exp_%s.root",count_i);
	}
//CARBON 21.98mm
if((string(count_i)).compare("0032_0001") == 0 || (string(count_i)).compare("0032_0002")== 0 || (string(count_i)).compare("0032_0003")== 0 || (string(count_i)).compare("0075_0001")== 0 || (string(count_i)).compare("0183_0001")== 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_2198_400amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 400;
	target_type = "c_2198";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_400_CTarget_2198_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0103_0001") == 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_2198_550amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 550;
	target_type = "c_2198";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_550_CTarget_2198_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0130_0001") == 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_2198_650amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 650;
	target_type = "c_2198";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_650_CTarget_2198_exp_%s.root",count_i);
	}
if((string(count_i)).compare("0170_0001") == 0){
	sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_2198_800amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);
	beam_energy = 800;
	target_type = "c_2198";
	sprintf(f_out_name,"/home/ge37liw/califa_cross_sec_800_CTarget_2198_exp_%s.root",count_i);
	}




//PLASTIC-------------
//sprintf(f_out_name,"/home/ge37liw/califa_opening_angle_ch2_run_%s.root",count_i);
////sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/p2p_data_fixed_angle/stitched_and_unpacked_main%s.root",count_i);
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/p2p_data_fixed_angle/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);

//EMPTY---------------
//sprintf(f_out_name,"/home/ge37liw/califa_opening_angle_empty_run.root"); //empty run -<<
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/empty_400amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i); //empty run-<<


//CARBON--------------21 mm
//sprintf(f_out_name,"/home/ge37liw/califa_opening_angle_carbon_run_%s.root",count_i);
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_2198_400amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i); //empty run-<<


//CARBON-------------- 5 mm
//sprintf(f_out_name,"/home/ge37liw/califa_opening_angle_thin_carbon_run_%s.root",count_i);
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_5_400amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i); //empty run-<<

//PLASTIC ------ 24.53 mm
//sprintf(f_out_name,"/home/ge37liw/califa_opening_angle_thick_ch2_run_%s.root",count_i);
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/plastic_target_2453_400amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);

//PLASTIC -------12.29mm, min bias
//sprintf(f_out_name,"/home/ge37liw/califa_opening_angle_min_bias_ch2_run_%s.root",count_i);
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/plastic_target_1229_400amev/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);


//sprintf(f_out_name,"/home/ge37liw/califa_opening_angle_ch2_run_%s.root",count_i);
//sprintf(f_out_name,"/home/ge37liw/califa_opening_angle_empty_run.root"); //empty run -<<
//sprintf(f_out_name,"/home/ge37liw/califa_opening_angle_simulation.root"); //simulation run -<<
//sprintf(f_out_name,"/home/ge37liw/califa_opening_angle_carbon_run_%s.root",count_i);
TFile * f = new TFile(f_out_name,"RECREATE");
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_2198_400amev/stitched_and_unpacked_main%s.root",count_i);
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/p2p_data_fixed_angle/stitched_and_unpacked_main%s.root",count_i);
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/empty_400amev/stitched_and_unpacked_main%s.root",count_i); //empty run-<<
//sprintf(fname,"/scratch5/ge37liw/updated_r3b/R3BRoot/macros/r3b/califa/sim_out.root"); //simulation run-<<

TChain* chain = new TChain("evt");
chain->Reset();
chain->Add(fname);


//BEGIN TCA

TClonesArray* CalifaHitData = new TClonesArray("R3BCalifaHitData",3);
R3BCalifaHitData** califahitdata;
TBranch *branchCalifaHitData = chain->GetBranch("CalifaHitData");
branchCalifaHitData->SetAddress(&CalifaHitData);

TClonesArray* CalifaCalData = new TClonesArray("R3BCalifaCrystalCalData",3);
R3BCalifaCrystalCalData** califacaltdata;
TBranch *branchCalifaCalData = chain->GetBranch("CalifaCrystalCalData");
branchCalifaCalData->SetAddress(&CalifaCalData);

//R3BEventHeader* DataCA = new R3BEventHeader();
//TBranch* branchData = chain->GetBranch("EventHeader.");
//branchData->SetAddress(&DataCA);

//end of TCA
Long64_t n_event_min_bias = 0;
Long64_t n_event_califa_p2p = 0;
Long64_t n_event_califa_p2p_inclusive = 0;
Long64_t n_evt_two_hits_sim = 0;
Long64_t n_evt_two_hits_sim_protons = 0;


//HISTOGRAMS
TH1F* h1_opang_p2p;
sprintf(hist_name, "Opening Angle of the two correlated protons");
h1_opang_p2p = new TH1F(hist_name,hist_name,52,22.15,152.15);
h1_opang_p2p->GetXaxis()->SetTitle("Opening Angle between two protons [degr]");
h1_opang_p2p->GetYaxis()->SetTitle("Counts");
h1_opang_p2p->GetXaxis()->CenterTitle(true);
h1_opang_p2p->GetYaxis()->CenterTitle(true);
h1_opang_p2p->GetYaxis()->SetLabelSize(0.045);
h1_opang_p2p->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_opang_p2p_nfns_positive;
sprintf(hist_name, "Opening Angle of the two correlated protons, if both nf and ns positive");
h1_opang_p2p_nfns_positive = new TH1F(hist_name,hist_name,52,22.15,152.15);
h1_opang_p2p_nfns_positive->GetXaxis()->SetTitle("Opening Angle between two protons [degr]");
h1_opang_p2p_nfns_positive->GetYaxis()->SetTitle("Counts");
h1_opang_p2p_nfns_positive->GetXaxis()->CenterTitle(true);
h1_opang_p2p_nfns_positive->GetYaxis()->CenterTitle(true);
h1_opang_p2p_nfns_positive->GetYaxis()->SetLabelSize(0.045);
h1_opang_p2p_nfns_positive->GetYaxis()->SetTitleSize(0.045);


TH1F* h1_min_bias_events;
sprintf(hist_name, "Number of min. bias events");
h1_min_bias_events = new TH1F(hist_name,hist_name,10,0,10);
h1_min_bias_events->GetXaxis()->SetTitle("Number of min. bias events");
h1_min_bias_events->GetYaxis()->SetTitle("Counts");
h1_min_bias_events->GetXaxis()->CenterTitle(true);
h1_min_bias_events->GetYaxis()->CenterTitle(true);
h1_min_bias_events->GetYaxis()->SetLabelSize(0.045);
h1_min_bias_events->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_p2p_events;
sprintf(hist_name, "Number of p2p  events with opening angle 80-90");
h1_p2p_events = new TH1F(hist_name,hist_name,10,0,10);
h1_p2p_events->GetXaxis()->SetTitle("Number p2p events");
h1_p2p_events->GetYaxis()->SetTitle("Counts");
h1_p2p_events->GetXaxis()->CenterTitle(true);
h1_p2p_events->GetYaxis()->CenterTitle(true);
h1_p2p_events->GetYaxis()->SetLabelSize(0.045);
h1_p2p_events->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_p2p_events_inclusive;
sprintf(hist_name, "Number of p2p  events inclusive");
h1_p2p_events_inclusive = new TH1F(hist_name,hist_name,10,0,10);
h1_p2p_events_inclusive->GetXaxis()->SetTitle("Number p2p events");
h1_p2p_events_inclusive->GetYaxis()->SetTitle("Counts");
h1_p2p_events_inclusive->GetXaxis()->CenterTitle(true);
h1_p2p_events_inclusive->GetYaxis()->CenterTitle(true);
h1_p2p_events_inclusive->GetYaxis()->SetLabelSize(0.045);
h1_p2p_events_inclusive->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_p2p_proton_energies;
sprintf(hist_name, "Proton energies for inclusive p2p events");
h1_p2p_proton_energies = new TH1F(hist_name,hist_name,700,0,700);
h1_p2p_proton_energies->GetXaxis()->SetTitle("Energy");
h1_p2p_proton_energies->GetYaxis()->SetTitle("Counts");
h1_p2p_proton_energies->GetXaxis()->CenterTitle(true);
h1_p2p_proton_energies->GetYaxis()->CenterTitle(true);
h1_p2p_proton_energies->GetYaxis()->SetLabelSize(0.045);
h1_p2p_proton_energies->GetYaxis()->SetTitleSize(0.045);


TH1F* h1_azimuthal_ang_diff;
sprintf(hist_name, "Azimuthal angular difference p2p events");
h1_azimuthal_ang_diff = new TH1F(hist_name,hist_name,160,100,260);
h1_azimuthal_ang_diff->GetXaxis()->SetTitle("Energy");
h1_azimuthal_ang_diff->GetYaxis()->SetTitle("Counts");
h1_azimuthal_ang_diff->GetXaxis()->CenterTitle(true);
h1_azimuthal_ang_diff->GetYaxis()->CenterTitle(true);
h1_azimuthal_ang_diff->GetYaxis()->SetLabelSize(0.045);
h1_azimuthal_ang_diff->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_num_califa_hits;
sprintf(hist_name, "Number of crystal hits for p2p events");
h1_num_califa_hits = new TH1F(hist_name,hist_name,100,0,100);
h1_num_califa_hits->GetXaxis()->SetTitle("Number of crystal hits");
h1_num_califa_hits->GetYaxis()->SetTitle("Counts");
h1_num_califa_hits->GetXaxis()->CenterTitle(true);
h1_num_califa_hits->GetYaxis()->CenterTitle(true);
h1_num_califa_hits->GetYaxis()->SetLabelSize(0.045);
h1_num_califa_hits->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_e_sum_p2p;
sprintf(hist_name, "Sum Energy for p2p events");
h1_e_sum_p2p = new TH1F(hist_name,hist_name,700,0,700);
h1_e_sum_p2p->GetXaxis()->SetTitle("Sum Energy in MeV");
h1_e_sum_p2p->GetYaxis()->SetTitle("Counts");
h1_e_sum_p2p->GetXaxis()->CenterTitle(true);
h1_e_sum_p2p->GetYaxis()->CenterTitle(true);
h1_e_sum_p2p->GetYaxis()->SetLabelSize(0.045);
h1_e_sum_p2p->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_nf_ns;
sprintf(hist_name, "Nf vs Ns for p2p events");
h2_nf_ns = new TH2F(hist_name,hist_name,500,0,200000,500,0,200000);
h2_nf_ns->GetXaxis()->SetTitle("Nf");
h2_nf_ns->GetYaxis()->SetTitle("Ns");
h2_nf_ns->GetXaxis()->CenterTitle(true);
h2_nf_ns->GetYaxis()->CenterTitle(true);
h2_nf_ns->GetYaxis()->SetLabelSize(0.045);
h2_nf_ns->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_nf_ns_carbon;
sprintf(hist_name, "Nf vs Ns for p2p events, carbon fragment");
h2_nf_ns_carbon = new TH2F(hist_name,hist_name,500,0,200000,500,0,200000);
h2_nf_ns_carbon->GetXaxis()->SetTitle("Nf");
h2_nf_ns_carbon->GetYaxis()->SetTitle("Ns");
h2_nf_ns_carbon->GetXaxis()->CenterTitle(true);
h2_nf_ns_carbon->GetYaxis()->CenterTitle(true);
h2_nf_ns_carbon->GetYaxis()->SetLabelSize(0.045);
h2_nf_ns_carbon->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_nf_ns_boron;
sprintf(hist_name, "Nf vs Ns for p2p events, boron fragment");
h2_nf_ns_boron = new TH2F(hist_name,hist_name,500,0,200000,500,0,200000);
h2_nf_ns_boron->GetXaxis()->SetTitle("Nf");
h2_nf_ns_boron->GetYaxis()->SetTitle("Ns");
h2_nf_ns_boron->GetXaxis()->CenterTitle(true);
h2_nf_ns_boron->GetYaxis()->CenterTitle(true);
h2_nf_ns_boron->GetYaxis()->SetLabelSize(0.045);
h2_nf_ns_boron->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_nf_ns_light;
sprintf(hist_name, "Nf vs Ns for p2p events, light fragment (charge < 4.5)");
h2_nf_ns_light = new TH2F(hist_name,hist_name,500,0,200000,500,0,200000);
h2_nf_ns_light->GetXaxis()->SetTitle("Nf");
h2_nf_ns_light->GetYaxis()->SetTitle("Ns");
h2_nf_ns_light->GetXaxis()->CenterTitle(true);
h2_nf_ns_light->GetYaxis()->CenterTitle(true);
h2_nf_ns_light->GetYaxis()->SetLabelSize(0.045);
h2_nf_ns_light->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_entries_califacal;
sprintf(hist_name, "entries califacal for p2p events");
h1_entries_califacal = new TH1F(hist_name,hist_name,300,0,300);
h1_entries_califacal->GetXaxis()->SetTitle("CalifaCal entries");
h1_entries_califacal->GetYaxis()->SetTitle("Counts");
h1_entries_califacal->GetXaxis()->CenterTitle(true);
h1_entries_califacal->GetYaxis()->CenterTitle(true);
h1_entries_califacal->GetYaxis()->SetLabelSize(0.045);
h1_entries_califacal->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_califacal_vs_opening_angle;
sprintf(hist_name, "Entries califacal vs opening angle for p2p events");
h2_califacal_vs_opening_angle = new TH2F(hist_name,hist_name,52,22.15,152.15,150,0,150);
h2_califacal_vs_opening_angle->GetXaxis()->SetTitle("Opening Angle");
h2_califacal_vs_opening_angle->GetYaxis()->SetTitle("Entries califacal");
h2_califacal_vs_opening_angle->GetXaxis()->CenterTitle(true);
h2_califacal_vs_opening_angle->GetYaxis()->CenterTitle(true);
h2_califacal_vs_opening_angle->GetYaxis()->SetLabelSize(0.045);
h2_califacal_vs_opening_angle->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_esum_vs_opening_angle;
sprintf(hist_name, "Esum vs opening angle for p2p events");
h2_esum_vs_opening_angle = new TH2F(hist_name,hist_name,52,22.15,152.15,700,0,700);
h2_esum_vs_opening_angle->GetXaxis()->SetTitle("Opening Angle");
h2_esum_vs_opening_angle->GetYaxis()->SetTitle("Esum");
h2_esum_vs_opening_angle->GetXaxis()->CenterTitle(true);
h2_esum_vs_opening_angle->GetYaxis()->CenterTitle(true);
h2_esum_vs_opening_angle->GetYaxis()->SetLabelSize(0.045);
h2_esum_vs_opening_angle->GetYaxis()->SetTitleSize(0.045);



TH1F* h1_charge_fragment_p2p;
sprintf(hist_name, "Charge for fragment for p2p events");
h1_charge_fragment_p2p = new TH1F(hist_name,hist_name,100,0,10);
h1_charge_fragment_p2p->GetXaxis()->SetTitle("Charge");
h1_charge_fragment_p2p->GetYaxis()->SetTitle("Counts");
h1_charge_fragment_p2p->GetXaxis()->CenterTitle(true);
h1_charge_fragment_p2p->GetYaxis()->CenterTitle(true);
h1_charge_fragment_p2p->GetYaxis()->SetLabelSize(0.045);
h1_charge_fragment_p2p->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_charge_vs_opening_angle;
sprintf(hist_name, "Charge vs opening angle for p2p events");
h2_charge_vs_opening_angle = new TH2F(hist_name,hist_name,100,0,10,52,22.15,152.15);
h2_charge_vs_opening_angle->GetXaxis()->SetTitle("Charge");
h2_charge_vs_opening_angle->GetYaxis()->SetTitle("Opening Angle");
h2_charge_vs_opening_angle->GetXaxis()->CenterTitle(true);
h2_charge_vs_opening_angle->GetYaxis()->CenterTitle(true);
h2_charge_vs_opening_angle->GetYaxis()->SetLabelSize(0.045);
h2_charge_vs_opening_angle->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_charge_vs_opening_angle_ns_nf_cut;
sprintf(hist_name, "Charge vs opening angle for p2p events for nf and ns positive");
h2_charge_vs_opening_angle_ns_nf_cut = new TH2F(hist_name,hist_name,100,0,10,52,22.15,152.15);
h2_charge_vs_opening_angle_ns_nf_cut->GetXaxis()->SetTitle("Charge");
h2_charge_vs_opening_angle_ns_nf_cut->GetYaxis()->SetTitle("Opening Angle");
h2_charge_vs_opening_angle_ns_nf_cut->GetXaxis()->CenterTitle(true);
h2_charge_vs_opening_angle_ns_nf_cut->GetYaxis()->CenterTitle(true);
h2_charge_vs_opening_angle_ns_nf_cut->GetYaxis()->SetLabelSize(0.045);
h2_charge_vs_opening_angle_ns_nf_cut->GetYaxis()->SetTitleSize(0.045);

TH1I* h1_time_diff_p2p;
sprintf(hist_name, "Time difference of CalifaHit - Event, for  p2p events");
h1_time_diff_p2p = new TH1I(hist_name,hist_name,400,-10000,10000);
h1_time_diff_p2p->GetXaxis()->SetTitle("Time diff [ns]");
h1_time_diff_p2p->GetYaxis()->SetTitle("Counts");
h1_time_diff_p2p->GetXaxis()->CenterTitle(true);
h1_time_diff_p2p->GetYaxis()->CenterTitle(true);
h1_time_diff_p2p->GetYaxis()->SetLabelSize(0.045);
h1_time_diff_p2p->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_time_diff_califa_vs_opang;
sprintf(hist_name, "Time difference between two califa hits vs opening angle for p2p events for nf and ns positive");
h2_time_diff_califa_vs_opang = new TH2F(hist_name,hist_name,400,0,4000,52,22.15,152.15);
h2_time_diff_califa_vs_opang->GetXaxis()->SetTitle("Time diff between two califa hits");
h2_time_diff_califa_vs_opang->GetYaxis()->SetTitle("Opening Angle");
h2_time_diff_califa_vs_opang->GetXaxis()->CenterTitle(true);
h2_time_diff_califa_vs_opang->GetYaxis()->CenterTitle(true);
h2_time_diff_califa_vs_opang->GetYaxis()->SetLabelSize(0.045);
h2_time_diff_califa_vs_opang->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_time_diff_califa_vs_charge;
sprintf(hist_name, "Time difference between two califa hits vs charge of fragment for p2p events for nf and ns positive");
h2_time_diff_califa_vs_charge = new TH2F(hist_name,hist_name,400,0,4000,100,0,10);
h2_time_diff_califa_vs_charge->GetXaxis()->SetTitle("Charge");
h2_time_diff_califa_vs_charge->GetYaxis()->SetTitle("Opening Angle");
h2_time_diff_califa_vs_charge->GetXaxis()->CenterTitle(true);
h2_time_diff_califa_vs_charge->GetYaxis()->CenterTitle(true);
h2_time_diff_califa_vs_charge->GetYaxis()->SetLabelSize(0.045);
h2_time_diff_califa_vs_charge->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_theta_vs_e;
sprintf(hist_name, "Theta vs kin. Energy of protonsfor p2p events for nf and ns positive");
h2_theta_vs_e = new TH2F(hist_name,hist_name,52,22.15,152.15,600,0,600);
h2_theta_vs_e->GetXaxis()->SetTitle("Theta Angle [degr]");
h2_theta_vs_e->GetYaxis()->SetTitle("Energy, MeV");
h2_theta_vs_e->GetXaxis()->CenterTitle(true);
h2_theta_vs_e->GetYaxis()->CenterTitle(true);
h2_theta_vs_e->GetYaxis()->SetLabelSize(0.045);
h2_theta_vs_e->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_p1x_vs_p2x;
sprintf(hist_name, "Px1 vs Px2 for p2p events for nf and ns positive");
h2_p1x_vs_p2x = new TH2F(hist_name,hist_name,350,-700,700,350,-700,700);
h2_p1x_vs_p2x->GetXaxis()->SetTitle("Px1");
h2_p1x_vs_p2x->GetYaxis()->SetTitle("Px2");
h2_p1x_vs_p2x->GetXaxis()->CenterTitle(true);
h2_p1x_vs_p2x->GetYaxis()->CenterTitle(true);
h2_p1x_vs_p2x->GetYaxis()->SetLabelSize(0.045);
h2_p1x_vs_p2x->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_p1y_vs_p2y;
sprintf(hist_name, "Py1 vs Py2 for p2p events for nf and ns positive");
h2_p1y_vs_p2y = new TH2F(hist_name,hist_name,350,-700,700,350,-700,700);
h2_p1y_vs_p2y->GetXaxis()->SetTitle("Py1");
h2_p1y_vs_p2y->GetYaxis()->SetTitle("Py2");
h2_p1y_vs_p2y->GetXaxis()->CenterTitle(true);
h2_p1y_vs_p2y->GetYaxis()->CenterTitle(true);
h2_p1y_vs_p2y->GetYaxis()->SetLabelSize(0.045);
h2_p1y_vs_p2y->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_p1z_vs_p2z;
sprintf(hist_name, "Pz1 vs Pz2 for p2p events for nf and ns positive");
h2_p1z_vs_p2z = new TH2F(hist_name,hist_name,280,0,1400,280,0,1400);
h2_p1z_vs_p2z->GetXaxis()->SetTitle("Pz1");
h2_p1z_vs_p2z->GetYaxis()->SetTitle("Pz2");
h2_p1z_vs_p2z->GetXaxis()->CenterTitle(true);
h2_p1z_vs_p2z->GetYaxis()->CenterTitle(true);
h2_p1z_vs_p2z->GetYaxis()->SetLabelSize(0.045);
h2_p1z_vs_p2z->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_e1_vs_e2;
sprintf(hist_name, "E1 vs E2 for p2p events for nf and ns positive");
h2_e1_vs_e2 = new TH2F(hist_name,hist_name,175,0,700,175,0,700);
h2_e1_vs_e2->GetXaxis()->SetTitle("E1");
h2_e1_vs_e2->GetYaxis()->SetTitle("E2");
h2_e1_vs_e2->GetXaxis()->CenterTitle(true);
h2_e1_vs_e2->GetYaxis()->CenterTitle(true);
h2_e1_vs_e2->GetYaxis()->SetLabelSize(0.045);
h2_e1_vs_e2->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_charge_vs_max_angle;
sprintf(hist_name, "Charge vs max angle for p2p and nf and  ns positive");
h2_charge_vs_max_angle = new TH2F(hist_name,hist_name,52,22.15,152.15,100,0,10);
h2_charge_vs_max_angle->GetXaxis()->SetTitle("Max theta_angle");
h2_charge_vs_max_angle->GetYaxis()->SetTitle("Charge");
h2_charge_vs_max_angle->GetXaxis()->CenterTitle(true);
h2_charge_vs_max_angle->GetYaxis()->CenterTitle(true);
h2_charge_vs_max_angle->GetYaxis()->SetLabelSize(0.045);
h2_charge_vs_max_angle->GetYaxis()->SetTitleSize(0.045);


TH1F* h1_charge_fragment_p2p_low_califa_mult;
sprintf(hist_name, "Charge for fragment for p2p events if califa mult < 5 for both protons");
h1_charge_fragment_p2p_low_califa_mult = new TH1F(hist_name,hist_name,100,0,10);
h1_charge_fragment_p2p_low_califa_mult->GetXaxis()->SetTitle("Charge");
h1_charge_fragment_p2p_low_califa_mult->GetYaxis()->SetTitle("Counts");
h1_charge_fragment_p2p_low_califa_mult->GetXaxis()->CenterTitle(true);
h1_charge_fragment_p2p_low_califa_mult->GetYaxis()->CenterTitle(true);
h1_charge_fragment_p2p_low_califa_mult->GetYaxis()->SetLabelSize(0.045);
h1_charge_fragment_p2p_low_califa_mult->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_charge_vs_opening_angle_ns_nf_cut_cal_low_mult;
sprintf(hist_name, "Charge vs opening angle for p2p events for nf and ns positive if califa mult < 5 for both protons");
h2_charge_vs_opening_angle_ns_nf_cut_cal_low_mult = new TH2F(hist_name,hist_name,100,0,10,52,22.15,152.15);
h2_charge_vs_opening_angle_ns_nf_cut_cal_low_mult->GetXaxis()->SetTitle("Charge");
h2_charge_vs_opening_angle_ns_nf_cut_cal_low_mult->GetYaxis()->SetTitle("Opening Angle");
h2_charge_vs_opening_angle_ns_nf_cut_cal_low_mult->GetXaxis()->CenterTitle(true);
h2_charge_vs_opening_angle_ns_nf_cut_cal_low_mult->GetYaxis()->CenterTitle(true);
h2_charge_vs_opening_angle_ns_nf_cut_cal_low_mult->GetYaxis()->SetLabelSize(0.045);
h2_charge_vs_opening_angle_ns_nf_cut_cal_low_mult->GetYaxis()->SetTitleSize(0.045);


TH2F* h2_charge_vs_esum;
sprintf(hist_name, "Charge vs energy Sum for p2p events for nf and ns positive");
h2_charge_vs_esum = new TH2F(hist_name,hist_name,100,0,10,700,0,700);
h2_charge_vs_esum->GetXaxis()->SetTitle("Charge");
h2_charge_vs_esum->GetYaxis()->SetTitle("Esum");
h2_charge_vs_esum->GetXaxis()->CenterTitle(true);
h2_charge_vs_esum->GetYaxis()->CenterTitle(true);
h2_charge_vs_esum->GetYaxis()->SetLabelSize(0.045);
h2_charge_vs_esum->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_theta_vs_phi;
sprintf(hist_name, "Theta vs Phi");
h2_theta_vs_phi = new TH2F(hist_name,hist_name,450,0,90,720,-180,180);
h2_theta_vs_phi->GetXaxis()->SetTitle("Theta [°]");
h2_theta_vs_phi->GetYaxis()->SetTitle("Phi [°]");
h2_theta_vs_phi->GetXaxis()->CenterTitle(true);
h2_theta_vs_phi->GetYaxis()->CenterTitle(true);
h2_theta_vs_phi->GetYaxis()->SetLabelSize(0.045);
h2_theta_vs_phi->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_theta_vs_phi_energy[10];
for (int i= 0; i < 10; i++){
sprintf(hist_name, "Theta vs Phi vs Energy, example %i",i+1);
h2_theta_vs_phi_energy[i] = new TH2F(hist_name,hist_name,52,22.15,152.15,60,-180,180);
h2_theta_vs_phi_energy[i]->GetXaxis()->SetTitle("Theta [°]");
h2_theta_vs_phi_energy[i]->GetYaxis()->SetTitle("Phi [°]");
h2_theta_vs_phi_energy[i]->GetXaxis()->CenterTitle(true);
h2_theta_vs_phi_energy[i]->GetYaxis()->CenterTitle(true);
h2_theta_vs_phi_energy[i]->GetYaxis()->SetLabelSize(0.045);
h2_theta_vs_phi_energy[i]->GetYaxis()->SetTitleSize(0.045);
}

//END OF HISTOGRAMS


if (simulation){

Long64_t nevents = chain->GetEntries(); //number of events in file with name "fname"
for(Long64_t i= 0; i < nevents;i++){
    Long64_t evtnr = i;

    if (i%100000==0)
        cout<<"Processing event "<<i<<endl;

	chain->GetEvent(i); //event_i_call
	//Int_t tpatbin = (DataCA->GetTpat() & (1 << 0));
	//if (tpatbin != 0){
	n_event_min_bias += 1;	
	Int_t entries_califa = CalifaHitData->GetEntries();
	Int_t entries_califa_cal = CalifaCalData->GetEntries();
	
	califahitdata = new R3BCalifaHitData*[entries_califa];
    vector<TLorentzVector> tl_vec_protons_lab;
    vector<TLorentzVector> tl_vec_protons_cms;
	vector<UInt_t> vec_numb_califa;
	vector<vector<double> > nf_ns;
	if(entries_califa > 1) n_evt_two_hits_sim +=1;
	for (Int_t j = 0; j < entries_califa; j++){
		califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);
		//Double_t energy_lab = (califahitdata[j]->GetEnergy())/1000.;
		Double_t energy_lab = (califahitdata[j]->GetEnergy());
		//if (energy_lab > 30){  //PROTONS
		if (energy_lab > 0.030){  //PROTONS
			UInt_t califa_numbers = califahitdata[j]->GetNbOfCrystalHits();
            vec_numb_califa.push_back(califa_numbers);
			TLorentzVector temp_tl_vec_proton;
			Double_t theta =  califahitdata[j]->GetTheta(); 
			Double_t phi = califahitdata[j]->GetPhi(); 
			Double_t total_momentum = sqrt(pow(energy_lab,2)+2*0.001*MASS_PROTON*energy_lab);
			Double_t px = total_momentum*sin(theta)*cos(phi);	
			Double_t py = total_momentum*sin(theta)*sin(phi);
			Double_t pz = total_momentum*cos(theta);
			Double_t total_energy = energy_lab + 0.001*MASS_PROTON;
			temp_tl_vec_proton.SetPxPyPzE(px,py,pz,total_energy);
			tl_vec_protons_lab.push_back(temp_tl_vec_proton);
			temp_tl_vec_proton.Boost(0,0,-beta_given);  //TODO: find boost vector
			tl_vec_protons_cms.push_back(temp_tl_vec_proton);
			vector<double> temp_nf_ns;
            temp_nf_ns.push_back(califahitdata[j]->GetNf());
            temp_nf_ns.push_back(califahitdata[j]->GetNs());
            nf_ns.push_back(temp_nf_ns);
			}
		}
	sort(tl_vec_protons_lab.begin(),tl_vec_protons_lab.end(),sortcol_lorentz);
	if (tl_vec_protons_lab.size() > 1 && abs((tl_vec_protons_lab[0]).Phi() - (tl_vec_protons_lab[1]).Phi()) < 3.8397243544 && abs((tl_vec_protons_lab[0]).Phi() - (tl_vec_protons_lab[1]).Phi()) > 2.4434609528) n_evt_two_hits_sim_protons +=1;
	if (tl_vec_protons_lab.size() == 2){
		Double_t diff_phi_angle_degr = (abs((tl_vec_protons_lab[0]).Phi() - (tl_vec_protons_lab[1]).Phi()))/PI*180.;
		//cout << "Azimuthal angular difference:\t" << diff_phi_angle_degr << endl;
		h1_azimuthal_ang_diff->Fill(diff_phi_angle_degr);
		}

	//if (entries_califa == 2 && tl_vec_protons_lab.size() == 2 && abs((tl_vec_protons_lab[0]).Phi() - (tl_vec_protons_lab[1]).Phi()) < 3.8397243544 && abs((tl_vec_protons_lab[0]).Phi() - (tl_vec_protons_lab[1]).Phi()) > 2.4434609528 && (tl_vec_protons_lab[0].E()-0.001*MASS_PROTON) > 0.070 && (tl_vec_protons_lab[0].E() - 0.001*MASS_PROTON) < 0.320 && (tl_vec_protons_lab[1].E() - 0.001*MASS_PROTON) > 0.070 && (tl_vec_protons_lab[1].E() - 0.001*MASS_PROTON) < 0.320){

	//now set just min limit for energy: 50 MeV, no upper limit
	if (entries_califa == 2 && tl_vec_protons_lab.size() == 2 && abs((tl_vec_protons_lab[0]).Phi() - (tl_vec_protons_lab[1]).Phi()) < 3.8397243544 && abs((tl_vec_protons_lab[0]).Phi() - (tl_vec_protons_lab[1]).Phi()) > 2.4434609528 && (tl_vec_protons_lab[0].E()-0.001*MASS_PROTON) > 0.070 && (tl_vec_protons_lab[1].E() - 0.001*MASS_PROTON) > 0.070 ){

			//get the opening angle
			Double_t cos_opening_angle_two_protons = sin(tl_vec_protons_lab[0].Theta())*sin(tl_vec_protons_lab[1].Theta())*cos(tl_vec_protons_lab[0].Phi()-tl_vec_protons_lab[1].Phi())+cos(tl_vec_protons_lab[0].Theta())*cos(tl_vec_protons_lab[1].Theta());
			Double_t opening_angle = acos(cos_opening_angle_two_protons)*180/PI;
			h1_opang_p2p->Fill(opening_angle);
			h1_p2p_proton_energies->Fill(1000*(tl_vec_protons_lab[0].E() -0.001*MASS_PROTON));
			h1_p2p_proton_energies->Fill(1000*(tl_vec_protons_lab[1].E() -0.001*MASS_PROTON));
			h1_num_califa_hits->Fill(vec_numb_califa[0]);
            h1_num_califa_hits->Fill(vec_numb_califa[1]);
			h2_nf_ns->Fill(1000000*nf_ns[0][0],1000000*nf_ns[0][1]);
            h2_nf_ns->Fill(1000000*nf_ns[1][0],1000000*nf_ns[1][1]);
			h1_e_sum_p2p->Fill(1000*((tl_vec_protons_lab[0].E()) - 0.001*MASS_PROTON+(tl_vec_protons_lab[1].E()) - 0.001*MASS_PROTON));
			h1_entries_califacal->Fill(entries_califa_cal);
			h2_califacal_vs_opening_angle->Fill(opening_angle,entries_califa_cal);
			h2_esum_vs_opening_angle->Fill(opening_angle,1000*((tl_vec_protons_lab[0].E()) - 0.001*MASS_PROTON+(tl_vec_protons_lab[1].E()) - 0.001*MASS_PROTON));

			//list the number of p2p-events
			if ( opening_angle > 80. && opening_angle < 90.){
				n_event_califa_p2p += 1;
				}
				//randomize the two protons
				double r = ((double) rand() / (RAND_MAX));
                Int_t first_index;
                first_index = (r > 0.5) ? 1 : 0;
                Int_t second_index = 1 - first_index;
                vector<TLorentzVector> tl_vec_protons_lab_rand;
                tl_vec_protons_lab_rand.push_back(tl_vec_protons_lab[first_index]);
                tl_vec_protons_lab_rand.push_back(tl_vec_protons_lab[second_index]);
                h2_theta_vs_e->Fill((tl_vec_protons_lab_rand[0].Theta())/PI*180,1000*((tl_vec_protons_lab_rand[0].E()) - 0.001*MASS_PROTON));
                h2_theta_vs_e->Fill((tl_vec_protons_lab_rand[1].Theta())/PI*180,1000*((tl_vec_protons_lab_rand[1].E()) - 0.001*MASS_PROTON));
                h2_p1x_vs_p2x->Fill(1000*tl_vec_protons_lab_rand[0].Px(),1000*tl_vec_protons_lab_rand[1].Px());
                h2_p1y_vs_p2y->Fill(1000*tl_vec_protons_lab_rand[0].Py(),1000*tl_vec_protons_lab_rand[1].Py());
                h2_p1z_vs_p2z->Fill(1000*tl_vec_protons_lab_rand[0].Pz(),1000*tl_vec_protons_lab_rand[1].Pz());
                h2_e1_vs_e2->Fill(1000*((tl_vec_protons_lab_rand[0].E()) - 0.001*MASS_PROTON),1000*((tl_vec_protons_lab_rand[1].E()) - 0.001*MASS_PROTON));	
				cout << "Px" << "   " << tl_vec_protons_lab_rand[0].Px() << "   " << " Py" << "   " << tl_vec_protons_lab_rand[0].Py() << "  " << "Pz" << "  " <<tl_vec_protons_lab_rand[0].Pz() << endl;
				cout << "mass:\t" << tl_vec_protons_lab_rand[0].M() << endl;
				cout << "energy:\t" << tl_vec_protons_lab_rand[0].E() << endl;
				cout << "momentum:\t" << sqrt(pow(tl_vec_protons_lab_rand[0].Px(),2)+pow(tl_vec_protons_lab_rand[0].Py(),2)+pow(tl_vec_protons_lab_rand[0].Pz(),2)) << endl;
			//inclusive number of p2p events, no cut on opening angle
			//n_event_califa_p2p_inclusive +=1;
			//update n_event_califa_p2p_inclusive: select only events with nf && ns positive
			//if(nf_ns[0][0]> 0 && nf_ns[0][1] > 0 && nf_ns[1][0]> 0 && nf_ns[1][1]> 0 && opening_angle < 90 ) n_event_califa_p2p_inclusive +=1;
			//no cut on opening angle
			if(nf_ns[0][0]> 0 && nf_ns[0][1] > 0 && nf_ns[1][0]> 0 && nf_ns[1][1]> 0) n_event_califa_p2p_inclusive +=1;
			}
		
		delete [] califahitdata;
		}
		//}
		//
		cout << "Events with at least 2 hits in CALIFA:\t" << n_evt_two_hits_sim << endl;
		cout << "Events with at least 2 hits more than 30 MEV in Califa:\t" << n_evt_two_hits_sim_protons << endl;
	int64_t beam_energy;
	//if ((string(count_i)).compare("sim_400") == 0) beam_energy = 400;
	//if ((string(count_i)).compare("sim_550") == 0) beam_energy = 550;
	//if ((string(count_i)).compare("sim_650") == 0) beam_energy = 650;
	//if ((string(count_i)).compare("sim_800") == 0) beam_energy = 800;

	if ((string(count_i)).find("sim_400_CH2Target_1229") != string::npos) beam_energy = 400,target_type = "ch2_1229";
	if ((string(count_i)).find("sim_550_CH2Target_1229") != string::npos) beam_energy = 550,target_type = "ch2_1229";
	if ((string(count_i)).find("sim_650_CH2Target_1229") != string::npos) beam_energy = 650,target_type = "ch2_1229";
	if ((string(count_i)).find("sim_800_CH2Target_1229") != string::npos) beam_energy = 800,target_type = "ch2_1229";
	if ((string(count_i)).find("sim_400_CH2Target_2453") != string::npos) beam_energy = 400,target_type = "ch2_2453";
	if ((string(count_i)).find("sim_550_CH2Target_2453") != string::npos) beam_energy = 550,target_type = "ch2_2453";
	if ((string(count_i)).find("sim_650_CH2Target_2453") != string::npos) beam_energy = 650,target_type = "ch2_2453";
	if ((string(count_i)).find("sim_800_CH2Target_2453") != string::npos) beam_energy = 800,target_type = "ch2_2453";
	if ((string(count_i)).find("sim_400_CTarget_54") != string::npos) beam_energy = 400,target_type = "c_54";
	if ((string(count_i)).find("sim_550_CTarget_54") != string::npos) beam_energy = 550,target_type = "c_54";
	if ((string(count_i)).find("sim_650_CTarget_54") != string::npos) beam_energy = 650,target_type = "c_54";
	if ((string(count_i)).find("sim_800_CTarget_54") != string::npos) beam_energy = 800,target_type = "c_54";
	if ((string(count_i)).find("sim_400_CTarget_1086") != string::npos) beam_energy = 400,target_type = "c_1086";
	if ((string(count_i)).find("sim_550_CTarget_1086") != string::npos) beam_energy = 550,target_type = "c_1086";
	if ((string(count_i)).find("sim_650_CTarget_1086") != string::npos) beam_energy = 650,target_type = "c_1086";
	if ((string(count_i)).find("sim_800_CTarget_1086") != string::npos) beam_energy = 800,target_type = "c_1086";
	if ((string(count_i)).find("sim_400_CTarget_2198") != string::npos) beam_energy = 400,target_type = "c_2198";
	if ((string(count_i)).find("sim_550_CTarget_2198") != string::npos) beam_energy = 550,target_type = "c_2198";
	if ((string(count_i)).find("sim_650_CTarget_2198") != string::npos) beam_energy = 650,target_type = "c_2198";
	if ((string(count_i)).find("sim_800_CTarget_2198") != string::npos) beam_energy = 800,target_type = "c_2198";
	
	myfile << "sim" << "\t" << target_type << "\t" <<  "NAN" << "\t" << beam_energy << "\t" << n_event_min_bias << "\t" << n_event_califa_p2p_inclusive << endl;
	myfile.close();
mystream << "sim" << "\t" << target_type << "\t" << "NAN" << "\t" << beam_energy << "\t" << n_event_min_bias << "\t" << n_event_califa_p2p_inclusive << "\n";
mystream.close();
	}
    
else{
R3BEventHeader* DataCA = new R3BEventHeader();
TBranch* branchData = chain->GetBranch("EventHeader.");
branchData->SetAddress(&DataCA);

TClonesArray* SofSciTcalData = new TClonesArray("R3BSofSciTcalData",2);
R3BSofSciTcalData** sofscitcaldata;
TBranch *branchSofSciTcalData = chain->GetBranch("SofSciTcalData");
branchSofSciTcalData->SetAddress(&SofSciTcalData);

TClonesArray* TwimHitData = new TClonesArray("R3BTwimHitData",2);
R3BTwimHitData** softwimhitdata;
TBranch *branchTwimHitData = chain->GetBranch("TwimHitData");
branchTwimHitData->SetAddress(&TwimHitData);

TClonesArray* SofTofWTcalData = new TClonesArray("R3BSofTofWTcalData",2);
R3BSofTofWTcalData** softofwtcaldata;
TBranch *branchSofTofWTcalData = chain->GetBranch("SofTofWTcalData");
branchSofTofWTcalData->SetAddress(&SofTofWTcalData);


Long64_t nevents = chain->GetEntries(); //number of events in file with name "fname"
//change here a little bit...
cout << "here I changed a little bit..." << endl;
cout << "size of iterating vector:\t" << vec_interating_loop.size() << endl;
if (vec_interating_loop.size() == 0){
    cout << "I do not have a file with interesting events, this will now take more time..." << endl;
        for(Long64_t i=0; i< nevents;i++){
				if (i%100000==0)
        			cout<<"Processing vector filling "<<i<<endl;
					vec_interating_loop.push_back(i);
//				chain->GetEvent(i); //event_i_call
//				Int_t entries_start = SofSciTcalData->GetEntriesFast();
//				if (entries_start ==2) vec_interating_loop.push_back(i);
	//			Int_t tpatbin;
	//			bool real_event=false;
	//			for (Int_t j=0; j < 16; j++){
	//				tpatbin = (DataCA->GetTpat() & (1 << j));
	//				if (tpatbin != 0 && j != 8 && j != 9) real_event=true;
	//				}
	//				if (real_event){
    //			            vec_interating_loop.push_back(i);
	//						}
							
    			                    }
                            }
for(Long64_t i : vec_interating_loop){
    Long64_t evtnr = i;
    if (i%100000==0)
        cout<<"Processing event "<<i<<endl;

	chain->GetEvent(i); //event_i_call
	//Int_t tpatbin;
	//bool real_event=false;
	//for (Int_t j=0; j < 16; j++){
	//	tpatbin = (DataCA->GetTpat() & (1 << j));
	//	if (tpatbin != 0 && j != 8 && j != 9) real_event=true;
	//	}
	//cout << "this is the eventtime:\t" << DataCA->GetTimeStamp() << endl;	
	softwimhitdata = new R3BTwimHitData*[1];
    softwimhitdata[0] = (R3BTwimHitData*)TwimHitData->At(0);
	
	Int_t entries_start = SofSciTcalData->GetEntriesFast();
	Int_t entries_tofw = SofTofWTcalData->GetEntriesFast();
	if (entries_start ==2){
	n_event_min_bias += 1;	
	//cout << "this is event:\t" << evtnr << endl;
	//cout << "this is tpat:\t" << DataCA->GetTpat() << endl;
	Int_t entries_califa = CalifaHitData->GetEntries();
	Int_t entries_califa_cal = CalifaCalData->GetEntries();
	
	califahitdata = new R3BCalifaHitData*[entries_califa];
    vector<TLorentzVector> tl_vec_protons_lab;
    vector<TLorentzVector> tl_vec_protons_cms;
	vector<UInt_t> vec_numb_califa;
	vector<vector<double> > nf_ns;
	vector<int64_t> califa_time;
	Int_t example_p2p_events = 0;

	for (Int_t j = 0; j < entries_califa; j++){
		califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);
		Double_t energy_lab = (califahitdata[j]->GetEnergy())/1000.;
		if (energy_lab > 30){  //PROTONS
			UInt_t califa_numbers = califahitdata[j]->GetNbOfCrystalHits();
			vec_numb_califa.push_back(califa_numbers);	
			TLorentzVector temp_tl_vec_proton;
			Double_t theta =  califahitdata[j]->GetTheta(); 
			Double_t phi = califahitdata[j]->GetPhi(); 
			h2_theta_vs_phi->Fill(theta/PI*180,phi/PI*180);
			Double_t total_momentum = sqrt(pow(energy_lab,2)+2*MASS_PROTON*energy_lab);
			Double_t px = total_momentum*sin(theta)*cos(phi);	
			Double_t py = total_momentum*sin(theta)*sin(phi);
			Double_t pz = total_momentum*cos(theta);
			Double_t total_energy = energy_lab + MASS_PROTON;
			temp_tl_vec_proton.SetPxPyPzE(px,py,pz,total_energy);
			tl_vec_protons_lab.push_back(temp_tl_vec_proton);
			temp_tl_vec_proton.Boost(0,0,-beta_given);  //TODO: find boost vector
			tl_vec_protons_cms.push_back(temp_tl_vec_proton);
			vector<double> temp_nf_ns;
			temp_nf_ns.push_back(califahitdata[j]->GetNf());
			temp_nf_ns.push_back(califahitdata[j]->GetNs());
			nf_ns.push_back(temp_nf_ns);
			califa_time.push_back(int64_t(califahitdata[j]->GetTime()-DataCA->GetTimeStamp())); 
			}
			sort(tl_vec_protons_lab.begin(),tl_vec_protons_lab.end(),sortcol_lorentz);
		}
	if (tl_vec_protons_lab.size() == 2){
        Double_t diff_phi_angle_degr = (abs((tl_vec_protons_lab[0]).Phi() - (tl_vec_protons_lab[1]).Phi()))/PI*180.;
        //cout << "Azimuthal angular difference:\t" << diff_phi_angle_degr << endl;
        h1_azimuthal_ang_diff->Fill(diff_phi_angle_degr);
        }
		//cut on the p2p events: you assume two protons, and azimuthal ang. difference = 180 +-15° && energy > 70 and smaller 320 MeV
	//if (entries_califa == 2 && tl_vec_protons_lab.size() == 2 && abs((tl_vec_protons_lab[0]).Phi() - (tl_vec_protons_lab[1]).Phi()) < 3.8397243544 && abs((tl_vec_protons_lab[0]).Phi() - (tl_vec_protons_lab[1]).Phi()) > 2.4434609528 && (tl_vec_protons_lab[0].E() - MASS_PROTON) > 70 && (tl_vec_protons_lab[0].E() - MASS_PROTON) < 320 && (tl_vec_protons_lab[1].E() - MASS_PROTON) > 70 && (tl_vec_protons_lab[1].E() - MASS_PROTON) < 320){

	//now set just min limit for energy: 50 MeV, no upper limit
	if (entries_califa == 2 && tl_vec_protons_lab.size() == 2 && abs((tl_vec_protons_lab[0]).Phi() - (tl_vec_protons_lab[1]).Phi()) < 3.8397243544 && abs((tl_vec_protons_lab[0]).Phi() - (tl_vec_protons_lab[1]).Phi()) > 2.4434609528 && (tl_vec_protons_lab[0].E() - MASS_PROTON) > 70 && (tl_vec_protons_lab[1].E() - MASS_PROTON) > 70){

			v_events.push_back(evtnr);
			//get the opening angle
			Double_t cos_opening_angle_two_protons = sin(tl_vec_protons_lab[0].Theta())*sin(tl_vec_protons_lab[1].Theta())*cos(tl_vec_protons_lab[0].Phi()-tl_vec_protons_lab[1].Phi())+cos(tl_vec_protons_lab[0].Theta())*cos(tl_vec_protons_lab[1].Theta());
			Double_t opening_angle = acos(cos_opening_angle_two_protons)*180/PI;
			h1_opang_p2p->Fill(opening_angle);
			h1_p2p_proton_energies->Fill((tl_vec_protons_lab[0].E()) - MASS_PROTON);
			h1_p2p_proton_energies->Fill((tl_vec_protons_lab[1].E()) - MASS_PROTON);
			h1_num_califa_hits->Fill(vec_numb_califa[0]);
			h1_num_califa_hits->Fill(vec_numb_califa[1]);

			h2_nf_ns->Fill(nf_ns[0][0],nf_ns[0][1]);
			h2_nf_ns->Fill(nf_ns[1][0],nf_ns[1][1]);
			h1_e_sum_p2p->Fill((tl_vec_protons_lab[0].E()) - MASS_PROTON+(tl_vec_protons_lab[1].E()) - MASS_PROTON);
			h1_entries_califacal->Fill(entries_califa_cal);
			h2_califacal_vs_opening_angle->Fill(opening_angle,entries_califa_cal);
			h2_esum_vs_opening_angle->Fill(opening_angle,(tl_vec_protons_lab[0].E()) - MASS_PROTON+(tl_vec_protons_lab[1].E()) - MASS_PROTON);
			if (opening_angle > 120) cout << "this is event with opening angle larger 120 degr:\t" << evtnr << endl;

			//list the number of p2p-events
			if ( opening_angle > 70. && opening_angle < 90.){
				n_event_califa_p2p += 1;
				//cout << "event p2p8090:\t" << evtnr << endl;
				}
			//inclusive number of p2p events, no cut on opening angle
			//n_event_califa_p2p_inclusive +=1;
			//update n_event_califa_p2p_inclusive: select only events with nf && ns positive
			if(nf_ns[0][0]> 0 && nf_ns[0][1] > 0 && nf_ns[1][0]> 0 && nf_ns[1][1]> 0 &&  opening_angle  < 90 && vec_numb_califa[0] < 5 && vec_numb_califa[1] < 5){
				//n_event_califa_p2p_inclusive +=1; 

				h1_time_diff_p2p->Fill(califa_time[0]);
				h1_time_diff_p2p->Fill(califa_time[1]);
				h2_time_diff_califa_vs_opang->Fill(abs(califa_time[1]-califa_time[0]),opening_angle);

				//randomize the two protons
				double r = ((double) rand() / (RAND_MAX));
				Int_t first_index;
				first_index = (r > 0.5) ? 1 : 0; 
				Int_t second_index = 1 - first_index;
				vector<TLorentzVector> tl_vec_protons_lab_rand;
				tl_vec_protons_lab_rand.push_back(tl_vec_protons_lab[first_index]); 
				tl_vec_protons_lab_rand.push_back(tl_vec_protons_lab[second_index]); 
				h2_theta_vs_e->Fill((tl_vec_protons_lab_rand[0].Theta())/PI*180,(tl_vec_protons_lab_rand[0].E()) - MASS_PROTON);
				h2_theta_vs_e->Fill((tl_vec_protons_lab_rand[1].Theta())/PI*180,(tl_vec_protons_lab_rand[1].E()) - MASS_PROTON);
				h2_p1x_vs_p2x->Fill(tl_vec_protons_lab_rand[0].Px(),tl_vec_protons_lab_rand[1].Px());
				h2_p1y_vs_p2y->Fill(tl_vec_protons_lab_rand[0].Py(),tl_vec_protons_lab_rand[1].Py());
				h2_p1z_vs_p2z->Fill(tl_vec_protons_lab_rand[0].Pz(),tl_vec_protons_lab_rand[1].Pz());
				h2_e1_vs_e2->Fill((tl_vec_protons_lab_rand[0].E()) - MASS_PROTON,(tl_vec_protons_lab_rand[1].E()) - MASS_PROTON);
				cout << "px1" <<  tl_vec_protons_lab_rand[0].Px() << endl;
				cout << "py1" <<  tl_vec_protons_lab_rand[0].Py() << endl;
				cout << "pz1" <<  tl_vec_protons_lab_rand[0].Pz() << endl;


				}
			//NOW PLOT ONLY OPENING ALGLE IF BOTH NF AND NS positive
			if(nf_ns[0][0]> 0 && nf_ns[0][1] > 0 && nf_ns[1][0]> 0 && nf_ns[1][1]> 0){
				h1_opang_p2p_nfns_positive->Fill(opening_angle);
				}

				//Double_t charge_val = softwimhitdata[0]->GetZcharge();
				Double_t charge_val = 5.;

				h1_charge_fragment_p2p->Fill(charge_val);	
				h2_charge_vs_opening_angle->Fill(charge_val,opening_angle);
				//cut on nf ns positive 
				if (nf_ns[0][0] > 0 && nf_ns[0][1] > 0 && nf_ns[1][0] > 0 && nf_ns[1][1] > 0){
					h2_charge_vs_opening_angle_ns_nf_cut->Fill(charge_val,opening_angle);	
					h2_time_diff_califa_vs_charge->Fill(abs(califa_time[1]-califa_time[0]),charge_val);
					Double_t max_angle_theta;
					max_angle_theta = (((tl_vec_protons_lab[0].Theta())/PI*180) > ((tl_vec_protons_lab[1].Theta())/PI*180)) ? ((tl_vec_protons_lab[0].Theta())/PI*180) : ((tl_vec_protons_lab[1].Theta())/PI*180); 
					h2_charge_vs_max_angle->Fill(max_angle_theta,charge_val);
					h2_charge_vs_esum->Fill(charge_val,(tl_vec_protons_lab[0].E()) - MASS_PROTON +(tl_vec_protons_lab[1].E()) - MASS_PROTON);
					if (vec_numb_califa[0] < 5 && vec_numb_califa[1] < 5){
						h1_charge_fragment_p2p_low_califa_mult->Fill(charge_val);
						h2_charge_vs_opening_angle_ns_nf_cut_cal_low_mult->Fill(charge_val,opening_angle);
						}
					//if (charge_val < 5.5 && opening_angle < 90  ){
					//no opening angle cut...
					if (charge_val < 5.5){
						n_event_califa_p2p_inclusive +=1;
						
						}	
					}
				//cut on charges and look at nf vs ns
				if (charge_val > 5.5){
					h2_nf_ns_carbon->Fill(nf_ns[0][0],nf_ns[0][1]);
					h2_nf_ns_carbon->Fill(nf_ns[1][0],nf_ns[1][1]);
					}
				if (charge_val <5.5 && charge_val > 4.5){
					h2_nf_ns_boron->Fill(nf_ns[0][0],nf_ns[0][1]);
					h2_nf_ns_boron->Fill(nf_ns[1][0],nf_ns[1][1]);
					}
				if (charge_val < 4.5){
					h2_nf_ns_light->Fill(nf_ns[0][0],nf_ns[0][1]);
					h2_nf_ns_light->Fill(nf_ns[1][0],nf_ns[1][1]);
					}
				}
			
		
		delete [] califahitdata;
		}	
		delete [] softwimhitdata;
		}


		//write important numbers to file -->empty runs!
		myfile << "exp" << "\t" << target_type << "\t" << (string(count_i)) << "\t" << beam_energy << "\t" << n_event_min_bias << "\t" << n_event_califa_p2p_inclusive << endl;
		mystream << "exp" << "\t" << target_type << "\t" << (string(count_i)) << "\t" << beam_energy << "\t" << n_event_min_bias << "\t" << n_event_califa_p2p_inclusive << "\n";
		mystream.close();
		myfile.close();

		}


h1_min_bias_events->SetBinContent(1,n_event_min_bias);
h1_p2p_events->SetBinContent(1,n_event_califa_p2p);
h1_p2p_events_inclusive->SetBinContent(1,n_event_califa_p2p_inclusive);
delete CalifaHitData;
if (!eventFile){
    cout << " I try to write to file..." << endl;
    t_events->Fill();
    f_events->Write();
    cout << "ending to write to file" << endl;
	cout << "size of vector:\t" <<v_events.size() << endl;
}
TList *l = new TList();
l->Add(h1_opang_p2p);
l->Add(h1_min_bias_events);
l->Add(h1_p2p_events);
l->Add(h1_p2p_events_inclusive);
l->Add(h1_p2p_proton_energies);
l->Add(h1_azimuthal_ang_diff);
l->Add(h1_num_califa_hits);
l->Add(h2_nf_ns);
l->Add(h1_e_sum_p2p);
l->Add(h1_entries_califacal);
l->Add(h2_califacal_vs_opening_angle);
l->Add(h1_opang_p2p_nfns_positive);
l->Add(h2_esum_vs_opening_angle);
l->Add(h1_charge_fragment_p2p);
l->Add(h2_charge_vs_opening_angle);
l->Add(h2_charge_vs_opening_angle_ns_nf_cut);
l->Add(h2_nf_ns_carbon);
l->Add(h2_nf_ns_boron);
l->Add(h2_nf_ns_light);
l->Add(h1_time_diff_p2p);
l->Add(h2_time_diff_califa_vs_opang);
l->Add(h2_time_diff_califa_vs_charge);
l->Add(h2_theta_vs_e);
l->Add(h2_p1x_vs_p2x);
l->Add(h2_p1y_vs_p2y);
l->Add(h2_p1z_vs_p2z);
l->Add(h2_e1_vs_e2);
l->Add(h2_charge_vs_max_angle);
l->Add(h1_charge_fragment_p2p_low_califa_mult);
l->Add(h2_charge_vs_opening_angle_ns_nf_cut_cal_low_mult);
l->Add(h2_charge_vs_esum);
l->Add(h2_theta_vs_phi);




l->Write("histlist", TObject::kSingleKey);
cout << "end of process, successful" <<endl;
cout << "------------STATISTICS--------" << endl;
cout << "min bias events:\t" << n_event_min_bias << endl;
cout << "p2p-events in CALIFA:\t" << n_event_califa_p2p << endl;
cout << "------------------------------" << endl;
}

