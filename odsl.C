#include <iostream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"
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
void odsl(char const count_i[50]){
char hist_name[500];
char f_out_name[500];
char f_events_name[500];
char fname[500];
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

//PLASTIC-------------
sprintf(f_out_name,"/home/ge37liw/califa_odsl_run_%s.root",count_i);
////////sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/p2p_data_fixed_angle/stitched_and_unpacked_main%s.root",count_i);
sprintf(fname,"/home/ge37liw/single_cluster_stitched_and_unpacked_main%s.root",count_i);
TFile * f = new TFile(f_out_name,"RECREATE");

TChain* chain = new TChain("evt");
chain->Reset();
chain->Add(fname);
//BEGIN TCA

TClonesArray* CalifaHitData = new TClonesArray("R3BCalifaHitData",3);
R3BCalifaHitData** califahitdata;
TBranch *branchCalifaHitData = chain->GetBranch("CalifaHitData");
branchCalifaHitData->SetAddress(&CalifaHitData);

R3BEventHeader* DataCA = new R3BEventHeader();
TBranch* branchData = chain->GetBranch("EventHeader.");
branchData->SetAddress(&DataCA);

TClonesArray* SofSciTcalData = new TClonesArray("R3BSofSciTcalData",2);
R3BSofSciTcalData** sofscitcaldata;
TBranch *branchSofSciTcalData = chain->GetBranch("SofSciTcalData");
branchSofSciTcalData->SetAddress(&SofSciTcalData);

TClonesArray* CalifaCalData = new TClonesArray("R3BCalifaCrystalCalData",3);
R3BCalifaCrystalCalData** califacaltdata;
TBranch *branchCalifaCalData = chain->GetBranch("CalifaCrystalCalData");
branchCalifaCalData->SetAddress(&CalifaCalData);

TClonesArray* SofTofWTcalData = new TClonesArray("R3BSofTofWTcalData",2);
R3BSofTofWTcalData** softofwtcaldata;
TBranch *branchSofTofWTcalData = chain->GetBranch("SofTofWTcalData");
branchSofTofWTcalData->SetAddress(&SofTofWTcalData);

TH2F* h2_theta_vs_phi;
sprintf(hist_name, "Theta vs Phi");
h2_theta_vs_phi = new TH2F(hist_name,hist_name,450,0,90,720,-180,180);
h2_theta_vs_phi->GetXaxis()->SetTitle("Theta [degree]");
h2_theta_vs_phi->GetYaxis()->SetTitle("Phi [degree]");
h2_theta_vs_phi->GetXaxis()->CenterTitle(true);
h2_theta_vs_phi->GetYaxis()->CenterTitle(true);
h2_theta_vs_phi->GetYaxis()->SetLabelSize(0.045);
h2_theta_vs_phi->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_califa_mult;
sprintf(hist_name, "Califa Multiplicity cal level");
h1_califa_mult = new TH1F(hist_name,hist_name,200,0,200);
h1_califa_mult->GetXaxis()->SetTitle("Multiplicity");
h1_califa_mult->GetYaxis()->SetTitle("Count");
h1_califa_mult->GetXaxis()->CenterTitle(true);
h1_califa_mult->GetYaxis()->CenterTitle(true);
h1_califa_mult->GetYaxis()->SetLabelSize(0.045);
h1_califa_mult->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_califa_mult_hit;
sprintf(hist_name, "Califa Multiplicity hit level");
h1_califa_mult_hit = new TH1F(hist_name,hist_name,200,0,200);
h1_califa_mult_hit->GetXaxis()->SetTitle("Multiplicity");
h1_califa_mult_hit->GetYaxis()->SetTitle("Count");
h1_califa_mult_hit->GetXaxis()->CenterTitle(true);
h1_califa_mult_hit->GetYaxis()->CenterTitle(true);
h1_califa_mult_hit->GetYaxis()->SetLabelSize(0.045);
h1_califa_mult_hit->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_energy_gamma;
sprintf(hist_name, "CALIFA energy spectrum in gamma range");
h1_energy_gamma = new TH1F(hist_name,hist_name,300,0,15);
h1_energy_gamma->GetXaxis()->SetTitle("Gamma Energy Spectrum");
h1_energy_gamma->GetYaxis()->SetTitle("Count");
h1_energy_gamma->GetXaxis()->CenterTitle(true);
h1_energy_gamma->GetYaxis()->CenterTitle(true);
h1_energy_gamma->GetYaxis()->SetLabelSize(0.045);
h1_energy_gamma->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_energy_proton;
sprintf(hist_name, "CALIFA energy spectrum in proton range");
h1_energy_proton = new TH1F(hist_name,hist_name,1400,0,700);
h1_energy_proton->GetXaxis()->SetTitle("High Energy Spectrum");
h1_energy_proton->GetYaxis()->SetTitle("Count");
h1_energy_proton->GetXaxis()->CenterTitle(true);
h1_energy_proton->GetYaxis()->CenterTitle(true);
h1_energy_proton->GetYaxis()->SetLabelSize(0.045);
h1_energy_proton->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_califa_time;
sprintf(hist_name, "CALIFA time");
h1_califa_time = new TH1F(hist_name,hist_name,800,-4000,4000);
h1_califa_time->GetXaxis()->SetTitle("Time diff main event trigger - califa hit");
h1_califa_time->GetYaxis()->SetTitle("Count");
h1_califa_time->GetXaxis()->CenterTitle(true);
h1_califa_time->GetYaxis()->CenterTitle(true);
h1_califa_time->GetYaxis()->SetLabelSize(0.045);
h1_califa_time->GetYaxis()->SetTitleSize(0.045);

TH2F* h2_theta_vs_phi_energy[10];
for (int i= 0; i < 10; i++){
sprintf(hist_name, "Theta vs Phi vs Energy, example %i",i+1);
h2_theta_vs_phi_energy[i] = new TH2F(hist_name,hist_name,52,22.15,152.15,60,-180,180);
h2_theta_vs_phi_energy[i]->GetXaxis()->SetTitle("Theta [degree]");
h2_theta_vs_phi_energy[i]->GetYaxis()->SetTitle("Phi [degree]");
h2_theta_vs_phi_energy[i]->GetXaxis()->CenterTitle(true);
h2_theta_vs_phi_energy[i]->GetYaxis()->CenterTitle(true);
h2_theta_vs_phi_energy[i]->GetYaxis()->SetLabelSize(0.045);
h2_theta_vs_phi_energy[i]->GetYaxis()->SetTitleSize(0.045);
}

Int_t number_examples_p2p = 0;
Long64_t nevents = chain->GetEntries(); //number of events in file with name "fname"
cout << "starting with the for loop" << endl;

if (vec_interating_loop.size() == 0){
	cout << "I do not have a file with interesting events, this will now take more time..." << endl;
	for(Long64_t i=0; i< nevents;i++){
		if (i%100000==0) cout<<"Processing vector filling "<<i<<endl;
		vec_interating_loop.push_back(i);
		}
	}

vector<int64_t> califa_time;
//TJ just to look at interesting events
vector<int64_t> vec_interesting_events {323988, 129485 , 402245, 584414,653444, 951429};
//for(Long64_t i : vec_interating_loop){
for(Long64_t i : vec_interesting_events){
	Long64_t evtnr = i;
	chain->GetEvent(i);
	Int_t entries_start = SofSciTcalData->GetEntriesFast();
	Int_t entries_tofw = SofTofWTcalData->GetEntriesFast();
	if (entries_start ==2 && entries_tofw > 1){
		Int_t entries_califa = CalifaHitData->GetEntries();
		Int_t entries_califa_cal = CalifaCalData->GetEntries();
		h1_califa_mult->Fill(entries_califa_cal);
		h1_califa_mult_hit->Fill(entries_califa);
		if (entries_califa > 100) v_events.push_back(evtnr);
		//if (entries_califa == entries_califa_cal){
			califahitdata = new R3BCalifaHitData*[entries_califa];
			for (Int_t j = 0; j < entries_califa; j++){
				califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);
				Double_t theta =  (califahitdata[j]->GetTheta())/PI*180.;
				Double_t phi = (califahitdata[j]->GetPhi())/PI*180;
				h1_califa_time->Fill(int64_t(califahitdata[j]->GetTime()-DataCA->GetTimeStamp()));
				//Double_t energy = (califahitdata[j]->GetEnergy())/1000.;
				Double_t energy = (califahitdata[j]->GetEnergy());
				h2_theta_vs_phi->Fill(theta,phi);
				if (energy < 15000){
					h1_energy_gamma->Fill(energy/1000.);
					}
				if (energy > 30000){
					h1_energy_proton->Fill(energy/1000.);
					}
				//if (evtnr == vec_interating_loop[number_examples_p2p] && number_examples_p2p < 10 && entries_califa > 5){
				if (evtnr == vec_interesting_events[number_examples_p2p] && number_examples_p2p < 10 && entries_califa > 5){
					cout << "This is the eventnr where I plot phi vs theta vs energy:\t" << evtnr << endl;
					for (Int_t j = 0; j < int(energy); j++){
						h2_theta_vs_phi_energy[number_examples_p2p]->Fill(theta,phi);
						}
					}
				}
			//if (evtnr == vec_interating_loop[number_examples_p2p] && number_examples_p2p < 10  && entries_califa > 5){
			if (evtnr == vec_interesting_events[number_examples_p2p] && number_examples_p2p < 10  && entries_califa > 5){

				if (entries_califa > 3) cout << "this is eventnr with multiplicit higher 3:\t" << evtnr << endl;
				if (entries_califa > 4)cout << "this is eventnr with multiplicit higher 4:\t" << evtnr << endl;
				if (entries_califa > 5)cout << "this is eventnr with multiplicit higher 5:\t" << evtnr << endl;
				number_examples_p2p +=1;
				}
				
			delete [] califahitdata; 
			//}
//		else {
//			cout << "For event:\t" << evtnr << "entries califa_cal is not equal entries_califa_hit" << endl;
//			 continue;
//			}
		}

//plot phi vs theta

//select 10 events with p2p and plot phi,theta, energy

}
if (!eventFile){
    cout << " I try to write to file..." << endl;
    t_events->Fill();
    f_events->Write();
    cout << "ending to write to file" << endl;
    cout << "size of vector:\t" <<v_events.size() << endl;
}

TList *l = new TList();
l->Add(h2_theta_vs_phi);
l->Add(h1_califa_mult);
l->Add(h1_califa_mult_hit);
l->Add(h1_energy_gamma);
l->Add(h1_energy_proton);
l->Add(h1_califa_time);
for (Int_t i= 0; i < 10; i++){
	l->Add(h2_theta_vs_phi_energy[i]);
}
l->Write("histlist",TObject::kSingleKey);
cout << "end of process, successful" << endl;

}




