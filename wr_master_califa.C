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
char f_out_name[500];

void wr_master_califa(char const count_i[50]){
	//name input file
	sprintf(fname,"/u/land/tobias_jenegger/test_s467_unpacker/stitched_and_unpacked_with_califa_wr_ts_main%s.root",count_i);
	//name output file
	sprintf(f_out_name,"/u/land/tobias_jenegger/test_s467_unpacker/histos_califa_wr_ts_main%s.root",count_i);
	//HISTOS
	TH1F* h1_wr_time_diff_califa_master;
	sprintf(hist_name, "WR time diff. CALIFA - MASTER");
	h1_wr_time_diff_califa_master = new TH1F(hist_name,hist_name,20000,-10000,10000);
	h1_wr_time_diff_califa_master->GetXaxis()->SetTitle("WR diff. ns");
	h1_wr_time_diff_califa_master->GetYaxis()->SetTitle("Counts");
	h1_wr_time_diff_califa_master->GetXaxis()->CenterTitle(true);
	h1_wr_time_diff_califa_master->GetYaxis()->CenterTitle(true);
	h1_wr_time_diff_califa_master->GetYaxis()->SetLabelSize(0.045);
	h1_wr_time_diff_califa_master->GetYaxis()->SetTitleSize(0.045);		
	//END OF HISTOS

	TChain* chain = new TChain("evt");
	chain->Reset();
	chain->Add(fname);

	//TClonesArrays
	TClonesArray* CalifaHitData = new TClonesArray("R3BCalifaHitData",3);
	R3BCalifaHitData** califahitdata;
	TBranch *branchCalifaHitData = chain->GetBranch("CalifaHitData");
	branchCalifaHitData->SetAddress(&CalifaHitData);
	//
	R3BEventHeader* DataCA = new R3BEventHeader();
	TBranch* branchData = chain->GetBranch("EventHeader.");
	branchData->SetAddress(&DataCA);
	//
	
	Int_t entries_califa_hit;	
	uint64_t wr_master;
	uint64_t wr_califa_hit;
	uint64_t wr_time_diff_master_califa;
	for(Long64_t i=0; i< nevents;i++){
	    Long64_t evtnr = i;
	
	    if (i%100000==0)
	        cout<<"Processing event "<<i<<endl;
		chain->GetEvent(i);	
		//get wr master timestamp
		wr_master = DataCA->GetTimeStamp();
		entries_califa_hit = CalifaHitData->GetEntries();
		califahitdata = new R3BCalifaHitData*[entries_califa_hit];
		for (Int_t j = 0.;j<entries_califa_hit;j++){
			califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);	
			wr_califa_hit = califahitdata[j]->GetTime();
			wr_time_diff_master_califa = wr_califa_hit - wr_master;
			h1_wr_time_diff_califa_master->Fill(double(wr_time_diff_master_califa));
			}
		delete [] califahitdata;
		}
	TFile * f = new TFile(f_out_name,"RECREATE");
	TList *l = new TList();
	l->Add(h1_wr_time_diff_califa_master);
	l->Write("histlist", TObject::kSingleKey);
	delete CalifaHitData;
	delete DataCA;
	cout << "end of process, successful" <<endl;
}
