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
#include "find_z.C"
#include "find_isotopes.C"
#include "all_params.C"
#include "find_psi_par.C"
#include <string>
using namespace std;

char fname[500];
bool sortcol( const vector<double>& v1,const vector<double>& v2 ) {
                                    return v1[1] > v2[1];
                                    }
bool sortcol_pseudo( const vector<double>& v1,const vector<double>& v2 ) {
                                    return v1[0] > v2[0];
                                    }
bool sortcol_lorentz( const TLorentzVector& v1,const TLorentzVector& v2 ) {
                                    return v1[3] > v2[3];
                                    }

TVector3 uniform_vec(){
    TRandom* event_gen_x = new TRandom();
    TRandom* event_gen_y = new TRandom();
    TRandom* event_gen_z = new TRandom();
	event_gen_x->SetSeed();
	event_gen_y->SetSeed();
	event_gen_z->SetSeed();
    Double_t x_comp = event_gen_x->Uniform(-1,1);
    Double_t y_comp = event_gen_y->Uniform(-1,1);
    Double_t z_comp = event_gen_z->Uniform(-1,1);
    TVector3 u_vec(x_comp,y_comp,z_comp);
    delete event_gen_x;
    delete event_gen_y;
    delete event_gen_z;
    return u_vec;
}


void measure_11b(char const count_i[50]){
cout << "thsi is count_i " << count_i << endl;
cout << "and this is count_i[50] " << count_i[50] << endl;
stringstream ss;
ss << count_i;
string my_runs;
ss >> my_runs;
all_params(my_runs);
char f_out_name[500];
char f_events_name[500];
char hist_name[500];
char f_parameters[500];
static TLorentzVector p_missing_old;
static TLorentzVector p_missing_new;
//DEFINITIONS OF THE HISTOGRAMS:

TH2F* h2_a_vs_q;
sprintf(hist_name, "A vs q ");
h2_a_vs_q = new TH2F(hist_name,hist_name,200,1,3,100,0,10);
h2_a_vs_q->GetXaxis()->SetTitle("A/q");
h2_a_vs_q->GetYaxis()->SetTitle("q");
h2_a_vs_q->GetXaxis()->CenterTitle(true);
h2_a_vs_q->GetYaxis()->CenterTitle(true);
h2_a_vs_q->GetYaxis()->SetLabelSize(0.045);
h2_a_vs_q->GetYaxis()->SetTitleSize(0.045);

//END OF IMPLEMENTATION OF HISTOS--------------------------------------------------
sprintf(f_events_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/vec_11bp2p_events_%s.root", count_i);
TFile* eventFile = TFile::Open(f_events_name,"READ");
vector<Long64_t> v_events;
vector<Long64_t> vec_interating_loop;
TFile * f_events;
TTree *t_events;
if (!eventFile){
	f_events = TFile::Open(f_events_name,"RECREATE");
	t_events = new TTree("tvec","Tree with vectors");
	t_events->Branch("v_events",&v_events);
	}
if (eventFile){
	TTreeReader myReader("tvec", eventFile);
	TTreeReaderValue<vector<Long64_t> >vec_good_events(myReader,"v_events");
	myReader.Next();
	vec_interating_loop = *vec_good_events;
	eventFile->Close();
	}





sprintf(f_out_name,"/home/ge37liw/a_vs_q_ch2_run_%s.root", count_i);

TFile * f = new TFile(f_out_name,"RECREATE");
sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/p2p_data_fixed_angle/with_time_info_nominal/stitched_and_unpacked_main%s.root",count_i);

sprintf(f_parameters,"/home/ge37liw/plots/histos/zero_with_miss/par_files_dir/par_run_%s.txt",count_i);
double t = find_z(f_parameters,fname);
cout << "this is the return value\t" << t << endl;
vector<double> s = find_isotopes(f_parameters,fname,t);
cout << "the isotope cuts are:\t" <<endl;\
cout << "10B:\t " << s[0] << endl;\
cout << "11B:\t" << s[1] << endl;\
cout << "11C:\t" << s[2] << endl;\
cout << "12C:\t" << s[3] << endl;

vector<double> u = find_psi_par(f_parameters,fname,t,s);
cout << "job done, now we can proceede with analysis" << endl;

//Finally the precise calculation with the parameters from Step1 -> Step3
cout << "finally precise calculation!---------------------------------" << endl;

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
TClonesArray* CalifaHitData = new TClonesArray("R3BCalifaHitData",3);
R3BCalifaHitData** califahitdata;
TBranch *branchCalifaHitData = chain->GetBranch("CalifaHitData");
branchCalifaHitData->SetAddress(&CalifaHitData);
//
R3BEventHeader* DataCA = new R3BEventHeader();
TBranch* branchData = chain->GetBranch("EventHeader.");
branchData->SetAddress(&DataCA);
//
//TClonesArray* WRData = new TClonesArray("R3BWRData",1);
//R3BWRData** wrmasterdata;
//TBranch *branchWRMasterData = chain->GetBranch("WRData");
//branchWRMasterData->SetAddress(&WRData);
//
//
Long64_t nevents = chain->GetEntries(); //number of events in file with name "fname"
for(Long64_t i=0; i< nevents;i++){
    Long64_t evtnr = i;

    if (i%100000==0)
        cout<<"Processing event "<<i<<endl;

	chain->GetEvent(i); //event_i_call
	entries_mw0 = Mwpc0HitData->GetEntries();
    entries_mw1 = Mwpc1HitData->GetEntries();
    entries_mw2 = Mwpc2HitData->GetEntries();
    entries_mw3 = Mwpc3HitData->GetEntries();
    entries_start = SofSciTcalData->GetEntriesFast();
    entries_tofw = SofTofWTcalData->GetEntriesFast();
    entries_califa = CalifaHitData->GetEntries();
    softwimhitdata = new R3BTwimHitData*[1];
    softwimhitdata[0] = (R3BTwimHitData*)TwimHitData->At(0);

    if (entries_start ==2 && entries_tofw == 2 && entries_mw0 == 1 && entries_mw1 == 1 && entries_mw2 == 1 && entries_mw3 == 1 && softwimhitdata[0]){
        charge_val = softwimhitdata[0]->GetZcharge();
        sofmwpc0hitdata = new R3BMwpcHitData*[1];
        sofmwpc1hitdata = new R3BMwpcHitData*[1];
        sofmwpc2hitdata = new R3BMwpcHitData*[1];
        sofmwpc3hitdata = new R3BMwpcHitData*[1];
        softofwtcaldata = new R3BSofTofWTcalData*[1];
        softofwmappeddata = new R3BSofTofWMappedData*[2];
        califahitdata = new R3BCalifaHitData*[entries_califa];
        sofmwpc0hitdata[0] = (R3BMwpcHitData*)Mwpc0HitData->At(0);
        sofmwpc1hitdata[0] = (R3BMwpcHitData*)Mwpc1HitData->At(0);
        sofmwpc2hitdata[0] = (R3BMwpcHitData*)Mwpc2HitData->At(0);
        sofmwpc3hitdata[0] = (R3BMwpcHitData*)Mwpc3HitData->At(0);
        softofwtcaldata[0] = (R3BSofTofWTcalData*)SofTofWTcalData->At(0);
        Double_t xMW0 = sofmwpc0hitdata[0]->GetX();
        Double_t xMW1 = sofmwpc1hitdata[0]->GetX();
        Double_t xMW2 = sofmwpc2hitdata[0]->GetX();
        Double_t xMW3 = sofmwpc3hitdata[0]->GetX();
		Double_t yMW0 = sofmwpc0hitdata[0]->GetY();
        Double_t yMW1 = sofmwpc1hitdata[0]->GetY();
        Double_t yMW2 = sofmwpc2hitdata[0]->GetY();
        Double_t yMW3 = sofmwpc3hitdata[0]->GetY();
if (xMW1 != -1000 && xMW2 != -1000 && xMW3 != -1000  && xMW0 != -1000 && yMW0 != -1000){
        Double_t psi_in = atan(((xMW2+xMW2_shift)-(xMW1+xMW1_shift))/(zM2-zM1));
        Double_t T_to_M2 = (zM2-zT)/(cos(psi_in));
        Double_t b1_cutoff = (xMW1+xMW1_shift)-tan(psi_in)*zM1;
        Double_t z_labMW3 = middle_zM3+cos(2*PI/5.)*(xMW3-xMW3_shift);
        Double_t x_labMW3 = middle_xM3+sin(2*PI/5.)*(xMW3-xMW3_shift);
        Double_t zB = (bGLAD_cutoff-b1_cutoff)/(tan(psi_in)-tan(PI/2.-alpha_G));
        Double_t xB = tan(psi_in)*zB+b1_cutoff;
        for (Int_t rec = 0; rec < 50; rec++){
            zC = (zGm_cutoff-b1_cutoff)/(tan(psi_in)-tan(PI/2.-alpha_G)) + 2*rec;
            xC = tan(psi_in)*zC+b1_cutoff;
            slope = (x_labMW3-xC)/(z_labMW3-zC);
            offset_slope = x_labMW3-slope*z_labMW3;
            psi_out_rec[rec] = -atan(slope);
            z_D = (D_cutoff-offset_slope)/(slope-tan(PI/2.-alpha_G));
            x_D = tan(PI/2.-alpha_G)*z_D +D_cutoff;
            l_diff[rec] = (sqrt((xC-xB)*(xC-xB)+(zC-zB)*(zC-zB))-sqrt((x_D-xC)*(x_D-xC)+(z_D-zC)*(z_D-zC)));
            z_pos_shift[rec] = zC;
            x_pos_shift[rec] = xC;
        }
        Int_t n = 50;
        TGraph *gr1 = new TGraph(n,psi_out_rec,l_diff);
        gr1->Fit("pol1","Q");
        TF1 *f3 =gr1->GetFunction("pol1");
        Double_t f3_slope = f3->GetParameter(1);
        Double_t f3_offset = f3->GetParameter(0);
        Double_t psi_out = -f3_offset/f3_slope;
        TGraph *gr2 = new TGraph(n,z_pos_shift,l_diff);
        gr2->Fit("pol1", "Q");
        TF1 *f4 = gr2->GetFunction("pol1");
        Double_t f4_slope = f4->GetParameter(1);
        Double_t f4_offset = f4->GetParameter(0);
        Double_t new_z_C = -f4_offset/f4_slope;

        TGraph *gr3 = new TGraph(n,x_pos_shift,l_diff);
        gr3->Fit("pol1","Q");
        TF1 *f5 = gr3->GetFunction("pol1");
        Double_t f5_slope = f5->GetParameter(1);
        Double_t f5_offset = f5->GetParameter(0);
        Double_t new_x_C = -f5_offset/f5_slope;
        Double_t z_diff = new_z_C-z_pos_shift[0];
        Double_t offset_new = new_x_C+tan(psi_out)*new_z_C;
        z_D = (D_cutoff-offset_new)/(-tan(psi_out)-tan(PI/2.-alpha_G));
        x_D = -tan(psi_out)*z_D + offset_new;
        Double_t zToFW = (cutoff_ToFW-offset_new)/(-tan(psi_out)-tan(2*PI/5.));
        Double_t xToFW = -tan(psi_out)*zToFW+offset_new;
        Double_t path_D_to_TOFW = (zToFW-z_D)/(cos(psi_out));
        Double_t BD = sqrt((x_D-xB)*(x_D-xB)+(z_D-zB)*(z_D-zB));

        Double_t m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
        Double_t rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
        Double_t w = 2*abs(asin(BD/(2*rho)));


        tot_length = start_to_target+ abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
        delete gr1;
        delete gr2;
        delete gr3;

        if (SofTofWTcalData && SofTofWTcalData->GetEntriesFast())
            {
            for (UShort_t i = 0; i < NbDets; i++)
                {
                for (UShort_t j = 0; j < NbChs; j++)
                    {
                    mult[i * NbChs + j] = 0;
                    }
                }
            for (Int_t ihit = 0; ihit < entries_tofw; ihit++)
                {
                R3BSofTofWMappedData* hitmapped = (R3BSofTofWMappedData*)SofTofWMappedData->At(ihit);
                if (!hitmapped)
                    continue;
                iDet = hitmapped->GetDetector() - 1;
                iCh = hitmapped->GetPmt() - 1;
                mult[iDet * NbChs + iCh]++;
                }
            for (Int_t ihit = 0; ihit < entries_tofw; ihit++)
                {
                R3BSofTofWTcalData* hittcal = (R3BSofTofWTcalData*)SofTofWTcalData->At(ihit);
                if (!hittcal)
                    continue;
                if (hittcal->GetPmt() == 3)
                    continue;
                iDet = hittcal->GetDetector() - 1;
                iCh = hittcal->GetPmt() - 1;
                iRawTimeNs[iDet * 2 + iCh] = hittcal->GetRawTimeNs();
                }
            for (Int_t s_hit = 0 ; s_hit < entries_start; s_hit++)
                {
                R3BSofSciTcalData* start_hits = (R3BSofSciTcalData*)SofSciTcalData->At(s_hit);
                if(!start_hits)
                    continue;
                iRawTimeStartNs[s_hit] = start_hits->GetRawTimeNs();

                }
                sofmwpc3hitdata[0]=(R3BMwpcHitData*)Mwpc3HitData->At(0);
            if ((sofmwpc3hitdata[0]->GetY() > -500.) && (sofmwpc3hitdata[0]->GetX() > -500.))  //cut on the mwpc3
            {
            Double_t raw_mwpc3 = sofmwpc3hitdata[0]->GetY();
            Double_t tofpos_ns = -1000;
            Double_t raw_tofpos_ns = -1000;
            Double_t raw_t_start = -1000000.;
            Double_t raw_time_of_flight = 0.;
            for (UShort_t i = 0; i < NbDets; i++)
                {
                if ((mult[i * NbChs] == 1) && (mult[i * NbChs + 1] == 1))
                    {
                    tofpos_ns = 0.5*(iRawTimeNs[i * NbChs + 1] - iRawTimeNs[i * NbChs]-offsets[i+1]);
                    raw_tofpos_ns = 0.5*(iRawTimeNs[i * NbChs + 1] - iRawTimeNs[i * NbChs]);
                    raw_t_start = 0.5*(iRawTimeStartNs[0]-iRawTimeStartNs[1]);
                    if (raw_t_start !=-1000000. && raw_tofpos_ns != -1000 && tofpos_ns != -1000)
                        {
                        Double_t tof_diff_up_down = 0.5*(iRawTimeNs[i * NbChs + 1] - iRawTimeNs[i * NbChs]);
                        raw_time_of_flight=(0.5*(iRawTimeNs[i*NbChs+1]+iRawTimeNs[i*NbChs]))-0.5*(iRawTimeStartNs[0]+iRawTimeStartNs[1]);
                        Double_t time_target_tof = raw_time_of_flight+tof_offs[i] -time_start_target;
                        Double_t path_from_target = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta = ((path_from_target/time_target_tof)*pow(10,6))/light_c;
                        Double_t mag_field = current*0.0006527728074785267;
                        Double_t a_q = (pow(10,-3)*((((mag_field*rho)/(beta))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
						h2_a_vs_q->Fill(a_q,charge_val);

}
}
}
}
}
}
 delete [] sofmwpc0hitdata;
 delete [] sofmwpc1hitdata;
 delete [] sofmwpc2hitdata;
 delete [] sofmwpc3hitdata;
 delete [] softofwtcaldata;
 delete [] softofwmappeddata;
 delete [] califahitdata;
}
delete [] softwimhitdata;
}

if (!eventFile){
	cout << " I try to write to file..." << endl;
	t_events->Fill();
	f_events->Write();
	cout << "ending to write to file" << endl;
}
cout << "almost done .. now only filling tlist...." << endl;
TList *l = new TList();
l->Add(h2_a_vs_q);

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
cout << "end of process, successful" <<endl;
}



