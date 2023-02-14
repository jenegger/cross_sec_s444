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


void clean_macro(char const count_i[50]){
cout << "thsi is count_i " << count_i << endl;
cout << "and this is count_i[50] " << count_i[50] << endl;
stringstream ss;
ss << count_i;
string my_runs;
ss >> my_runs;
all_params(my_runs);
char f_out_name[500];
char hist_name[500];
char f_parameters[500];
static TLorentzVector p_missing_old;
static TLorentzVector p_missing_new;



//OUTPUT FILENAME

//sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/plus_48/foo_%s.root", count_i);
//sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/minus_16/foo_%s.root", count_i);
//sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/minus_24/foo_%s.root", count_i);
sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/plus_24/foo_%s.root", count_i);
//sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/foo_nominal_%s.root", count_i);
//for carbon
//sprintf(f_out_name,"/scratch5/ge37liw/macros_s444/analysis_p2p/carbon_target_foo_%s.root", count_i);
//TJ for testing reasons
//sprintf(f_out_name,"/home/ge37liw/plots/histos/zero_with_miss/fast_macro/minus2_4_phi_test_automatic_src_califa_%s.root", count_i);
//sprintf(f_out_name,"/home/ge37liw/plots/histos/zero_with_miss/fast_macro/plus2_4_phi_test_automatic_src_califa_%s.root", count_i);

TFile * f = new TFile(f_out_name,"RECREATE");



//INPUT FILENAME

//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/stitched_and_unpacked_main%s.root",count_i);
//now using the feature of random angles
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/p2p_data_rand_angle/stitched_and_unpacked_main%s.root",count_i);
//this is with fixed angles
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/p2p_data_fixed_angle/stitched_and_unpacked_main%s.root",count_i);
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/p2p_data_fixed_angle/plus_48/stitched_and_unpacked_main%s.root",count_i);
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/p2p_data_fixed_angle/minus_16/stitched_and_unpacked_main%s.root",count_i);
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/p2p_data_fixed_angle/minus_24/stitched_and_unpacked_main%s.root",count_i);
sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/p2p_data_fixed_angle/plus_24/stitched_and_unpacked_main%s.root",count_i);

//for carbon target
//sprintf(fname,"/scratch5/ge37liw/unpacked_s444_data/carbon_target/stitched_and_unpacked_main%s.root",count_i);
//TJ:for test reasons...
//sprintf(fname,"/scratch8/ge37liw/workingspace/data/root_files/all_ts_unpack/sweep_target/ts_cone_cluster_%s_minus2_4cm.root",count_i);
//sprintf(fname,"/scratch8/ge37liw/workingspace/data/root_files/all_ts_unpack/sweep_target/ts_cone_cluster_%s_plus2_4cm.root",count_i);
//END TJ
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

						//it has to beconsidered that the beam does not go parallel to the z-direction, but has deflections in x and y. Therefore also for later analysis and boosting the right trajectory 
						//before the 12C hits the target has to be found. I do this by using the info of the target point and the x-y info of the MW0
						TVector3 v_mw0(xMW0+xMW0_shift,yMW0+yMW0_shift,-2520); //vector of the MW0
						TVector3 v_target(0,0,-684.5); //vector of target position
						TVector3 v_diff_mw0_target;
						v_diff_mw0_target = v_mw0 - v_target;
						Double_t c12_arzimuth = v_diff_mw0_target.Phi(); // get arzimuthal angle of the incoming 12C beam (careful, it's negative..., or?)
						Double_t c12_theta = v_diff_mw0_target.Theta(); // get polar angle of the incoming 12C beam (careful, it's negative ..., or?)
						Double_t beta_x = abs(sin(c12_theta)*cos(c12_arzimuth));
						Double_t beta_y = abs(sin(c12_theta)*sin(c12_arzimuth));
						Double_t beta_z = abs(cos(c12_theta));
						TVector3 v_beta_12c (beta_x*beta_given,beta_y*beta_given,beta_z*beta_given);
						//compute position direction of fragment after target using info of mw1 and mw2
						TVector3 v_mw1(xMW1+xMW1_shift,yMW1+yMW1_shift,zM1);
						TVector3 v_mw2(xMW2+xMW2_shift,yMW2+yMW2_shift,zM2);
						TVector3 v_diff_mw2_mw1;
						v_diff_mw2_mw1 = v_mw2 - v_mw1;
						Double_t fragment_theta = v_diff_mw2_mw1.Theta();
						Double_t fragment_arzimuth = v_diff_mw2_mw1.Phi();
						Double_t fragment_x = sin(fragment_theta)*cos(fragment_arzimuth);
						Double_t fragment_y = sin(fragment_theta)*sin(fragment_arzimuth);
						Double_t fragment_z = cos(fragment_theta);
						if (charge_val< t+0.8 && charge_val> t){ //Z = 6
							//cout << "evtnr:\t" << evtnr << endl;
							h2_xmw2_vs_xmw3->Fill(xMW2,xMW3);
							if(xMW2 > 400) cout << "large mw2_x value:\t" << xMW2 << endl;
							}
                        if (charge_val< t && charge_val> t-1.) //Z = 5 before it was charge_val< charge_cut && charge_val> charge_cut-0.8, too straight ....
                        {
//############################################################################################################
//------------BEGIN OF 11B ANALYSIS---------------------------------------------------------------------------
//############################################################################################################
                        if (a_q > s[0] && a_q < (s[0]+2*s[3]) && entries_califa >= 2){  //11B
						//	cout << "inside 11b eventnumber:\t" << evtnr << endl;
						psi_in = u[3] - u[4]*xMW3 - u[5];
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
                        Double_t gamma_precise = 1/(sqrt(1-beta_precise*beta_precise));
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
						h1_beta_precise_11b->Fill(beta_precise);
						// checking MWPCs Y-orientation:
						h2_corr_y_mws->Fill((yMW2+yMW2_shift)-(yMW1+yMW1_shift),(yMW3+yMW3_shift)-(yMW1+yMW1_shift));
						h2_corr_y_mws_minusmw3->Fill((yMW2+yMW2_shift)-(yMW1+yMW1_shift),-(yMW3+yMW3_shift)-(yMW1+yMW1_shift));
						//TJ checking position mw2_x and mw3_y
						h1_x_position_mw2_11b->Fill(xMW2+xMW2_shift);	
						h1_y_position_mw3_11b->Fill(yMW3+yMW3_shift);
						//compute position direction of fragment after the target using the target pos. and the MW1/2 (x-pos) and MWPC3 (y-pos.)
//--> dangerous!! set fragment y to -y, TODO: change it back after usage
						fragment_y = -(yMW3+yMW3_shift)/(pathlength_precise-760.4); //substract path between MW3 and TOFW
						fragment_x =  v_diff_mw2_mw1.X()/v_diff_mw2_mw1.Mag();
						fragment_z = v_diff_mw2_mw1.Z()/v_diff_mw2_mw1.Mag();	
						//beginning after isotope has been selected:
						//TLorentzVector of 11B in the lab and cms 
						Double_t tot_mom11B = 10255.09731131*(beta_precise/(sqrt(1-(beta_precise*beta_precise)))); //here now right mass with mass defect
						//TODO: this energy calculation is wrong, be careful and recheck it..
						//Double_t energy_11B = sqrt(pow(tot_mom11B,2)+pow(11*938.272,2));
						Double_t energy_11B = (1./(sqrt(1-(beta_precise*beta_precise))))*10255.09731131;  //here now right mass with mass defect
						TLorentzVector v_11B_lab(tot_mom11B*fragment_x,tot_mom11B*fragment_y,tot_mom11B*fragment_z,energy_11B);
						TLorentzVector v_11B_cms = v_11B_lab;
						v_11B_cms.Boost(-v_beta_12c);
						h1_px_11B_lab_test->Fill(v_11B_lab.Px());
						h1_py_11B_lab_test->Fill(v_11B_lab.Py());
						//Lorentz vector of the 12C in rest 
						TLorentzVector p_12C_cms(0,0,0,11177.92337806);

						//TLorentzVector of the target in lab and cms frame
						TLorentzVector v_target_lab(0,0,0,938.272);
						TLorentzVector v_target_cms = v_target_lab;
						v_target_cms.Boost(-v_beta_12c);
						//fill entries of CALIFA in vectors, no matter if they are protons or gammas..
						//variales needed:
						vector<vector<double> > v_CALIFA_evts_lab_energy; //here I insert vectors v = (E_lab,theta_lab,phi_lab)
						vector<vector<double> > v_CALIFA_evts_cms_energy; //here I inserv vectors v = (E_cms,theta_lab,phi_lab)
						//for latest update I change to: v_CALIFA_evts_cms_energy(E_cms,theta_lab,phi_lab,Nf,Ns) 
						vector<TLorentzVector> v_lorentz_CALIFA_evts_lab; //here I insert Lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_lab)
						vector<TLorentzVector> v_lorentz_CALIFA_evts_cms; //here I insert Lorentz vector v = (p_x_cms,p_y_cms,p_z_cms,E_cms)
						vector<TLorentzVector> v_lorentz_PROTONS_lab; //here I inserrt lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_lab) where E_lab > 30MeV
						vector<TLorentzVector> v_lorentz_PROTONS_cms; //here I inserrt lorentz vector v = (p_x_cms,p_y_cms,p_z_cms,E_cms) 
						vector<TLorentzVector> v_lorentz_GAMMAS_lab; //here I inserrt lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_lab) where E_lab < 10MeV for now... (change to 30MeV??)
						vector<TLorentzVector> v_lorentz_GAMMAS_cms; //here I inserrt lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_cms) 
						
						for (Int_t j = 0.;j<entries_califa;j++){
							vector<double> v_temp(3);
							vector<double> v_temp_cms(5);
							califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);
							v_temp[0] = (califahitdata[j]->GetEnergy())/1000.;
							//randomize angles to remove artifacts...
							//if (califahitdata[j]->GetTheta() > 0.45){
							//	Int_t rand_angl = rand() % 4500 + 1;
							//	v_temp[1] = califahitdata[j]->GetTheta() - 0.0225 + 0.00001*rand_angl;
							//	}
							//else {
							//	v_temp[1] = califahitdata[j]->GetTheta();
							//	}
							//
							v_temp[1] = califahitdata[j]->GetTheta();
							v_temp[2] = califahitdata[j]->GetPhi();
							v_temp_cms[0] = v_temp[0]*(1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(v_temp[1])); //Doppler correction
							v_temp_cms[1] = v_temp[1];
						    v_temp_cms[2] = v_temp[2];
							v_temp_cms[3] = (califahitdata[j]->GetNf());
							v_temp_cms[4] = (califahitdata[j]->GetNs());
							v_CALIFA_evts_lab_energy.push_back(v_temp);
							v_CALIFA_evts_cms_energy.push_back(v_temp_cms);
							v_temp.clear();
							v_temp_cms.clear();
							
						}
						
						
						for (int CALIFA_events = 0; CALIFA_events < entries_califa; CALIFA_events++){
						
						if(v_CALIFA_evts_lab_energy[CALIFA_events][0] < 15) { //gammas
						//create lorentz vector in the laboratory system
						TLorentzVector temp_lorentz_lab;
						TLorentzVector temp_lorentz_cms;
						temp_lorentz_lab[0] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[1] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[2] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][1]); //coordinate transform
						temp_lorentz_lab[3] = v_CALIFA_evts_lab_energy[CALIFA_events][0];
						v_lorentz_CALIFA_evts_lab.push_back(temp_lorentz_lab);
						v_lorentz_GAMMAS_lab.push_back(temp_lorentz_lab);
						temp_lorentz_lab.Boost(-v_beta_12c);
						v_lorentz_CALIFA_evts_cms.push_back(temp_lorentz_lab);
						v_lorentz_GAMMAS_cms.push_back(temp_lorentz_lab);
						temp_lorentz_lab = TLorentzVector();
						
						}
						else {//protons
						
						TLorentzVector temp_lorentz_lab;
						temp_lorentz_lab[0] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[1] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[2] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*cos(v_CALIFA_evts_lab_energy[CALIFA_events][1]); //coordinate transform
						temp_lorentz_lab[3] = v_CALIFA_evts_lab_energy[CALIFA_events][0]+938.272;
						v_lorentz_CALIFA_evts_lab.push_back(temp_lorentz_lab);
						v_lorentz_PROTONS_lab.push_back(temp_lorentz_lab);
						temp_lorentz_lab.Boost(-v_beta_12c);
						v_lorentz_CALIFA_evts_cms.push_back(temp_lorentz_lab);
						v_lorentz_PROTONS_cms.push_back(temp_lorentz_lab);
						temp_lorentz_lab = TLorentzVector();	
						}
						}
						//now sort the vectors according to the Energy.......
						//3d vectors
						sort(v_CALIFA_evts_cms_energy.begin(),v_CALIFA_evts_cms_energy.end(),sortcol_pseudo);
						sort(v_CALIFA_evts_lab_energy.begin(),v_CALIFA_evts_lab_energy.end(),sortcol_pseudo);
						//lorentz vectors
						sort(v_lorentz_CALIFA_evts_lab.begin(),v_lorentz_CALIFA_evts_lab.end(),sortcol_lorentz);
						sort(v_lorentz_CALIFA_evts_cms.begin(),v_lorentz_CALIFA_evts_cms.end(),sortcol_lorentz);
						sort(v_lorentz_GAMMAS_lab.begin(),v_lorentz_GAMMAS_lab.end(),sortcol_lorentz);
						sort(v_lorentz_GAMMAS_cms.begin(),v_lorentz_GAMMAS_cms.end(),sortcol_lorentz);
						sort(v_lorentz_PROTONS_lab.begin(),v_lorentz_PROTONS_lab.end(),sortcol_lorentz);
						sort(v_lorentz_PROTONS_cms.begin(),v_lorentz_PROTONS_cms.end(),sortcol_lorentz);

						//make restriction of having exactly two protons (with more than 30 MeV in lab frame...)
						if (v_CALIFA_evts_lab_energy.size() >= 2 && v_CALIFA_evts_lab_energy[0][0] > 30 && v_CALIFA_evts_lab_energy[1][0] > 30){
							if (v_CALIFA_evts_lab_energy.size() == 2 || (v_CALIFA_evts_lab_energy.size() > 2 && v_CALIFA_evts_lab_energy[2][0] < 30)){

								//TJ check
								//cout << "make cut on the two protons to be larger 30MeV, eventnr:\t" << evtnr << endl;
								//cout << "corresponding theta_sum:\t" << (v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1])  << endl;
								//cout << "corresponding phidiff:\t" << abs(v_CALIFA_evts_cms_energy[0][2]-v_CALIFA_evts_cms_energy[1][2]) << endl;
								//real opening angle calculation:
								Double_t theta_angle_proton1 = v_CALIFA_evts_lab_energy[0][1];
								Double_t theta_angle_proton2 =  v_CALIFA_evts_lab_energy[1][1];
								Double_t phi_angle_proton1 = v_CALIFA_evts_lab_energy[0][2]; 
								Double_t phi_angle_proton2 = v_CALIFA_evts_lab_energy[1][2];
								Double_t cos_opening_angle_two_protons = sin(theta_angle_proton1)*sin(theta_angle_proton2)*cos(phi_angle_proton2-phi_angle_proton1)+cos(theta_angle_proton1)*cos(theta_angle_proton2);
								Double_t opening_angle = acos(cos_opening_angle_two_protons);
								//cout << "this is the opening angle:\t" << opening_angle/PI*180. << endl;
								int rand_num;
                                rand_num = rand() % 2;
								//now plot randomized theta1 theta2 correlations
								Double_t rand_angle1 = v_CALIFA_evts_lab_energy[rand_num][1];
								Double_t rand_angle2 = v_CALIFA_evts_lab_energy[1-rand_num][1];
								h1_opening_angle_p1_p2->Fill((opening_angle/PI)*180.);	
								h2_theta1_vs_theta2_rand->Fill(rand_angle1*180./PI,rand_angle2*180./PI);
								h1_theta_sum_p1_p2->Fill(theta_angle_proton1*180./PI,theta_angle_proton2*180./PI);
								//cout << "phi_1:\t" << phi_angle_proton1*180/PI << endl;
								//cout << "phi_2:\t" << phi_angle_proton2*180/PI << endl;
								//plot theta_messel vs theta wixhausen
								//proton 1 is messel, proton 2 is wixhausen
								if (abs(phi_angle_proton1*180/PI) < 90. && abs(phi_angle_proton2*180/PI) > 90.){
									h2_theta1_vs_theta2_mes_wix->Fill(theta_angle_proton1*180./PI,theta_angle_proton2*180./PI);		
									}
								//proton 1 is wixhausen, proton 2 is messel
								if (abs(phi_angle_proton2*180/PI) < 90. && abs(phi_angle_proton1*180/PI) > 90.){
									h2_theta1_vs_theta2_mes_wix->Fill(theta_angle_proton2*180./PI,theta_angle_proton1*180./PI);
									}

								//make a cut on the good p2p events:
								if ((v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1]) > 1.308996939 && (v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1]) < 1.4835298642 \
									&& abs(v_CALIFA_evts_cms_energy[0][2]-v_CALIFA_evts_cms_energy[1][2]) < 3.6651914292 && abs(v_CALIFA_evts_cms_energy[0][2]-v_CALIFA_evts_cms_energy[1][2]) > 2.617993878){
									h2_theta1_vs_theta2_rand_good->Fill(rand_angle1*180./PI,rand_angle2*180./PI);
								//cout << "good p2p event:\t" << evtnr << endl;
								//cout << "proton lab one momentum:\t" << sqrt((v_lorentz_PROTONS_lab[0]).Px()*(v_lorentz_PROTONS_lab[0]).Px()+(v_lorentz_PROTONS_lab[0]).Py()*(v_lorentz_PROTONS_lab[0]).Py() +\
								//(v_lorentz_PROTONS_lab[0]).Pz()*(v_lorentz_PROTONS_lab[0]).Pz()) << endl;
								//cout << "proton lab two momentum:\t" << sqrt((v_lorentz_PROTONS_lab[1]).Px()*(v_lorentz_PROTONS_lab[1]).Px()+(v_lorentz_PROTONS_lab[1]).Py()*(v_lorentz_PROTONS_lab[1]).Py() +\
								//(v_lorentz_PROTONS_lab[1]).Pz()*(v_lorentz_PROTONS_lab[1]).Pz()) << endl;
								//cout << "momentum of 11B in 12CMS:\t" << sqrt(v_11B_cms.Px()*v_11B_cms.Px() + v_11B_cms.Py()*v_11B_cms.Py() + v_11B_cms.Pz()*v_11B_cms.Pz()) << endl;
								//cout << "z_momentum of 11B in 12C CMS:\t" << v_11B_cms.Pz() << endl;
								//cout << "z_mom of two protons lab sys:\t" << (v_lorentz_PROTONS_lab[0]).Pz()+(v_lorentz_PROTONS_lab[1]).Pz() << endl;
								//cout << "z_mom of two protons in CMS system:\t" << (v_lorentz_PROTONS_cms[0]).Pz()+(v_lorentz_PROTONS_cms[1]).Pz() << endl;
								
								//check momentum correlations
								TLorentzVector p_knocked_proton = v_lorentz_PROTONS_cms[0] +v_lorentz_PROTONS_cms[1]-v_target_cms;
								h2_pi_x_y_cms->Fill(p_knocked_proton.Px(),p_knocked_proton.Py());		
								h2_pi_x_z_cms->Fill(p_knocked_proton.Pz(),p_knocked_proton.Px());		

								Double_t value_11B_cms = sqrt(pow(v_11B_cms.Px(),2)+pow(v_11B_cms.Py(),2)+pow(v_11B_cms.Pz(),2));
								h2_p11bval_vs_mw2x->Fill(xMW2+xMW2_shift,value_11B_cms);
								h2_p11bval_vs_mw2y->Fill(yMW2+yMW2_shift,value_11B_cms);
								h2_p11b_x_y_cms->Fill(v_11B_cms.Px(),v_11B_cms.Py());
								h2_p11b_x_z_cms->Fill(v_11B_cms.Pz(),v_11B_cms.Px());

								//provisional...
								TLorentzVector p_missing_prov = v_lorentz_PROTONS_cms[0] +v_lorentz_PROTONS_cms[1]-v_target_cms;
								Double_t angle_pi_11B_p2p = p_missing_prov.Angle(v_11B_cms.Vect());
								h1_cos_angle_pi_11B_p2p->Fill(cos(angle_pi_11B_p2p));	
								TVector3 p_miss_prov_v = p_missing_prov.Vect();
								TVector3 p_11b_prov_v = v_11B_cms.Vect();
								TVector3 p_11b_lab_v = v_11B_lab.Vect();
								TLorentzVector p_missing_lab = v_lorentz_PROTONS_lab[0] +v_lorentz_PROTONS_lab[1] + v_target_lab;
								TVector3 p_miss_prov_lab_v = p_missing_lab.Vect();
								h2_p_missing_vs_p_11B_p2p->Fill(p_miss_prov_v.Mag(),p_11b_prov_v.Mag());	
								h2_p_missing_vs_p_11B_p2p_x->Fill(p_miss_prov_v.X(),p_11b_prov_v.X());
								h2_p_missing_vs_p_11B_p2p_y->Fill(p_miss_prov_v.Y(),p_11b_prov_v.Y());
								h2_p_missing_vs_p_11B_p2p_z->Fill(p_miss_prov_v.Z(),p_11b_prov_v.Z());
								h2_pz_11B_vs_protons_lab->Fill(p_11b_lab_v.Z(),p_missing_lab.Z());
								h2_energy_vs_angle_lab->Fill((v_lorentz_PROTONS_lab[0].Theta())*180./PI,v_lorentz_PROTONS_lab[0].E() -938.27208816);
                                h2_energy_vs_angle_lab->Fill((v_lorentz_PROTONS_lab[1].Theta())*180./PI,v_lorentz_PROTONS_lab[1].E() -938.27208816);
								if (v_lorentz_PROTONS_lab[0].Theta() > 0.6108652382 && v_lorentz_PROTONS_lab[1].Theta() > 0.6108652382){
									//cout << "I am in the p2p large angle sector, with eventnr:\t" << evtnr << endl;
                                    h2_p_missing_vs_p_11B_p2p_z_large_angle->Fill(p_miss_prov_v.Z(),p_11b_prov_v.Z());
                                    h1_e_sum_protons_large_angle->Fill(v_lorentz_PROTONS_lab[0].E()+v_lorentz_PROTONS_lab[1].E() -2*938.27208816);
                                    }
                                if (v_lorentz_PROTONS_lab[0].Theta() < 0.6108652382 || v_lorentz_PROTONS_lab[1].Theta() < 0.6108652382){
                                    h2_p_missing_vs_p_11B_p2p_z_small_angle->Fill(p_miss_prov_v.Z(),p_11b_prov_v.Z());
                                    h1_e_sum_protons_small_angle->Fill(v_lorentz_PROTONS_lab[0].E()+v_lorentz_PROTONS_lab[1].E() -2*938.27208816);
                                    }								

								//Cut on the gamma peaks
								if (v_lorentz_GAMMAS_cms.size()){
								if ((v_lorentz_GAMMAS_cms[0][3]>1.9 && v_lorentz_GAMMAS_cms[0][3] < 2.3) || (v_lorentz_GAMMAS_cms[0][3]>4 && v_lorentz_GAMMAS_cms[0][3] < 5.2)){
									h2_p_missing_vs_p_11B_p2p_x_gamma_peaks->Fill(p_miss_prov_v.X(),p_11b_prov_v.X());
									h2_p_missing_vs_p_11B_p2p_y_gamma_peaks->Fill(p_miss_prov_v.Y(),p_11b_prov_v.Y());
									h2_p_missing_vs_p_11B_p2p_z_gamma_peaks->Fill(p_miss_prov_v.Z(),p_11b_prov_v.Z());
									h2_p_missing_vs_p_11B_p2p_gamma_peaks->Fill(p_miss_prov_v.Mag(),p_11b_prov_v.Mag());

									
									}
									}
					
								//plot p1_z + p2_z in lab frame
								h1_p1_p2_mom_z_lab->Fill((v_lorentz_PROTONS_lab[0]).Pz()+(v_lorentz_PROTONS_lab[1]).Pz()); //should peak at around 950 MeV/c
								//plot p1_z + p2_z in 12C CMS frame
								h1_p1_p2_mom_z_cms->Fill((v_lorentz_PROTONS_cms[0]).Pz()+(v_lorentz_PROTONS_cms[1]).Pz()); //should peak at around  -950 MeV/c
								//plot p1_x + p2_x in lab frame
								h1_p1_p2_mom_x_lab->Fill((v_lorentz_PROTONS_lab[0]).Px()+(v_lorentz_PROTONS_lab[1]).Px()); //should peak at round 0 MeV/c
								//plot p1_y + p2_y in lab frame
								h1_p1_p2_mom_y_lab->Fill((v_lorentz_PROTONS_lab[0]).Py()+(v_lorentz_PROTONS_lab[1]).Py()); //should peak at round 0 MeV/c
								//plot 11B p_z in lab frame
								h1_pz_11B_lab->Fill(v_11B_lab.Pz()); //should peak at around 10255 MeV/c
								//plot 11B p_z in 12C CMS frame
								h1_pz_11B_cms->Fill(v_11B_cms.Pz()); //should peak at around 0 MeV/c
								//plot 11B p_x in lab frame
								h1_px_11B_lab->Fill(v_11B_lab.Px()); //should peak at around 0 MeV/c
								//plot 11B p_y in lab frame
								h1_py_11B_lab->Fill(v_11B_lab.Py()); //should peak at around 0 MeV/c	
								//plot 11B px vs p1_x + p2_x in lab frame
								h2_px_11B_vs_protons_lab->Fill(v_11B_lab.Px(),(v_lorentz_PROTONS_lab[0]).Px()+(v_lorentz_PROTONS_lab[1]).Px()); //should be linear, no offset, slope ~1
								//plot 11B py vs p1_y + p2_y in lab frame
								h2_py_11B_vs_protons_lab->Fill(v_11B_lab.Py(),(v_lorentz_PROTONS_lab[0]).Py()+(v_lorentz_PROTONS_lab[1]).Py()); //should be linear, no offset, slope ~1
							
								TLorentzVector p_missing = v_lorentz_PROTONS_cms[0] +v_lorentz_PROTONS_cms[1]-v_target_cms;
								//calculate separation energy according to Exclusive quasi-free proton knockout from oxygen
								//isotopes at intermediate energies
								Double_t separation_energy = (1-gamma_precise)*938.272-gamma_precise*((v_lorentz_PROTONS_lab[0]).E()+(v_lorentz_PROTONS_lab[1]).E() -2*938.272) \
								+ beta_precise*gamma_precise*((v_lorentz_PROTONS_lab[0]).Pz()+(v_lorentz_PROTONS_lab[1]).Pz()) - (v_11B_cms.Px()*v_11B_cms.Px()+v_11B_cms.Py()*v_11B_cms.Py()\
								+v_11B_cms.Pz()*v_11B_cms.Pz())/(2*10255.09731131);
								Double_t separation_energy_valerii = (gamma_precise -1)*938.272 + gamma_precise*((v_lorentz_PROTONS_lab[0]).E()+(v_lorentz_PROTONS_lab[1]).E() -2*938.272) \
								-beta_precise*gamma_precise*((v_lorentz_PROTONS_lab[0]).Pz()+(v_lorentz_PROTONS_lab[1]).Pz())+(pow(p_missing.Px(),2)+pow(p_missing.Py(),2)+pow(p_missing.Pz(),2))/(2*10255.09731131);
								h1_sep_energy_11B_p2p_valerii->Fill(separation_energy_valerii);	
								//cout << "separation energy:\t" << separation_energy << endl;
								//make separation energy for different beam energies .... to test...
								//400 AMeV
								TLorentzVector v_11B_400amev = v_11B_lab;
								Double_t beta_400amev = 0.714549; 
								Double_t gamma_400amev = 1.42942;
								TVector3 v_beta_400amev(0,0,beta_400amev);
								v_11B_400amev.Boost(-v_beta_400amev); 	
								Double_t separation_energy_400amev = (1-gamma_400amev)*938.272-gamma_400amev*((v_lorentz_PROTONS_lab[0]).E()+(v_lorentz_PROTONS_lab[1]).E() -2*938.272) \
                                + beta_400amev*gamma_400amev*((v_lorentz_PROTONS_lab[0]).Pz()+(v_lorentz_PROTONS_lab[1]).Pz()) - (v_11B_400amev.Px()*v_11B_400amev.Px()+v_11B_400amev.Py()*v_11B_400amev.Py()\
                                +v_11B_400amev.Pz()*v_11B_400amev.Pz())/(2*10255.09731131);	
								h1_sep_energy_11B_p2p_400->Fill(separation_energy_400amev);
								//384 AMeV
								TLorentzVector v_11B_384amev = v_11B_lab;
								Double_t beta_384amev = 0.706118; 
								Double_t gamma_384amev = 1.41224;
								TVector3 v_beta_384amev(0,0,beta_384amev);
								v_11B_384amev.Boost(-v_beta_384amev); 	
								Double_t separation_energy_384amev = (1-gamma_384amev)*938.272-gamma_384amev*((v_lorentz_PROTONS_lab[0]).E()+(v_lorentz_PROTONS_lab[1]).E() -2*938.272) \
                                + beta_384amev*gamma_384amev*((v_lorentz_PROTONS_lab[0]).Pz()+(v_lorentz_PROTONS_lab[1]).Pz()) - (v_11B_384amev.Px()*v_11B_384amev.Px()+v_11B_384amev.Py()*v_11B_384amev.Py()\
                                +v_11B_384amev.Pz()*v_11B_384amev.Pz())/(2*10255.09731131);	
								h1_sep_energy_11B_p2p_384->Fill(separation_energy_384amev);


								//375 AMeV
								TLorentzVector v_11B_375amev = v_11B_lab;
								Double_t beta_375amev = 0.701192; 
								Double_t gamma_375amev = 1.40258;
								TVector3 v_beta_375amev(0,0,beta_375amev);
								v_11B_375amev.Boost(-v_beta_375amev); 	
								Double_t separation_energy_375amev = (1-gamma_375amev)*938.272-gamma_375amev*((v_lorentz_PROTONS_lab[0]).E()+(v_lorentz_PROTONS_lab[1]).E() -2*938.272) \
                                + beta_375amev*gamma_375amev*((v_lorentz_PROTONS_lab[0]).Pz()+(v_lorentz_PROTONS_lab[1]).Pz()) - (v_11B_375amev.Px()*v_11B_375amev.Px()+v_11B_375amev.Py()*v_11B_375amev.Py()\
                                +v_11B_375amev.Pz()*v_11B_375amev.Pz())/(2*10255.09731131);	
								h1_sep_energy_11B_p2p_375->Fill(separation_energy_375amev);

								//350 AMeV
								
								TLorentzVector v_11B_350amev = v_11B_lab;
								Double_t beta_350amev = 0.686763; 
								Double_t gamma_350amev = 1.37574;
								TVector3 v_beta_350amev(0,0,beta_350amev);
								v_11B_350amev.Boost(-v_beta_350amev); 	
								Double_t separation_energy_350amev = (1-gamma_350amev)*938.272-gamma_350amev*((v_lorentz_PROTONS_lab[0]).E()+(v_lorentz_PROTONS_lab[1]).E() -2*938.272) \
                                + beta_350amev*gamma_350amev*((v_lorentz_PROTONS_lab[0]).Pz()+(v_lorentz_PROTONS_lab[1]).Pz()) - (v_11B_350amev.Px()*v_11B_350amev.Px()+v_11B_350amev.Py()*v_11B_350amev.Py()\
                                +v_11B_350amev.Pz()*v_11B_350amev.Pz())/(2*10255.09731131);	
								h1_sep_energy_11B_p2p_350->Fill(separation_energy_350amev);
			
								//300 AMeV
								TLorentzVector v_11B_300amev = v_11B_lab;
								Double_t beta_300amev = 0.654117; 
								Double_t gamma_300amev = 1.32206;
								TVector3 v_beta_300amev(0,0,beta_300amev);
								v_11B_300amev.Boost(-v_beta_300amev); 	
								Double_t separation_energy_300amev = (1-gamma_300amev)*938.272-gamma_300amev*((v_lorentz_PROTONS_lab[0]).E()+(v_lorentz_PROTONS_lab[1]).E() -2*938.272) \
                                + beta_300amev*gamma_300amev*((v_lorentz_PROTONS_lab[0]).Pz()+(v_lorentz_PROTONS_lab[1]).Pz()) - (v_11B_300amev.Px()*v_11B_300amev.Px()+v_11B_300amev.Py()*v_11B_300amev.Py()\
                                +v_11B_300amev.Pz()*v_11B_300amev.Pz())/(2*10255.09731131);	
								h1_sep_energy_11B_p2p_300->Fill(separation_energy_300amev);

								//425 AMeV
								TLorentzVector v_11B_425amev = v_11B_lab;
								Double_t beta_425amev = 0.726948; 
								Double_t gamma_425amev = 1.45626;
								TVector3 v_beta_425amev(0,0,beta_425amev);
								v_11B_425amev.Boost(-v_beta_425amev); 	
								Double_t separation_energy_425amev = (1-gamma_425amev)*938.272-gamma_425amev*((v_lorentz_PROTONS_lab[0]).E()+(v_lorentz_PROTONS_lab[1]).E() -2*938.272) \
                                + beta_425amev*gamma_425amev*((v_lorentz_PROTONS_lab[0]).Pz()+(v_lorentz_PROTONS_lab[1]).Pz()) - (v_11B_425amev.Px()*v_11B_425amev.Px()+v_11B_425amev.Py()*v_11B_425amev.Py()\
                                +v_11B_425amev.Pz()*v_11B_425amev.Pz())/(2*10255.09731131);	
								
								h1_sep_energy_11B_p2p_425->Fill(separation_energy_425amev);
								h1_sep_energy_11B_p2p->Fill(separation_energy);
								Double_t own_sep_energy = (1-gamma_precise)*938.272-gamma_precise*((v_lorentz_PROTONS_lab[0]).E()+(v_lorentz_PROTONS_lab[1]).E() -2*938.272) \
                                + beta_precise*gamma_precise*((v_lorentz_PROTONS_lab[0]).Pz()+(v_lorentz_PROTONS_lab[1]).Pz()) + (pow(p_missing.Px(),2)+pow(p_missing.Py(),2)+pow(p_missing.Pz(),2))/(2*938.272);
								h1_sep_energy_11B_p2p_own->Fill(own_sep_energy);
								//now I use my separation energy formula with beam energy = 384
								TLorentzVector p_missing_384 = v_lorentz_PROTONS_lab[0] + v_lorentz_PROTONS_lab[1] - v_target_lab;
								p_missing_384.Boost(-v_beta_384amev);
								 Double_t own_sep_energy_384mev = (1-gamma_384amev)*938.272-gamma_384amev*((v_lorentz_PROTONS_lab[0]).E()+(v_lorentz_PROTONS_lab[1]).E() -2*938.272) \
                                + beta_384amev*gamma_384amev*((v_lorentz_PROTONS_lab[0]).Pz()+(v_lorentz_PROTONS_lab[1]).Pz()) + (pow(p_missing_384.Px(),2)+pow(p_missing_384.Py(),2)+pow(p_missing_384.Pz(),2))/(2*938.272);
								h1_sep_energy_11B_p2p_own_384_mev->Fill(own_sep_energy_384mev);	
								}
								//TJ-END check
								Double_t biggest_proton_e = v_lorentz_PROTONS_cms[0][3];
								Double_t big_proton_e = v_lorentz_PROTONS_cms[1][3];
								TLorentzVector biggest_proton_lab = v_lorentz_PROTONS_cms[0];
								TLorentzVector big_proton_lab = v_lorentz_PROTONS_cms[1];
								biggest_proton_lab.Boost(v_beta_12c);
								big_proton_lab.Boost(v_beta_12c);	
								h2_E_sum_vs_theta_sum->Fill(biggest_proton_e+big_proton_e-2*938,(biggest_proton_lab.Theta()+big_proton_lab.Theta())*(180./PI));
								//Definition of four-vector p_missing: in cms frame
								TLorentzVector p_missing = v_lorentz_PROTONS_cms[0] +v_lorentz_PROTONS_cms[1]-v_target_cms;
								//Definition of missing Energy = m_proton - e_miss (where e_miss is equal to the energy component of p_missing)
								Double_t missing_energy = 938.272 - p_missing.E();
								//Reconstructing missing mass as for src 10B (where the reconstructed mass was equal to the neutron mass)
								//here I expect a peak at around 0...
								TLorentzVector p_12C_cms(0,0,0,11177.92337806);
								TLorentzVector tl_reco = p_12C_cms + v_target_cms -  v_lorentz_PROTONS_cms[0] - v_lorentz_PROTONS_cms[1] - v_11B_cms;		
								//construct momentum of p_missing:
								Double_t test_mom = sqrt(p_missing.Px()*p_missing.Px()+p_missing.Py()*p_missing.Py()+p_missing.Pz()*p_missing.Pz());
								Double_t test_mom_11b = sqrt(v_11B_cms.Px()*v_11B_cms.Px() + v_11B_cms.Py()*v_11B_cms.Py() + v_11B_cms.Pz()*v_11B_cms.Pz());
								Double_t M_missing = tl_reco.Mag();	
								h1_missing_mass_11B->Fill(M_missing);
								//TJ here I insert and compare the momenta of the participants
								//TLorentzVector t_output =  v_lorentz_PROTONS_cms[0] + v_lorentz_PROTONS_cms[1] + v_11B_cms;
								//TLorentzVector t_input = p_12C_cms + v_target_cms;
								//Double_t input_momentum = sqrt(t_input.Px()*t_input.Px()+t_input.Py()*t_input.Py()+t_input.Pz()*t_input.Pz());
								//Double_t output_momentum = sqrt(t_output.Px()*t_output.Px()+t_output.Py()*t_output.Py()+t_output.Pz()*t_output.Pz());
								//cout << "input momentum:\t" << input_momentum << endl;
								//cout << "output momentum:\t" << output_momentum << endl;
								//cout << "combined energy of the two protons" << v_CALIFA_evts_lab_energy[0][0]+v_CALIFA_evts_lab_energy[1][0] << endl;
								if ((v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1]) > 1.308996939 && (v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1]) < 1.4835298642 \
									&& abs(v_CALIFA_evts_cms_energy[0][2]-v_CALIFA_evts_cms_energy[1][2]) < 3.6651914292 && abs(v_CALIFA_evts_cms_energy[0][2]-v_CALIFA_evts_cms_energy[1][2]) > 2.617993878){
									h1_missing_mass_11B_peak_theta->Fill(M_missing);  	
								}
								//TJ--END
								h1_missing_energy_11B_p2p->Fill(missing_energy);
								h2_e_miss_vs_theta_sum_p2p_11B->Fill(missing_energy,((v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1])/PI)*180.);
								//Filling histo: missing Energy(cms) vs theta1+theta2
								//Fillin histo: missing momentum (cms) where missing_mom is equal to the magnitude of the 3-vector part of p_missing
								//TJ_NOTE: check that you do not get artifacts from clustering... therefore arzimuthal angle > 25degrees
								if(abs(v_CALIFA_evts_lab_energy[0][2] - v_CALIFA_evts_lab_energy[1][2]) > 0.436332313){
									h1_p_missing_11B->Fill(sqrt(pow(p_missing[0],2)+pow(p_missing[1],2)+pow(p_missing[2],2))); //mom_miss = sqrt(p_x^2+p_y^2+p_z^2)
								}
								//Filling missing momentum histo  for x and z component:
								h1_p_missing_11B_z->Fill(p_missing[2]);
								h1_p_missing_11B_x->Fill(p_missing[0]);	
								//As we do just have the x,z info for the 11B fragment (MWPCs) we also need to project the p1 and p2 protons to the x,z plane for making \
								//some angular distribution analysis.... :
								TLorentzVector proton_1_x_z = v_lorentz_PROTONS_cms[0];
								//now setting the y contribution to 0:
								proton_1_x_z.SetPy(0.);
								//same for proton 2:
								TLorentzVector proton_2_x_z = v_lorentz_PROTONS_cms[1];
								proton_2_x_z.SetPy(0.);
								//p_missing in the CMS x,z plane...:
								TLorentzVector p_missing_x_z = proton_1_x_z + proton_2_x_z - v_target_cms;
								
								//Angle between p_missing_x_z (3-vec) and 11B (3-vec) in the x-z plane
								Double_t angle_pi_11B = p_missing_x_z.Angle(v_11B_cms.Vect());
								//filling histogram showing angle between p_i and 11B in CMS
								h1_angle_pi_11B->Fill(angle_pi_11B*(180/PI));
								h1_angle_pi_11B_cos->Fill(cos(angle_pi_11B));

								//making cut on the QE (p,2p) reaction: theta1 + theta2 < 90° and phi1+phi2 = 180° +- 40°
								if (v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1] < PI/2. && abs(v_CALIFA_evts_cms_energy[0][2]-v_CALIFA_evts_cms_energy[1][2]) < 3.83972 \
									&& abs(v_CALIFA_evts_cms_energy[0][2]-v_CALIFA_evts_cms_energy[1][2]) > 2.44346){
									//Short Analysis for Time Stitching etc...
									//Int_t entries_wr = DataCA->GetEntries();
									//if (entries_califa == 2 && entries_wr){
									if (entries_califa == 2){
										uint64_t wr_ts = DataCA->GetTimeStamp();	
								//		//wrmasterdata = new R3BEventHeader*[entries_wr];
								//		//wrmasterdata[0] = (R3BEventHeader*)DataCA->At(0);
								//		//get time difference between frst WRCALIFA and WRMaster
								//		h1_wr_first_wr_master->Fill((califahitdata[0]->GetTime())-(wrmasterdata[0]->GetTimeStamp()));
								//		//get time difference between second WRCALIFA and WRMaster
								//		h1_wr_second_wr_master->Fill((califahitdata[1]->GetTime())-(wrmasterdata[0]->GetTimeStamp()));

										h1_wr_first_wr_master->Fill((califahitdata[0]->GetTime())-wr_ts);
										h1_wr_second_wr_master->Fill((califahitdata[1]->GetTime())- wr_ts);

										}
									else
										//cout << "NOT ENOUGH WR entries for eventnr.:\t" << evtnr << "   " << "something went wrong..." << endl;
									



									//--------------------------------------

									h2_missing_energy_vs_t1_t2_11B_p2p->Fill(missing_energy,(v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1])*(180/PI));
									h1_missing_energy_11B_p2p_cut->Fill(missing_energy);	
									h1_p_missing_11B_z_cut->Fill(p_missing[2]);
									h1_p_missing_11B_y_cut->Fill(p_missing[1]);
									h1_p_missing_11B_x_cut->Fill(p_missing[0]);
									//same but with the 11B momentum in the 12C frame:
									h1_mom11B_z_cut->Fill(v_11B_cms.Pz());
									h1_mom11B_x_cut->Fill(v_11B_cms.Px());
									h1_mom11B_y_cut->Fill(v_11B_cms.Py());
									h1_mom11B_cut->Fill(sqrt(pow(v_11B_cms.Pz(),2)+pow(v_11B_cms.Px(),2)+pow(v_11B_cms.Py(),2)));
									h1_p_missing_11B_cut->Fill(sqrt(pow(p_missing[0],2)+pow(p_missing[1],2)+pow(p_missing[2],2)));
									//TJ insert here
									int tpat = DataCA->GetTpat();

									Int_t tpatbin;
									for (Int_t i = 0; i < 16; i++){
										tpatbin = (DataCA->GetTpat() & (1 << i));
										if (tpatbin != 0){
											h2_p_missing_11B_vs_tpat->Fill(sqrt(pow(p_missing[0],2)+pow(p_missing[1],2)+pow(p_missing[2],2)),(double) (i+1));
											}
										}
									//simulate neutron evaporation after p2p scattering
									TVector3 uni_vec = uniform_vec();
									while(uni_vec.Mag() > 1){
									    uni_vec = uniform_vec();
									    }
									TRandom* neutron_gen = new TRandom();
									neutron_gen->SetSeed();
									//Double_t neutron_mom_val = neutron_gen->Gaus(240,20); around 30MeV energy, gaussian distribution, sigma 30 MeV/c
									//Double_t neutron_mom_val = neutron_gen->Uniform(140,240); //uniform distribution above the neutron separation energy at around 10 MeV
									//now I try to do it with a landau distribution as from paper: https://reader.elsevier.com/reader/sd/pii/0029558262901050?token=E8B0A60F8BE1FF2F6EFD8D406BF8FC2D68ED2418D6CCC2BED4FF93C76866E0F2CB6227906785ECE39039C88448121E97&originRegion=eu-west-1&originCreation=20220720091318
									//
									Double_t neutron_mom_val = -1000;
									while (neutron_mom_val < 0){
										neutron_mom_val = neutron_gen->Landau(45,0.8);
										}
									TVector3 neutron_vec(neutron_mom_val*sin(uni_vec.Theta())*cos(uni_vec.Phi()),neutron_mom_val*sin(uni_vec.Theta())*sin(uni_vec.Phi()),neutron_mom_val*cos(uni_vec.Theta()));
									TVector3 neutron_evap_plus_p_miss(neutron_vec[0]+p_missing[0],neutron_vec[1]+p_missing[1],neutron_vec[2]+p_missing[2]);
									//Fill Histogram
									h1_p_missing_neutron_evap_sim->Fill(neutron_evap_plus_p_miss.Mag());
									//delete TRandom
									delete neutron_gen;
									
									//simpulate p,ppn scattering, actually not really correct way to do it...
									p_missing_new = p_missing;
									if (p_missing_old[3] != 0){
										//adding up momentum vectors
										TVector3 added_miss_mom_vec(p_missing_new[0]+p_missing_old[0],p_missing_new[1]+p_missing_old[1],p_missing_new[2]+p_missing_old[2]);
										//Fill the histogram with added_miss_mom_vec.Mag()
										h1_p_missing_double_scattering->Fill(added_miss_mom_vec.Mag());
										}
									p_missing_old = p_missing_new;

									//TJ end here insertion	
									h2_p_missing_vs_p_11B_cut->Fill(sqrt(pow(p_missing[0],2)+pow(p_missing[1],2)+pow(p_missing[2],2)),sqrt(pow(v_11B_cms[0],2)+pow(v_11B_cms[1],2)+pow(v_11B_cms[2],2)));
									h1_missing_mass_11B_p2p_cut->Fill(p_missing.M());
									h1_sep_12cp2p_cut->Fill(938.27208816-p_missing.M());
									//cout << "This is event:\t" << evtnr << "\t with missing mass:\t" << p_missing.M() << endl;
									//analyze if both protons hit IPhos or Barrel region...
									Double_t theta_0 = v_CALIFA_evts_lab_energy[0][1];
									Double_t theta_1 = v_CALIFA_evts_lab_energy[1][1];
									if (theta_0 < 0.75049157836 && theta_1 < 0.75049157836){
										//cout << (*angle_comb)[0] << endl;
										//cout << angle_comb[0] << endl;
										h1_reco_mass_11B_p2p_cut_two_iPhos->Fill(p_missing.M());
										
										h2_reco_mass_angle_comb_cut->Fill(p_missing.M(),0);
										}
									else if ((theta_0 < 0.75049157836 && theta_1 > 0.75049157836) || (theta_0 > 0.75049157836 && theta_1 < 0.75049157836)){
										//cout << (*angle_comb)[1] << endl;
										//cout << angle_comb[1] << endl;
										h1_reco_mass_11B_p2p_cut_iPhos_barrel->Fill(p_missing.M());
										h2_reco_mass_angle_comb_cut->Fill(p_missing.M(),1);
										if (theta_0 > theta_1){
											h2_reco_mass_11B_iphos_angle->Fill(p_missing.M(),(theta_1*180.)/PI);
											h2_reco_mass_11B_barrel_angle->Fill(p_missing.M(),(theta_0*180.)/PI);
											}
										else {
											h2_reco_mass_11B_iphos_angle->Fill(p_missing.M(),(theta_0*180.)/PI);
                                            h2_reco_mass_11B_barrel_angle->Fill(p_missing.M(),(theta_1*180.)/PI);
		
											}
										}
									else if (theta_0 > 0.75049157836 && theta_1 > 0.75049157836){
										//cout << (*angle_comb)[2] << endl;
										//cout << angle_comb[2] << endl;
										h1_reco_mass_11B_p2p_cut_two_barrel->Fill(p_missing.M());
										h2_reco_mass_angle_comb_cut->Fill(p_missing.M(),2);
										}
									//----------------------------------------------------
									h1_excit_energy_11B_cut->Fill((p_12C_cms + v_target_cms - v_lorentz_PROTONS_cms[0]- v_lorentz_PROTONS_cms[1]).M() - 10255.09731131);
									//here I compare the excitation energy and missing energy formula....
									h2_missing_e_2meth_11B_p2p_cut->Fill(missing_energy,(p_12C_cms + v_target_cms - v_lorentz_PROTONS_cms[0]- v_lorentz_PROTONS_cms[1]).M() - 10255.09731131);
					
									int rand_num;
                                	rand_num = rand() % 2 ;

                                	//computing P_y as Valerii does in his thesis (page 66)
                                	Double_t P_y_knocked_out = sqrt(pow(v_lorentz_PROTONS_lab[rand_num][0],2)+pow(v_lorentz_PROTONS_lab[rand_num][1],2)+pow(v_lorentz_PROTONS_lab[rand_num][2],2))*sin(v_lorentz_PROTONS_lab[rand_num].Theta())*sin(v_lorentz_PROTONS_lab[rand_num].Phi()-v_lorentz_PROTONS_lab[1-rand_num].Phi());
									Double_t P_fragment = v_11B_lab.Py(); 
									h2_corr_y_proton_fragm->Fill(P_fragment,P_y_knocked_out);
									h2_corr_x_proton_fragm->Fill(v_11B_lab.Px(),P_y_knocked_out);
									//not fully convinced from the above formula... what I get as formula P_y = Q_k*sin(theta_k)*sin(phi_k) - Q_i*sin(theta_i)*sin(phi_i) 
									//where k is the knocked out proton in the 12C and i is the target proton (at rest in laboratory frame)
									Double_t P_y_knocked_out_tj = sqrt(pow(v_lorentz_PROTONS_lab[rand_num][0],2)+pow(v_lorentz_PROTONS_lab[rand_num][1],2)+pow(v_lorentz_PROTONS_lab[rand_num][2],2))*sin(v_lorentz_PROTONS_lab[rand_num].Theta())*sin(v_lorentz_PROTONS_lab[rand_num].Phi()) -sqrt(pow(v_lorentz_PROTONS_lab[1-rand_num][0],2)+pow(v_lorentz_PROTONS_lab[1-rand_num][1],2)+pow(v_lorentz_PROTONS_lab[1-rand_num][2],2))*sin(v_lorentz_PROTONS_lab[1-rand_num].Theta())*sin(v_lorentz_PROTONS_lab[1-rand_num].Phi());
									h2_corr_y_proton_fragm_tj->Fill(P_fragment,P_y_knocked_out_tj);	
								//by setting the y-value as above for the particles to 0 (see: proton_2_x_z, p_missing_x_z,...) 
								//we get an anguar distortion. Therefore I want to use here p_missing instead and together with 
								//the angular cut compute cos(p_i and p_11B) again:
								Double_t angle_pi_11B_3d = p_missing.Angle(v_11B_cms.Vect());
								h1_angle_pi_11B_cos_cut->Fill(cos(angle_pi_11B_3d));	
								//in the paper they calculate the angle in the x-y plane, so I will do...
								TVector3 p3_missing_x_y(p_missing.Px(),p_missing.Py(),0);
								TVector3 p3_11b_x_y(v_11B_cms.Px(),v_11B_cms.Py(),0);
								h1_angle_pi_11B_cos_cut_plane->Fill(cos(p3_missing_x_y.Angle(p3_11b_x_y)));
								if (cos(p3_missing_x_y.Angle(p3_11b_x_y)) < 0){
									h2_corr_y_proton_fragm_cos_less->Fill(v_11B_lab.Py(),P_y_knocked_out);
									}
								if (cos(p3_missing_x_y.Angle(p3_11b_x_y)) > 0){
									h2_corr_y_proton_fragm_cos_more->Fill(v_11B_lab.Py(),P_y_knocked_out);
									}
								//make plot h1_angle_pi_11B_cos_cut_plane mroe restictive: 77 < theta < 85 and E_miss > -60:
								if ((v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1]) > 1.343903524 && (v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1]) < 1.4835298642 && missing_energy > -60.){
									h1_angle_pi_11B_cos_hcut_plane->Fill(cos(p3_missing_x_y.Angle(p3_11b_x_y)));	
									} 

								//another test:
								//check arzimuth angle of p1+p2, if positivem,give b11 an angle in y-direction of -0.03 rad, otherwise vice versa...
								Double_t p1_plus_p2_phi = (v_lorentz_PROTONS_cms[0] +v_lorentz_PROTONS_cms[1]).Phi();
									if (p1_plus_p2_phi > 0){
										TLorentzVector fake_y(0,tot_mom11B*sin(-0.03),0,0);
										TLorentzVector v_11B_cms_fake_y = (v_11B_lab + fake_y);	
										v_11B_cms_fake_y.Boost(-v_beta_12c);
										Double_t fake_angle = p_missing.Angle(v_11B_cms_fake_y.Vect());
										h1_angle_pi_11B_cos_cut_fake_y->Fill(cos(fake_angle));			
										h2_p_missing_vs_p_11B_cut_fake_y->Fill(sqrt(pow(p_missing[0],2)+pow(p_missing[1],2)+pow(p_missing[2],2)),sqrt(pow(v_11B_cms_fake_y[0],2)+pow(v_11B_cms_fake_y[1],2)+pow(v_11B_cms_fake_y[2],2)));
										}
									else if (p1_plus_p2_phi < 0){
										TLorentzVector fake_y(0,tot_mom11B*sin(0.03),0,0);
										TLorentzVector v_11B_cms_fake_y = (v_11B_lab + fake_y);
										v_11B_cms_fake_y.Boost(-v_beta_12c);
										Double_t fake_angle = p_missing.Angle(v_11B_cms_fake_y.Vect());
										h1_angle_pi_11B_cos_cut_fake_y->Fill(cos(fake_angle));	
										}
								
								
									}
								
								//-------------------------------------------------------------------------------		
								
							if (v_lorentz_GAMMAS_cms.size()){
							//	//h1_gamma_energyE_max_val_11B->Fill(v_lorentz_GAMMAS_cms[0][3]);
							//	h1_gamma_energyE_max_val_11B_lab->Fill(v_lorentz_GAMMAS_lab[0][3]);
							//	//checking multiplicity...
							//	h2_gamma_energy_11B_vs_mult->Fill(v_lorentz_GAMMAS_cms[0][3],v_lorentz_GAMMAS_cms.size());
							//	//check angle distribution (in lab frame) vs largest energy in cm frame
							//	TLorentzVector highest_gamma_temp = v_lorentz_GAMMAS_cms[0];
							//	highest_gamma_temp.Boost(-v_beta_12c);
							//	h2_gamma_energy_11_vs_angle->Fill(v_lorentz_GAMMAS_cms[0][3],highest_gamma_temp.Theta()*(180./PI));
							//	//summing up the total energy deposited in CALIFA. As angle take the one with highest energy...
							//	Double_t sum_energy_gamma = highest_gamma_temp[3];
							//	for (int g = 1; g < v_lorentz_GAMMAS_cms.size();++g){
							//		TLorentzVector temp_gamma_v_lorentz = v_lorentz_GAMMAS_cms[g];
							//		v_lorentz_GAMMAS_cms[g].Boost(-v_beta_12c);
							//		sum_energy_gamma += v_lorentz_GAMMAS_cms[g][3];
							//	}
							//	Double_t doppl_corr_summ_gamma = sum_energy_gamma*(1/sqrt(1-beta_given*beta_given))*(1-beta_given*cos(highest_gamma_temp.Theta()));
							//	h1_gamma_energyE_sum_11B->Fill(doppl_corr_summ_gamma);
									
								//making cut on the QE (p,2p) reaction: theta1 + theta2 < 90° and phi1+phi2 = 180° +- 40°
								if (v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1] < PI/2. && abs(v_CALIFA_evts_cms_energy[0][2]-v_CALIFA_evts_cms_energy[1][2]) < 3.83972  && abs(v_CALIFA_evts_cms_energy[0][2]-v_CALIFA_evts_cms_energy[1][2]) > 2.44346){

									//TJ insertion
									//

								if (v_lorentz_GAMMAS_cms.size()){	
								vector< vector<double > > vect_gamma_energy_angle;
								for (Int_t k = 0; k < v_lorentz_GAMMAS_cms.size();k++){
									Double_t energy_gamma_cms = v_lorentz_GAMMAS_cms[k].E();
									v_lorentz_GAMMAS_cms[k].Boost(v_beta_12c);
									Double_t theta_gamma_lab = (v_lorentz_GAMMAS_cms[k].Theta())*180/PI;
									vector<Double_t> temp_vec{energy_gamma_cms,theta_gamma_lab};
									vect_gamma_energy_angle.push_back(temp_vec);
									//Fill histograms with all gammas:
									h1_gamma_energyE_11B->Fill(energy_gamma_cms);
									h2_gamma_energyE_vs_angle_11B->Fill(theta_gamma_lab,energy_gamma_cms);
									temp_vec.clear();
								}
								sort(vect_gamma_energy_angle.begin(),vect_gamma_energy_angle.end(),sortcol_pseudo);
								h1_gamma_energyE_max_val_11B->Fill(vect_gamma_energy_angle[0][0]);
								h1_gamma_energyE_max_val_11B_fine_bin->Fill(vect_gamma_energy_angle[0][0]);
								h2_gamma_energyE_max_vs_angle_11B->Fill(vect_gamma_energy_angle[0][1],vect_gamma_energy_angle[0][0]);
								}


									//TJ END OF INSERTION
										
									h1_gamma_energyE_max_val_11B_cut->Fill(v_lorentz_GAMMAS_cms[0][3]);	
									h2_gamma_energy_11B_vs_mult_cut->Fill(v_lorentz_GAMMAS_cms[0][3],v_lorentz_GAMMAS_cms.size());
									//h1_gamma_energyE_sum_11B_cut->Fill(doppl_corr_summ_gamma);
									//boosting back to 11B frame instead as incorrectly to the 12C frame...
									TLorentzVector tl_highest_g_11B = v_lorentz_GAMMAS_lab[0];
									tl_highest_g_11B.Boost(-sin(psi_in)*beta_precise,0,-cos(psi_in)*beta_precise);
									h1_gamma_energyE_max_val_11Bframe_cut->Fill(tl_highest_g_11B[3]);
									//cut on multiplicity < 3
									if (v_lorentz_GAMMAS_cms.size() < 3){
										h1_gamma_energyE_max_val_11B_m_cut->Fill(v_lorentz_GAMMAS_cms[0][3]);
									}
									//TJ start look at gamma spectrum with high x and y momenta of p_i and p11b in 12C frame
									if (abs(p_missing.Py()) > 100 && abs(v_11B_cms.Py()) > 100){
										h1_gamma_spec_high_y_comp->Fill(v_lorentz_GAMMAS_cms[0][3]);	
										}
									if (abs(p_missing.Px()) > 100 && abs(v_11B_cms.Px()) > 100){
										h1_gamma_spec_high_x_comp->Fill(v_lorentz_GAMMAS_cms[0][3]);
										}
									if (sqrt(pow(p_missing[0],2)+pow(p_missing[1],2)+pow(p_missing[2],2)) > 500){
										h1_gamma_spec_high_missing->Fill(v_lorentz_GAMMAS_cms[0][3]);
										}
									if (sqrt(pow(p_missing[0],2)+pow(p_missing[1],2)+pow(p_missing[2],2)) > 400){
										h1_gamma_spec_high_missing_400->Fill(v_lorentz_GAMMAS_cms[0][3]);
										}
									if (sqrt(pow(p_missing[0],2)+pow(p_missing[1],2)+pow(p_missing[2],2)) > 300){
										h1_gamma_spec_high_missing_300->Fill(v_lorentz_GAMMAS_cms[0][3]);
										}

									//TJ END
								}
								//--------------------------------------------------------------------------------------

								}	
								
							if (v_lorentz_PROTONS_lab.size() > 2 && v_CALIFA_evts_lab_energy[2][0] > 30){
								
								//cout << "more than 3 protons, really? " << endl;
								}
								}
							}
						
					}
                        }

//############################################################################################################
////------------END OF 11B ANALYSIS-------------------------------------------------------------------------
////##########################################################################################################

//############################################################################################################
//--------------BEGIN OF 10B ANALYSIS-------------------------------------------------------------------------
//############################################################################################################
                        if (a_q < s[0] && a_q > (s[0]-0.15) && entries_califa >= 2 && charge_val< t && charge_val> t-1.){  //10B
                        //psi_in = mean_psi_out_10b-slope_10b*xMW3 - angle_offs_10b;
						psi_in = u[0]-u[1]*xMW3 -u[2];
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
						Double_t gamma_precise = 1/(sqrt(1-beta_precise*beta_precise));
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;

						//compute position direction of fragment after the target using the target pos. and the MW1/2 (x-pos) and MWPC3 (y-pos.)
						////--> dangerous!! set fragment y to -y, TODO: change it back after usage
						fragment_y = -(yMW3+yMW3_shift)/(pathlength_precise-760.4); //substract distance MW3 to TOFW
						fragment_x =  v_diff_mw2_mw1.X()/v_diff_mw2_mw1.Mag();
						fragment_z = v_diff_mw2_mw1.Z()/v_diff_mw2_mw1.Mag();

						//TJ check position mw2_x and mw3_y
						h1_x_position_mw2_10b->Fill(xMW2+xMW2_shift);
						h1_y_position_mw3_10b->Fill(yMW3+yMW3_shift);
						//beginning after isotope has been selected:
						//TLorentzVector of 10B in the lab and cms 
						Double_t tot_mom10B = 9326.98688128*(beta_precise/(sqrt(1-(beta_precise*beta_precise)))); //now with correct mass defect
						//TODO: check this formula, it doesn't seem to be correct, as I use the lorentz factor twice...
						//Double_t energy_10B = sqrt(pow(tot_mom10B,2)+pow(10*938.272,2));
						Double_t energy_10B = (1./(sqrt(1-(beta_precise*beta_precise))))*9326.98688128; //now with correct mass defect
						TLorentzVector v_10B_lab(tot_mom10B*fragment_x,tot_mom10B*fragment_y,tot_mom10B*fragment_z,energy_10B);
						TLorentzVector v_10B_cms = v_10B_lab;
						v_10B_cms.Boost(-v_beta_12c);
						h1_px_10B_lab_test->Fill(v_10B_lab.Px());
						h1_py_10B_lab_test->Fill(v_10B_lab.Py());
						//cout << "this is the charge of 10b:\t" << charge_val << endl;
						//TLorentzVector of the target in lab and cms frame
						TLorentzVector v_target_lab(0,0,0,938.272);
						TLorentzVector v_target_cms = v_target_lab;
						v_target_cms.Boost(-v_beta_12c);
						//fill entries of CALIFA in vectors, no matter if they are protons or gammas..
						//variales needed:
						vector<vector<double> > v_CALIFA_evts_lab_energy; //here I insert vectors v = (E_lab,theta_lab,phi_lab)
						vector<vector<double> > v_CALIFA_evts_cms_energy; //here I inserv vectors v = (E_cms,theta_lab,phi_lab)
						//for latest update I change to: v_CALIFA_evts_cms_energy(E_cms,theta_lab,phi_lab,Nf,Ns)
						vector<TLorentzVector> v_lorentz_CALIFA_evts_lab; //here I insert Lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_lab)
						vector<TLorentzVector> v_lorentz_CALIFA_evts_cms; //here I insert Lorentz vector v = (p_x_cms,p_y_cms,p_z_cms,E_cms)
						vector<TLorentzVector> v_lorentz_PROTONS_lab; //here I inserrt lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_lab) where E_lab > 30MeV
						vector<TLorentzVector> v_lorentz_PROTONS_cms; //here I inserrt lorentz vector v = (p_x_cms,p_y_cms,p_z_cms,E_cms) 
						vector<TLorentzVector> v_lorentz_GAMMAS_lab; //here I inserrt lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_lab) where E_lab < 30MeV
						vector<TLorentzVector> v_lorentz_GAMMAS_cms; //here I inserrt lorentz vector v = (p_x_lab,p_y_lab,p_z_lab,E_cms) 
						TVector3 v_boost(0,0,beta_precise*light_c);
						
						for (Int_t j = 0.;j<entries_califa;j++){
							vector<double> v_temp(3);
							vector<double> v_temp_cms(5);
							califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);
							v_temp[0] = (califahitdata[j]->GetEnergy())/1000.;
							//randomize angles to remove artifacts...	
							//if (califahitdata[j]->GetTheta() > 0.45){
                            //    Int_t rand_angl = rand() % 4500 + 1;
                            //    v_temp[1] = califahitdata[j]->GetTheta() - 0.0225 + 0.00001*rand_angl;
                            //    }
                            //else {
                            //    v_temp[1] = califahitdata[j]->GetTheta();
                            //    }	

							v_temp[1] = califahitdata[j]->GetTheta();
							v_temp[2] = califahitdata[j]->GetPhi();
							v_temp_cms[0] = v_temp[0]*(1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(v_temp[1])); //Doppler correction
							v_temp_cms[1] = v_temp[1];
						    v_temp_cms[2] = v_temp[2];
							v_temp_cms[3] = (califahitdata[j]->GetNf());
							v_temp_cms[4] = (califahitdata[j]->GetNs());
							v_CALIFA_evts_lab_energy.push_back(v_temp);
							v_CALIFA_evts_cms_energy.push_back(v_temp_cms);
							v_temp.clear();
							v_temp_cms.clear();
							
						}
						
						
						for (int CALIFA_events = 0; CALIFA_events < entries_califa; CALIFA_events++){
						
						if(v_CALIFA_evts_lab_energy[CALIFA_events][0] < 30) { //gammas
						//create lorentz vector in the laboratory system
						TLorentzVector temp_lorentz_lab;
						TLorentzVector temp_lorentz_cms;
						temp_lorentz_lab[0] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[1] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[2] = (v_CALIFA_evts_lab_energy[CALIFA_events][0])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][1]); //coordinate transform
						temp_lorentz_lab[3] = v_CALIFA_evts_lab_energy[CALIFA_events][0];
						v_lorentz_CALIFA_evts_lab.push_back(temp_lorentz_lab);
						v_lorentz_GAMMAS_lab.push_back(temp_lorentz_lab);
						//temp_lorentz_lab.Boost(v_boost);
						temp_lorentz_lab.Boost(-v_beta_12c);
						v_lorentz_CALIFA_evts_cms.push_back(temp_lorentz_lab);
						v_lorentz_GAMMAS_cms.push_back(temp_lorentz_lab);
						temp_lorentz_lab = TLorentzVector();
						
						}
						else {//protons
						
						TLorentzVector temp_lorentz_lab;
						temp_lorentz_lab[0] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*cos(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[1] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*sin(v_CALIFA_evts_lab_energy[CALIFA_events][1])*sin(v_CALIFA_evts_lab_energy[CALIFA_events][2]); //coordinate transform
						temp_lorentz_lab[2] = sqrt(pow(((v_CALIFA_evts_lab_energy[CALIFA_events][0])+938.272),2)-pow(938.272,2))*cos(v_CALIFA_evts_lab_energy[CALIFA_events][1]); //coordinate transform
						temp_lorentz_lab[3] = v_CALIFA_evts_lab_energy[CALIFA_events][0]+938.272;
						v_lorentz_CALIFA_evts_lab.push_back(temp_lorentz_lab);
						v_lorentz_PROTONS_lab.push_back(temp_lorentz_lab);
						temp_lorentz_lab.Boost(-v_beta_12c);
						v_lorentz_CALIFA_evts_cms.push_back(temp_lorentz_lab);
						v_lorentz_PROTONS_cms.push_back(temp_lorentz_lab);
						temp_lorentz_lab = TLorentzVector();	
						}
						}
						//now sort the vectors according to the Energy.......
						//3d vectors
						sort(v_CALIFA_evts_cms_energy.begin(),v_CALIFA_evts_cms_energy.end(),sortcol_pseudo);
						sort(v_CALIFA_evts_lab_energy.begin(),v_CALIFA_evts_lab_energy.end(),sortcol_pseudo);
						//lorentz vectors
						sort(v_lorentz_CALIFA_evts_lab.begin(),v_lorentz_CALIFA_evts_lab.end(),sortcol_lorentz);
						sort(v_lorentz_CALIFA_evts_cms.begin(),v_lorentz_CALIFA_evts_cms.end(),sortcol_lorentz);
						sort(v_lorentz_GAMMAS_lab.begin(),v_lorentz_GAMMAS_lab.end(),sortcol_lorentz);
						sort(v_lorentz_GAMMAS_cms.begin(),v_lorentz_GAMMAS_cms.end(),sortcol_lorentz);
						sort(v_lorentz_PROTONS_lab.begin(),v_lorentz_PROTONS_lab.end(),sortcol_lorentz);
						sort(v_lorentz_PROTONS_cms.begin(),v_lorentz_PROTONS_cms.end(),sortcol_lorentz);
						int tpat = DataCA->GetTpat();	
						//make restriction of having two protons (= more than 30 MeV in lab frame...) and psi difference greater than  25 degrees
						if (v_CALIFA_evts_lab_energy.size() >= 2 && v_CALIFA_evts_lab_energy[0][0] > 30 && v_CALIFA_evts_lab_energy[1][0] > 30 && abs(v_CALIFA_evts_lab_energy[0][2] - v_CALIFA_evts_lab_energy[1][2]) > 0.436332313){
								
								h1_mom_abs_10B->Fill(sqrt(v_10B_cms.Px()*v_10B_cms.Px()+v_10B_cms.Py()*v_10B_cms.Py()+v_10B_cms.Pz()*v_10B_cms.Pz()));
								h1_mom_x_10B->Fill(v_10B_cms.Px());
								h1_mom_y_10B->Fill(v_10B_cms.Py());
								h1_mom_z_10B->Fill(v_10B_cms.Pz());
								//Definition of 12C four-vector in its cms frame:
								TLorentzVector p_12C_cms(0,0,0,11177.92337806);
								//Definition of four-vector p_missing: in cms frame
								TLorentzVector p_missing = p_12C_cms + v_target_cms -  v_lorentz_PROTONS_cms[0] - v_lorentz_PROTONS_cms[1] - v_10B_cms;
								//Definition of M_missing is the reaction missing mass (see Nature paper of Valerii)
								Double_t M_missing = p_missing.Mag();
								h1_missing_mass_10B->Fill(M_missing);	
								//Lorentz vector corresponding to the initial projectile proton in the 12C frame
								TLorentzVector p_proton_proj = v_lorentz_PROTONS_cms[0] +v_lorentz_PROTONS_cms[1]-v_target_cms;
								Double_t missing_energy = 938.272 - p_proton_proj.E();
								h2_e_miss_vs_theta_sum_p2p_10B->Fill(missing_energy,((v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1])/PI)*180.);
								//Filling histo: missing momentum (cms) where missing_mom is equal to the magnitude of the 3-vector part of p_proton_proj. It's actually the momentum of the projectile proton before scattering...
								Double_t p_proton_proj_abs_val = sqrt(pow(p_proton_proj[0],2)+pow(p_proton_proj[1],2)+pow(p_proton_proj[2],2));
								h1_p_missing_10B->Fill(p_proton_proj_abs_val); //mom_miss = sqrt(p_x^2+p_y^2+p_z^2)

								//TJ start here with new analysis path...
								h2_px_10B_vs_protons_lab->Fill(v_10B_cms.Px(),p_proton_proj[0]);
								h2_py_10B_vs_protons_lab->Fill(v_10B_cms.Py(),p_proton_proj[1]);	
								

								//now with the real lab-system four-vectors
								h2_px_10B_vs_protons_lab_real->Fill(v_10B_lab.Px(),(v_lorentz_PROTONS_lab[0]).Px()+(v_lorentz_PROTONS_lab[1]).Px());
								h2_py_10B_vs_protons_lab_real->Fill(v_10B_lab.Py(),(v_lorentz_PROTONS_lab[0]).Py()+(v_lorentz_PROTONS_lab[1]).Py());
								//make tvector in x-y plane to get the angle between 10B and pi
								TVector3 p_10B_cms_xy(v_10B_cms.Px(),v_10B_cms.Py(),0.);
								TVector3 p_i_cms_xy(p_proton_proj[0],p_proton_proj[1],0);
								Double_t angle_pi_10b_cms_xy = p_10B_cms_xy.Angle(p_i_cms_xy);
								h1_cos_angle_10b_pi_xy_cms->Fill(cos(angle_pi_10b_cms_xy));	
								Int_t tpatbin;
								for (Int_t i = 0; i < 16; i++){
									tpatbin = (DataCA->GetTpat() & (1 << i));
									if (tpatbin != 0){
										h2_p_missing_10B_vs_tpat->Fill(p_proton_proj_abs_val,(double) (i+1));
										if (i+1 == 4){  //cut on all events wit tpat == 4 , ie. spill on + sofstart + califa + neuland (mainly neutron evaporation.... ?)
											h1_cos_angle_10b_pi_xy_cms_tpat4->Fill(cos(angle_pi_10b_cms_xy));
											}
										if (i+1 == 2){ // cut on all events with tpat == 2, ie. spill on + sofstart + califa
											h1_cos_angle_10b_pi_xy_cms_tpat2->Fill(cos(angle_pi_10b_cms_xy));
											}
										
										}
									}
								h2_mom_p_i_vs_cos->Fill(cos(angle_pi_10b_cms_xy),p_proton_proj_abs_val);
								h2_mom_10b_vs_cos->Fill(cos(angle_pi_10b_cms_xy),sqrt(pow(v_10B_cms.Px(),2) + pow(v_10B_cms.Py(),2) + pow(v_10B_cms.Pz(),2)));
								h2_p_miss_vs_10B_mom->Fill((sqrt(v_10B_cms.Px()*v_10B_cms.Px()+v_10B_cms.Py()*v_10B_cms.Py()+v_10B_cms.Pz()*v_10B_cms.Pz())),p_proton_proj_abs_val);
								h2_p_miss_vs_10B_mom_transv->Fill(sqrt(v_10B_cms.Px()*v_10B_cms.Px()+v_10B_cms.Py()*v_10B_cms.Py()),sqrt(pow(p_proton_proj[0],2)+pow(p_proton_proj[1],2)));
								//now really the usage of the 4-vectors in the lab system
								h1_pz_10B_lab->Fill(v_10B_lab.Pz());
								h1_px_10B_lab->Fill(v_10B_lab.Px());
								h1_py_10B_lab->Fill(v_10B_lab.Py());
								
								if (p_proton_proj_abs_val > 700 && (sqrt(v_10B_cms.Px()*v_10B_cms.Px()+v_10B_cms.Py()*v_10B_cms.Py()+v_10B_cms.Pz()*v_10B_cms.Pz())) < 150){
									h1_10b_p_i_cos_cuts->Fill(cos(angle_pi_10b_cms_xy));
									}
								if (p_proton_proj_abs_val > 950 && (sqrt(v_10B_cms.Px()*v_10B_cms.Px()+v_10B_cms.Py()*v_10B_cms.Py()+v_10B_cms.Pz()*v_10B_cms.Pz())) < 150){
									h1_10b_p_i_cos_cuts_from_sim->Fill(cos(angle_pi_10b_cms_xy));	
									}

								//TJ END with new analysis path...
								
								
								//Filling histo missing_mom (CMS) of the projectile proton before QFS vs the theta angles of the two outgoing protons
								h2_missing_energy_vs_t1_t2_10B->Fill(sqrt(p_proton_proj*p_proton_proj),(v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1])*(180/PI));
								
								//now we want to check the anglular distributions between the proton and neutron in the SRC\
								//therefore we need the CMS momentum of the interacting proton (= 3-vec of p_i, see calculation here...) and the CMS momentum of the \
								//do the cosine of the proton neutron pair and 10B in the x-y plane:
								TVector3 p3_missing_x_y(p_missing.Px(),p_missing.Py(),0.); //neutron vector
								TVector3 p3_p_i_x_y(p_proton_proj.Px(),p_proton_proj.Py(),0.);//initial proton vector 
								TVector3 p3_10B_x_y(v_10B_cms.Px(),v_10B_cms.Py(),0.);
								h1_angle_pi_p_n_xy_plane->Fill(cos(p3_missing_x_y.Angle(p3_p_i_x_y)));//histo p_i vs p_n in xy plane
								//histo of 10B vs src pair in xy plane
								h1_angle_src_vs_10B_cos_xy_plane->Fill(cos(((p3_p_i_x_y-p3_missing_x_y)*0.5).Angle(p3_10B_x_y)));
							    //histo of 10B vs src_pair in 3d
							    h1_angle_src_vs_10B_cos_3d->Fill(cos(((p_missing.Vect()-p_proton_proj.Vect())*0.5).Angle(v_10B_cms.Vect())));
								//use same restriction as in paper + phi_diff > 100 degree 
								if((p_proton_proj.Vect()).Mag() > 350. && ((v_CALIFA_evts_cms_energy[0][1]+v_CALIFA_evts_cms_energy[1][1])/PI)*180. > 63 && missing_energy > -110 && missing_energy < 240){
									h1_angle_src_vs_10B_cos_3d_paper_cut->Fill(cos(((p_missing.Vect()-p_proton_proj.Vect())*0.5).Angle(v_10B_cms.Vect())));	
								}

								//neutron from the SRC (= 3-vector of p_missing); put again everything into x-z plane!
								TLorentzVector p_i_x_z = v_lorentz_PROTONS_cms[0] + v_lorentz_PROTONS_cms[1] - v_target_cms;
								TLorentzVector p_missing_x_z = p_missing;
								p_missing_x_z.SetPy(0.);
								//do the x-z projection
								p_i_x_z.SetPy(0.);
								//also x-z projection for the 10B fragment
								v_10B_cms.SetPy(0.);
								//Angle between p_i (projectile proton(SRC)) and the according neutron 
								Double_t angle_pi_pn_src = p_i_x_z.Angle(p_missing_x_z.Vect());
								//Angle between src pair and 10B fragment:-> now with updated version p_rel
								Double_t angle_src_pair_10B_frag = ((p_i_x_z-p_missing_x_z)*0.5).Angle(v_10B_cms.Vect());
								//Fill histogram with angular distribution between p_i and p_n in x-z plane:
								h1_angle_pi_pn_10B->Fill(angle_pi_pn_src*(180/PI));
								h1_angle_pi_pn_10B_cos->Fill(cos(angle_pi_pn_src));
								h1_angle_src_vs_10B_cos->Fill(cos(angle_src_pair_10B_frag));
								
								TVector3 sum_mom_inner_part_cms = p_proton_proj.Vect() + p_missing.Vect() + v_10B_cms.Vect(); //sum of the momentum vectors
								//of the 12C inner patricles (should be around 0 from momentum conservation)
								h2_sum_mom_vs_cos_p_i_p_n_xy->Fill(cos(p3_missing_x_y.Angle(p3_p_i_x_y)),sum_mom_inner_part_cms.Mag());
								//Histogram theta1 vs theta2 randomizing which is particle 1 and 2
							    int rand_num;
								rand_num = rand() % 2 ;
								
								//looking at the angular distributions of the two protons for different p_missing regions
								if (p_proton_proj_abs_val < 250){
									h2_theta1_theta2_10b_lower250->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
									h2_n_f_n_s_10b_lower250->Fill((v_CALIFA_evts_cms_energy[rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[rand_num][4]/1000.));
									h2_n_f_n_s_10b_lower250->Fill((v_CALIFA_evts_cms_energy[1-rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[1-rand_num][4]/1000.));
									if (tpat == 4){
										h2_theta1_theta2_10b_lower250_neu->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
										}
									}
								if (p_proton_proj_abs_val > 250 && p_proton_proj_abs_val < 650){
									h2_theta1_theta2_10b_lower650->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
									h2_n_f_n_s_10b_lower650->Fill((v_CALIFA_evts_cms_energy[rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[rand_num][4]/1000.));
									h2_n_f_n_s_10b_lower650->Fill((v_CALIFA_evts_cms_energy[1-rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[1-rand_num][4]/1000.));
									if (tpat == 4){
										h2_theta1_theta2_10b_lower650_neu->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
										}

									}
								if (p_proton_proj_abs_val > 650){
									h2_theta1_theta2_10b_higher650->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);	
									h1_angle_pi_p_n_xy_plane_mom_cut->Fill(cos(p3_missing_x_y.Angle(p3_p_i_x_y)));
									h1_angle_src_vs_10B_cos_xy_plane_mom_cut->Fill(cos(((p3_p_i_x_y-p3_missing_x_y)*0.5).Angle(p3_10B_x_y)));
									h2_n_f_n_s_10b_higher650->Fill((v_CALIFA_evts_cms_energy[rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[rand_num][4]/1000.));
									h2_n_f_n_s_10b_higher650->Fill((v_CALIFA_evts_cms_energy[1-rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[1-rand_num][4]/1000.));
									if (tpat == 4){
										h2_theta1_theta2_10b_higher650_neu->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
										}
									}
								
								h2_theta1_theta2_10b->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
								h2_n_f_n_s_10b->Fill((v_CALIFA_evts_cms_energy[rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[rand_num][4]/1000.));
								h2_n_f_n_s_10b->Fill((v_CALIFA_evts_cms_energy[1-rand_num][3])/1000.,(v_CALIFA_evts_cms_energy[1-rand_num][4]/1000.));
								if (tpat == 4){
									h2_theta1_theta2_10b_neu->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
									}
								//making mass cut on the reconstucted neutron mass M_missing 
								if (M_missing > 850 && M_missing < 1100){
									h2_theta1_theta2_10b_m_cut->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
									h1_p_missing_10B_m_cut->Fill(sqrt(pow(p_proton_proj[0],2)+pow(p_proton_proj[1],2)+pow(p_proton_proj[2],2)));
									//here making also cut on the reconstructed mass of the initial proton:
									if(p_proton_proj.Mag() > 844 && p_proton_proj.Mag() < 960){
										h2_theta1_theta2_10b_double_cut->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);
										} 
									}
								
								//here I calculate the mass the mass of the initial proton: p_missing = p1 + p2 -p_tg -> sqrt(pmissing^2) \
								//using also the y component
								Double_t M_mass_rec = p_proton_proj.Mag();
								h1_missing_mass_10B_rec->Fill(M_mass_rec);
								//use M_mass_rec as cut variable to check the momentum distribution...for eg 10B:
								if (M_mass_rec > 830 && M_mass_rec < 970){
									h1_p_missing_10B_p_i_cut->Fill(sqrt(pow(p_proton_proj[0],2)+pow(p_proton_proj[1],2)+pow(p_proton_proj[2],2)));
									}
								//ok, now vice versa, I use the momentum of the initial proton as cut variable. To see scr-pairs I cut on mom_p_i > 900
								if ((sqrt(pow(p_proton_proj[0],2)+pow(p_proton_proj[1],2)+pow(p_proton_proj[2],2))) > 900){
									h2_theta1_theta2_10b_miss_mom_cut->Fill(((v_CALIFA_evts_cms_energy[rand_num][1])/PI)*180,((v_CALIFA_evts_cms_energy[1-rand_num][1])/PI)*180);	
									h1_missing_mass_10B_rec_miss_cut->Fill(M_mass_rec);
									h1_missing_mass_10B_neutron_miss_cut->Fill(M_missing);
									}

							//gamma analysis for 12C(p,2p)10B reaction
							if (v_lorentz_GAMMAS_cms.size()){
                                h1_gamma_energyE_max_val_10B->Fill(v_lorentz_GAMMAS_cms[0][3]);
                                           }
				
								
							if (v_lorentz_PROTONS_lab.size() > 2){
								
								//cout << "more than 3 protons, really? 10B " << endl;
								}
							}
						

                        }
//############################################################################################################
//--------------END OF 10B ANALYSIS---------------------------------------------------------------------------
//############################################################################################################

                        
                        if (charge_val< t+0.8 && charge_val> t) //Z = 6
                        {
                        if (a_q > s[1] && a_q < (s[1] +2*s[5])){  //12C
                        //psi_in = mean_psi_out_12c-slope_12c*xMW3 - angle_offs_12c;
						psi_in = u[9] -u[10]*xMW3 - u[11];
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
                        if (time_target_tof < 35.13){
                        h1_tof_detnr_strange_12c->Fill(i+1);
                        h1_mw3_fy_strange_12c->Fill(raw_mwpc3);
                        h2_tof_path_strange12c->Fill(pathlength_precise,time_target_tof);
                        h2_radius_path_strange12c->Fill(rho,time_target_tof);
                        h1_tofwall_posy_fast[i]->Fill(tof_diff_up_down);
                        }
                        if (time_target_tof> 35.35 && time_target_tof<35.61){
                        h1_tof_detnr_12c->Fill(i+1);
                        h1_mw3_fy_12c->Fill(raw_mwpc3);
                        h2_tof_path_12c->Fill(pathlength_precise,time_target_tof);
                        h2_radius_path_12c->Fill(rho,time_target_tof);
                        h1_tofwall_posy_medium[i]->Fill(tof_diff_up_down);
                        }
                        if (a_q_precise > 1.997){
                        h1_tof_detnr_strange_12c_large_aq->Fill(i+1);
                        h1_mw3_fy_strange_12c_large_aq->Fill(raw_mwpc3);
                        h2_tof_path_strange_12c_large_aq->Fill(pathlength_precise,time_target_tof);
                        h2_radius_path_strange_12c_large_aq->Fill(rho,time_target_tof);
                        h1_tofwall_posy_slow[i]->Fill(tof_diff_up_down);
                        }
                        h2_z_vs_a_q->Fill(a_q_precise,charge_val);
                        h2_tof_vs_aq_fix_g_12c->Fill(a_q_precise,time_target_tof);
                        h2_mw2x_vs_tof_12c->Fill(xMW2,time_target_tof);
                        h2_mw3x_vs_tof_12c->Fill(xMW3,time_target_tof);
                        h2_charge_vs_tof_12c->Fill(charge_val,time_target_tof);
                        h2_psi_in_vs_tof_12c->Fill(psi_in,time_target_tof);
                        h2_detnr_vs_tof_12c->Fill(i+1,time_target_tof);
                        h2_radius_vs_beta_12c->Fill(beta_precise,rho);
                        h2_radius_vs_tof_12c->Fill(rho,time_target_tof);
						//check beta_precise
						h1_beta_precise_12c->Fill(beta_precise);

						//look at the gamma spectrum for 12C, can I see the 4.4MeV 2+ excited state?
						vector<vector<double> > v_CALIFA_evts_lab_energy; //here I insert vectors v = (E_lab,theta_lab,phi_lab)
						vector<vector<double> > v_CALIFA_evts_cms_energy; //here I inserv vectors v = (E_cms,theta_lab,phi_lab)
						for (Int_t j = 0.;j<entries_califa;j++){
							vector<double> v_temp(3);
							vector<double> v_temp_cms(5);
							califahitdata[j] = (R3BCalifaHitData*)CalifaHitData->At(j);
							v_temp[0] = (califahitdata[j]->GetEnergy())/1000.;
							v_temp[1] = califahitdata[j]->GetTheta();
							v_temp[2] = califahitdata[j]->GetPhi();
							v_temp_cms[0] = v_temp[0]*(1/(sqrt(1-beta_precise*beta_precise)))*(1-beta_precise*cos(v_temp[1])); //Doppler correction
							v_temp_cms[1] = v_temp[1];
							v_temp_cms[2] = v_temp[2];
							v_temp_cms[3] = (califahitdata[j]->GetNf());
							v_temp_cms[4] = (califahitdata[j]->GetNs());
							v_CALIFA_evts_lab_energy.push_back(v_temp);
							v_CALIFA_evts_cms_energy.push_back(v_temp_cms);
							v_temp.clear();
							v_temp_cms.clear();
							}
						for (int CALIFA_events = 0; CALIFA_events < entries_califa; CALIFA_events++){
							h1_gamma_spec_lab->Fill(v_CALIFA_evts_lab_energy[CALIFA_events][0]);	
							h1_gamma_spec_cms->Fill(v_CALIFA_evts_cms_energy[CALIFA_events][0]);
							if (v_CALIFA_evts_lab_energy[CALIFA_events][1] < 0.75){
								h1_gamma_spec_lab_iphos->Fill(v_CALIFA_evts_lab_energy[CALIFA_events][0]);
								h1_gamma_spec_cms_iphos->Fill(v_CALIFA_evts_cms_energy[CALIFA_events][0]);
								}
							if (v_CALIFA_evts_lab_energy[CALIFA_events][1] > 0.75 && v_CALIFA_evts_lab_energy[CALIFA_events][1] < 1.047){
								h1_gamma_spec_lab_to60->Fill(v_CALIFA_evts_lab_energy[CALIFA_events][0]);
								h1_gamma_spec_cms_to60->Fill(v_CALIFA_evts_cms_energy[CALIFA_events][0]);
								}
							if (v_CALIFA_evts_lab_energy[CALIFA_events][1] > 1.047){
								h1_gamma_spec_lab_larger60->Fill(v_CALIFA_evts_lab_energy[CALIFA_events][0]);
								h1_gamma_spec_cms_larger60->Fill(v_CALIFA_evts_cms_energy[CALIFA_events][0]);	
								}
							}
						//cout << "entries_wr_califa:\t" << entries_wr_califa << "    " << "entries_califa" << entries_califa <<  "   "   << "eventnr:\t" << evtnr << endl;	
						//TJ end of looking at gamma spectrum from 12C
                        }

                        if (a_q < s[1]-0.06 && a_q > (s[1]-1.65*s[4]) && entries_califa >= 2){  //11C before it was a_q < aq_cut_z_6 && a_q > (aq_cut_z_6-2*width_11c), too bad for 11c...
                        //psi_in = mean_psi_out_11c-slope_11c*xMW3 - angle_offs_11c;
						psi_in = u[6] - u[7]*xMW3 - u[8];
                        m = (cos(psi_out)-cos(psi_in))/(sin(psi_out)+sin(psi_in));
                        rho =  (Leff/(2*sin((psi_in+psi_out)/2)))*sqrt(pow(((m+tan(alpha_G))/(1-m*tan(alpha_G))),2)+1);
                        w = 2*abs(asin(BD/(2*rho)));
                        Double_t pathlength_precise = abs((zB-zT)/(cos(psi_in)))+rho*w+abs((zToFW-z_D)/cos(psi_out));
                        Double_t beta_precise = ((pathlength_precise/time_target_tof)*pow(10,6))/light_c;
                        Double_t a_q_precise = (pow(10,-3)*((((mag_field*rho)/(beta_precise))*1.602176634*pow(10,-19))/(1.66053906660*pow(10,-27)*light_c)))/gamma_given;
                        h2_z_vs_a_q->Fill(a_q_precise,charge_val);
                        //cout << a_q_precise << "and   " << time_target_tof << endl;
                        h2_tof_vs_aq_fix_g_11c->Fill(a_q_precise,time_target_tof);
                        }

                        }
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



