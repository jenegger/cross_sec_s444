void get_cross_sec_macro(char const count_i[50]){

string input_setting = string(count_i);
if (input_setting.compare("plastic") == 0){
TFile* file = TFile::Open("/home/ge37liw/califa_combined_ch2_runs.root","READ");

TList *list = (TList*)file->Get("histlist");
TH1F* h1_inclusive = (TH1F*)list->FindObject("Number of p2p  events inclusive");

Long64_t num_min_bias = 13768654;
Long64_t num_reactions = h1_inclusive->Integral();

Double_t empty_reactions = 29.*num_min_bias/2058403.;
Double_t corr_factor_sim = pow(10,5)/9282.;
Double_t numb_reactions_corr = (num_reactions - empty_reactions)*corr_factor_sim;
Double_t T = 0.94*1.229*6.02214076*pow(10,23)/14.02658;
Double_t cross_section = (numb_reactions_corr/(num_min_bias*T))*pow(10,24);

cout << "**********CROSS SECTION MEASUREMENT***********" << endl; 
cout << "**********TARGET------------PLASTIC, 12.29 mm,density = 0.94 g/cm**3***********" << endl; 
cout << "Number of CH2 targets per cm**2:\t" << T << endl;
cout << "number of min bias events:\t" << num_min_bias << endl;
cout << "number of reactions, no correction:\t" << num_reactions << endl;
cout << "corrected reactions" << "\t" << numb_reactions_corr << endl;
cout << "12C(p,2p)X inclusive:\t" << endl;
cout << 1000*cross_section << "  mbarn" << endl;
cout << "*******************************************************************************" << endl;
}
else if (input_setting.compare("plastic_thick") == 0){
TFile* file = TFile::Open("/home/ge37liw/califa_combined_thick_ch2_runs.root","READ");

TList *list = (TList*)file->Get("histlist");
TH1F* h1_inclusive = (TH1F*)list->FindObject("Number of p2p  events inclusive");

Long64_t num_min_bias = 4414267;
Long64_t num_reactions = h1_inclusive->Integral();

Double_t empty_reactions = 29.*num_min_bias/2058403.;
Double_t corr_factor_sim = pow(10,5)/9282.;
Double_t numb_reactions_corr = (num_reactions - empty_reactions)*corr_factor_sim;
Double_t T = 0.94*2.453*6.02214076*pow(10,23)/14.02658;
Double_t cross_section = (numb_reactions_corr/(num_min_bias*T))*pow(10,24);

cout << "**********CROSS SECTION MEASUREMENT***********" << endl; 
cout << "**********TARGET------------PLASTIC, 24.53 mm,density = 0.94 g/cm**3***********" << endl; 
cout << "Number of CH2 targets per cm**2:\t" << T << endl;
cout << "number of min bias events:\t" << num_min_bias << endl;
cout << "number of reactions, no correction:\t" << num_reactions << endl;
cout << "corrected reactions" << "\t" << numb_reactions_corr << endl;
cout << "12C(p,2p)X inclusive:\t" << endl;
cout << 1000*cross_section << "  mbarn" << endl;
cout << "*******************************************************************************" << endl;
}
else if (input_setting.compare("carbon_thick") == 0){
	TFile* file = TFile::Open("/home/ge37liw/califa_opening_angle_carbon_run_0183_0001.root","READ");
	TList *list = (TList*)file->Get("histlist");
	TH1F* h1_inclusive = (TH1F*)list->FindObject("Number of p2p  events inclusive");
	
	Long64_t num_min_bias = 2604021;
	Long64_t num_reactions = h1_inclusive->Integral();
	
	Double_t empty_reactions = 29.*num_min_bias/2058403.;
	Double_t corr_factor_sim = pow(10,5)/9282.;
	Double_t numb_reactions_corr = (num_reactions - empty_reactions)*corr_factor_sim;
	Double_t T = 1.84*2.198*6.02214076*pow(10,23)/12.011;
	Double_t cross_section = (numb_reactions_corr/(num_min_bias*T))*pow(10,24);
	
	cout << "**********CROSS SECTION MEASUREMENT***********" << endl; 
	cout << "**********TARGET ------------CARBON, 21.98 mm,density = 1.84 g/cm**3***********" << endl; 
	cout << "Number of carbon targets per cm**2:\t" << T << endl;
	cout << "number of min bias events:\t" << num_min_bias << endl;
	cout << "number of reactions, no correction:\t" << num_reactions << endl;
	cout << "corrected reactions" << "\t" << numb_reactions_corr << endl;
	cout << "12C(p,2p)X inclusive:\t" << endl;
	cout << 1000*cross_section << "  mbarn" << endl;
	cout << "*******************************************************************************" << endl;
	}

else if (input_setting.compare("carbon_thin") == 0){
	TFile* file = TFile::Open("/home/ge37liw/califa_opening_angle_thin_carbon_run_0066_0001.root","READ");
	TList *list = (TList*)file->Get("histlist");
	TH1F* h1_inclusive = (TH1F*)list->FindObject("Number of p2p  events inclusive");
	
	Long64_t num_min_bias = 1897157;
	Long64_t num_reactions = h1_inclusive->Integral();
	
	Double_t empty_reactions = 29.*num_min_bias/2058403.;
	Double_t corr_factor_sim = pow(10,5)/9282.;
	Double_t numb_reactions_corr = (num_reactions - empty_reactions)*corr_factor_sim;
	Double_t T = 1.84*0.54*6.02214076*pow(10,23)/12.011;
	Double_t cross_section = (numb_reactions_corr/(num_min_bias*T))*pow(10,24);
	
	cout << "**********CROSS SECTION MEASUREMENT***********" << endl; 
	cout << "**********TARGET ------------CARBON, 5.4mm,density = 1.84 g/cm**3***********" << endl; 
	cout << "Number of carbon targets per cm**2:\t" << T << endl;
	cout << "number of min bias events:\t" << num_min_bias << endl;
	cout << "number of reactions, no correction:\t" << num_reactions << endl;
	cout << "corrected reactions" << "\t" << numb_reactions_corr << endl;
	cout << "12C(p,2p)X inclusive:\t" << endl;
	cout << 1000*cross_section << "  mbarn" << endl;
	cout << "****************************************************************************" << endl;
	}
else if (input_setting.compare("single_plastic_runs") == 0){
	vector<string> vec_runs{"0082_0001","0083_0001","0086_0001","0088_0001","0089_0001","0089_0002","0089_0003","0089_0004"};
	cout << "**********CROSS SECTION MEASUREMENT*****TARGET------------PLASTIC, 12.29 mm,density = 0.94 g/cm**3****" << endl;
	for (string i : vec_runs){
		char f_name[500];
		sprintf(f_name,"/home/ge37liw/califa_opening_angle_ch2_run_%s.root",i.c_str());
		TFile* eventfile = TFile::Open(f_name,"READ");
		TList *list = (TList*)eventfile->Get("histlist");
		TH1F* h1_inclusive = (TH1F*)list->FindObject("Number of p2p  events inclusive");
		Long64_t num_reactions = h1_inclusive->Integral();
		TH1F* h1_min_bias = (TH1F*)list->FindObject("Number of min. bias events");		
		Long64_t num_min_bias = h1_min_bias->Integral();		

		Double_t empty_reactions = 29.*num_min_bias/2058403.;
		Double_t corr_factor_sim = pow(10,5)/9282.;
		Double_t numb_reactions_corr = (num_reactions - empty_reactions)*corr_factor_sim;
		Double_t T = 0.94*1.229*6.02214076*pow(10,23)/14.02658;
		Double_t cross_section = (numb_reactions_corr/(num_min_bias*T))*pow(10,24);
		cout << "-----THIS IS RUN" << i << "-----------------------" << endl;	
		cout << "Number of CH2 targets per cm**2:\t" << T << endl;
		cout << "number of min bias events:\t" << num_min_bias << endl;
		cout << "number of reactions, no correction:\t" << num_reactions << endl;
		cout << "corrected reactions" << "\t" << numb_reactions_corr << endl;
		cout << "12C(p,2p)X inclusive:\t" << endl;
		cout << 1000*cross_section << "  mbarn" << endl;
		cout << "*******************************************************************************" << endl;
		}
	}
else if (input_setting.compare("single_plastic_min_bias_runs") == 0){
	vector<string> vec_runs{"0064_0001","0180_0001","0188_0001"};
	cout << "**********CROSS SECTION MEASUREMENT*****TARGET------------PLASTIC, 12.29 mm,density = 0.94 g/cm**3****" << endl;
	for (string i : vec_runs){
		char f_name[500];
		sprintf(f_name,"/home/ge37liw/califa_opening_angle_min_bias_ch2_run_%s.root",i.c_str());
		TFile* eventfile = TFile::Open(f_name,"READ");
		TList *list = (TList*)eventfile->Get("histlist");
		TH1F* h1_inclusive = (TH1F*)list->FindObject("Number of p2p  events inclusive");
		Long64_t num_reactions = h1_inclusive->Integral();
		TH1F* h1_min_bias = (TH1F*)list->FindObject("Number of min. bias events");		
		Long64_t num_min_bias = h1_min_bias->Integral();		

		Double_t empty_reactions = 29.*num_min_bias/2058403.;
		Double_t corr_factor_sim = pow(10,5)/9282.;
		Double_t numb_reactions_corr = (num_reactions - empty_reactions)*corr_factor_sim;
		Double_t T = 0.94*1.229*6.02214076*pow(10,23)/14.02658;
		Double_t cross_section = (numb_reactions_corr/(num_min_bias*T))*pow(10,24);
		cout << "-----THIS IS RUN" << i << "-----------------------" << endl;	
		cout << "Number of CH2 targets per cm**2:\t" << T << endl;
		cout << "number of min bias events:\t" << num_min_bias << endl;
		cout << "number of reactions, no correction:\t" << num_reactions << endl;
		cout << "corrected reactions" << "\t" << numb_reactions_corr << endl;
		cout << "12C(p,2p)X inclusive:\t" << endl;
		cout << 1000*cross_section << "  mbarn" << endl;
		cout << "*******************************************************************************" << endl;
		}
	}

else cout << "you have to choose as input \"carbon_thin\", \"carbon_thick\", \"plastic\", \"plastic_thick\", \"single_plastic_min_bias_runs\"  or \"single_plastic_runs\" " << endl; 
}
