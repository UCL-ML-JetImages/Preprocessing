#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include "TH2.h"
#include <cmath>
using namespace std;

// MACRO TO TRANSLATE, ROTATE, AND WEIGHT PLOTS
// ADJUST LINES 17, 106, 241, 242 WHEN SWAPPING BETWEEN SIGNAL AND BACKGROUND ROOT FILES 
// ADJUST LINES 181, 237 TO CHANGE WEIGHTS
// ADJUST LINES 234-236 FOR JET RESCALING



int background()
{

// Open files
	TFile *file = TFile::Open("QCD_background.root");	// switch file between signal and background here
	TTree *tree = (TTree*)file -> Get("FlatSubstructureJetTree");
	TCanvas *can1 = new TCanvas ("can1", "Eta Phi Plot", 600, 400);
	
// Define histogram
	TH2F *histo1 = new TH2F("histo1", "	",50,-1,1,50,-1,1);	// this seems to be a histogram struct or class 
                                                                //check documentation for input arguments
// Initialize
	float clus_N;
	float clus_E;		
	float clus_P;	
	float fjet_eta;
	float fjet_phi;
	float fjet_E;
	float fjet_pt;
	float max;
	float weight;
	float clus_dist;	
	float theta;
	float clus_Nfin;
	float clus_Pfin;
	float clus_1_rawE;
	float clus_2_LC;
	float clus_3_ISO;
	float clus_4_LAT;
	float clus_5_LON;
	float clus_6_SEC_LAM;
	float clus_7_SEC_R;
	float clus_8_CEN_LAM;
	float clus_9_CEN_MAG;
	float clus_10_ENG_POS;
	float clus_11_EM_PRO;
	float clus_12_ENG_FR;
	float clus_13_FIR;

    // Declaration & initialization of pointer to vector float and set to 0
	vector <float> *clus_eta 	    	= 0;
	vector <float> *clus_energy 		= 0;		
	vector <float> *clus_phi 	    	= 0;
	vector <float> *clus_rawE 		    = 0;
	vector <float> *clus_LCEoverEME 	= 0;
	vector <float> *clus_ISOLATION 		= 0;	
	vector <float> *clus_LATERAL		= 0;	
	vector <float> *clus_LONGITUDINAL 	= 0;	
	vector <float> *clus_SECOND_LAMBDA 	= 0;
	vector <float> *clus_SECOND_R 		= 0;			
	vector <float> *clus_CENTER_LAMBDA 	= 0;			
	vector <float> *clus_CENTER_MAG 	= 0;
	vector <float> *clus_ENG_POS	 	= 0;				
	vector <float> *clus_EM_PROBABILITY = 0;			
	vector <float> *clus_ENG_FRAC_MAX 	= 0;			
	vector <float> *clus_FIRST_ENG_DENS = 0;
					
	int fjet_fatjet_dRmatched_particle_flavor;

	const float PI = 3.1415927;
	const int W_boson = 24;

    //not quite sure what the X -> Y means    
    //this is probably how the root files are data structured
// Set variables
	tree -> SetBranchAddress("clus_eta", 		&clus_eta);
	tree -> SetBranchAddress("clus_phi", 		&clus_phi);
	tree -> SetBranchAddress("clus_E", 		&clus_energy);
	tree -> SetBranchAddress("clus_rawE",		&clus_rawE);
	tree -> SetBranchAddress("clus_LCEoverEME",	&clus_LCEoverEME);
	tree -> SetBranchAddress("clus_ISOLATION",	&clus_ISOLATION);
	tree -> SetBranchAddress("clus_LATERAL",	&clus_LATERAL);
	tree -> SetBranchAddress("clus_LONGITUDINAL", 	&clus_LONGITUDINAL);
	tree -> SetBranchAddress("clus_SECOND_LAMBDA",	&clus_SECOND_LAMBDA);
	tree -> SetBranchAddress("clus_SECOND_R",	&clus_SECOND_R);
	tree -> SetBranchAddress("clus_CENTER_LAMBDA", 	&clus_CENTER_LAMBDA);
	tree -> SetBranchAddress("clus_CENTER_MAG", 	&clus_CENTER_MAG);
	tree -> SetBranchAddress("clus_ENG_POS",	&clus_ENG_POS);
	tree -> SetBranchAddress("clus_EM_PROBABILITY", &clus_EM_PROBABILITY);
	tree -> SetBranchAddress("clus_ENG_FRAC_MAX",	&clus_ENG_FRAC_MAX);
	tree -> SetBranchAddress("clus_FIRST_ENG_DENS",	&clus_FIRST_ENG_DENS);
	tree -> SetBranchAddress("fjet_pt", 		&fjet_pt);
	tree -> SetBranchAddress("fjet_eta", 		&fjet_eta);
	tree -> SetBranchAddress("fjet_phi", 		&fjet_phi);	
	tree -> SetBranchAddress("fjet_E", 		&fjet_E);	
	tree -> SetBranchAddress("fjet_fatjet_dRmatched_particle_flavor", &fjet_fatjet_dRmatched_particle_flavor);		
	int tree_size = tree -> GetEntries(); // ?
	

// Fill histo1 with clusters	
	//tree_size=2; // ADJUST NUMBER OF JETS TO LOOK AT
	int clus_hits;

	for (int i = 0; i<tree_size; i++) {
		tree->GetEntry(i);
		
		// select only W bosons jets
		//while ((fjet_fatjet_dRmatched_particle_flavor == W_boson) || (fjet_fatjet_dRmatched_particle_flavor == -W_boson)) {	//PGID, comment out for background file	
				
			// determine which cluster has the largest energy
			for (int k = 0; k<clus_eta->size(); k++) {
				if (k == 0)
					max = clus_energy->at(k);
				if (clus_energy->at(k) > max)
					max = clus_energy->at(k);
				}		

			// run through each cluster for all clusters in a jet	
			for (int j = 0; j<clus_eta->size(); j++) {
		 
				// for only the cluster with the largest energy...
				while (clus_energy->at(j) == max) {
					
					// get initial cluster	
					clus_N = clus_eta->at(j);
					clus_P = clus_phi->at(j);
					clus_E = clus_energy->at(j);				
			
					// first shift to center at jet (0,0)
					clus_N = clus_N - fjet_eta;
					//cout << "Clus_N: " << clus_N << endl;		
					clus_P = clus_P - fjet_phi;
					//cout << "Clus_P: " << clus_P << endl;				
					
					// correct bounds for overshift	(azimuthal periodicity)
					if (clus_P > PI) 
						float clus_P = clus_P - 2*PI;	
			
					else if (clus_P < -PI) 
						float clus_P = clus_P + 2*PI ;				

					else 
						float clus_P = clus_P;			
					
					// find angle between cluster and the y-axis
					float clus_dist = sqrt( pow((clus_N-0),2) + pow((clus_P-0),2) );
					
					// correct for global theta
						// quadrant 1
						if (clus_N > 0 && clus_P > 0) {
							float theta = asin(clus_N/clus_dist);
							//cout << "theta #1: " << theta << endl; 
							} 
					
						// quadrant 2
						else if (clus_N < 0 && clus_P > 0) {
							clus_N = abs(clus_N);
							theta = asin(clus_N/clus_dist);
							theta = 2*PI - theta;
							//cout << "theta #2: " << theta << endl; 
							}
						
						// quadrant 3
						else if (clus_N < 0 && clus_P < 0) {
							clus_N = abs(clus_N);
							theta = asin(clus_N/clus_dist);
							theta = PI + theta;
							//cout << "theta #3: " << theta << endl; 
							} 
						
						// quadrant 4 
						else if (clus_N > 0 && clus_P < 0) {
							theta = asin(clus_N/clus_dist);
							theta = PI - theta;
							//cout << "theta #4: " << theta << endl; 
							} 	
						
						else
						//cout << "No Quadrant Selected" << endl; 			
					
					// rotate this cluster to the y-axis
					clus_N = 0;
					clus_P = clus_dist; 
					histo1 -> Fill(clus_N, clus_P, clus_E);
		
					break;
					}
					
			}		
			
			// run through each cluster for all clusters in a jet again	
			for (int q = 0; q<clus_eta->size(); q++) {					
					
				// rotate other clusters by theta 
				//while (clus_energy->at(q) != max) {
				
					// get initial cluster	
					
					clus_N 		= clus_eta->at(q);
					clus_P 		= clus_phi->at(q);
					float clus_E	 	= clus_energy->at(q);
					float clus_1_rawE 	= clus_rawE		->at(q);
					float clus_2_LC		= clus_LCEoverEME	->at(q); // var 1 (use in report)
					float clus_3_ISO	= clus_ISOLATION	->at(q); // var 2 (use in report)	
					float clus_4_LAT	= clus_LATERAL		->at(q);
					float clus_5_LON	= clus_LONGITUDINAL	->at(q); // var 3 (use in report)
					float clus_6_SEC_LAM	= clus_SECOND_LAMBDA	->at(q);
					float clus_7_SEC_R	= clus_SECOND_R		->at(q);
					float clus_8_CEN_LAM	= clus_CENTER_LAMBDA	->at(q);
					float clus_9_CEN_MAG	= clus_CENTER_MAG	->at(q);
					float clus_10_ENG_POS	= clus_ENG_POS		->at(q);
					float clus_11_EM_PRO	= clus_EM_PROBABILITY	->at(q); // var 4 (use in report)
					float clus_12_ENG_FR	= clus_ENG_FRAC_MAX	->at(q);
					float clus_13_FIR	= clus_FIRST_ENG_DENS	->at(q);				
			
					float clus_theta = 2*atan(exp(-clus_N));
					float clus_et = abs(clus_E*cos(clus_theta)); //	transverse energy, var 5 (use in report)	
					
					// first shift to centralize at jet (0,0)
					clus_N = clus_N - fjet_eta;				
					clus_P = clus_P - fjet_phi;				
					
					// correct bounds for overshift	
					if (clus_P > PI) 
						float clus_P = clus_P - 2*PI;	
			
					else if (clus_P < -PI) 
						float clus_P = clus_P + 2*PI ;				

					else 
						float clus_P = clus_P;					
		
					// rotate the clusters
					
					//float sf = 300000/fjet_pt; // scale
					float sf = 1.0;
					float clus_Nfin = (clus_N*cos(theta) - clus_P*sin(theta))*sf;
					float clus_Pfin = (clus_N*sin(theta) + clus_P*cos(theta))*sf;
		
					
					histo1 -> Fill(clus_Nfin, clus_Pfin, clus_11_EM_PRO);	// switch out third value to weight differently				
					clus_hits++;
					}
					
		//break;	//PGID, comment out for background file		
		//}       //PGID, comment out for background file	
												
	}		
	
// Plot histo1
	histo1 -> GetYaxis() -> SetTitle("Cluster Azimuthal Angle (Phi)");
	histo1 -> GetYaxis() -> SetTitleOffset(2);
	//histo1 -> GetYaxis() -> SetRangeUser(-0.6,0.7);
	histo1 -> GetXaxis() -> SetTitle("Cluster Pseudorapidity (Eta)");
	histo1 -> GetXaxis() -> SetTitleOffset(2);
//	histo1 -> GetXaxis() -> SetRangeUser(-0.7,0.7);
	histo1 -> GetZaxis() -> SetTitle("clus_ET");
	histo1 -> GetZaxis() -> SetTitleOffset(1.35);
   	gStyle->SetOptStat(0); // stat box display off
   	histo1 -> Draw("COLZ");
   	can1 -> SaveAs ("LSB_background_50by50_Clus_EMPROnoscale.png");
	return 0;
}			
