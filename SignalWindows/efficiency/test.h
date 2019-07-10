//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul  9 15:09:37 2019 by ROOT version 6.16/00
// from TTree L1PiXTRKTree/L1PiXTRKTree
// found on file: resultsv2.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

// Fixed size dimensions of array or collections stored in the TTree if any.
#include <iostream>
#include <fstream>
#include <string>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include "./roi_v2.h"
#include "./withEM_v2.h"
#include "./withoutEM_v2.h"

using namespace std;

class test {

	private:
		map<TString, TH1*> maphist;

	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Fixed size dimensions of array or collections stored in the TTree if any.

		// Declaration of leaf types
		vector<float>   *propgenElPartPhi;
		vector<float>   *propgenElPartEta;
		vector<float>   *propgenElPartPt;
		Int_t           bRecHitN;
		vector<int>     *bRecHitLayer;
		vector<float>   *bRecHitGx;
		vector<float>   *bRecHitGy;
		vector<float>   *bRecHitGz;
		Int_t           fRecHitN;
		vector<int>     *fRecHitDisk;
		vector<float>   *fRecHitGx;
		vector<float>   *fRecHitGy;
		vector<float>   *fRecHitGz;
		Int_t           EgN;
		vector<float>   *egCrysClusterEt;
		vector<float>   *egCrysClusterEta;
		vector<float>   *egCrysClusterPhi;
		vector<float>   *egCrysClusterGx;
		vector<float>   *egCrysClusterGy;
		vector<float>   *egCrysClusterGz;

		// List of branches
		TBranch        *b_propgenElPartPhi;   //!
		TBranch        *b_propgenElPartEta;   //!
		TBranch        *b_propgenElPartPt;   //!
		TBranch        *b_bRecHitN;   //!
		TBranch        *b_bRecHitLayer;   //!
		TBranch        *b_bRecHitGx;   //!
		TBranch        *b_bRecHitGy;   //!
		TBranch        *b_bRecHitGz;   //!
		TBranch        *b_fRecHitN;   //!
		TBranch        *b_fRecHitDisk;   //!
		TBranch        *b_fRecHitGx;   //!
		TBranch        *b_fRecHitGy;   //!
		TBranch        *b_fRecHitGz;   //!
		TBranch        *b_EgN;   //!
		TBranch        *b_egCrysClusterEt;   //!
		TBranch        *b_egCrysClusterEta;   //!
		TBranch        *b_egCrysClusterPhi;   //!
		TBranch        *b_egCrysClusterGx;   //!
		TBranch        *b_egCrysClusterGy;   //!
		TBranch        *b_egCrysClusterGz;   //!

		test(TTree *tree=0);
		virtual ~test();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual void     Loop();
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);

		int Ele, Pos;
		int skip;
		float EGN;
		int eta_region;

		int withoutEM_count_Ele, withEM_count_Ele;
		bool PixTrkPassed; bool TrkIsoPassed;
		int pass_count_wo4thPix, pass_count_wo3thPix, pass_count_wo2thPix, pass_count_wo1thPix;
		int woEM_pass_Ele_count_wo4thPix, woEM_pass_Ele_count_wo3thPix, woEM_pass_Ele_count_wo2thPix, woEM_pass_Ele_count_wo1thPix;
		int wEM_pass_Ele_count_wo4thPix, wEM_pass_Ele_count_wo3thPix, wEM_pass_Ele_count_wo2thPix, wEM_pass_Ele_count_wo1thPix;

		double all_cut_pass_eg;
		int all_cut_pass_Ele, withoutEM_pass_Ele, withEM_pass_Ele;
		int all_cut_pass_Pos;
		int fourth_layer_missing;
		int third_layer_missing;
		int second_layer_missing;
		int first_layer_missing;

		int bit1;
		int bit2;
		int trigger_bit_width_;
		int trigger_bit_width_iso_;
		int pix_comb_;

		bool debug;

		double  L1_Dphi_cut1, L1_Dphi_cut2;
		double  L2_Dphi_cut1, L2_Dphi_cut2;
		double  L3_Dphi_cut1, L3_Dphi_cut2;
		double  L4_Dphi_cut1, L4_Dphi_cut2;
		double  D1_Dphi_cut1, D1_Dphi_cut2;
		double  D2_Dphi_cut1, D2_Dphi_cut2;
		double  D3_Dphi_cut1, D3_Dphi_cut2;

		double dPhi012;
		double dPhi013;
		double dPhi014;
		double dPhi023;
		double dPhi024;
		double dPhi034;

		double  L012_DPhi_cut1, L012_DPhi_cut2;

		double  L013_DPhi_cut1, L013_DPhi_cut2;

		double  L014_DPhi_cut1, L014_DPhi_cut2;

		double  L023_DPhi_cut1, L023_DPhi_cut2;

		double  L024_DPhi_cut1, L024_DPhi_cut2;

		double  L034_DPhi_cut1, L034_DPhi_cut2;

		double  L123_DPhi_cut1, L123_DPhi_cut2;
		double  L123_DEta_cut1, L123_DEta_cut2;

		double  L124_DPhi_cut1, L124_DPhi_cut2;
		double  L124_DEta_cut1, L124_DEta_cut2;

		double  L134_DPhi_cut1, L134_DPhi_cut2;
		double  L134_DEta_cut1, L134_DEta_cut2;

		double  L234_DPhi_cut1, L234_DPhi_cut2;
		double  L234_DEta_cut1, L234_DEta_cut2;

		double  L12_eta_upper, L13_eta_upper, L14_eta_upper, L23_eta_upper, L24_eta_upper, L34_eta_upper;
		double  L12_phi_upper, L13_phi_upper, L14_phi_upper, L23_phi_upper, L24_phi_upper, L34_phi_upper;
		double  L12_eta_bellow, L13_eta_bellow, L14_eta_bellow, L23_eta_bellow, L24_eta_bellow, L34_eta_bellow;
		double  L12_phi_bellow, L13_phi_bellow, L14_phi_bellow, L23_phi_bellow, L24_phi_bellow, L34_phi_bellow;
		double  L12_R_bellow, L13_R_bellow, L14_R_bellow, L23_R_bellow, L24_R_bellow, L34_R_bellow;

		double dPhi;
		double dEta;
		double dPhi_1, dPhi_2, dPhi_3;
		double dEta_1, dEta_2, dEta_3;
		TVector3 first_temp, second_temp;
		int _pass_Ele, _pass_Pos;

		int L012_pass_Ele, L012_pass_Pos; 
		int L013_pass_Ele, L013_pass_Pos;
		int L014_pass_Ele, L014_pass_Pos;
		int L023_pass_Ele, L023_pass_Pos;
		int L024_pass_Ele, L024_pass_Pos;
		int L034_pass_Ele, L034_pass_Pos;
		int L123_pass_Ele, L123_pass_Pos;
		int L124_pass_Ele, L124_pass_Pos;
		int L134_pass_Ele, L134_pass_Pos;
		int L234_pass_Ele, L234_pass_Pos;

		int L12_EM_Ele, L12_EM_Pos; 
		int L13_EM_Ele, L13_EM_Pos;
		int L14_EM_Ele, L14_EM_Pos;
		int L23_EM_Ele, L23_EM_Pos;
		int L24_EM_Ele, L24_EM_Pos;
		int L34_EM_Ele, L34_EM_Pos;

		std::vector<TVector3> first_layer_hits;
		std::vector<TVector3> second_layer_hits;
		std::vector<TVector3> third_layer_hits;
		std::vector<TVector3> fourth_layer_hits;

		double r; // r for radius of pixel tracker layer
		int layers[5];  // initialize as 0, layers contain # of hits on each pixel layer

		TVector3 emvector;
		float EgEt;
		float EgEta;
		float EgPhi;

		std::vector<int> first_layer_hits_Ele_or_Pos;
		std::vector<int> second_layer_hits_Ele_or_Pos;
		std::vector<int> third_layer_hits_Ele_or_Pos;
		std::vector<int> fourth_layer_hits_Ele_or_Pos;
		std::vector<int> hitted_layers;

		void MakeHistograms(TString hname, int nbins, float xmin, float xmax);
		TH1* GetHist(TString hname);
		void FillHist(TString histname, float value, float w, float xmin, float xmax, int nbins);

		void StorePixelHit( int region);
		double StandaloneDPhi( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit );
		double StandaloneDEta( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit );
		double EMmatchingDEta(TVector3& first_layer, TVector3& second_layer, TVector3& egvector);
		double EMmatchingDPhi(TVector3& first_layer, TVector3& second_layer, TVector3& egvector);
		int Signal_window_check( double upper, double value, double lower, int Ele_Pos);
		void FillCutFlow(TString cut, float weight);
		void SetROI(int region);
		void SetSingalBoundary(int region, double eg_dphi, double eg_deta, double sa_dphi, double sa_deta);
		void TriggeringWith_1st2ndPixel(int nthFirstHit, int nthSecondHit);
		void TriggeringWith_1st3rdPixel(int nthFirstHit, int nthSecondHit);
		void TriggeringWith_2nd3rdPixel(int nthFirstHit, int nthSecondHit);
		void TriggeringWithout_4thPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);
		void TriggeringWithout_3rdPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);
		void TriggeringWithout_2ndPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);
		void TriggeringWithout_1stPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);

		void TriggeringWith_1st2ndPixel_v2(int nthFirstHit, int nthSecondHit);
		void TriggeringWith_1st3rdPixel_v2(int nthFirstHit, int nthSecondHit);
		void TriggeringWith_1st4thPixel_v2(int nthFirstHit, int nthSecondHit);
		void TriggeringWith_2nd3rdPixel_v2(int nthFirstHit, int nthSecondHit);
		void TriggeringWith_2nd4thPixel_v2(int nthFirstHit, int nthSecondHit);
		void TriggeringWith_3rd4thPixel_v2(int nthFirstHit, int nthSecondHit);

		inline float deltaPhi(float phi1, float phi2) { 
			float result = phi1 - phi2;
			while (result > float(M_PI)) result -= float(2*M_PI);
			while (result <= -float(M_PI)) result += float(2*M_PI);
			return result;
		}

		TFile *outfile;
		TTree* pixtrk_tree;

		int count_Entry; 
		int pass_egobjects_check;
		int ntnEg2; 
		int event_denominator; 
		int event_nominator; 

		int nPix123_segments;
		int nPix124_segments;
		int nPix134_segments;
		int nPix234_segments;

		vector<float> ntEgEt; 
		vector<float> ntEgEta; 
		vector<float> ntEgPhi; 

		vector<float> ntL1TkEgEt; 
		vector<float> ntL1TkEgEta; 
		vector<float> ntL1TkEgPhi; 

		vector<float> ntL1TkLooseEgEt; 
		vector<float> ntL1TkLooseEgEta; 
		vector<float> ntL1TkLooseEgPhi; 

		float matchedEgEt; 
		float matchedEgEta; 
		float matchedEgPhi; 
		int   fired;

		vector<int> PiXTRKbit;
		vector<int> pix_comb;
		vector<int> trigger_bit_width;
		vector<int> trigger_bit_width_iso;

		vector<bool> ntCl_match; 
		vector<bool> isTrack_match; 
		vector<float> chi2; 
		vector<float> track_dr;
		vector<bool> withoutEM_match; 
		vector<bool> withEM_match; 

		vector<int> ntfirstPix;  
		vector<int> ntsecondPix; 
		vector<int> ntthirdPix;  
		vector<int> ntfourthPix; 

		float nt_lastSimtkpt;
		float nt_initialSimtkpt;

		float nt_genPhi;
		float nt_genEta;
		float nt_genPt;

};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	if (tree == 0) {
		TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("results.root");
		if (!f || !f->IsOpen()) {
			f = new TFile("results.root");
		}
		TDirectory * dir = (TDirectory*)f->Get("results.root:/l1PiXTRKTree");
		dir->GetObject("L1PiXTRKTree",tree);

	}
	Init(tree);

	Ele = 1, Pos = 2;
	skip = 0;

	outfile = new TFile("eff.root","recreate");
	pixtrk_tree = new TTree("t","t");

	count_Entry = 1;
	pixtrk_tree->Branch("totalEvent", &count_Entry, "count_Entry/I");   
	pixtrk_tree->Branch("totalEgN", &EgN, "EgN/F");
	pixtrk_tree->Branch("ntnEg2", &ntnEg2, "ntnEg2/I");

	pixtrk_tree->Branch("ntEgEt",&ntEgEt);
	pixtrk_tree->Branch("ntEgEta",&ntEgEta);
	pixtrk_tree->Branch("ntEgPhi",&ntEgPhi);

	pixtrk_tree->Branch("ntL1TkEgEt",&ntL1TkEgEt);
	pixtrk_tree->Branch("ntL1TkEgEta",&ntL1TkEgEta);
	pixtrk_tree->Branch("ntL1TkEgPhi",&ntL1TkEgPhi);

	pixtrk_tree->Branch("ntL1TkLooseEgEt",&ntL1TkLooseEgEt);
	pixtrk_tree->Branch("ntL1TkLooseEgEta",&ntL1TkLooseEgEta);
	pixtrk_tree->Branch("ntL1TkLooseEgPhi",&ntL1TkLooseEgPhi);

	pixtrk_tree->Branch("PiXTRKbit",&PiXTRKbit);
	pixtrk_tree->Branch("trigger_bit_width",&trigger_bit_width);
	pixtrk_tree->Branch("trigger_bit_width_iso",&trigger_bit_width_iso);
	pixtrk_tree->Branch("pix_comb",&pix_comb);

	pixtrk_tree->Branch("ntCl_match",&ntCl_match);
	pixtrk_tree->Branch("withoutEM_match",&withoutEM_match);
	pixtrk_tree->Branch("withEM_match",&withEM_match);

	pixtrk_tree->Branch("nt_genPhi",&nt_genPhi,"nt_genPhi/F");
	pixtrk_tree->Branch("nt_genEta",&nt_genEta,"nt_genEta/F");
	pixtrk_tree->Branch("nt_genPt",&nt_genPt,"nt_genPt/F");

	pixtrk_tree->Branch("nt_lastSimtkpt",&nt_lastSimtkpt,"nt_lastSimtkpt/F");
	pixtrk_tree->Branch("nt_initialSimtkpt",&nt_initialSimtkpt,"nt_initialSimtkpt/F");

	pixtrk_tree->Branch("fired",&fired,"fired/I");

}

test::~test()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void test::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set object pointer
	propgenElPartPhi = 0;
	propgenElPartEta = 0;
	propgenElPartPt = 0;
	bRecHitLayer = 0;
	bRecHitGx = 0;
	bRecHitGy = 0;
	bRecHitGz = 0;
	fRecHitDisk = 0;
	fRecHitGx = 0;
	fRecHitGy = 0;
	fRecHitGz = 0;
	egCrysClusterEt = 0;
	egCrysClusterEta = 0;
	egCrysClusterPhi = 0;
	egCrysClusterGx = 0;
	egCrysClusterGy = 0;
	egCrysClusterGz = 0;
	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("propgenElPartPhi", &propgenElPartPhi, &b_propgenElPartPhi);
	fChain->SetBranchAddress("propgenElPartEta", &propgenElPartEta, &b_propgenElPartEta);
	fChain->SetBranchAddress("propgenElPartPt", &propgenElPartPt, &b_propgenElPartPt);
	fChain->SetBranchAddress("bRecHitN", &bRecHitN, &b_bRecHitN);
	fChain->SetBranchAddress("bRecHitLayer", &bRecHitLayer, &b_bRecHitLayer);
	fChain->SetBranchAddress("bRecHitGx", &bRecHitGx, &b_bRecHitGx);
	fChain->SetBranchAddress("bRecHitGy", &bRecHitGy, &b_bRecHitGy);
	fChain->SetBranchAddress("bRecHitGz", &bRecHitGz, &b_bRecHitGz);
	fChain->SetBranchAddress("fRecHitN", &fRecHitN, &b_fRecHitN);
	fChain->SetBranchAddress("fRecHitDisk", &fRecHitDisk, &b_fRecHitDisk);
	fChain->SetBranchAddress("fRecHitGx", &fRecHitGx, &b_fRecHitGx);
	fChain->SetBranchAddress("fRecHitGy", &fRecHitGy, &b_fRecHitGy);
	fChain->SetBranchAddress("fRecHitGz", &fRecHitGz, &b_fRecHitGz);
	fChain->SetBranchAddress("EgN", &EgN, &b_EgN);
	fChain->SetBranchAddress("egCrysClusterEt", &egCrysClusterEt, &b_egCrysClusterEt);
	fChain->SetBranchAddress("egCrysClusterEta", &egCrysClusterEta, &b_egCrysClusterEta);
	fChain->SetBranchAddress("egCrysClusterPhi", &egCrysClusterPhi, &b_egCrysClusterPhi);
	fChain->SetBranchAddress("egCrysClusterGx", &egCrysClusterGx, &b_egCrysClusterGx);
	fChain->SetBranchAddress("egCrysClusterGy", &egCrysClusterGy, &b_egCrysClusterGy);
	fChain->SetBranchAddress("egCrysClusterGz", &egCrysClusterGz, &b_egCrysClusterGz);
	Notify();
}

Bool_t test::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void test::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t test::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}

void test::MakeHistograms(TString hname, int nbins, float xmin, float xmax){

	maphist[hname] =  new TH1F(hname.Data(),hname.Data(),nbins,xmin,xmax);

}

TH1* test::GetHist(TString hname){

	TH1* h = NULL;
	std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
	if(mapit != maphist.end()) return mapit->second;

	return h;

}

void test::FillHist(TString histname, float value, float w, float xmin, float xmax, int nbins){

	if(GetHist(histname)) GetHist(histname)->Fill(value, w);
	else{
		//     cout << "Making histogram..." << endl;
		MakeHistograms(histname, nbins, xmin, xmax);
		if(GetHist(histname)) GetHist(histname)->Fill(value, w);
	}
}

double test::StandaloneDPhi( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit ){
	if( first_hit == 0 ){

		TVector3 temp;

		if( second_hit == 1 ) temp.SetXYZ( first_layer_hits[which_second_hit].X(), first_layer_hits[which_second_hit].Y(), first_layer_hits[which_second_hit].Z() );
		if( second_hit == 2 ) temp.SetXYZ( second_layer_hits[which_second_hit].X(), second_layer_hits[which_second_hit].Y(), second_layer_hits[which_second_hit].Z() );
		if( second_hit == 3 ) temp.SetXYZ( third_layer_hits[which_second_hit].X(), third_layer_hits[which_second_hit].Y(), third_layer_hits[which_second_hit].Z() );


		if( third_hit == 2 ) return deltaPhi( (second_layer_hits[which_third_hit] - temp).Phi(), temp.Phi());
		if( third_hit == 3 ) return deltaPhi( (third_layer_hits[which_third_hit] - temp).Phi(),  temp.Phi());
		if( third_hit == 4 ) return deltaPhi( (fourth_layer_hits[which_third_hit] - temp).Phi(), temp.Phi());

	}
	if( first_hit != 0 ){
		TVector3 temp_first_layer;
		TVector3 temp_second_layer;
		TVector3 temp_third_layer;

		if( first_hit == 1 ) temp_first_layer = first_layer_hits[which_first_hit];
		if( first_hit == 2 ) temp_first_layer = second_layer_hits[which_first_hit];

		if( second_hit == 2 ) temp_second_layer = second_layer_hits[which_second_hit];
		if( second_hit == 3 ) temp_second_layer = third_layer_hits[which_second_hit];


		if( third_hit == 3 ) temp_third_layer = third_layer_hits[which_third_hit];
		if( third_hit == 4 ) temp_third_layer = fourth_layer_hits[which_third_hit];

		return deltaPhi( (temp_third_layer - temp_second_layer).Phi(), (temp_second_layer - temp_first_layer).Phi());
	}
	return 0.;
}

double test::StandaloneDEta( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit ){

	TVector3 temp_first_layer;
	TVector3 temp_second_layer;
	TVector3 temp_third_layer;

	if( first_hit == 1 ) temp_first_layer = first_layer_hits[which_first_hit];
	if( first_hit == 2 ) temp_first_layer = second_layer_hits[which_first_hit];

	if( second_hit == 2 ) temp_second_layer = second_layer_hits[which_second_hit];
	if( second_hit == 3 ) temp_second_layer = third_layer_hits[which_second_hit];


	if( third_hit == 3 ) temp_third_layer = third_layer_hits[which_third_hit];
	if( third_hit == 4 ) temp_third_layer = fourth_layer_hits[which_third_hit];

	return (temp_third_layer - temp_second_layer).PseudoRapidity() - (temp_second_layer - temp_first_layer).PseudoRapidity();
}

double test::EMmatchingDEta(TVector3& first_layer, TVector3& second_layer, TVector3& egvector){

	TVector3 pixelVector = second_layer - first_layer;
	TVector3 EM_pixelVector = egvector - second_layer;
	return EM_pixelVector.Eta() - pixelVector.Eta();

}

double test::EMmatchingDPhi(TVector3& first_layer, TVector3& second_layer, TVector3& egvector){


	TVector3 pixelVector = second_layer - first_layer;
	TVector3 EM_pixelVector = egvector - second_layer;

	return deltaPhi( EM_pixelVector.Phi(), pixelVector.Phi() );
}

int test::Signal_window_check( double upper, double value, double lower, int Ele_Pos){

	if( Ele_Pos == 1 ){ // 1 is Electron
		if( value <= upper && value >= lower){
			return true;
		}
		else
			return false;
	}
	if( Ele_Pos == 2 ){ // 2 is Positron
		if( value >= -upper && value <= -lower){
			return true;
		}
		else
			return false;
	}
	return 0;
}

void test::FillCutFlow(TString cut, float weight){


	if(GetHist("cutflow")) {
		GetHist("cutflow")->Fill(cut,weight);

	}
	else{
		test::MakeHistograms("cutflow", 6,0.,6.);

		GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
		GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"MinEtCut");
		GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"EtaCut");
		GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"DRCut");
		GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"PtErrCut");
		GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"EvtCut");


	}
}
void test::StorePixelHit(int region){

	for(int a=0; a<bRecHitN; a++){
		int Dphi_Ele_pass = 0;
		int Dphi_Pos_pass = 0;
		double Dphi = 0.;
		int el_or_po = 0; // electron = 1, positron = 2, both = 3
		TVector3 current_hit;
		current_hit.SetXYZ( bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a) );
		Dphi = deltaPhi( current_hit.Phi(), EgPhi);


		if( region < 5 ){
			if( bRecHitLayer->at(a) == 1 ){ // First layer
				if( Dphi < L1_Dphi_cut1 && Dphi > L1_Dphi_cut2 ){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L1_Dphi_cut1 && Dphi < -L1_Dphi_cut2 ){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[1]++;
					first_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
					first_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
		}
		if( region == 1 || region == 2 || region == 3 ){
			if( bRecHitLayer->at(a) == 2 ){ // Second layer
				if( Dphi < L2_Dphi_cut1 && Dphi > L2_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L2_Dphi_cut1 && Dphi < -L2_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[2]++;
					second_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
					second_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
		}
		if( region == 1 || region == 2 ){
			if( bRecHitLayer->at(a) == 3 ){ // Third layer 
				if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[3]++;
					third_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
					third_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
		}
		if( region == 1 ){
			if( bRecHitLayer->at(a) == 4 ){ // Fourth layer
				if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[4]++;
					fourth_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
					fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
		}
	}
	for(int a=0; a<fRecHitN; a++){
		int Dphi_Ele_pass = 0;
		int Dphi_Pos_pass = 0;
		double Dphi = 0.;
		int el_or_po = 0; // electron = 1, positron = 2, both = 3
		TVector3 current_hit;
		current_hit.SetXYZ( fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) );
		Dphi = deltaPhi(current_hit.Phi(), EgPhi);

		if( region == 5 ){

			if( fRecHitDisk->at(a) == 1 ){ // first disk
				if( Dphi < L1_Dphi_cut1 && Dphi > L1_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L1_Dphi_cut1 && Dphi < -L1_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[1]++;
					first_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					first_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}

			if( fRecHitDisk->at(a) == 2 ){ // second disk
				if( Dphi < L2_Dphi_cut1 && Dphi > L2_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L2_Dphi_cut1 && Dphi < -L2_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[2]++;
					second_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					second_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
			if( fRecHitDisk->at(a) == 3 ){ // third disk
				if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[3]++;
					third_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					third_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
			if( fRecHitDisk->at(a) == 4 ){ // fourth disk
				if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[4]++;
					fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
		}

		if( region == 6 ){
			if( fRecHitDisk->at(a) == 2 ){ // second disk
				if( Dphi < L1_Dphi_cut1 && Dphi > L1_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L1_Dphi_cut1 && Dphi < -L1_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[1]++;
					first_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					first_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}

			if( fRecHitDisk->at(a) == 3 ){ // third disk
				if( Dphi < L2_Dphi_cut1 && Dphi > L2_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L2_Dphi_cut1 && Dphi < -L2_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[2]++;
					second_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					second_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
			if( fRecHitDisk->at(a) == 4 ){ // fourth disk
				if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[3]++;
					third_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					third_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}

			if( fRecHitDisk->at(a) == 5 ){ // fifth disk
				if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[4]++;
					fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}

		}

		if( region == 4 ){
			if( fRecHitDisk->at(a) == 1 ){ // first disk
				if( Dphi < L1_Dphi_cut1 && Dphi > L1_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L1_Dphi_cut1 && Dphi < -L1_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[2]++;
					second_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					second_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
			if( fRecHitDisk->at(a) == 2 ){ // second disk
				if( Dphi < L2_Dphi_cut1 && Dphi > L2_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L2_Dphi_cut1 && Dphi < -L2_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[3]++;
					third_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					third_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
			if( fRecHitDisk->at(a) == 3 ){ // third disk
				if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[4]++;
					fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
		}

		if( region == 3 ){
			if( fRecHitDisk->at(a) == 1 ){ // first disk
				if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[3]++;
					third_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					third_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
			if( fRecHitDisk->at(a) == 2 ){ // second disk
				if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[4]++;
					fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
		}
		if( region == 2 ){
			if( fRecHitDisk->at(a) == 1 ){ // first disk
				if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
					Dphi_Ele_pass = 1; el_or_po = 1;
				}
				if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
					Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
				}
				if( Dphi_Ele_pass || Dphi_Pos_pass ){
					layers[4]++;
					fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
					fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
				}
			}
		}

	}
}

void test::SetROI(int region){


	float upper_width = 0.055;
	float lower_width = 0.055;

	L1_Dphi_cut1 = ROI_func(region, 0, EgEt);
	L1_Dphi_cut2 = ROI_func(region, 0, EgEt);

	L1_Dphi_cut1 = L1_Dphi_cut1 + upper_width;
	L1_Dphi_cut2 = L1_Dphi_cut2 - lower_width;

	L2_Dphi_cut1 = ROI_func(region, 1, EgEt);
	L2_Dphi_cut2 = ROI_func(region, 1, EgEt);

	L2_Dphi_cut1 = L2_Dphi_cut1 + upper_width;
	L2_Dphi_cut2 = L2_Dphi_cut2 - lower_width;

	L3_Dphi_cut1 = ROI_func(region, 2, EgEt);
	L3_Dphi_cut2 = ROI_func(region, 2, EgEt);

	L3_Dphi_cut1 = L3_Dphi_cut1 + upper_width;
	L3_Dphi_cut2 = L3_Dphi_cut2 - lower_width;

	L4_Dphi_cut1 = ROI_func(region, 3, EgEt);
	L4_Dphi_cut2 = ROI_func(region, 3, EgEt);

	L4_Dphi_cut1 = L4_Dphi_cut1 + upper_width;
	L4_Dphi_cut2 = L4_Dphi_cut2 - lower_width;

}

void test::SetSingalBoundary(int region, double eg_dphi, double eg_deta, double sa_dphi, double sa_deta){

	float EG_pixel_dphi_upper_width = eg_dphi;
	float EG_pixel_dphi_lower_width = eg_dphi;

	// pixel-EG 
	L12_phi_upper =  SW_func2_dphi_v2(region, 0, EgEt);
	L12_phi_bellow = SW_func2_dphi_v2(region, 0, EgEt);

	L12_phi_upper = L12_phi_upper   + EG_pixel_dphi_upper_width;
	L12_phi_bellow = L12_phi_bellow - EG_pixel_dphi_lower_width;

	L13_phi_upper =  SW_func2_dphi_v2(region, 1, EgEt);
	L13_phi_bellow = SW_func2_dphi_v2(region, 1, EgEt);

	L13_phi_upper =  L13_phi_upper  + EG_pixel_dphi_upper_width;
	L13_phi_bellow = L13_phi_bellow - EG_pixel_dphi_lower_width;

	L14_phi_upper =  SW_func2_dphi_v2(region, 2, EgEt);
	L14_phi_bellow = SW_func2_dphi_v2(region, 2, EgEt);

	L14_phi_upper =  L14_phi_upper  + EG_pixel_dphi_upper_width;
	L14_phi_bellow = L14_phi_bellow - EG_pixel_dphi_lower_width;

	L23_phi_upper =  SW_func2_dphi_v2(region, 3, EgEt);
	L23_phi_bellow = SW_func2_dphi_v2(region, 3, EgEt);

	L23_phi_upper =  L23_phi_upper  + EG_pixel_dphi_upper_width;
	L23_phi_bellow = L23_phi_bellow - EG_pixel_dphi_lower_width;

	L24_phi_upper =  SW_func2_dphi_v2(region, 4, EgEt);
	L24_phi_bellow = SW_func2_dphi_v2(region, 4, EgEt);

	L24_phi_upper =  L24_phi_upper  + EG_pixel_dphi_upper_width;
	L24_phi_bellow = L24_phi_bellow - EG_pixel_dphi_lower_width;

	L34_phi_upper =  SW_func2_dphi_v2(region, 5, EgEt);
	L34_phi_bellow = SW_func2_dphi_v2(region, 5, EgEt);

	L34_phi_upper = L34_phi_upper   + EG_pixel_dphi_upper_width;
	L34_phi_bellow = L34_phi_bellow - EG_pixel_dphi_lower_width;

	float EG_pixel_deta_upper_width = eg_deta;
	float EG_pixel_deta_lower_width = eg_deta;

	L12_eta_upper =  SW_func2_deta_v2(EgEt);
	L12_eta_bellow = SW_func2_deta_v2(EgEt);

	L12_eta_upper =  L12_eta_upper  + EG_pixel_deta_upper_width;
	L12_eta_bellow = L12_eta_bellow - EG_pixel_deta_lower_width;

	L13_eta_upper =  SW_func2_deta_v2(EgEt);
	L13_eta_bellow = SW_func2_deta_v2(EgEt);

	L13_eta_upper =  L13_eta_upper  + EG_pixel_deta_upper_width;
	L13_eta_bellow = L13_eta_bellow - EG_pixel_deta_lower_width;

	L14_eta_upper =  SW_func2_deta_v2(EgEt);
	L14_eta_bellow = SW_func2_deta_v2(EgEt);

	L14_eta_upper =  L14_eta_upper  + EG_pixel_deta_upper_width;
	L14_eta_bellow = L14_eta_bellow - EG_pixel_deta_lower_width;

	L23_eta_upper =  SW_func2_deta_v2(EgEt);
	L23_eta_bellow = SW_func2_deta_v2(EgEt);

	L23_eta_upper =  L23_eta_upper  + EG_pixel_deta_upper_width;
	L23_eta_bellow = L23_eta_bellow - EG_pixel_deta_lower_width;

	L24_eta_upper =  SW_func2_deta_v2(EgEt);
	L24_eta_bellow = SW_func2_deta_v2(EgEt);

	L24_eta_upper =  L24_eta_upper  + EG_pixel_deta_upper_width;
	L24_eta_bellow = L24_eta_bellow - EG_pixel_deta_lower_width;

	L34_eta_upper =  SW_func2_deta_v2(EgEt);
	L34_eta_bellow = SW_func2_deta_v2(EgEt);

	L34_eta_upper =  L34_eta_upper  + EG_pixel_deta_upper_width;
	L34_eta_bellow = L34_eta_bellow - EG_pixel_deta_lower_width;

	// pixel-pixel 
	float pixel_pixel_dphi_upper_width = sa_dphi;
	float pixel_pixel_dphi_lower_width = sa_dphi;

	L012_DPhi_cut1 = SW_func1_dphi_v2(region, 0, EgEt);
	L012_DPhi_cut2 = SW_func1_dphi_v2(region, 0, EgEt);

	L012_DPhi_cut1 = L012_DPhi_cut1 + pixel_pixel_dphi_upper_width;
	L012_DPhi_cut2 = L012_DPhi_cut2 - pixel_pixel_dphi_lower_width;

	L013_DPhi_cut1 = SW_func1_dphi_v2(region, 1, EgEt);
	L013_DPhi_cut2 = SW_func1_dphi_v2(region, 1, EgEt);

	L013_DPhi_cut1 = L013_DPhi_cut1 + pixel_pixel_dphi_upper_width;
	L013_DPhi_cut2 = L013_DPhi_cut2 - pixel_pixel_dphi_lower_width;

	L014_DPhi_cut1 = SW_func1_dphi_v2(region, 2, EgEt);
	L014_DPhi_cut2 = SW_func1_dphi_v2(region, 2, EgEt);

	L014_DPhi_cut1 = L014_DPhi_cut1 + pixel_pixel_dphi_upper_width;
	L014_DPhi_cut2 = L014_DPhi_cut2 - pixel_pixel_dphi_lower_width;

	L023_DPhi_cut1 = SW_func1_dphi_v2(region, 3, EgEt);
	L023_DPhi_cut2 = SW_func1_dphi_v2(region, 3, EgEt);

	L023_DPhi_cut1 = L023_DPhi_cut1 + pixel_pixel_dphi_upper_width;
	L023_DPhi_cut2 = L023_DPhi_cut2 - pixel_pixel_dphi_lower_width;

	L024_DPhi_cut1 = SW_func1_dphi_v2(region, 4, EgEt);
	L024_DPhi_cut2 = SW_func1_dphi_v2(region, 4, EgEt);

	L024_DPhi_cut1 = L024_DPhi_cut1 + pixel_pixel_dphi_upper_width;
	L024_DPhi_cut2 = L024_DPhi_cut2 - pixel_pixel_dphi_lower_width;

	L034_DPhi_cut1 = SW_func1_dphi_v2(region, 5, EgEt);
	L034_DPhi_cut2 = SW_func1_dphi_v2(region, 5, EgEt);

	L034_DPhi_cut1 = L034_DPhi_cut1 + pixel_pixel_dphi_upper_width;
	L034_DPhi_cut2 = L034_DPhi_cut2 - pixel_pixel_dphi_lower_width;

	L123_DPhi_cut1 = SW_func1_dphi_v2(region, 6, EgEt);
	L123_DPhi_cut2 = SW_func1_dphi_v2(region, 6, EgEt);

	L123_DPhi_cut1 = L123_DPhi_cut1 + pixel_pixel_dphi_upper_width;
	L123_DPhi_cut2 = L123_DPhi_cut2 - pixel_pixel_dphi_lower_width;

	L124_DPhi_cut1 = SW_func1_dphi_v2(region, 7, EgEt);
	L124_DPhi_cut2 = SW_func1_dphi_v2(region, 7 ,EgEt);

	L124_DPhi_cut1 = L124_DPhi_cut1 + pixel_pixel_dphi_upper_width;
	L124_DPhi_cut2 = L124_DPhi_cut2 - pixel_pixel_dphi_lower_width;

	L134_DPhi_cut1 = SW_func1_dphi_v2(region, 8, EgEt);
	L134_DPhi_cut2 = SW_func1_dphi_v2(region, 8, EgEt);

	L134_DPhi_cut1 = L134_DPhi_cut1 + pixel_pixel_dphi_upper_width;
	L134_DPhi_cut2 = L134_DPhi_cut2 - pixel_pixel_dphi_lower_width;

	L234_DPhi_cut1 = SW_func1_dphi_v2(region, 9, EgEt);
	L234_DPhi_cut2 = SW_func1_dphi_v2(region, 9, EgEt);

	L234_DPhi_cut1 = L234_DPhi_cut1 + pixel_pixel_dphi_upper_width;
	L234_DPhi_cut2 = L234_DPhi_cut2 - pixel_pixel_dphi_lower_width;

	float pixel_pixel_deta_upper_width = sa_deta;
	float pixel_pixel_deta_lower_width = sa_deta;

	L123_DEta_cut1 = SW_func1_deta_v2(EgEt);
	L123_DEta_cut2 = SW_func1_deta_v2(EgEt);

	L123_DEta_cut1 = L123_DEta_cut1 + pixel_pixel_deta_upper_width;
	L123_DEta_cut2 = L123_DEta_cut2 - pixel_pixel_deta_lower_width;

	L124_DEta_cut1 = SW_func1_deta_v2(EgEt);
	L124_DEta_cut2 = SW_func1_deta_v2(EgEt);

	L124_DEta_cut1 = L124_DEta_cut1 + pixel_pixel_deta_upper_width;
	L124_DEta_cut2 = L124_DEta_cut2 - pixel_pixel_deta_lower_width;

	L134_DEta_cut1 = SW_func1_deta_v2(EgEt);
	L134_DEta_cut2 = SW_func1_deta_v2(EgEt);

	L134_DEta_cut1 = L134_DEta_cut1 + pixel_pixel_deta_upper_width;
	L134_DEta_cut2 = L134_DEta_cut2 - pixel_pixel_deta_lower_width;

	L234_DEta_cut1 = SW_func1_deta_v2(EgEt);
	L234_DEta_cut2 = SW_func1_deta_v2(EgEt);

	L234_DEta_cut1 = L234_DEta_cut1 + pixel_pixel_deta_upper_width;
	L234_DEta_cut2 = L234_DEta_cut2 - pixel_pixel_deta_lower_width;

}
void test::TriggeringWithout_4thPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

	dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit);
	dPhi013 = StandaloneDPhi( 0, 1, 3, 0, nthFirstHit, nthThirdHit );
	dPhi023 = StandaloneDPhi( 0, 2, 3, 0, nthSecondHit, nthThirdHit);

	dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);
	dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);

	dPhi_2 = EMmatchingDPhi(first_layer_hits[nthFirstHit], third_layer_hits[nthThirdHit], emvector);
	dEta_2 = EMmatchingDEta(first_layer_hits[nthFirstHit], third_layer_hits[nthThirdHit], emvector);

	dPhi_3 = EMmatchingDPhi(second_layer_hits[nthSecondHit], third_layer_hits[nthThirdHit], emvector);
	dEta_3 = EMmatchingDEta(second_layer_hits[nthSecondHit], third_layer_hits[nthThirdHit], emvector);

	L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
	L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );
	L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
	L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );
	L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
	L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );

	if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
		L12_EM_Ele = 1;
	if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
		L12_EM_Pos = 1;

	if( Signal_window_check(L13_phi_upper, dPhi_2, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_2, L13_eta_bellow, Ele)  )
		L13_EM_Ele = 1;
	if( Signal_window_check(L13_phi_upper, dPhi_2, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_2, L13_eta_bellow, Ele)  )
		L13_EM_Pos = 1;

	if( Signal_window_check(L23_phi_upper, dPhi_3, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_3, L23_eta_bellow, Ele)  )
		L23_EM_Ele = 1;
	if( Signal_window_check(L23_phi_upper, dPhi_3, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_3, L23_eta_bellow, Ele)  )
		L23_EM_Pos = 1;

	if( Signal_window_check(L123_DPhi_cut1, dPhi, L123_DPhi_cut2, Ele ) && Signal_window_check(L123_DEta_cut1, dEta, L123_DEta_cut2, Ele )  )
		L123_pass_Ele = 1;
	if( Signal_window_check(L123_DPhi_cut1, dPhi, L123_DPhi_cut2, Pos ) && Signal_window_check(L123_DEta_cut1, dEta, L123_DEta_cut2, Ele )  ) 
		L123_pass_Pos = 1;

}

void test::TriggeringWithout_3rdPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

	dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit);
	dPhi014 = StandaloneDPhi( 0, 1, 4, 0, nthFirstHit, nthThirdHit );
	dPhi024 = StandaloneDPhi( 0, 2, 4, 0, nthSecondHit, nthThirdHit);

	dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);
	dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);

	dPhi_2 = EMmatchingDPhi(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
	dEta_2 = EMmatchingDEta(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);

	dPhi_3 = EMmatchingDPhi(second_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
	dEta_3 = EMmatchingDEta(second_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);

	L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
	L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );
	L014_pass_Ele = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Ele );
	L014_pass_Pos = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Pos );
	L024_pass_Ele = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Ele );
	L024_pass_Pos = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Pos );

	if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  ) 
		L12_EM_Ele = 1;
	if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
		L12_EM_Pos = 1;

	if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Ele) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele)  )
		L14_EM_Ele = 1;
	if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Pos) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele)  )
		L14_EM_Pos = 1;

	if( Signal_window_check(L24_phi_upper, dPhi_3, L24_phi_bellow, Ele) && Signal_window_check(L24_eta_upper, dEta_3, L24_eta_bellow, Ele)  ) 
		L24_EM_Ele = 1;
	if( Signal_window_check(L24_phi_upper, dPhi_3, L24_phi_bellow, Pos) && Signal_window_check(L24_eta_upper, dEta_3, L24_eta_bellow, Ele)  )
		L24_EM_Pos = 1;

	if( Signal_window_check(L124_DPhi_cut1, dPhi, L124_DPhi_cut2, Ele ) && Signal_window_check(L124_DEta_cut1, dEta, L124_DEta_cut2, Ele )  )
		L124_pass_Ele = 1;
	if( Signal_window_check(L124_DPhi_cut1, dPhi, L124_DPhi_cut2, Pos ) && Signal_window_check(L124_DEta_cut1, dEta, L124_DEta_cut2, Ele )  )
		L124_pass_Pos = 1;
}

void test::TriggeringWithout_2ndPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

	dPhi013 = StandaloneDPhi( 0, 1, 3, 0, nthFirstHit, nthSecondHit);
	dPhi014 = StandaloneDPhi( 0, 1, 4, 0, nthFirstHit, nthThirdHit );
	dPhi034 = StandaloneDPhi( 0, 3, 4, 0, nthSecondHit, nthThirdHit);


	dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
	dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);

	dPhi_2 = EMmatchingDPhi(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
	dEta_2 = EMmatchingDEta(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);

	dPhi_3 = EMmatchingDPhi(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
	dEta_3 = EMmatchingDEta(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);

	L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
	L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );
	L014_pass_Ele = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Ele );
	L014_pass_Pos = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Pos );
	L034_pass_Ele = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Ele );
	L034_pass_Pos = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Pos );

	if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele)  ) 
		L13_EM_Ele = 1;
	if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele)  )
		L13_EM_Pos = 1;

	if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Ele) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele)  )
		L14_EM_Ele = 1;
	if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Pos) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele)  )
		L14_EM_Pos = 1;

	if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Ele) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele)  )
		L34_EM_Ele = 1;
	if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Pos) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele)  ) 
		L34_EM_Pos = 1;

	if( Signal_window_check(L134_DPhi_cut1, dPhi, L134_DPhi_cut2, Ele ) && Signal_window_check(L134_DEta_cut1, dEta, L134_DEta_cut2, Ele )  )
		L134_pass_Ele = 1;
	if( Signal_window_check(L134_DPhi_cut1, dPhi, L134_DPhi_cut2, Pos ) && Signal_window_check(L134_DEta_cut1, dEta, L134_DEta_cut2, Ele )  ) 
		L134_pass_Pos = 1;
}

void test::TriggeringWithout_1stPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

	dPhi023 = StandaloneDPhi( 0, 2, 3, 0, nthFirstHit, nthSecondHit);
	dPhi024 = StandaloneDPhi( 0, 2, 4, 0, nthFirstHit, nthThirdHit );
	dPhi034 = StandaloneDPhi( 0, 3, 4, 0, nthSecondHit, nthThirdHit);

	dPhi_1 = EMmatchingDPhi(second_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
	dEta_1 = EMmatchingDEta(second_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);

	dPhi_2 = EMmatchingDPhi(second_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
	dEta_2 = EMmatchingDEta(second_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);

	dPhi_3 = EMmatchingDPhi(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
	dEta_3 = EMmatchingDEta(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);

	L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
	L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );
	L024_pass_Ele = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Ele );
	L024_pass_Pos = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Pos );
	L034_pass_Ele = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Ele );
	L034_pass_Pos = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Pos );

	if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele)  )
		L23_EM_Ele = 1;
	if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele)  )
		L23_EM_Pos = 1;

	if( Signal_window_check(L24_phi_upper, dPhi_2, L24_phi_bellow, Ele) && Signal_window_check(L24_eta_upper, dEta_2, L24_eta_bellow, Ele)  )
		L24_EM_Ele = 1;
	if( Signal_window_check(L24_phi_upper, dPhi_2, L24_phi_bellow, Pos) && Signal_window_check(L24_eta_upper, dEta_2, L24_eta_bellow, Ele)  )
		L24_EM_Pos = 1;

	if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Ele) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele)  )
		L34_EM_Ele = 1;
	if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Pos) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele)  )
		L34_EM_Pos = 1;

	if( Signal_window_check(L234_DPhi_cut1, dPhi, L234_DPhi_cut2, Ele ) && Signal_window_check(L234_DEta_cut1, dEta, L234_DEta_cut2, Ele )  )
		L234_pass_Ele = 1;
	if( Signal_window_check(L234_DPhi_cut1, dPhi, L234_DPhi_cut2, Pos ) && Signal_window_check(L234_DEta_cut1, dEta, L234_DEta_cut2, Ele )  )
		L234_pass_Pos = 1;
}
void test::TriggeringWith_1st2ndPixel(int nthFirstHit, int nthSecondHit){
	dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
	L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
	L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );

	if( L012_pass_Ele &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(second_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || second_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
	if( L012_pass_Pos &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(second_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || second_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_1st3rdPixel(int nthFirstHit, int nthSecondHit){
	dPhi013 = StandaloneDPhi( 0, 1, 3, 0, nthFirstHit, nthSecondHit );
	L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
	L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );

	if( L013_pass_Ele &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
	if( L013_pass_Pos &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;

}
void test::TriggeringWith_2nd3rdPixel(int nthFirstHit, int nthSecondHit){
	dPhi023 = StandaloneDPhi( 0, 2, 3, 0, nthFirstHit, nthSecondHit );
	L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
	L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );

	if( L023_pass_Ele &&
			(second_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || second_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
	if( L023_pass_Pos &&
			(second_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || second_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_1st2ndPixel_v2(int nthFirstHit, int nthSecondHit){
	dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
	L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
	L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );

	dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);
	dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);

	if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
		L12_EM_Ele = 1;
	if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
		L12_EM_Pos = 1;

	if( L012_pass_Ele && L12_EM_Ele &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(second_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || second_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
	if( L012_pass_Pos && L12_EM_Pos &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(second_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || second_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_1st3rdPixel_v2(int nthFirstHit, int nthSecondHit){
	dPhi013 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
	L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
	L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );

	dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
	dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);

	if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele)  )
		L13_EM_Ele = 1;
	if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele)  )
		L13_EM_Pos = 1;

	if( L013_pass_Ele && L13_EM_Ele &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
	if( L013_pass_Pos && L13_EM_Pos &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_1st4thPixel_v2(int nthFirstHit, int nthSecondHit){
	dPhi014 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
	L014_pass_Ele = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Ele );
	L014_pass_Pos = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Pos );

	dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], fourth_layer_hits[nthSecondHit], emvector);
	dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], fourth_layer_hits[nthSecondHit], emvector);

	if( Signal_window_check(L14_phi_upper, dPhi_1, L14_phi_bellow, Ele) && Signal_window_check(L14_eta_upper, dEta_1, L14_eta_bellow, Ele)  )
		L14_EM_Ele = 1;
	if( Signal_window_check(L14_phi_upper, dPhi_1, L14_phi_bellow, Pos) && Signal_window_check(L14_eta_upper, dEta_1, L14_eta_bellow, Ele)  )
		L14_EM_Pos = 1;

	if( L014_pass_Ele && L14_EM_Ele &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(fourth_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || fourth_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
	if( L014_pass_Pos && L14_EM_Pos &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(fourth_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || fourth_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_2nd3rdPixel_v2(int nthFirstHit, int nthSecondHit){
	dPhi023 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
	L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
	L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );

	dPhi_1 = EMmatchingDPhi(second_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
	dEta_1 = EMmatchingDEta(second_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);

	if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele)  )
		L23_EM_Ele = 1;
	if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele)  )
		L23_EM_Pos = 1;

	if( L023_pass_Ele && L23_EM_Ele &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
	if( L023_pass_Pos && L23_EM_Pos &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_2nd4thPixel_v2(int nthFirstHit, int nthSecondHit){
	dPhi024 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
	L024_pass_Ele = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Ele );
	L024_pass_Pos = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Pos );

	dPhi_1 = EMmatchingDPhi(second_layer_hits[nthFirstHit], fourth_layer_hits[nthSecondHit], emvector);
	dEta_1 = EMmatchingDEta(second_layer_hits[nthFirstHit], fourth_layer_hits[nthSecondHit], emvector);

	if( Signal_window_check(L24_phi_upper, dPhi_1, L24_phi_bellow, Ele) && Signal_window_check(L24_eta_upper, dEta_1, L24_eta_bellow, Ele)  )
		L24_EM_Ele = 1;
	if( Signal_window_check(L24_phi_upper, dPhi_1, L24_phi_bellow, Pos) && Signal_window_check(L24_eta_upper, dEta_1, L24_eta_bellow, Ele)  )
		L24_EM_Pos = 1;

	if( L024_pass_Ele && L24_EM_Ele &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
	if( L024_pass_Pos && L24_EM_Pos &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_3rd4thPixel_v2(int nthFirstHit, int nthSecondHit){
	dPhi034 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
	L034_pass_Ele = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Ele );
	L034_pass_Pos = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Pos );

	dPhi_1 = EMmatchingDPhi(third_layer_hits[nthFirstHit], fourth_layer_hits[nthSecondHit], emvector);
	dEta_1 = EMmatchingDEta(third_layer_hits[nthFirstHit], fourth_layer_hits[nthSecondHit], emvector);

	if( Signal_window_check(L34_phi_upper, dPhi_1, L34_phi_bellow, Ele) && Signal_window_check(L34_eta_upper, dEta_1, L34_eta_bellow, Ele)  )
		L34_EM_Ele = 1;
	if( Signal_window_check(L34_phi_upper, dPhi_1, L34_phi_bellow, Pos) && Signal_window_check(L34_eta_upper, dEta_1, L34_eta_bellow, Ele)  )
		L34_EM_Pos = 1;

	if( L034_pass_Ele && L34_EM_Ele &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
	if( L034_pass_Pos && L34_EM_Pos &&
			(first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
			(third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

#endif // #ifdef test_cxx
