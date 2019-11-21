//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov  8 02:00:58 2019 by ROOT version 6.18/00
// from TTree L1PiXTRKTree/L1PiXTRKTree
// found on file: result.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class test {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   vector<float>   *propgenElPartPhi;
   vector<float>   *propgenElPartEta;
   vector<float>   *propgenElPartPt;
   vector<int>     *propgenElPartPdgID;
   vector<float>   *propgenElPartX;
   vector<float>   *propgenElPartY;
   vector<float>   *propgenElPartZ;
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
   vector<float>   *egCrysClusterEt;
   vector<float>   *egCrysClusterEta;
   vector<float>   *egCrysClusterPhi;
   vector<float>   *egCrysClusterGx;
   vector<float>   *egCrysClusterGy;
   vector<float>   *egCrysClusterGz;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_propgenElPartPhi;   //!
   TBranch        *b_propgenElPartEta;   //!
   TBranch        *b_propgenElPartPt;   //!
   TBranch        *b_propgenElPartPdgID;   //!
   TBranch        *b_propgenElPartX;   //!
   TBranch        *b_propgenElPartY;   //!
   TBranch        *b_propgenElPartZ;   //!
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

   void StorePixelHit( int region);

   vector<float> v1R;
   vector<float> v2R;
   vector<float> v3R;
   vector<float> v4R;
   vector<float> v1z;
   vector<float> v2z;
   vector<float> v3z;
   vector<float> v4z;

   std::vector<TVector3> first_layer_hits;
   std::vector<TVector3> second_layer_hits;
   std::vector<TVector3> third_layer_hits;
   std::vector<TVector3> fourth_layer_hits;

   inline float deltaPhi(float phi1, float phi2){
	   float result = phi1 - phi2;
	   while(result > float(M_PI)) result -= float(2*M_PI);
	   while(result <= -float(M_PI)) result += float(2*M_PI);
	   return result;
   }

};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("result.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("result.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("result.root:/l1PiXTRKTree");
      dir->GetObject("L1PiXTRKTree",tree);

   }
   Init(tree);
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
   propgenElPartPdgID = 0;
   propgenElPartX = 0;
   propgenElPartY = 0;
   propgenElPartZ = 0;
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

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("propgenElPartPhi", &propgenElPartPhi, &b_propgenElPartPhi);
   fChain->SetBranchAddress("propgenElPartEta", &propgenElPartEta, &b_propgenElPartEta);
   fChain->SetBranchAddress("propgenElPartPt", &propgenElPartPt, &b_propgenElPartPt);
   fChain->SetBranchAddress("propgenElPartPdgID", &propgenElPartPdgID, &b_propgenElPartPdgID);
   fChain->SetBranchAddress("propgenElPartX", &propgenElPartX, &b_propgenElPartX);
   fChain->SetBranchAddress("propgenElPartY", &propgenElPartY, &b_propgenElPartY);
   fChain->SetBranchAddress("propgenElPartZ", &propgenElPartZ, &b_propgenElPartZ);
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

void test::StorePixelHit(int region){

	for(int a=0; a<bRecHitN; a++){
		float R = sqrt(pow(bRecHitGx->at(a),2)+pow(bRecHitGy->at(a),2));
		float Z = bRecHitGz->at(a);

		if( region < 5 ){
			if( bRecHitLayer->at(a) == 1 ){ // First layer
				first_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
				v1R.push_back(R); v1z.push_back(Z); 
			}
		}
		if( region == 1 || region == 2 || region == 3 ){
			if( bRecHitLayer->at(a) == 2 ){ // Second layer
				second_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
				v2R.push_back(R); v2z.push_back(Z); 
			}
		}
		if( region == 1 || region == 2 ){
			if( bRecHitLayer->at(a) == 3 ){ // Third layer 
				third_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
				v3R.push_back(R); v3z.push_back(Z); 
			}
		}
		if( region == 1 ){
			if( bRecHitLayer->at(a) == 4 ){ // Fourth layer
				fourth_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
				v4R.push_back(R); v4z.push_back(Z); 
			}
		}
	} // bRecHit Loop

	for(int a=0; a<fRecHitN; a++){
		float R = sqrt(pow(fRecHitGx->at(a),2)+pow(fRecHitGy->at(a),2));
		float Z = fRecHitGz->at(a);

		if( region == 5 ){
			if( fRecHitDisk->at(a) == 1 ){ // first disk
				first_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v1R.push_back(R); v1z.push_back(Z); 
			}

			if( fRecHitDisk->at(a) == 2 ){ // second disk
				second_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v2R.push_back(R); v2z.push_back(Z); 
			}
			if( fRecHitDisk->at(a) == 3 ){ // third disk
				third_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v3R.push_back(R); v3z.push_back(Z); 
			}
			if( fRecHitDisk->at(a) == 4 ){ // fourth disk
				fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v4R.push_back(R); v4z.push_back(Z); 
			}
		}

		if( region == 6 ){
			if( fRecHitDisk->at(a) == 2 ){ // second disk
				first_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v1R.push_back(R); v1z.push_back(Z); 
			}
			if( fRecHitDisk->at(a) == 3 ){ // third disk
				second_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v2R.push_back(R); v2z.push_back(Z); 
			}
			if( fRecHitDisk->at(a) == 4 ){ // fourth disk
				third_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v3R.push_back(R); v3z.push_back(Z); 
			}
			if( fRecHitDisk->at(a) == 5 ){ // fifth disk
				fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v4R.push_back(R); v4z.push_back(Z); 
			}
		}

		if( region == 4 ){
			if( fRecHitDisk->at(a) == 1 ){ // first disk
				second_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v2R.push_back(R); v2z.push_back(Z); 
			}
			if( fRecHitDisk->at(a) == 2 ){ // second disk
				third_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v3R.push_back(R); v3z.push_back(Z); 
			}
			if( fRecHitDisk->at(a) == 3 ){ // third disk
				fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v4R.push_back(R); v4z.push_back(Z); 
			}
		}

		if( region == 3 ){
			if( fRecHitDisk->at(a) == 1 ){ // first disk
				third_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v3R.push_back(R); v3z.push_back(Z); 
			}
			if( fRecHitDisk->at(a) == 2 ){ // second disk
				fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v4R.push_back(R); v4z.push_back(Z); 
			}
		}
		if( region == 2 ){
			if( fRecHitDisk->at(a) == 1 ){ // first disk
				fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
				v4R.push_back(R); v4z.push_back(Z); 
			}
		}
	} // fRecHit Loop
}

#endif // #ifdef test_cxx
