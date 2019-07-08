//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jul  7 22:39:31 2019 by ROOT version 6.16/00
// from TTree L1PiXTRKTree/L1PiXTRKTree
// found on file: results.root
//////////////////////////////////////////////////////////

#ifndef SW_h
#define SW_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class SW {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *genPartPhi;
   vector<float>   *genPartEta;
   vector<float>   *genPartPt;
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
   TBranch        *b_genPartPhi;   //!
   TBranch        *b_genPartEta;   //!
   TBranch        *b_genPartPt;   //!
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

   SW(TTree *tree=0);
   virtual ~SW();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   TVector3 emvector;
   float EgEt;
   float EgEta;
   float EgPhi;

   int layers[5];
   //
   std::vector<TVector3> first_layer_hits;
   std::vector<TVector3> second_layer_hits;
   std::vector<TVector3> third_layer_hits;
   std::vector<TVector3> fourth_layer_hits;

   void StorePixelHit(int region);

   inline float deltaPhi(float phi1, float phi2) {
     float result = phi1 - phi2;
     while (result > float(M_PI)) result -= float(2*M_PI);
     while (result <= -float(M_PI)) result += float(2*M_PI);
     return result;
   }


   TFile *file;
   TTree* sw_tree;

   // output
   vector<float> ntptErr;
   vector<float> ntEgEt;
   vector<float> ntEgEta;
   vector<float> ntEgPhi;

   // pixel-EG dphi
   vector<float> ntPix1EGdphi;
   vector<float> ntPix2EGdphi;
   vector<float> ntPix3EGdphi;
   vector<float> ntPix4EGdphi;

   // pixel segment-EG dphi
   vector<float> ntPix12EGdphi;
   vector<float> ntPix13EGdphi;
   vector<float> ntPix14EGdphi;
   vector<float> ntPix23EGdphi;
   vector<float> ntPix24EGdphi;
   vector<float> ntPix34EGdphi;

   // pixel segment-EG deta
   vector<float> ntPix12EGdeta;
   vector<float> ntPix13EGdeta;
   vector<float> ntPix14EGdeta;
   vector<float> ntPix23EGdeta;
   vector<float> ntPix24EGdeta;
   vector<float> ntPix34EGdeta;

   // pixel segment-pixel segment
   vector<float> ntPix012dphi;
   vector<float> ntPix013dphi;
   vector<float> ntPix014dphi;
   vector<float> ntPix023dphi;
   vector<float> ntPix024dphi;
   vector<float> ntPix034dphi;
   vector<float> ntPix123dphi;
   vector<float> ntPix124dphi;
   vector<float> ntPix134dphi;
   vector<float> ntPix234dphi;
   vector<float> ntPix123deta;
   vector<float> ntPix124deta;
   vector<float> ntPix134deta;
   vector<float> ntPix234deta;
};

#endif

#ifdef SW_cxx
SW::SW(TTree *tree) : fChain(0) 
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

   file = new TFile("sw.root","recreate");
   sw_tree = new TTree("t","t");

   sw_tree->Branch("ntEgEt",&ntEgEt);
   sw_tree->Branch("ntEgEta",&ntEgEta);
   sw_tree->Branch("ntEgPhi",&ntEgPhi);
   sw_tree->Branch("ntptErr",&ntptErr);

   sw_tree->Branch("ntPix1EGdphi",&ntPix1EGdphi);
   sw_tree->Branch("ntPix2EGdphi",&ntPix2EGdphi);
   sw_tree->Branch("ntPix3EGdphi",&ntPix3EGdphi);
   sw_tree->Branch("ntPix4EGdphi",&ntPix4EGdphi);

   sw_tree->Branch("ntPix12EGdphi",&ntPix12EGdphi);
   sw_tree->Branch("ntPix13EGdphi",&ntPix13EGdphi);
   sw_tree->Branch("ntPix14EGdphi",&ntPix14EGdphi);
   sw_tree->Branch("ntPix23EGdphi",&ntPix23EGdphi);
   sw_tree->Branch("ntPix24EGdphi",&ntPix24EGdphi);
   sw_tree->Branch("ntPix34EGdphi",&ntPix34EGdphi);

   sw_tree->Branch("ntPix12EGdeta",&ntPix12EGdeta);
   sw_tree->Branch("ntPix13EGdeta",&ntPix13EGdeta);
   sw_tree->Branch("ntPix14EGdeta",&ntPix14EGdeta);
   sw_tree->Branch("ntPix23EGdeta",&ntPix23EGdeta);
   sw_tree->Branch("ntPix24EGdeta",&ntPix24EGdeta);
   sw_tree->Branch("ntPix34EGdeta",&ntPix34EGdeta);

   sw_tree->Branch("ntPix012dphi",&ntPix012dphi);
   sw_tree->Branch("ntPix013dphi",&ntPix013dphi);
   sw_tree->Branch("ntPix014dphi",&ntPix014dphi);
   sw_tree->Branch("ntPix023dphi",&ntPix023dphi);
   sw_tree->Branch("ntPix024dphi",&ntPix024dphi);
   sw_tree->Branch("ntPix034dphi",&ntPix034dphi);
   sw_tree->Branch("ntPix123dphi",&ntPix123dphi);
   sw_tree->Branch("ntPix124dphi",&ntPix124dphi);
   sw_tree->Branch("ntPix134dphi",&ntPix134dphi);
   sw_tree->Branch("ntPix234dphi",&ntPix234dphi);
   sw_tree->Branch("ntPix123deta",&ntPix123deta);
   sw_tree->Branch("ntPix124deta",&ntPix124deta);
   sw_tree->Branch("ntPix134deta",&ntPix134deta);
   sw_tree->Branch("ntPix234deta",&ntPix234deta);

}

SW::~SW()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SW::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SW::LoadTree(Long64_t entry)
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

void SW::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   genPartPhi = 0;
   genPartEta = 0;
   genPartPt = 0;
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

   fChain->SetBranchAddress("genPartPhi", &genPartPhi, &b_genPartPhi);
   fChain->SetBranchAddress("genPartEta", &genPartEta, &b_genPartEta);
   fChain->SetBranchAddress("genPartPt", &genPartPt, &b_genPartPt);
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

void SW::StorePixelHit(int region){

 // barrel pixels
 for(int a=0; a<bRecHitN; a++){

    double Dphi = 0.;
    TVector3 current_hit;
    current_hit.SetXYZ( bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a) );
    Dphi = deltaPhi( current_hit.Phi(), EgPhi);
    if(fabs(Dphi) > 0.1) continue;

    if( region == 1 ){
      if( bRecHitLayer->at(a) == 1 ){ // First layer
           layers[1]++;
           first_layer_hits.push_back(current_hit);
      }
      if( bRecHitLayer->at(a) == 2 ){ // First layer
           layers[2]++;
           second_layer_hits.push_back(current_hit);
      }
      if( bRecHitLayer->at(a) == 3 ){ // First layer
           layers[3]++;
           third_layer_hits.push_back(current_hit);
      }
      if( bRecHitLayer->at(a) == 4 ){ // First layer
           layers[4]++;
           fourth_layer_hits.push_back(current_hit);
      }
    } // |eta| < 0.8  

    if( region == 2 ){
      if( bRecHitLayer->at(a) == 1 ){ // First layer
           layers[1]++;
           first_layer_hits.push_back(current_hit);
      }   
      if( bRecHitLayer->at(a) == 2 ){ // First layer
           layers[2]++;
           second_layer_hits.push_back(current_hit);
      }   
      if( bRecHitLayer->at(a) == 3 ){ // First layer
           layers[3]++;
           third_layer_hits.push_back(current_hit);
      }   
    } // 0.8 < |eta| < 1.4 

    if( region == 3 ){
      if( bRecHitLayer->at(a) == 1 ){ // First layer
           layers[1]++;
           first_layer_hits.push_back(current_hit);
      }   
      if( bRecHitLayer->at(a) == 2 ){ // First layer
           layers[2]++;
           second_layer_hits.push_back(current_hit);
      }   
    } // 1.4 < |eta| < 1.7 

    if( region == 4 ){
      if( bRecHitLayer->at(a) == 1 ){ // First layer
           layers[1]++;
           first_layer_hits.push_back(current_hit);
      }   
    } // 1.7 < |eta| < 2.1
 }

 // disk pixels
 for(int a=0; a<fRecHitN; a++){

    double Dphi = 0.;
    TVector3 current_hit;
    current_hit.SetXYZ( fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) );
    Dphi = deltaPhi( current_hit.Phi(), EgPhi);
    if(fabs(Dphi) > 0.1) continue;

    if( region == 2 ){
      if( fRecHitDisk->at(a) == 1 ){ // First disk
           layers[4]++;
           fourth_layer_hits.push_back(current_hit);
      }   
    } // 0.8 < |eta| < 1.4 

    if( region == 3 ){
      if( fRecHitDisk->at(a) == 1 ){ // First disk
           layers[3]++;
           third_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 2 ){ // First disk
           layers[4]++;
           fourth_layer_hits.push_back(current_hit);
      }
    } // 1.4 < |eta| < 1.7

    if( region == 4 ){
      if( fRecHitDisk->at(a) == 1 ){ // First disk
           layers[2]++;
           second_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 2 ){ // Second disk
           layers[3]++;
           third_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 3 ){ // Third disk
           layers[4]++;
           fourth_layer_hits.push_back(current_hit);
      }
    } // 1.7 < |eta| < 2.1

    if( region == 5 ){
      if( fRecHitDisk->at(a) == 1 ){ // First disk
           layers[1]++;
           first_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 2 ){ // Second disk
           layers[2]++;
           second_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 3 ){ // Third disk
           layers[3]++;
           third_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 4 ){ // Fourth disk
           layers[4]++;
           fourth_layer_hits.push_back(current_hit);
      }
    } // 2.7 < |eta| < 3.0

    if( region == 6 ){
      if( fRecHitDisk->at(a) == 2 ){ // First disk
           layers[1]++;
           first_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 3 ){ // Second disk
           layers[2]++;
           second_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 4 ){ // Third disk
           layers[3]++;
           third_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 5 ){ // Fourth disk
           layers[4]++;
           fourth_layer_hits.push_back(current_hit);
      }
    } // 2.7 < |eta| < 3.0

 }

}

Bool_t SW::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SW::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SW::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SW_cxx
