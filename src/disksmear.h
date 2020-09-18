//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov 28 14:40:42 2019 by ROOT version 6.06/00
// from TTree L1PiXTRKTree/L1PiXTRKTree
// found on file: result.root
//////////////////////////////////////////////////////////

#ifndef disksmear_h
#define disksmear_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class disksmear {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           tmPevent;
   vector<float>   *tmPpropgenElPartPhi;
   vector<float>   *tmPpropgenElPartEta;
   vector<float>   *tmPpropgenElPartPt;
   vector<float>   *tmPpropgenElPartX;
   vector<float>   *tmPpropgenElPartY;
   vector<float>   *tmPpropgenElPartZ;
   Int_t           tmPbRecHitN;
   vector<int>     *tmPbRecHitLayer;
   vector<float>   *tmPbRecHitGx;
   vector<float>   *tmPbRecHitGy;
   vector<float>   *tmPbRecHitGz;
   Int_t           tmPfRecHitN;
   vector<int>     *tmPfRecHitDisk;
   vector<float>   *tmPfRecHitGx;
   vector<float>   *tmPfRecHitGy;
   vector<float>   *tmPfRecHitGz;
   vector<float>   *tmPegCrysClusterEem;
   vector<float>   *tmPegCrysClusterEhad;
   vector<float>   *tmPegCrysClusterEt;
   vector<float>   *tmPegCrysClusterEta;
   vector<float>   *tmPegCrysClusterPhi;
   vector<float>   *tmPegCrysClusterGx;
   vector<float>   *tmPegCrysClusterGy;
   vector<float>   *tmPegCrysClusterGz;

   // List of branches
   TBranch        *b_tmPevent;   //!
   TBranch        *b_tmPpropgenElPartPhi;   //!
   TBranch        *b_tmPpropgenElPartEta;   //!
   TBranch        *b_tmPpropgenElPartPt;   //!
   TBranch        *b_tmPpropgenElPartX;   //!
   TBranch        *b_tmPpropgenElPartY;   //!
   TBranch        *b_tmPpropgenElPartZ;   //!
   TBranch        *b_tmPbRecHitN;   //!
   TBranch        *b_tmPbRecHitLayer;   //!
   TBranch        *b_tmPbRecHitGx;   //!
   TBranch        *b_tmPbRecHitGy;   //!
   TBranch        *b_tmPbRecHitGz;   //!
   TBranch        *b_tmPfRecHitN;   //!
   TBranch        *b_tmPfRecHitDisk;   //!
   TBranch        *b_tmPfRecHitGx;   //!
   TBranch        *b_tmPfRecHitGy;   //!
   TBranch        *b_tmPfRecHitGz;   //!
   TBranch        *b_tmPegCrysClusterEem;   //!
   TBranch        *b_tmPegCrysClusterEhad;   //!
   TBranch        *b_tmPegCrysClusterEt;   //!
   TBranch        *b_tmPegCrysClusterEta;   //!
   TBranch        *b_tmPegCrysClusterPhi;   //!
   TBranch        *b_tmPegCrysClusterGx;   //!
   TBranch        *b_tmPegCrysClusterGy;   //!
   TBranch        *b_tmPegCrysClusterGz;   //!

   disksmear(TTree *tree=0);
   virtual ~disksmear();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   Int_t           event;
   vector<float>   propgenElPartPhi;
   vector<float>   propgenElPartEta;
   vector<float>   propgenElPartPt;
   vector<float>   propgenElPartX;
   vector<float>   propgenElPartY;
   vector<float>   propgenElPartZ;
   Int_t           bRecHitN;
   vector<int>     bRecHitLayer;
   vector<float>   bRecHitGx;
   vector<float>   bRecHitGy;
   vector<float>   bRecHitGz;
   Int_t           fRecHitN;
   vector<int>     fRecHitDisk;
   vector<float>   fRecHitGx;
   vector<float>   fRecHitGy;
   vector<float>   fRecHitGz;
   vector<float>   egCrysClusterEem;
   vector<float>   egCrysClusterEhad;
   vector<float>   egCrysClusterEt;
   vector<float>   egCrysClusterEta;
   vector<float>   egCrysClusterPhi;
   vector<float>   egCrysClusterGx;
   vector<float>   egCrysClusterGy;
   vector<float>   egCrysClusterGz;
   
   TFile *result;
   TTree *mytree;
};

#endif

#ifdef disksmear_cxx
disksmear::disksmear(TTree *tree) : fChain(0) 
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

   result = new TFile("results.root","RECREATE");
   result->mkdir("l1PiXTRKTree");
   result->cd("l1PiXTRKTree");

   mytree = new TTree("L1PiXTRKTree","L1PiXTRKTree");

   //mytree->Branch("event", &event, "event/I"); 
   mytree->Branch("propgenElPartPhi", &propgenElPartPhi);
   mytree->Branch("propgenElPartEta", &propgenElPartEta);
   mytree->Branch("propgenElPartPt", &propgenElPartPt);
   mytree->Branch("propgenElPartX", &propgenElPartX);
   mytree->Branch("propgenElPartY", &propgenElPartY);
   mytree->Branch("propgenElPartZ", &propgenElPartZ);

   mytree->Branch("bRecHitN", &bRecHitN, "bRecHitN/I");
   mytree->Branch("bRecHitLayer", &bRecHitLayer);
   mytree->Branch("bRecHitGx", &bRecHitGx);
   mytree->Branch("bRecHitGy", &bRecHitGy);
   mytree->Branch("bRecHitGz", &bRecHitGz);

   mytree->Branch("fRecHitN", &fRecHitN, "fRecHitN/I");
   mytree->Branch("fRecHitDisk", &fRecHitDisk);
   mytree->Branch("fRecHitGx", &fRecHitGx);
   mytree->Branch("fRecHitGy", &fRecHitGy);
   mytree->Branch("fRecHitGz", &fRecHitGz);

   mytree->Branch("egCrysClusterEem", &egCrysClusterEem);
   mytree->Branch("egCrysClusterEhad", &egCrysClusterEhad);
   mytree->Branch("egCrysClusterEt", &egCrysClusterEt);
   mytree->Branch("egCrysClusterEta", &egCrysClusterEta);
   mytree->Branch("egCrysClusterPhi", &egCrysClusterPhi);
   mytree->Branch("egCrysClusterGx",&egCrysClusterGx);
   mytree->Branch("egCrysClusterGy",&egCrysClusterGy);
   mytree->Branch("egCrysClusterGz",&egCrysClusterGz);
}

disksmear::~disksmear()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t disksmear::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t disksmear::LoadTree(Long64_t entry)
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

void disksmear::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   tmPpropgenElPartPhi = 0;
   tmPpropgenElPartEta = 0;
   tmPpropgenElPartPt = 0;
   tmPpropgenElPartX = 0;
   tmPpropgenElPartY = 0;
   tmPpropgenElPartZ = 0;
   tmPbRecHitLayer = 0;
   tmPbRecHitGx = 0;
   tmPbRecHitGy = 0;
   tmPbRecHitGz = 0;
   tmPfRecHitDisk = 0;
   tmPfRecHitGx = 0;
   tmPfRecHitGy = 0;
   tmPfRecHitGz = 0;
   tmPegCrysClusterEem = 0;
   tmPegCrysClusterEhad = 0;
   tmPegCrysClusterEt = 0;
   tmPegCrysClusterEta = 0;
   tmPegCrysClusterPhi = 0;
   tmPegCrysClusterGx = 0;
   tmPegCrysClusterGy = 0;
   tmPegCrysClusterGz = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("tmPevent", &tmPevent, &b_tmPevent);
   fChain->SetBranchAddress("tmPpropgenElPartPhi", &tmPpropgenElPartPhi, &b_tmPpropgenElPartPhi);
   fChain->SetBranchAddress("tmPpropgenElPartEta", &tmPpropgenElPartEta, &b_tmPpropgenElPartEta);
   fChain->SetBranchAddress("tmPpropgenElPartPt", &tmPpropgenElPartPt, &b_tmPpropgenElPartPt);
   fChain->SetBranchAddress("tmPpropgenElPartX", &tmPpropgenElPartX, &b_tmPpropgenElPartX);
   fChain->SetBranchAddress("tmPpropgenElPartY", &tmPpropgenElPartY, &b_tmPpropgenElPartY);
   fChain->SetBranchAddress("tmPpropgenElPartZ", &tmPpropgenElPartZ, &b_tmPpropgenElPartZ);
   fChain->SetBranchAddress("tmPbRecHitN", &tmPbRecHitN, &b_tmPbRecHitN);
   fChain->SetBranchAddress("tmPbRecHitLayer", &tmPbRecHitLayer, &b_tmPbRecHitLayer);
   fChain->SetBranchAddress("tmPbRecHitGx", &tmPbRecHitGx, &b_tmPbRecHitGx);
   fChain->SetBranchAddress("tmPbRecHitGy", &tmPbRecHitGy, &b_tmPbRecHitGy);
   fChain->SetBranchAddress("tmPbRecHitGz", &tmPbRecHitGz, &b_tmPbRecHitGz);
   fChain->SetBranchAddress("tmPfRecHitN", &tmPfRecHitN, &b_tmPfRecHitN);
   fChain->SetBranchAddress("tmPfRecHitDisk", &tmPfRecHitDisk, &b_tmPfRecHitDisk);
   fChain->SetBranchAddress("tmPfRecHitGx", &tmPfRecHitGx, &b_tmPfRecHitGx);
   fChain->SetBranchAddress("tmPfRecHitGy", &tmPfRecHitGy, &b_tmPfRecHitGy);
   fChain->SetBranchAddress("tmPfRecHitGz", &tmPfRecHitGz, &b_tmPfRecHitGz);
   fChain->SetBranchAddress("tmPegCrysClusterEem", &tmPegCrysClusterEem, &b_tmPegCrysClusterEem);
   fChain->SetBranchAddress("tmPegCrysClusterEhad", &tmPegCrysClusterEhad, &b_tmPegCrysClusterEhad);
   fChain->SetBranchAddress("tmPegCrysClusterEt", &tmPegCrysClusterEt, &b_tmPegCrysClusterEt);
   fChain->SetBranchAddress("tmPegCrysClusterEta", &tmPegCrysClusterEta, &b_tmPegCrysClusterEta);
   fChain->SetBranchAddress("tmPegCrysClusterPhi", &tmPegCrysClusterPhi, &b_tmPegCrysClusterPhi);
   fChain->SetBranchAddress("tmPegCrysClusterGx", &tmPegCrysClusterGx, &b_tmPegCrysClusterGx);
   fChain->SetBranchAddress("tmPegCrysClusterGy", &tmPegCrysClusterGy, &b_tmPegCrysClusterGy);
   fChain->SetBranchAddress("tmPegCrysClusterGz", &tmPegCrysClusterGz, &b_tmPegCrysClusterGz);
   Notify();
}

Bool_t disksmear::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void disksmear::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t disksmear::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef disksmear_cxx
