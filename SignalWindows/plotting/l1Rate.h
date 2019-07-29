//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 12 09:06:21 2018 by ROOT version 6.10/05
// from TTree t/t
// found on file: Tree_MinBias_PU200.root
//////////////////////////////////////////////////////////

#ifndef l1Rate_h
#define l1Rate_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class l1Rate {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           totalEvent;
   Float_t         totalEgN;
   Int_t           ntnEg2;
   Int_t           L1TkEleN;
   Int_t           L1TkEleIsoN;
   vector<float>   *ntEgEt;
   vector<float>   *ntL1TkEleEt;
   vector<float>   *ntL1TkEleIsoEt;
   vector<float>   *ntEgEta;
   vector<bool>    *ntCl_match;
   vector<bool>    *ntCl_iso_match;
   vector<float>   *ntL1TkEgEt;
   vector<float>   *ntL1TkEgEta;
   vector<float>   *ntL1TkEgPhi;
   vector<float>   *ntL1TkLooseEgEt;
   vector<float>   *ntL1TkLooseEgEta;
   vector<float>   *ntL1TkLooseEgPhi;

   // List of branches
   TBranch        *b_count_Entry;   //!
   TBranch        *b_EgN;   //!
   TBranch        *b_ntnEg2;   //!
   TBranch        *b_event_denominator;   //!
   TBranch        *b_event_nominator;   //!
   TBranch        *b_L1TkEleN;   //!
   TBranch        *b_L1TkEleIsoN;   //!
   TBranch        *b_ntEgEt;   //!
   TBranch        *b_ntL1TkEleEt;   //!
   TBranch        *b_ntL1TkEleIsoEt;   //!
   TBranch        *b_ntEgEta;   //!
   TBranch        *b_ntCl_match;   //!
   TBranch        *b_ntCl_iso_match;   //!
   TBranch        *b_ntL1TkEgEt;   //!
   TBranch        *b_ntL1TkEgEta;   //!
   TBranch        *b_ntL1TkEgPhi;   //!
   TBranch        *b_ntL1TkLooseEgEt;   //!
   TBranch        *b_ntL1TkLooseEgEta;   //!
   TBranch        *b_ntL1TkLooseEgPhi;   //!

   l1Rate(TTree *tree=0);
   virtual ~l1Rate();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef l1Rate_cxx
l1Rate::l1Rate(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("result_all.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("result_all.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
}

l1Rate::~l1Rate()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t l1Rate::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t l1Rate::LoadTree(Long64_t entry)
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

void l1Rate::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ntEgEt = 0;
   ntL1TkEleEt = 0;
   ntL1TkEleIsoEt = 0;
   ntEgEta = 0;
   ntCl_match = 0;
   ntCl_iso_match = 0;
   ntL1TkEgEt = 0;
   ntL1TkEgEta = 0;
   ntL1TkEgPhi = 0;
   ntL1TkLooseEgEt = 0;
   ntL1TkLooseEgEta = 0;
   ntL1TkLooseEgPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("totalEvent", &totalEvent, &b_count_Entry);
   fChain->SetBranchAddress("totalEgN", &totalEgN, &b_EgN);
   fChain->SetBranchAddress("ntnEg2", &ntnEg2, &b_ntnEg2);
   fChain->SetBranchAddress("L1TkEleN", &L1TkEleN, &b_L1TkEleN);
   fChain->SetBranchAddress("L1TkEleIsoN", &L1TkEleIsoN, &b_L1TkEleIsoN);
   fChain->SetBranchAddress("ntEgEt", &ntEgEt, &b_ntEgEt);
   fChain->SetBranchAddress("ntL1TkEleEt", &ntL1TkEleEt, &b_ntL1TkEleEt);
   fChain->SetBranchAddress("ntL1TkEleIsoEt", &ntL1TkEleIsoEt, &b_ntL1TkEleIsoEt);
   fChain->SetBranchAddress("ntEgEta", &ntEgEta, &b_ntEgEta);
   fChain->SetBranchAddress("ntCl_match", &ntCl_match, &b_ntCl_match);
   fChain->SetBranchAddress("ntCl_iso_match", &ntCl_iso_match, &b_ntCl_iso_match);
   fChain->SetBranchAddress("ntL1TkEgEt", &ntL1TkEgEt, &b_ntL1TkEgEt);
   fChain->SetBranchAddress("ntL1TkEgEta", &ntL1TkEgEta, &b_ntL1TkEgEta);
   fChain->SetBranchAddress("ntL1TkEgPhi", &ntL1TkEgPhi, &b_ntL1TkEgPhi);
   fChain->SetBranchAddress("ntL1TkLooseEgEt", &ntL1TkLooseEgEt, &b_ntL1TkLooseEgEt);
   fChain->SetBranchAddress("ntL1TkLooseEgEta", &ntL1TkLooseEgEta, &b_ntL1TkLooseEgEta);
   fChain->SetBranchAddress("ntL1TkLooseEgPhi", &ntL1TkLooseEgPhi, &b_ntL1TkLooseEgPhi);
   Notify();
}

Bool_t l1Rate::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void l1Rate::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t l1Rate::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef l1Rate_cxx
