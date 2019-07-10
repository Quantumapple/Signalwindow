//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul  9 15:45:14 2019 by ROOT version 6.16/00
// from TTree t/t
// found on file: eff.root
//////////////////////////////////////////////////////////

#ifndef plot_h
#define plot_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class plot {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           totalEvent;
   Float_t         totalEgN;
   Int_t           ntnEg2;
   vector<float>   *ntEgEt;
   vector<float>   *ntEgEta;
   vector<float>   *ntEgPhi;
   vector<float>   *ntL1TkEgEt;
   vector<float>   *ntL1TkEgEta;
   vector<float>   *ntL1TkEgPhi;
   vector<float>   *ntL1TkLooseEgEt;
   vector<float>   *ntL1TkLooseEgEta;
   vector<float>   *ntL1TkLooseEgPhi;
   vector<int>     *PiXTRKbit;
   vector<int>     *trigger_bit_width;
   vector<int>     *trigger_bit_width_iso;
   vector<int>     *pix_comb;
   vector<bool>    *ntCl_match;
   vector<bool>    *withoutEM_match;
   vector<bool>    *withEM_match;
   Float_t         nt_genPhi;
   Float_t         nt_genEta;
   Float_t         nt_genPt;
   Float_t         nt_lastSimtkpt;
   Float_t         nt_initialSimtkpt;
   Int_t           fired;

   // List of branches
   TBranch        *b_count_Entry;   //!
   TBranch        *b_EgN;   //!
   TBranch        *b_ntnEg2;   //!
   TBranch        *b_ntEgEt;   //!
   TBranch        *b_ntEgEta;   //!
   TBranch        *b_ntEgPhi;   //!
   TBranch        *b_ntL1TkEgEt;   //!
   TBranch        *b_ntL1TkEgEta;   //!
   TBranch        *b_ntL1TkEgPhi;   //!
   TBranch        *b_ntL1TkLooseEgEt;   //!
   TBranch        *b_ntL1TkLooseEgEta;   //!
   TBranch        *b_ntL1TkLooseEgPhi;   //!
   TBranch        *b_PiXTRKbit;   //!
   TBranch        *b_trigger_bit_width;   //!
   TBranch        *b_trigger_bit_width_iso;   //!
   TBranch        *b_pix_comb;   //!
   TBranch        *b_ntCl_match;   //!
   TBranch        *b_withoutEM_match;   //!
   TBranch        *b_withEM_match;   //!
   TBranch        *b_nt_genPhi;   //!
   TBranch        *b_nt_genEta;   //!
   TBranch        *b_nt_genPt;   //!
   TBranch        *b_nt_lastSimtkpt;   //!
   TBranch        *b_nt_initialSimtkpt;   //!
   TBranch        *b_fired;   //!

   plot(TTree *tree=0);
   virtual ~plot();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef plot_cxx
plot::plot(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("eff.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("eff.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
}

plot::~plot()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t plot::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t plot::LoadTree(Long64_t entry)
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

void plot::Init(TTree *tree)
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
   ntEgEta = 0;
   ntEgPhi = 0;
   ntL1TkEgEt = 0;
   ntL1TkEgEta = 0;
   ntL1TkEgPhi = 0;
   ntL1TkLooseEgEt = 0;
   ntL1TkLooseEgEta = 0;
   ntL1TkLooseEgPhi = 0;
   PiXTRKbit = 0;
   trigger_bit_width = 0;
   trigger_bit_width_iso = 0;
   pix_comb = 0;
   ntCl_match = 0;
   withoutEM_match = 0;
   withEM_match = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("totalEvent", &totalEvent, &b_count_Entry);
   fChain->SetBranchAddress("totalEgN", &totalEgN, &b_EgN);
   fChain->SetBranchAddress("ntnEg2", &ntnEg2, &b_ntnEg2);
   fChain->SetBranchAddress("ntEgEt", &ntEgEt, &b_ntEgEt);
   fChain->SetBranchAddress("ntEgEta", &ntEgEta, &b_ntEgEta);
   fChain->SetBranchAddress("ntEgPhi", &ntEgPhi, &b_ntEgPhi);
   fChain->SetBranchAddress("ntL1TkEgEt", &ntL1TkEgEt, &b_ntL1TkEgEt);
   fChain->SetBranchAddress("ntL1TkEgEta", &ntL1TkEgEta, &b_ntL1TkEgEta);
   fChain->SetBranchAddress("ntL1TkEgPhi", &ntL1TkEgPhi, &b_ntL1TkEgPhi);
   fChain->SetBranchAddress("ntL1TkLooseEgEt", &ntL1TkLooseEgEt, &b_ntL1TkLooseEgEt);
   fChain->SetBranchAddress("ntL1TkLooseEgEta", &ntL1TkLooseEgEta, &b_ntL1TkLooseEgEta);
   fChain->SetBranchAddress("ntL1TkLooseEgPhi", &ntL1TkLooseEgPhi, &b_ntL1TkLooseEgPhi);
   fChain->SetBranchAddress("PiXTRKbit", &PiXTRKbit, &b_PiXTRKbit);
   fChain->SetBranchAddress("trigger_bit_width", &trigger_bit_width, &b_trigger_bit_width);
   fChain->SetBranchAddress("trigger_bit_width_iso", &trigger_bit_width_iso, &b_trigger_bit_width_iso);
   fChain->SetBranchAddress("pix_comb", &pix_comb, &b_pix_comb);
   fChain->SetBranchAddress("ntCl_match", &ntCl_match, &b_ntCl_match);
   fChain->SetBranchAddress("withoutEM_match", &withoutEM_match, &b_withoutEM_match);
   fChain->SetBranchAddress("withEM_match", &withEM_match, &b_withEM_match);
   fChain->SetBranchAddress("nt_genPhi", &nt_genPhi, &b_nt_genPhi);
   fChain->SetBranchAddress("nt_genEta", &nt_genEta, &b_nt_genEta);
   fChain->SetBranchAddress("nt_genPt", &nt_genPt, &b_nt_genPt);
   fChain->SetBranchAddress("nt_lastSimtkpt", &nt_lastSimtkpt, &b_nt_lastSimtkpt);
   fChain->SetBranchAddress("nt_initialSimtkpt", &nt_initialSimtkpt, &b_nt_initialSimtkpt);
   fChain->SetBranchAddress("fired", &fired, &b_fired);
   Notify();
}

Bool_t plot::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void plot::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t plot::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef plot_cxx
