#define disksmear_cxx
#include "disksmear.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void disksmear::Loop()
{
//   In a ROOT session, you can do:
//      root> .L disksmear.C
//      root> disksmear t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   //for (Long64_t jentry=0; jentry<1;jentry++) {
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //if( jentry%10000 == 0 ) cout << "Event: " << jentry << endl;

      propgenElPartEta.clear();
      propgenElPartPhi.clear();
      propgenElPartPt.clear();
      propgenElPartX.clear();
      propgenElPartY.clear();
      propgenElPartZ.clear();

      egCrysClusterEta.clear();
      egCrysClusterPhi.clear();
      egCrysClusterEem.clear();
      egCrysClusterEhad.clear();
      egCrysClusterEt.clear();
      egCrysClusterGx.clear();
      egCrysClusterGy.clear();
      egCrysClusterGz.clear();

      bRecHitN = 0;
      bRecHitGx.clear();
      bRecHitGy.clear();
      bRecHitGz.clear();
      bRecHitLayer.clear();

      fRecHitN = 0;
      fRecHitGx.clear();
      fRecHitGy.clear();
      fRecHitGz.clear();
      fRecHitDisk.clear();

      for(int i=0; i < tmPpropgenElPartEta->size(); i++){
          float genEta = tmPpropgenElPartEta->at(i);
          float genPhi = tmPpropgenElPartPhi->at(i);
          float genPt  = tmPpropgenElPartPt->at(i);
          float genX   = tmPpropgenElPartX->at(i);
          float genY   = tmPpropgenElPartY->at(i);
          float genZ   = tmPpropgenElPartZ->at(i);

          propgenElPartEta.push_back(genEta);
          propgenElPartPhi.push_back(genPhi);
          propgenElPartPt.push_back(genPt);
          propgenElPartX.push_back(genX);
          propgenElPartY.push_back(genY);
          propgenElPartZ.push_back(genZ);
      }

      for(int i=0; i < tmPegCrysClusterEt->size(); i++){
          float egEta = tmPegCrysClusterEta->at(i);
          float egPhi = tmPegCrysClusterPhi->at(i);
          float egEem  = tmPegCrysClusterEem->at(i);
          float egEhad  = tmPegCrysClusterEhad->at(i);
          float egEt  = tmPegCrysClusterEt->at(i);
          float egX   = tmPegCrysClusterGx->at(i);
          float egY   = tmPegCrysClusterGy->at(i);
          float egZ   = tmPegCrysClusterGz->at(i);

          egCrysClusterEta.push_back(egEta);
          egCrysClusterPhi.push_back(egPhi);
          egCrysClusterEem.push_back(egEem);
          egCrysClusterEhad.push_back(egEhad);
          egCrysClusterEt.push_back(egEt);
          egCrysClusterGx.push_back(egX);
          egCrysClusterGy.push_back(egY);
          egCrysClusterGz.push_back(egZ);
      }


      for(int i=0; i < tmPbRecHitGx->size(); i++){
          float bx = tmPbRecHitGx->at(i);
          float by = tmPbRecHitGy->at(i);
          float bz = tmPbRecHitGz->at(i);
          int   bl = tmPbRecHitLayer->at(i);

          bRecHitGx.push_back(bx);
          bRecHitGy.push_back(by);
          bRecHitGz.push_back(bz);
          bRecHitLayer.push_back(bl);
      }
      bRecHitN = bRecHitGx.size();

      for(int i=0; i < tmPfRecHitGx->size(); i++){

          float fx = tmPfRecHitGx->at(i);
          float fy = tmPfRecHitGy->at(i);
          float fz = tmPfRecHitGz->at(i);
          int   fd = tmPfRecHitDisk->at(i);
          float fz_smear, fx_smear, fy_smear;
          TLorentzVector vpix;

          //TVector3 pix(fx, fy, fz);
          vpix.SetVect(TVector3(fx,fy,fz));
          float pix_Eta = vpix.Eta();
          //cout << "pix_Eta : " << pix_Eta << endl;
          int eta_region = 0;
          if( fabs(pix_Eta) <= 0.8 ) eta_region =1; 
          if( fabs(pix_Eta) <= 1.4 && fabs(pix_Eta) > 0.8 ) eta_region =2;
          if( fabs(pix_Eta) <= 1.7 && fabs(pix_Eta) > 1.4 ) eta_region =3;
          if( fabs(pix_Eta) <= 2.1 && fabs(pix_Eta) > 1.7 ) eta_region =4;
          if( fabs(pix_Eta) <= 2.7 && fabs(pix_Eta) > 2.1 ) eta_region =5;
          if( fabs(pix_Eta) <= 3.0 && fabs(pix_Eta) > 2.7 ) eta_region =6;

          // x_smear = y_smear
          float x_smear = 0.00072;

          float z1_smear = 0.01;
          float z2_smear = 0.01;
          float z3_smear = 0.01;
          float z4_smear = 0.01;
          float z5_smear = 0.01;

          if( eta_region == 1 ){
              z1_smear = 0.01;

              x_smear = 0.00072;
          }    
          if( eta_region == 2 ){
              z1_smear = 0.0011;
              x_smear = 0.00072;

          }    
          if( eta_region == 3 ){
              z1_smear = 0.0044;
              z2_smear = 0.0055;

              x_smear = 0.00072;
          }    
          if( eta_region == 4 ){
              z1_smear = 0.0111;
              z2_smear = 0.0121;
              z3_smear = 0.0131;

              x_smear = 0.00072;
          }    
          if( eta_region == 5 ){
              z1_smear = 0.0091;
              z2_smear = 0.0101;
              z3_smear = 0.0117;
              z4_smear = 0.0125;

              x_smear = 0.00072;
          }    
          if( eta_region == 6 ){
              z2_smear = 0.0162; 
              z3_smear = 0.0162;
              z4_smear = 0.0162;
              z5_smear = 0.0162;

              x_smear = 0.00072;
          }

          fx_smear = gRandom->Gaus(fx,x_smear);
          fy_smear = gRandom->Gaus(fy,x_smear);
          if( fd == 1 ) fz_smear = gRandom->Gaus(fz,z1_smear);
          if( fd == 2 ) fz_smear = gRandom->Gaus(fz,z2_smear);
          if( fd == 3 ) fz_smear = gRandom->Gaus(fz,z3_smear);
          if( fd == 4 ) fz_smear = gRandom->Gaus(fz,z4_smear);
          if( fd == 5 ) fz_smear = gRandom->Gaus(fz,z5_smear);

          fRecHitGx.push_back(fx_smear);
          fRecHitGy.push_back(fy_smear);
          fRecHitGz.push_back(fz_smear);
          fRecHitDisk.push_back(fd);
      }
      fRecHitN = fRecHitGx.size();

      mytree->Fill();
   }
   result->Write();
}
