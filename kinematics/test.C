#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void test::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   TH1F *h1 = new TH1F("dz1st2nd","1st-2nd",1000,-0.1,0.1);
   TH1F *h2 = new TH1F("dz1st3rd","1st-3rd",1000,-0.1,0.1);
   TH1F *h3 = new TH1F("dz1st4th","1st-4th",1000,-0.1,0.1); 
   TH1F *h4 = new TH1F("dz2nd3rd","2nd-3rd",1000,-0.1,0.1);
   TH1F *h5 = new TH1F("dz2nd4th","2nd-4th",1000,-0.1,0.1); 
   TH1F *h6 = new TH1F("dz3rd4th","3rd-4th",1000,-0.1,0.1);

   TH1F *h1z = new TH1F("dz1st","1st",1000,-2*1e-6,2*1e-6); 
   TH1F *h2z = new TH1F("dz2nd","2nd",1000,-2*1e-6,2*1e-6); 
   TH1F *h3z = new TH1F("dz3rd","3rd",1000,-2*1e-6,2*1e-6); 
   TH1F *h4z = new TH1F("dz4th","4th",1000,-2*1e-6,2*1e-6); 
   
   TH1F *a1 = new TH1F("dEcl1st2nd-1st3rd","1st2nd - 1st3rd #Delta#eta",1000,-0.05,0.05);
   TH1F *a2 = new TH1F("dEcl1st2nd-1st4th","1st2nd - 1st4th #Delta#eta",1000,-0.05,0.05);
   TH1F *a3 = new TH1F("dEcl1st2nd-2nd3rd","1st2nd - 2nd3rd #Delta#eta",1000,-0.05,0.05);
   TH1F *a4 = new TH1F("dEcl1st2nd-2nd4th","1st2nd - 2nd4th #Delta#eta",1000,-0.05,0.05);
   
   TH1F *b1 = new TH1F("dEcl1st3rd-1st4th","1st3rd - 1st4th #Delta#eta",1000,-0.05,0.05);
   TH1F *b2 = new TH1F("dEcl1st3rd-2nd3rd","1st3rd - 2nd3rd #Delta#eta",1000,-0.05,0.05);
   TH1F *b3 = new TH1F("dEcl1st3rd-3rd4th","1st3rd - 3rd4th #Delta#eta",1000,-0.05,0.05);
   
   TH1F *c1 = new TH1F("dEcl1st4th-2nd4th","1st4th - 2nd4th #Delta#eta",1000,-0.05,0.05);
   TH1F *c2 = new TH1F("dEcl1st4th-3rd4th","1st4th - 3rd4th #Delta#eta",1000,-0.05,0.05);
   
   TH1F *d1 = new TH1F("dEcl2nd3rd-2nd4th","2nd3rd - 2nd4th #Delta#eta",1000,-0.05,0.05);
   TH1F *d2 = new TH1F("dEcl2nd3rd-3rd4th","2nd3rd - 3rd4th #Delta#eta",1000,-0.05,0.05);
   
   TH1F *e1 = new TH1F("dEcl2nd4th-3rd4th","2nd4th - 3rd4th #Delta#eta",1000,-0.05,0.05);

   TH1F *f1 = new TH1F("dEpv1st2nd","PV1st-PV2nd #Delta#eta",1000,-0.05,0.05);
   TH1F *f2 = new TH1F("dEpv1st3rd","PV1st-PV3rd #Delta#eta",1000,-0.05,0.05);
   TH1F *f3 = new TH1F("dEpv1st4th","PV1st-PV4th #Delta#eta",1000,-0.05,0.05); 
   TH1F *f4 = new TH1F("dEpv2nd3rd","PV2nd-PV3rd #Delta#eta",1000,-0.05,0.05);
   TH1F *f5 = new TH1F("dEpv2nd4th","PV2nd-PV4th #Delta#eta",1000,-0.05,0.05); 
   TH1F *f6 = new TH1F("dEpv3rd4th","PV3rd-PV4th #Delta#eta",1000,-0.05,0.05); 

   TH1F *g1 = new TH1F("dPhiDiff1", "[PV1st-1st2nd]-[1st2nd-2nd3rd]", 1000, -0.5, 0.5);
   TH1F *g2 = new TH1F("dPhiDiff2", "[PV1st-1st2nd]-[1st2nd-2nd4th]", 1000, -0.5, 0.5);
   TH1F *g3 = new TH1F("dPhiDiff3", "[PV1st-1st3rd]-[1st3rd-3rd4th]", 1000, -0.5, 0.5);
   TH1F *g4 = new TH1F("dPhiDiff4", "[PV2nd-2nd3rd]-[2nd3rd-3rd4th]", 1000, -0.5, 0.5);

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      v1R.clear();
      v2R.clear();
      v3R.clear();
      v4R.clear();
      v1z.clear();
      v2z.clear();
      v3z.clear();
      v4z.clear();
      first_layer_hits.clear();
      second_layer_hits.clear();
      third_layer_hits.clear();
      fourth_layer_hits.clear();

      int region = 0;
      float genEta = fabs(propgenElPartEta->at(0));
      if( genEta < 0.8 ) region = 1; 
      if( genEta > 0.8 && genEta < 1.4 ) region = 2; 
      if( genEta > 1.4 && genEta < 1.7 ) region = 3; 
      if( genEta > 1.7 && genEta < 2.1 ) region = 4; 
      if( genEta > 2.1 && genEta < 2.7 ) region = 5; 
      if( genEta > 2.7 && genEta < 3.0 ) region = 6;
      if( region != 1 ) continue;

      //cout << "Event: " << jentry << ", region: " << region << endl;
      //cout << bRecHitN << endl;
      StorePixelHit(region);

      int size1 = v1R.size();
      int size2 = v2R.size();
      int size3 = v3R.size();
      int size4 = v4R.size();
      float genX = propgenElPartX->at(0)/10.;
      float genY = propgenElPartY->at(0)/10.;
      float genZ = propgenElPartZ->at(0)/10.;
     
      // genZ - vertex (theta method)
      for(int j = 0; j < first_layer_hits.size(); j++)
      {
          TVector3 v1pv(first_layer_hits[j].X()-genX, first_layer_hits[j].Y()-genY, first_layer_hits[j].Z()-genZ);
          float slope = tan(v1pv.Theta());
          float R = sqrt(pow(first_layer_hits[j].X(),2)+pow(first_layer_hits[j].Y(),2));
          float Z = first_layer_hits[j].Z();
          float vertex = Z - R/slope;
          h1z->Fill(genZ-vertex);
      }
      
      for(int j = 0; j < second_layer_hits.size(); j++)
      {
          TVector3 v2pv(second_layer_hits[j].X()-genX, second_layer_hits[j].Y()-genY, second_layer_hits[j].Z()-genZ);
          float slope = tan(v2pv.Theta());
          float R = sqrt(pow(second_layer_hits[j].X(),2)+pow(second_layer_hits[j].Y(),2));
          float Z = second_layer_hits[j].Z();
          float vertex = Z - R/slope;
          h2z->Fill(genZ-vertex);
      }
      
      for(int j = 0; j < third_layer_hits.size(); j++)
      {
          TVector3 v3pv(third_layer_hits[j].X()-genX, third_layer_hits[j].Y()-genY, third_layer_hits[j].Z()-genZ);
          float slope = tan(v3pv.Theta());
          float R = sqrt(pow(third_layer_hits[j].X(),2)+pow(third_layer_hits[j].Y(),2));
          float Z = third_layer_hits[j].Z();
          float vertex = Z - R/slope;
          h3z->Fill(genZ-vertex);
      }
      
      for(int j = 0; j < fourth_layer_hits.size(); j++)
      {
          TVector3 v4pv(fourth_layer_hits[j].X()-genX, fourth_layer_hits[j].Y()-genY, fourth_layer_hits[j].Z()-genZ);
          float slope = tan(v4pv.Theta());
          float R = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
          float Z = fourth_layer_hits[j].Z();
          float vertex = Z - R/slope;
          h4z->Fill(genZ-vertex);
      }

      // genZ - vertex (two pixel clusters method)
      for(int j = 0; j < size1; j++)
      {
          for(int k = 0; k < size2; k++)
          {
              float nom = (v1z.at(j)*v2R.at(k) - v1R.at(j)*v2z.at(k));
              float denom = v2R.at(k)-v1R.at(j);
              h1->Fill(genZ-nom/denom);
          }
      }
      
      for(int j = 0; j < size1; j++)
      {
          for(int k = 0; k < size3; k++)
          {
              float nom = (v1z.at(j)*v3R.at(k) - v1R.at(j)*v3z.at(k));
              float denom = v3R.at(k)-v1R.at(j);
              h2->Fill(genZ-nom/denom);
          }
      }
      
      for(int j = 0; j < size1; j++)
      {
          for(int k = 0; k < size4; k++)
          {
              float nom = (v1z.at(j)*v4R.at(k) - v1R.at(j)*v4z.at(k));
              float denom = v4R.at(k)-v1R.at(j);
              h3->Fill(genZ-nom/denom);
          }
      }
      
      for(int j = 0; j < size2; j++)
      {
          for(int k = 0; k < size3; k++)
          {
              float nom = (v2z.at(j)*v3R.at(k) - v2R.at(j)*v3z.at(k));
              float denom = v3R.at(k)-v2R.at(j);
              h4->Fill(genZ-nom/denom);
          }
      }
      
      for(int j = 0; j < size2; j++)
      {
          for(int k = 0; k < size4; k++)
          {
              float nom = (v2z.at(j)*v4R.at(k) - v2R.at(j)*v4z.at(k));
              float denom = v4R.at(k)-v2R.at(j);
              h5->Fill(genZ-nom/denom);
          }
      }
      
      for(int j = 0; j < size3; j++)
      {
          for(int k = 0; k < size4; k++)
          {
              float nom = (v3z.at(j)*v4R.at(k) - v3R.at(j)*v4z.at(k));
              float denom = v4R.at(k)-v3R.at(j);
              h6->Fill(genZ-nom/denom);
          }
      }

      // delta Eta and delta Phi difference
      TVector3 v12;
      TVector3 v13;
      TVector3 v23;
      
      // 123 combi
      for(int j = 0; j < first_layer_hits.size(); j++)
      {
          TVector3 v1(first_layer_hits[j].X(), first_layer_hits[j].Y(), first_layer_hits[j].Z());
          TVector3 v1pv(first_layer_hits[j].X()-genX, first_layer_hits[j].Y()-genY, first_layer_hits[j].Z()-genZ);
          for(int k = 0; k < second_layer_hits.size(); k++)
          {
              TVector3 v2(second_layer_hits[k].X(), second_layer_hits[k].Y(), second_layer_hits[k].Z());
              TVector3 v2pv(second_layer_hits[k].X()-genX, second_layer_hits[k].Y()-genY, second_layer_hits[k].Z()-genZ);
              v12 = v2 - v1;
              float dphi1 = deltaPhi(v12.Phi(), v1pv.Phi());
              for(int l = 0; l < third_layer_hits.size(); l++)
              {
                  TVector3 v3(third_layer_hits[k].X(), third_layer_hits[k].Y(), third_layer_hits[k].Z());
                  TVector3 v3pv(third_layer_hits[l].X()-genX, third_layer_hits[l].Y()-genY, third_layer_hits[l].Z()-genZ);
                  v13 = v3 - v1;
                  v23 = v3 - v2;
                  float dphi2 = deltaPhi(v23.Phi(), v12.Phi());
                  a1->Fill(v13.Eta()-v12.Eta());
                  a3->Fill(v23.Eta()-v12.Eta());
                  b2->Fill(v23.Eta()-v13.Eta());
                  f1->Fill(v2pv.Eta()-v1pv.Eta());
                  f2->Fill(v3pv.Eta()-v1pv.Eta());
                  f4->Fill(v3pv.Eta()-v2pv.Eta());
                  g1->Fill(dphi2-dphi1);
              }
          }
      }
    
      
      TVector3 v14;
      TVector3 v24;
      
      // 124 combi
      for(int j = 0; j < first_layer_hits.size(); j++)
      {
          TVector3 v1(first_layer_hits[j].X(), first_layer_hits[j].Y(), first_layer_hits[j].Z());
          TVector3 v1pv(first_layer_hits[j].X()-genX, first_layer_hits[j].Y()-genY, first_layer_hits[j].Z()-genZ);
          for(int k = 0; k < second_layer_hits.size(); k++)
          {
              TVector3 v2(second_layer_hits[k].X(), second_layer_hits[k].Y(), second_layer_hits[k].Z());
              TVector3 v2pv(second_layer_hits[k].X()-genX, second_layer_hits[k].Y()-genY, second_layer_hits[k].Z()-genZ);
              v12 = v2 - v1;
              float dphi1 = deltaPhi(v12.Phi(), v1pv.Phi());
              for(int l = 0; l < fourth_layer_hits.size(); l++)
              {
                  TVector3 v4(fourth_layer_hits[l].X(), fourth_layer_hits[l].Y(), fourth_layer_hits[l].Z());
                  TVector3 v4pv(fourth_layer_hits[l].X()-genX, fourth_layer_hits[l].Y()-genY, fourth_layer_hits[l].Z()-genZ);
                  v14 = v4 - v1;
                  v24 = v4 - v2;
                  float dphi2 = deltaPhi(v24.Phi(), v12.Phi());
                  a2->Fill(v14.Eta()-v12.Eta()); 
                  a4->Fill(v24.Eta()-v12.Eta());  
                  c1->Fill(v24.Eta()-v14.Eta());
                  f3->Fill(v4pv.Eta()-v1pv.Eta());
                  f5->Fill(v4pv.Eta()-v2pv.Eta());
                  g2->Fill(dphi2-dphi1);
              }
          }
      }

      
      TVector3 v34;

      // 134 combi
      for(int j = 0; j < first_layer_hits.size(); j++)
      {
          TVector3 v1(first_layer_hits[j].X(), first_layer_hits[j].Y(), first_layer_hits[j].Z());
          TVector3 v1pv(first_layer_hits[j].X()-genX, first_layer_hits[j].Y()-genY, first_layer_hits[j].Z()-genZ);
          for(int k = 0; k < third_layer_hits.size(); k++)
          {
              TVector3 v3(third_layer_hits[k].X(), third_layer_hits[k].Y(), third_layer_hits[k].Z());
              TVector3 v3pv(third_layer_hits[k].X()-genX, third_layer_hits[k].Y()-genY, third_layer_hits[k].Z()-genZ);
              v13 = v3 - v1;
              float dphi1 = deltaPhi(v13.Phi(), v1pv.Phi());
              for(int l = 0; l < fourth_layer_hits.size(); l++)
              {
                  TVector3 v4(fourth_layer_hits[l].X(), fourth_layer_hits[l].Y(), fourth_layer_hits[l].Z());
                  TVector3 v4pv(fourth_layer_hits[l].X()-genX, fourth_layer_hits[l].Y()-genY, fourth_layer_hits[l].Z()-genZ);
                  v14 = v4 - v1;
                  v34 = v4 - v3;
                  float dphi2 = deltaPhi(v34.Phi(), v13.Phi());
                  b1->Fill(v14.Eta()-v13.Eta()); 
                  b3->Fill(v34.Eta()-v13.Eta());  
                  c2->Fill(v34.Eta()-v14.Eta());
                  f6->Fill(v4pv.Eta()-v3pv.Eta());
                  g3->Fill(dphi2-dphi1);
              }
          }
      }

      // 234 combi
      for(int j = 0; j < second_layer_hits.size(); j++)
      {
          TVector3 v2(second_layer_hits[j].X(), second_layer_hits[j].Y(), second_layer_hits[j].Z());
          TVector3 v2pv(second_layer_hits[j].X()-genX, second_layer_hits[j].Y()-genY, second_layer_hits[j].Z()-genZ);
          for(int k = 0; k < third_layer_hits.size(); k++)
          {
              TVector3 v3(third_layer_hits[k].X(), third_layer_hits[k].Y(), third_layer_hits[k].Z());
              v23 = v3 - v2;
              float dphi1 = deltaPhi(v23.Phi(), v2pv.Phi());
              for(int l = 0; l < fourth_layer_hits.size(); l++)
              {
                  TVector3 v4(fourth_layer_hits[l].X(), fourth_layer_hits[l].Y(), fourth_layer_hits[l].Z());
                  v24 = v4 - v2;
                  v34 = v4 - v3;
                  float dphi2 = deltaPhi(v34.Phi(), v23.Phi());
                  d1->Fill(v24.Eta()-v23.Eta()); 
                  d2->Fill(v34.Eta()-v23.Eta());  
                  e1->Fill(v34.Eta()-v24.Eta());
                  g4->Fill(dphi2-dphi1);
              }
          }
      }
   
      

   } // event loop


   TFile *output = new TFile("region1.root","recreate");
   h1->Write();
   h2->Write();
   h3->Write();
   h4->Write();
   h5->Write();
   h6->Write();
   h1z->Write();
   h2z->Write();
   h3z->Write();
   h4z->Write();
   a1->Write();
   a2->Write();
   a3->Write();
   a4->Write();
   b1->Write();
   b2->Write();
   b3->Write();
   c1->Write();
   c2->Write();
   d1->Write();
   d2->Write();
   e1->Write();
   f1->Write();
   f2->Write();
   f3->Write(); 
   f4->Write();
   f5->Write(); 
   f6->Write(); 
   g1->Write();
   g2->Write();
   g3->Write();
   g4->Write();

   output->Close();
}



