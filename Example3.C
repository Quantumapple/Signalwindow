/*
This macro shows how to access the particle-level reference for reconstructed objects.
It is also shown how to loop over the jet constituents.

root -l examples/Example3.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

//------------------------------------------------------------------------------

struct TestPlots
{
    TH1 *fL1XY;
};

//------------------------------------------------------------------------------

// Check ExRootResult.h
//void BookHistograms(ExRootResult *result, TestPlots *plots)
void BookHistograms(ExRootResult *result)
{
  TLegend *legend;
  TPaveText *comment;

  /*
  plots->fL1XY = result->AddHist2D(
          "1stBarrelXY", "1st pixel layer in xy plane", 
          "X [mm]", "Y [mm]", 
          400, -40, 40,
          400, -40, 40,
          0, 0
          );
  */
}

//------------------------------------------------------------------------------

int GenRandom()
{
    srand((unsigned int)time(0));

    int nRandom = rand() % 8 + 1;

    return nRandom;
}



//void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
void AnalyseEvents(ExRootTreeReader *treeReader, TFile *result, TTree *tt)
{
    int bRecHitN = 0;
    std::vector<int> bRecHitLayer;
    std::vector<float> bRecHitGx; 
    std::vector<float> bRecHitGy; 
    std::vector<float> bRecHitGz;
    
    int fRecHitN = 0;
    std::vector<int> fRecHitDisk;
    std::vector<float> fRecHitGx; 
    std::vector<float> fRecHitGy; 
    std::vector<float> fRecHitGz; 

    std::vector<float> propgenElPartPhi;
    std::vector<float> propgenElPartEta;
    std::vector<float> propgenElPartPt;

    int EgN = 0;
    std::vector<float> egCrysClusterEt;
    std::vector<float> egCrysClusterEta;
    std::vector<float> egCrysClusterPhi;
    std::vector<float> egCrysClusterGx;
    std::vector<float> egCrysClusterGy;
    std::vector<float> egCrysClusterGz;
    
    tt->Branch("propgenElPartPhi", &propgenElPartPhi); 
    tt->Branch("propgenElPartEta", &propgenElPartEta); 
    tt->Branch("propgenElPartPt", &propgenElPartPt); 
    
    tt->Branch("bRecHitN", &bRecHitN, "bRecHitN/I");
    tt->Branch("bRecHitLayer", &bRecHitLayer);
    tt->Branch("bRecHitGx", &bRecHitGx);
    tt->Branch("bRecHitGy", &bRecHitGy);
    tt->Branch("bRecHitGz", &bRecHitGz);
    
    tt->Branch("fRecHitN", &fRecHitN, "fRecHitN/I");
    tt->Branch("fRecHitDisk", &fRecHitDisk);
    tt->Branch("fRecHitGx", &fRecHitGx);
    tt->Branch("fRecHitGy", &fRecHitGy);
    tt->Branch("fRecHitGz", &fRecHitGz);

    tt->Branch("EgN", &EgN, "EgN/I");
    tt->Branch("egCrysClusterEt", &egCrysClusterEt);
    tt->Branch("egCrysClusterEta", &egCrysClusterEta);
    tt->Branch("egCrysClusterPhi", &egCrysClusterPhi);
    tt->Branch("egCrysClusterGx",&egCrysClusterGx);
    tt->Branch("egCrysClusterGy",&egCrysClusterGy);
    tt->Branch("egCrysClusterGz",&egCrysClusterGz);

    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchTower = treeReader->UseBranch("ECalTower");
    TClonesArray *branchTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchPixelL1 = treeReader->UseBranch("TrackPixel1");
    TClonesArray *branchPixelL2 = treeReader->UseBranch("TrackPixel2");
    TClonesArray *branchPixelL3 = treeReader->UseBranch("TrackPixel3");
    TClonesArray *branchPixelL4 = treeReader->UseBranch("TrackPixel4");
    TClonesArray *branchPixelD1 = treeReader->UseBranch("DiskPixel1");
    TClonesArray *branchPixelD2 = treeReader->UseBranch("DiskPixel2");
    TClonesArray *branchPixelD3 = treeReader->UseBranch("DiskPixel3");
    TClonesArray *branchPixelD4 = treeReader->UseBranch("DiskPixel4");
    TClonesArray *branchPixelD5 = treeReader->UseBranch("DiskPixel5");

    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl;

    GenParticle *particle;

    Track *trackL1;
    Track *trackL2;
    Track *trackL3;
    Track *trackL4;
    
    Track *trackD1;
    Track *trackD2;
    Track *trackD3;
    Track *trackD4;
    Track *trackD5;
    
    Tower *tower;
    Track *eCalTrack;

    TLorentzVector egVector;

    Long64_t entry;

    bool kFlag = false;

    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {
        bRecHitN = 0;
        bRecHitLayer.clear();
        bRecHitGx.clear();
        bRecHitGy.clear();
        bRecHitGz.clear();

        fRecHitN = 0;
        fRecHitDisk.clear();
        fRecHitGx.clear(); 
        fRecHitGy.clear(); 
        fRecHitGz.clear(); 

        propgenElPartPhi.clear();
        propgenElPartEta.clear();
        propgenElPartPt.clear();

        EgN = 0;
        egCrysClusterEta.clear();
        egCrysClusterPhi.clear();
        egCrysClusterEt.clear();
        egCrysClusterGx.clear();
        egCrysClusterGy.clear();
        egCrysClusterGz.clear();

        //cout << "Event loop start" << endl;
        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);
        kFlag = false;

        // Loop over Genparticle
        for(unsigned int k = 0; k < branchParticle->GetEntriesFast(); k++)
        {
            particle = (GenParticle*) branchParticle->At(k);
            propgenElPartPhi.push_back(particle->Phi);
            propgenElPartEta.push_back(particle->Eta);
            propgenElPartPt.push_back(particle->PT);
        }

        int TowerN = branchTower->GetEntriesFast();
        for(unsigned int k = 0; k < TowerN; k++)
        {
            eCalTower = (Tower*) branchTower->At(k);
            double eCalX, eCalY, eCalZ, eCalEt, eCalEta, eCalPhi;
            double j11 = 0, k11, l11;

            //ECal tower phi
            eCalPhi = eCalTower->Phi;

            //ECal tower eta
            eCalEta = eCalTower->Eta;

            //ECal tower Et
            eCalEt = eCalTower->ET;

            // 1.29m = 129cm
            //ECal position (unit: cm)
            eCalX = 129*cos(eCalPhi);   // x = rcos(phi)
            eCalY = 129*sin(eCalPhi);   // y = rsin(phi)
            double theta = 2.*atan(exp(-eCalEta));
            eCalZ = 129*cos(theta);     // z = rcos(theta)

            egCrysClusterEta.push_back(eCalEta);
            egCrysClusterPhi.push_back(eCalPhi);
            egCrysClusterEt.push_back(eCalEt);
            egCrysClusterGx.push_back(eCalX);
            egCrysClusterGy.push_back(eCalY);
            egCrysClusterGz.push_back(eCalZ);
        }
    
        EgN = egCrysClusterEta.size();

        // Loop over 1st pixel layer
        for(unsigned int k = 0; k < branchPixelL1->GetEntriesFast(); k++)
        {
            trackL1 = (Track*) branchPixelL1->At(k);
            if( fabs(trackL1->ZOuter) < 200 )
            {
                if( sqrt(pow(trackL1->XOuter,2)+pow(trackL1->YOuter,2)) < 29.9 ) continue; 
                bRecHitLayer.push_back(1);
                bRecHitGx.push_back(trackL1->XOuter/10.);
                bRecHitGy.push_back(trackL1->YOuter/10.);
                bRecHitGz.push_back(trackL1->ZOuter/10.);
            }
        }

        // Loop over 2nd pixel layer
        for(unsigned int k = 0; k < branchPixelL2->GetEntriesFast(); k++)
        {
            trackL2 = (Track*) branchPixelL2->At(k);
            if( fabs(trackL2->ZOuter) < 200 )
            {
                if( sqrt(pow(trackL2->XOuter,2)+pow(trackL2->YOuter,2)) < 69.9 ) continue; 
                bRecHitLayer.push_back(2);
                bRecHitGx.push_back(trackL2->XOuter/10.);
                bRecHitGy.push_back(trackL2->YOuter/10.);
                bRecHitGz.push_back(trackL2->ZOuter/10.);
            }
        }

        // Loop over 3rd pixel layer
        for(unsigned int k = 0; k < branchPixelL3->GetEntriesFast(); k++)
        {
            trackL3 = (Track*) branchPixelL3->At(k);
            if( fabs(trackL3->ZOuter) < 200 )
            {
                if( sqrt(pow(trackL3->XOuter,2)+pow(trackL3->YOuter,2)) < 117.9 ) continue; 
                bRecHitLayer.push_back(3);
                bRecHitGx.push_back(trackL3->XOuter/10.);
                bRecHitGy.push_back(trackL3->YOuter/10.);
                bRecHitGz.push_back(trackL3->ZOuter/10.);
            }
        }

        // Loop over 4th pixel layer
        for(unsigned int k = 0; k < branchPixelL4->GetEntriesFast(); k++)
        {
            trackL4 = (Track*) branchPixelL4->At(k);
            if( fabs(trackL4->ZOuter) < 200 )
            {
                if( sqrt(pow(trackL4->XOuter,2)+pow(trackL4->YOuter,2)) < 157.9 ) continue; 
                bRecHitLayer.push_back(4);
                bRecHitGx.push_back(trackL4->XOuter/10.);
                bRecHitGy.push_back(trackL4->YOuter/10.);
                bRecHitGz.push_back(trackL4->ZOuter/10.);
            }
        }
        bRecHitN = bRecHitGx.size();

        // Loop over 1st pixel Disk
        for(unsigned int k = 0; k < branchPixelD1->GetEntriesFast(); k++)
        {
            Int_t random1 = GenRandom();
            trackD1 = (Track*) branchPixelD1->At(k);
            float R = sqrt(pow(trackD1->XOuter,2)+pow(trackD1->YOuter,2));
            if( fabs(trackD1->ZOuter) > 248 && R > 28.0 )
                //if( fabs(trackD1->ZOuter) > 248 )
            {
                fRecHitDisk.push_back(1);
                fRecHitGx.push_back(trackD1->XOuter/10.);
                fRecHitGy.push_back(trackD1->YOuter/10.);

                fRecHitGz.push_back(trackD1->ZOuter/10.);
            }
        }

        // Loop over 2nd pixel Disk
        for(unsigned int k = 0; k < branchPixelD2->GetEntriesFast(); k++)
        {
            trackD2 = (Track*) branchPixelD2->At(k);
            float R = sqrt(pow(trackD2->XOuter,2)+pow(trackD2->YOuter,2));
            if( fabs(trackD2->ZOuter) > 315 && R > 28.0 )
            {
                fRecHitDisk.push_back(2);
                fRecHitGx.push_back(trackD2->XOuter/10.);
                fRecHitGy.push_back(trackD2->YOuter/10.);
                fRecHitGz.push_back(trackD2->ZOuter/10.);
            }
        }

        // Loop over 3rd pixel Disk
        for(unsigned int k = 0; k < branchPixelD3->GetEntriesFast(); k++)
        {
            trackD3 = (Track*) branchPixelD3->At(k);
            float R = sqrt(pow(trackD3->XOuter,2)+pow(trackD3->YOuter,2));
            if( fabs(trackD3->ZOuter) > 407 && R > 28.0 )
            {
                fRecHitDisk.push_back(3);
                fRecHitGx.push_back(trackD3->XOuter/10.);
                fRecHitGy.push_back(trackD3->YOuter/10.);
                fRecHitGz.push_back(trackD3->ZOuter/10.);
            }
        }

        // Loop over 4th pixel Disk
        for(unsigned int k = 0; k < branchPixelD4->GetEntriesFast(); k++)
        {
            trackD4 = (Track*) branchPixelD4->At(k);
            float R = sqrt(pow(trackD4->XOuter,2)+pow(trackD4->YOuter,2));
            if( fabs(trackD4->ZOuter) > 521 && R > 28.0 )
            {
                fRecHitDisk.push_back(4);
                fRecHitGx.push_back(trackD4->XOuter/10.);
                fRecHitGy.push_back(trackD4->YOuter/10.);
                fRecHitGz.push_back(trackD4->ZOuter/10.);
            }
        }

        // Loop over 5th pixel Disk
        for(unsigned int k = 0; k < branchPixelD5->GetEntriesFast(); k++)
        {
            trackD5 = (Track*) branchPixelD5->At(k);
            float R = sqrt(pow(trackD5->XOuter,2)+pow(trackD5->YOuter,2));
            if( fabs(trackD5->ZOuter) > 667 && R > 28.0 )
            {
                fRecHitDisk.push_back(5);
                fRecHitGx.push_back(trackD5->XOuter/10.);
                fRecHitGy.push_back(trackD5->YOuter/10.);
                fRecHitGz.push_back(trackD5->ZOuter/10.);
            }
        }

        fRecHitN = fRecHitGx.size();


        tt->Fill();
    } // event loop
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
    result->Print("png");
}

//------------------------------------------------------------------------------

void Example3(const char *inputFile)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  //ExRootResult *result = new ExRootResult();
  
  //TFile *result = new TFile("results.root","RECREATE");
  //TFile *result = new TFile("MinBiasResults.root","RECREATE");
  TFile *result = new TFile("temp.root","RECREATE");
  result->mkdir("l1PiXTRKTree");
  result->cd("l1PiXTRKTree");
  
  TTree *tree = new TTree("L1PiXTRKTree","L1PiXTRKTree");

  AnalyseEvents(treeReader, result, tree);
  result->Write();

  cout << "** Exiting..." << endl;

  //delete plots;
  //delete tree;
  //delete result;
  delete treeReader;
  delete chain;

  result->Close();
}

//------------------------------------------------------------------------------
