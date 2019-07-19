// C++ libraries
#include <iostream>
#include <algorithm> //std::random_shuffle
#include <vector>
#include <ctime>
#include <cstdlib>

using namespace std;

// ROOT libraries
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TFrame.h"
#include "TSystem.h"

//------------------------------------------------------------------------------

void AnalyseEvents(TTree *fTree1, TTree *fTree2, TFile *result, TTree *tt)
{

    TBranch *b_sigGenPt;
    TBranch *b_sigGenPhi;
    TBranch *b_sigGenEta;

    TBranch *b_sigEgEt;
    TBranch *b_sigEgEta;
    TBranch *b_sigEgPhi;
    TBranch *b_sigEgGx;
    TBranch *b_sigEgGy;
    TBranch *b_sigEgGz;

    TBranch *b_sigbLa;
    TBranch *b_sigbGx;
    TBranch *b_sigbGy;
    TBranch *b_sigbGz;

    TBranch *b_sigfDi;
    TBranch *b_sigfGx;
    TBranch *b_sigfGy;
    TBranch *b_sigfGz;

    std::vector<float> *sigGenPt  = 0;
    std::vector<float> *sigGenPhi = 0;
    std::vector<float> *sigGenEta = 0;

    std::vector<float> *sigEgEt  = 0;
    std::vector<float> *sigEgEta = 0;
    std::vector<float> *sigEgPhi = 0;
    std::vector<float> *sigEgGx  = 0;
    std::vector<float> *sigEgGy  = 0;
    std::vector<float> *sigEgGz  = 0;

    std::vector<int> *sigbLa   = 0;
    std::vector<float> *sigbGx = 0;
    std::vector<float> *sigbGy = 0;
    std::vector<float> *sigbGz = 0;

    std::vector<int> *sigfDi   = 0;
    std::vector<float> *sigfGx = 0;
    std::vector<float> *sigfGy = 0;
    std::vector<float> *sigfGz = 0;

    fTree1->SetBranchAddress("propgenElPartPt", &sigGenPt, &b_sigGenPt);
    fTree1->SetBranchAddress("propgenElPartPhi", &sigGenPhi, &b_sigGenPhi);
    fTree1->SetBranchAddress("propgenElPartEta", &sigGenEta, &b_sigGenEta);

    fTree1->SetBranchAddress("egCrysClusterEt", &sigEgEt, &b_sigEgEt);
    fTree1->SetBranchAddress("egCrysClusterEta", &sigEgEta, &b_sigEgEta);
    fTree1->SetBranchAddress("egCrysClusterPhi", &sigEgPhi, &b_sigEgPhi);
    fTree1->SetBranchAddress("egCrysClusterGx", &sigEgGx, &b_sigEgGx);
    fTree1->SetBranchAddress("egCrysClusterGy", &sigEgGy, &b_sigEgGy);
    fTree1->SetBranchAddress("egCrysClusterGz", &sigEgGz, &b_sigEgGz);

    fTree1->SetBranchAddress("bRecHitLayer", &sigbLa, &b_sigbLa);
    fTree1->SetBranchAddress("bRecHitGx", &sigbGx, &b_sigbGx);
    fTree1->SetBranchAddress("bRecHitGy", &sigbGy, &b_sigbGy);
    fTree1->SetBranchAddress("bRecHitGz", &sigbGz, &b_sigbGz);

    fTree1->SetBranchAddress("fRecHitDisk", &sigfDi, &b_sigfDi);
    fTree1->SetBranchAddress("fRecHitGx", &sigfGx, &b_sigfGx);
    fTree1->SetBranchAddress("fRecHitGy", &sigfGy, &b_sigfGy);
    fTree1->SetBranchAddress("fRecHitGz", &sigfGz, &b_sigfGz);

    TBranch *b_puEgEt;
    TBranch *b_puEgEta;
    TBranch *b_puEgPhi;
    TBranch *b_puEgGx;
    TBranch *b_puEgGy;
    TBranch *b_puEgGz;

    TBranch *b_pubLa;
    TBranch *b_pubGx;
    TBranch *b_pubGy;
    TBranch *b_pubGz;

    TBranch *b_pufDi;
    TBranch *b_pufGx;
    TBranch *b_pufGy;
    TBranch *b_pufGz;

    std::vector<float> *puEgEt  = 0;
    std::vector<float> *puEgEta = 0;
    std::vector<float> *puEgPhi = 0;
    std::vector<float> *puEgGx  = 0;
    std::vector<float> *puEgGy  = 0;
    std::vector<float> *puEgGz  = 0;

    std::vector<int> *pubLa   = 0;
    std::vector<float> *pubGx = 0;
    std::vector<float> *pubGy = 0;
    std::vector<float> *pubGz = 0;

    std::vector<int> *pufDi   = 0;
    std::vector<float> *pufGx = 0;
    std::vector<float> *pufGy = 0;
    std::vector<float> *pufGz = 0;

    fTree2->SetBranchAddress("egCrysClusterEt", &puEgEt, &b_puEgEt);
    fTree2->SetBranchAddress("egCrysClusterEta", &puEgEta, &b_puEgEta);
    fTree2->SetBranchAddress("egCrysClusterPhi", &puEgPhi, &b_puEgPhi);
    fTree2->SetBranchAddress("egCrysClusterGx", &puEgGx, &b_puEgGx);
    fTree2->SetBranchAddress("egCrysClusterGy", &puEgGy, &b_puEgGy);
    fTree2->SetBranchAddress("egCrysClusterGz", &puEgGz, &b_puEgGz);

    fTree2->SetBranchAddress("bRecHitLayer", &pubLa, &b_pubLa);
    fTree2->SetBranchAddress("bRecHitGx", &pubGx, &b_pubGx);
    fTree2->SetBranchAddress("bRecHitGy", &pubGy, &b_pubGy);
    fTree2->SetBranchAddress("bRecHitGz", &pubGz, &b_pubGz);

    fTree2->SetBranchAddress("fRecHitDisk", &pufDi, &b_pufDi);
    fTree2->SetBranchAddress("fRecHitGx", &pufGx, &b_pufGx);
    fTree2->SetBranchAddress("fRecHitGy", &pufGy, &b_pufGy);
    fTree2->SetBranchAddress("fRecHitGz", &pufGz, &b_pufGz);

    bool debug = false;

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

    tt->Branch("egCrysClusterEt", &egCrysClusterEt);
    tt->Branch("egCrysClusterEta", &egCrysClusterEta);
    tt->Branch("egCrysClusterPhi", &egCrysClusterPhi);
    tt->Branch("egCrysClusterGx",&egCrysClusterGx);
    tt->Branch("egCrysClusterGy",&egCrysClusterGy);
    tt->Branch("egCrysClusterGz",&egCrysClusterGz);

    Long64_t allEntries1 = fTree1->GetEntries();
    Long64_t allEntries2 = fTree2->GetEntries();

    if( debug ) cout << "Number of events in signal: " << allEntries1 << ", in QCD: " << allEntries2 << endl;
    cout << "** Chain contains " << allEntries1 << " events" << endl;

    // Loop events
    for(unsigned int i = 0; i < allEntries1; i++)
    {
        // vector clearing
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

        egCrysClusterEta.clear();
        egCrysClusterPhi.clear();
        egCrysClusterEt.clear();
        egCrysClusterGx.clear();
        egCrysClusterGy.clear();
        egCrysClusterGz.clear();

        // Load trees and branches
        Long64_t jentry1 = fTree1->LoadTree(i);
        Long64_t jentry2 = fTree2->LoadTree(i);

        // single electron 
        b_sigGenPt->GetEntry(jentry1);
        b_sigGenPhi->GetEntry(jentry1);
        b_sigGenEta->GetEntry(jentry1);

        b_sigEgEt->GetEntry(jentry1);
        b_sigEgEta->GetEntry(jentry1);
        b_sigEgPhi->GetEntry(jentry1);
        b_sigEgGx->GetEntry(jentry1);
        b_sigEgGy->GetEntry(jentry1);
        b_sigEgGz->GetEntry(jentry1);

        b_sigbLa->GetEntry(jentry1);
        b_sigbGx->GetEntry(jentry1);
        b_sigbGy->GetEntry(jentry1);
        b_sigbGz->GetEntry(jentry1);

        b_sigfDi->GetEntry(jentry1);
        b_sigfGx->GetEntry(jentry1);
        b_sigfGy->GetEntry(jentry1);
        b_sigfGz->GetEntry(jentry1);

        // QCD 
        b_puEgEt->GetEntry(jentry2);
        b_puEgEta->GetEntry(jentry2);
        b_puEgPhi->GetEntry(jentry2);
        b_puEgGx->GetEntry(jentry2);
        b_puEgGy->GetEntry(jentry2);
        b_puEgGz->GetEntry(jentry2);

        b_pubLa->GetEntry(jentry2);
        b_pubGx->GetEntry(jentry2);
        b_pubGy->GetEntry(jentry2);
        b_pubGz->GetEntry(jentry2);

        b_pufDi->GetEntry(jentry2);
        b_pufGx->GetEntry(jentry2);
        b_pufGy->GetEntry(jentry2);
        b_pufGz->GetEntry(jentry2);

        cout << "Event: " << i+1 << endl;
        //for(unsigned int j = 0; j < sigEgPhi->size(); j++)
        //{
        //    cout << "Signal: " << sigEgPhi->at(j) << ' ';
        //}
        //cout << endl;
        //for(unsigned int j = 0; j < puEgPhi->size(); j++)
        //{
        //    cout << "QCD: " << puEgPhi->at(j) << ' ';
        //}
        //cout << endl

        // Loop over Genparticle
        for(unsigned int k = 0; k < sigGenPt->size(); k++)
        {
            propgenElPartPhi.push_back(sigGenPhi->at(k));
            propgenElPartEta.push_back(sigGenEta->at(k));
            propgenElPartPt.push_back(sigGenPt->at(k));
        }

        // Loop over L1 egamma
        int EgN = sigEgEt->size();
        for(unsigned int k = 0; k < EgN; k++)
        {
            egCrysClusterEta.push_back(sigEgEta->at(k));
            egCrysClusterPhi.push_back(sigEgPhi->at(k));
            egCrysClusterEt.push_back(sigEgEt->at(k)); 
            egCrysClusterGx.push_back(sigEgGx->at(k)); 
            egCrysClusterGy.push_back(sigEgGy->at(k)); 
            egCrysClusterGz.push_back(sigEgGz->at(k)); 
        }

        EgN = puEgEt->size();
        for(unsigned int k = 0; k < EgN; k++)
        {
            egCrysClusterEta.push_back(puEgEta->at(k));
            egCrysClusterPhi.push_back(puEgPhi->at(k));
            egCrysClusterEt.push_back(puEgEt->at(k)); 
            egCrysClusterGx.push_back(puEgGx->at(k)); 
            egCrysClusterGy.push_back(puEgGy->at(k)); 
            egCrysClusterGz.push_back(puEgGz->at(k)); 
        }

        // Loop over pixel barrels
        int bRecHitN = sigbLa->size();
        for(unsigned int k = 0; k < bRecHitN; k++)
        {
            bRecHitLayer.push_back(sigbLa->at(k));
            bRecHitGx.push_back(sigbGx->at(k));
            bRecHitGy.push_back(sigbGy->at(k));
            bRecHitGz.push_back(sigbGz->at(k));
        }
        bRecHitN = pubLa->size();
        for(unsigned int k = 0; k < bRecHitN; k++)
        {
            bRecHitLayer.push_back(pubLa->at(k));
            bRecHitGx.push_back(pubGx->at(k));
            bRecHitGy.push_back(pubGy->at(k));
            bRecHitGz.push_back(pubGz->at(k));
        }

        // Loop over pixel disks
        int fRecHitN = sigfDi->size();
        for(unsigned int k = 0; k < fRecHitN; k++)
        {
            fRecHitDisk.push_back(sigfDi->at(k));
            fRecHitGx.push_back(sigfGx->at(k));
            fRecHitGy.push_back(sigfGy->at(k));
            fRecHitGz.push_back(sigfGz->at(k));
        }
        fRecHitN = pufDi->size();
        for(unsigned int k = 0; k < fRecHitN; k++)
        {
            fRecHitDisk.push_back(pufDi->at(k));
            fRecHitGx.push_back(pufGx->at(k));
            fRecHitGy.push_back(pufGy->at(k));
            fRecHitGz.push_back(pufGz->at(k));
        }
        fRecHitN = fRecHitGx.size();

        //cout << endl << "Before mixing" << endl;
        //for(std::vector<float>::iterator it = egCrysClusterPhi.begin(); it != egCrysClusterPhi.end(); ++it)
        //{
        //    std::cout << ' ' << *it;
        //}
        //cout << endl;

        // random shuffle data in vector
        // use built-in random generator:
        std::random_shuffle( egCrysClusterPhi.begin(), egCrysClusterPhi.end() );
        std::random_shuffle( egCrysClusterEta.begin(), egCrysClusterEta.end() );
        std::random_shuffle( egCrysClusterEt.begin(), egCrysClusterEt.end() );
        std::random_shuffle( egCrysClusterGx.begin(), egCrysClusterGx.end() );
        std::random_shuffle( egCrysClusterGy.begin(), egCrysClusterGy.end() );
        std::random_shuffle( egCrysClusterGz.begin(), egCrysClusterGz.end() );

        //cout << endl << "After mixing" << endl;
        //for(std::vector<float>::iterator it = egCrysClusterPhi.begin(); it != egCrysClusterPhi.end(); ++it)
        //{
        //    std::cout << ' ' << *it;
        //}
        //cout << endl;

        tt->Fill();

    } // event loop

    fTree1->ResetBranchAddresses();

}


//------------------------------------------------------------------------------

void mix()
{

    TFile *file1 = new TFile("SingleElNOPU.root","OPEN");
    TFile *file2 = new TFile("qcd200PU.root","OPEN");

    TTree *fTree1 = (TTree*)file1->Get("l1PiXTRKTree/L1PiXTRKTree");
    TTree *fTree2 = (TTree*)file2->Get("l1PiXTRKTree/L1PiXTRKTree");

    //file1->GetObject("l1PiXTRKTree/L1PiXTRKTree", fTree1);
    //file2->GetObject("l1PiXTRKTree/L1PiXTRKTree", fTree2);

    //TChain *chain1 = new TChain("l1PiXTRKTree/L1PiXTRKTree");
    //chain1->Add("SingleElNOPU.root");
    //TChain *chain2 = new TChain("l1PiXTRKTree/L1PiXTRKTree");
    //chain2->Add("qcd200PU.root");

    TFile *output = new TFile("temp.root","RECREATE");
    output->mkdir("l1PiXTRKTree");
    output->cd("l1PiXTRKTree");

    TTree *tt = new TTree("L1PiXTRKTree","L1PiXTRKTree");

    auto nevent1 = fTree1->GetEntries();
    auto nevent2 = fTree2->GetEntries();
    cout << nevent1 << ", " << nevent2 << endl;

    //Long64_t nentries1 = fTree1->GetEntries();
    //Long64_t nentries2 = fTree2->GetEntries();

    if( nevent1 != nevent2 ) {
        cout << "You should match number entries between two files" << endl;
        exit(0);
    }

    AnalyseEvents(fTree1, fTree2, output, tt);
    output->Write();
}
