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

//---------------------  Function definition  ----------------------------------

float ratiofunc(int type, int count);
void AnalyseEvents(TTree *fTree1, TTree *fTree2, TFile *result, TTree *tt);

//------------------------------------------------------------------------------

void mix13_r56()
{

    TFile *file1 = new TFile("input.root","OPEN");
    TFile *file2 = new TFile("../SE_NoPU.root","OPEN");

    TTree *fTree1 = (TTree*)file1->Get("l1PiXTRKTree/L1PiXTRKTree");
    TTree *fTree2 = (TTree*)file2->Get("l1PiXTRKTree/L1PiXTRKTree");

    TFile *output = new TFile("merged_0.root","RECREATE");
    output->mkdir("l1PiXTRKTree");
    output->cd("l1PiXTRKTree");

    TTree *tt = new TTree("L1PiXTRKTree","L1PiXTRKTree");

    auto nevent1 = fTree1->GetEntries();
    auto nevent2 = fTree2->GetEntries();
    cout << nevent1 << ", " << nevent2 << endl;

    //Long64_t nentries1 = fTree1->GetEntries();
    //Long64_t nentries2 = fTree2->GetEntries();

    //if( nevent1 != nevent2 ) {
    //    cout << "You should match number entries between two files" << endl;
    //    exit(0);
    //}

    AnalyseEvents(fTree1, fTree2, output, tt);
    output->Write();
}

//------------------------------------------------------------------------------

float ratiofunc(int type, int count)
{
    double out = 0.;

    TF1 *f1 = new TF1("f1","[0]*exp([1]*pow(x,[2]))+[3]+[4]*x+[5]*x*x+[6]*pow(x,3)+[7]*pow(x,4)",0.003,1.);
    TF1 *f2 = new TF1("f2","[0]*exp([1]*pow(x,[2]))+[3]+[4]*x+[5]*x*x+[6]*pow(x,3)+[7]*pow(x,4)+[8]*pow(x,5)+[9]*pow(x,6)",0.0022,0.95);
    TF1 *f3 = new TF1("f3","[0]*exp([1]*pow(x,[2]))+[3]+[4]*x+[5]*x*x",0.0013,1.);
    TF1 *f4 = new TF1("f4","[0]*exp([1]*pow(x,[2]))+[3]+[4]*x+[5]*x*x+[6]*pow(x,3)",0.001,1.);
    TF1 *f5 = new TF1("f5","[0]*exp([1]*pow(x,[2]))+[3]+[4]*x+[5]*x*x+[6]*pow(x,3)+[7]*pow(x,4)+[8]*pow(x,5)",0.003,1.);
    TF1 *f6 = new TF1("f6","[0]*exp([1]*pow(x,[2]))+[3]+[4]*x+[5]*x*x+[6]*pow(x,3)+[7]*pow(x,4)+[8]*pow(x,5)",0.0035,0.97);
    
    f1->SetParameters(0.0678573,0.609147,-0.294585,-0.228479,0.386804,-0.676396,0.607686,-0.214137);
    f2->SetParameters(0.235596,1.98923,-0.0388373,-2.22096,2.40337,-7.17586,15.467,-21.9822,17.8493,-6.11948);
    f3->SetParameters(0.0119969, 0.286603, -0.470787, 0.00172564, 0.0233817, -0.0255747); 
    f4->SetParameters(0.842313,-27.2558,0.779984,0.0101008,-0.0328682,0.0618374,-0.029615);
    f5->SetParameters(0.00846556,3.84928,-0.0322101,-0.533008,0.34178,-0.289075,-0.0497208,0.330497,-0.187506,0);
    f6->SetParameters(0.000331564,3.7808,-0.149153,-0.0571487,0.144501,-0.207412,0.172996,-0.0825695,0.0145481,0);

    if( type == 1 && count == 0 ){
        out = fabs(f1->GetRandom(0.003, 0.32));
    }
    if( type == 1 && count == 1 ){
        out = fabs(f2->GetRandom(0.0022, 0.32));
    }
    if( type == 1 && count > 1 ){
        out = fabs(f3->GetRandom(0.0013, 0.22));
    }
    if( type == 2 && count == 0 ){
        out = fabs(f4->GetRandom());
    }
    if( type == 2 && count == 1 ){
        out = fabs(f5->GetRandom());
    }
    if( type == 2 && count > 1 ){
        out = fabs(f6->GetRandom());
    }
    
    return out;
}


void AnalyseEvents(TTree *fTree1, TTree *fTree2, TFile *result, TTree *tt)
{
    gRandom->SetSeed(0);

    TBranch *b_puGenPt;
    TBranch *b_puGenPhi;
    TBranch *b_puGenEta;
    TBranch *b_puGenX;
    TBranch *b_puGenY;
    TBranch *b_puGenZ;

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

    std::vector<float> *puGenPt  = 0;
    std::vector<float> *puGenPhi = 0;
    std::vector<float> *puGenEta = 0;
    std::vector<float> *puGenX = 0;
    std::vector<float> *puGenY = 0;
    std::vector<float> *puGenZ = 0;

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

    fTree1->SetBranchAddress("propgenElPartPt", &puGenPt, &b_puGenPt);
    fTree1->SetBranchAddress("propgenElPartPhi", &puGenPhi, &b_puGenPhi);
    fTree1->SetBranchAddress("propgenElPartEta", &puGenEta, &b_puGenEta);
    fTree1->SetBranchAddress("propgenElPartX", &puGenX, &b_puGenX);
    fTree1->SetBranchAddress("propgenElPartY", &puGenY, &b_puGenY);
    fTree1->SetBranchAddress("propgenElPartZ", &puGenZ, &b_puGenZ);

    fTree1->SetBranchAddress("egCrysClusterEt", &puEgEt, &b_puEgEt);
    fTree1->SetBranchAddress("egCrysClusterEta", &puEgEta, &b_puEgEta);
    fTree1->SetBranchAddress("egCrysClusterPhi", &puEgPhi, &b_puEgPhi);
    fTree1->SetBranchAddress("egCrysClusterGx", &puEgGx, &b_puEgGx);
    fTree1->SetBranchAddress("egCrysClusterGy", &puEgGy, &b_puEgGy);
    fTree1->SetBranchAddress("egCrysClusterGz", &puEgGz, &b_puEgGz);

    fTree1->SetBranchAddress("bRecHitLayer", &pubLa, &b_pubLa);
    fTree1->SetBranchAddress("bRecHitGx", &pubGx, &b_pubGx);
    fTree1->SetBranchAddress("bRecHitGy", &pubGy, &b_pubGy);
    fTree1->SetBranchAddress("bRecHitGz", &pubGz, &b_pubGz);

    fTree1->SetBranchAddress("fRecHitDisk", &pufDi, &b_pufDi);
    fTree1->SetBranchAddress("fRecHitGx", &pufGx, &b_pufGx);
    fTree1->SetBranchAddress("fRecHitGy", &pufGy, &b_pufGy);
    fTree1->SetBranchAddress("fRecHitGz", &pufGz, &b_pufGz);

    TBranch *b_sigGenPt;
    TBranch *b_sigGenPhi;
    TBranch *b_sigGenEta;
    TBranch *b_sigGenX;
    TBranch *b_sigGenY;
    TBranch *b_sigGenZ;

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
    std::vector<float> *sigGenX = 0;
    std::vector<float> *sigGenY = 0;
    std::vector<float> *sigGenZ = 0;

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

    fTree2->SetBranchAddress("propgenElPartPt", &sigGenPt, &b_sigGenPt);
    fTree2->SetBranchAddress("propgenElPartPhi", &sigGenPhi, &b_sigGenPhi);
    fTree2->SetBranchAddress("propgenElPartEta", &sigGenEta, &b_sigGenEta);
    fTree2->SetBranchAddress("propgenElPartX", &sigGenX, &b_sigGenX);
    fTree2->SetBranchAddress("propgenElPartY", &sigGenY, &b_sigGenY);
    fTree2->SetBranchAddress("propgenElPartZ", &sigGenZ, &b_sigGenZ);

    fTree2->SetBranchAddress("egCrysClusterEt", &sigEgEt, &b_sigEgEt);
    fTree2->SetBranchAddress("egCrysClusterEta", &sigEgEta, &b_sigEgEta);
    fTree2->SetBranchAddress("egCrysClusterPhi", &sigEgPhi, &b_sigEgPhi);
    fTree2->SetBranchAddress("egCrysClusterGx", &sigEgGx, &b_sigEgGx);
    fTree2->SetBranchAddress("egCrysClusterGy", &sigEgGy, &b_sigEgGy);
    fTree2->SetBranchAddress("egCrysClusterGz", &sigEgGz, &b_sigEgGz);

    fTree2->SetBranchAddress("bRecHitLayer", &sigbLa, &b_sigbLa);
    fTree2->SetBranchAddress("bRecHitGx", &sigbGx, &b_sigbGx);
    fTree2->SetBranchAddress("bRecHitGy", &sigbGy, &b_sigbGy);
    fTree2->SetBranchAddress("bRecHitGz", &sigbGz, &b_sigbGz);

    fTree2->SetBranchAddress("fRecHitDisk", &sigfDi, &b_sigfDi);
    fTree2->SetBranchAddress("fRecHitGx", &sigfGx, &b_sigfGx);
    fTree2->SetBranchAddress("fRecHitGy", &sigfGy, &b_sigfGy);
    fTree2->SetBranchAddress("fRecHitGz", &sigfGz, &b_sigfGz);


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
    std::vector<float> propgenElPartX;
    std::vector<float> propgenElPartY;
    std::vector<float> propgenElPartZ;

    std::vector<float> egCrysClusterEt;
    std::vector<float> egCrysClusterEta;
    std::vector<float> egCrysClusterPhi;
    std::vector<float> egCrysClusterGx;
    std::vector<float> egCrysClusterGy;
    std::vector<float> egCrysClusterGz;

    tt->Branch("propgenElPartPhi", &propgenElPartPhi); 
    tt->Branch("propgenElPartEta", &propgenElPartEta); 
    tt->Branch("propgenElPartPt", &propgenElPartPt); 
    tt->Branch("propgenElPartX", &propgenElPartX); 
    tt->Branch("propgenElPartY", &propgenElPartY); 
    tt->Branch("propgenElPartZ", &propgenElPartZ); 

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

    //int mergedCheck[6] = {0};
    int mergedCheck[100000] = {0};

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
        propgenElPartX.clear();
        propgenElPartY.clear();
        propgenElPartZ.clear();

        egCrysClusterEta.clear();
        egCrysClusterPhi.clear();
        egCrysClusterEt.clear();
        egCrysClusterGx.clear();
        egCrysClusterGy.clear();
        egCrysClusterGz.clear();

        // Load trees and branches
        Long64_t jentry1 = fTree1->LoadTree(i);
        Long64_t jentry2 = fTree2->LoadTree(i);

        // single electron 200PU 
        b_puGenPt->GetEntry(jentry1);
        b_puGenPhi->GetEntry(jentry1);
        b_puGenEta->GetEntry(jentry1);
        b_puGenX->GetEntry(jentry1);
        b_puGenY->GetEntry(jentry1);
        b_puGenZ->GetEntry(jentry1);

        b_puEgEt->GetEntry(jentry1);
        b_puEgEta->GetEntry(jentry1);
        b_puEgPhi->GetEntry(jentry1);
        b_puEgGx->GetEntry(jentry1);
        b_puEgGy->GetEntry(jentry1);
        b_puEgGz->GetEntry(jentry1);

        b_pubLa->GetEntry(jentry1);
        b_pubGx->GetEntry(jentry1);
        b_pubGy->GetEntry(jentry1);
        b_pubGz->GetEntry(jentry1);

        b_pufDi->GetEntry(jentry1);
        b_pufGx->GetEntry(jentry1);
        b_pufGy->GetEntry(jentry1);
        b_pufGz->GetEntry(jentry1);

        cout << endl;

        // Loop over Genparticle
        for(unsigned int k = 0; k < puGenPt->size(); k++)
        {
            propgenElPartPhi.push_back(puGenPhi->at(k));
            propgenElPartEta.push_back(puGenEta->at(k));
            propgenElPartPt.push_back(puGenPt->at(k));
            propgenElPartX.push_back(puGenX->at(k));
            propgenElPartY.push_back(puGenY->at(k));
            propgenElPartZ.push_back(puGenZ->at(k));
        }

        float dr_cut = 0.1;
        float closest_dr = 9999.;
        int closest_eg = 0;
        // Loop over L1 egamma
        int EgN = puEgEt->size();
        for(unsigned int k = 0; k < EgN; k++)
        {
            egCrysClusterEta.push_back(puEgEta->at(k));
            egCrysClusterPhi.push_back(puEgPhi->at(k));
            egCrysClusterEt.push_back(puEgEt->at(k)); 
            egCrysClusterGx.push_back(puEgGx->at(k)); 
            egCrysClusterGy.push_back(puEgGy->at(k)); 
            egCrysClusterGz.push_back(puEgGz->at(k));

            float dPhi = puGenPhi->at(0) - puEgPhi->at(k);
            if( dPhi >= float(M_PI)) dPhi -= float(2*M_PI);
            if( dPhi < -float(M_PI)) dPhi += float(2*M_PI);

            float current_dr = sqrt(pow(dPhi,2)+pow(puGenEta->at(0)-puEgEta->at(k),2));
            if( puEgEt->at(k) < 10 ) continue;
            if( current_dr < closest_dr ){
                closest_dr = current_dr;
                closest_eg = k;
            }
        }

        // Loop over pixel barrels
        int bN = pubLa->size();
        for(unsigned int k = 0; k < bN; k++)
        {
            bRecHitLayer.push_back(pubLa->at(k));
            bRecHitGx.push_back(pubGx->at(k));
            bRecHitGy.push_back(pubGy->at(k));
            bRecHitGz.push_back(pubGz->at(k));
        }

        // Loop over pixel disks
        int fN = pufDi->size();
        for(unsigned int k = 0; k < fN; k++)
        {
            fRecHitDisk.push_back(pufDi->at(k));
            fRecHitGx.push_back(pufGx->at(k));
            fRecHitGy.push_back(pufGy->at(k));
            fRecHitGz.push_back(pufGz->at(k));
        }

        float dR = 9999.;
        float dZ = 9999.;
        float dphi, deta;

        //int cut[6] = {1000, 37, 3, 2, 1, 0}; //cut4
        int N_merged = 0;
        int Type = 0;
        //int Nflag[5] = {0};
        //if( mergedCheck[1] < 35 ) Nflag[1] = 1;
        //if( mergedCheck[2] < 3 ) Nflag[2] = 1;
        //if( mergedCheck[3] < 2 ) Nflag[3] = 1;
        //if( mergedCheck[4] < 1 ) Nflag[4] = 1;

        if(closest_dr < dr_cut && closest_dr != 9999. ) {    
            if( fabs(puEgEta->at(closest_eg)) < 1.7 ) Type=1;
            if( fabs(puEgEta->at(closest_eg)) > 1.7 ) Type=2;
            float ratio = 0;
            float ratio0 = ratiofunc(Type, 0);
            float ratio1 = ratiofunc(Type, 1);
            float ratio2 = ratiofunc(Type, 2);

            float etaflag = fabs(puEgEta->at(closest_eg));
            if( etaflag > 2.1 )
            {
                for(unsigned int j = 0; j < allEntries2; j++) {
                    Long64_t jentry3 = fTree2->LoadTree(j);

                    //if( Type == 1 && N_merged == 1 ) break;
                    //if( N_merged > 1000 ) break;
                    //ratio = ratio0;
                    if( N_merged < 1 ) ratio = ratio0;
                    if( N_merged == 1 ) ratio = ratio1;
                    if( N_merged > 1 ) ratio = ratio2;

                    // single electron NoPU
                    b_sigGenPt->GetEntry(jentry3);
                    b_sigGenPhi->GetEntry(jentry3);
                    b_sigGenEta->GetEntry(jentry3);
                    b_sigGenX->GetEntry(jentry3);
                    b_sigGenY->GetEntry(jentry3);
                    b_sigGenZ->GetEntry(jentry3);

                    b_sigEgEt->GetEntry(jentry3);
                    b_sigEgEta->GetEntry(jentry3);
                    b_sigEgPhi->GetEntry(jentry3);
                    b_sigEgGx->GetEntry(jentry3);
                    b_sigEgGy->GetEntry(jentry3);
                    b_sigEgGz->GetEntry(jentry3);

                    b_sigbLa->GetEntry(jentry3);
                    b_sigbGx->GetEntry(jentry3);
                    b_sigbGy->GetEntry(jentry3);
                    b_sigbGz->GetEntry(jentry3);

                    b_sigfDi->GetEntry(jentry3);
                    b_sigfGx->GetEntry(jentry3);
                    b_sigfGy->GetEntry(jentry3);
                    b_sigfGz->GetEntry(jentry3);

                    dphi = puGenPhi->at(0) - sigGenPhi->at(0);
                    if( dphi >= float(M_PI)) dphi -= float(2*M_PI);
                    if( dphi < -float(M_PI)) dphi += float(2*M_PI);

                    deta = puGenEta->at(0) - sigGenEta->at(0);
                    dR = sqrt(pow(dphi,2)+pow(deta,2));
                    dZ = puGenZ->at(0) - sigGenZ->at(0);
                    float siggenpt = sigGenPt->at(0);
                    float pugenpt = puGenPt->at(0);

                    if( fabs(dR) < 0.3 && (siggenpt < pugenpt*ratio*1.20 && siggenpt > pugenpt*ratio*0.80) ) {
                        if( fabs(dZ) < 0.3 )
                        {
                            N_merged++;

                            for(unsigned int k = 0; k < sigGenPt->size(); k++)
                            {
                                propgenElPartPhi.push_back(sigGenPhi->at(k));
                                propgenElPartEta.push_back(sigGenEta->at(k));
                                propgenElPartPt.push_back(sigGenPt->at(k));
                                propgenElPartX.push_back(sigGenX->at(k));
                                propgenElPartY.push_back(sigGenY->at(k));
                                propgenElPartZ.push_back(sigGenZ->at(k));
                            }

                            EgN = sigEgEt->size();
                            for(unsigned int k = 0; k < EgN; k++)
                            {
                                egCrysClusterEta.push_back(sigEgEta->at(k));
                                egCrysClusterPhi.push_back(sigEgPhi->at(k));
                                egCrysClusterEt.push_back(sigEgEt->at(k)); 
                                egCrysClusterGx.push_back(sigEgGx->at(k)); 
                                egCrysClusterGy.push_back(sigEgGy->at(k)); 
                                egCrysClusterGz.push_back(sigEgGz->at(k)); 
                            }

                            bN = sigbLa->size();
                            for(unsigned int k = 0; k < bN; k++)
                            {
                                bRecHitLayer.push_back(sigbLa->at(k));
                                bRecHitGx.push_back(sigbGx->at(k));
                                bRecHitGy.push_back(sigbGy->at(k));
                                bRecHitGz.push_back(sigbGz->at(k));
                            }

                            fN = sigfDi->size();
                            for(unsigned int k = 0; k < fN; k++)
                            {
                                fRecHitDisk.push_back(sigfDi->at(k));
                                fRecHitGx.push_back(sigfGx->at(k));
                                fRecHitGy.push_back(sigfGy->at(k));
                                fRecHitGz.push_back(sigfGz->at(k));
                            }
                            Double_t random_flag = gRandom->Uniform(0, 1);
                            cout << "  Region 5/6 " << N_merged << " happened" << endl;
                            if( random_flag < 0.7 ) break;
                        }
                    } // endcap region 
                } // NoPU loop
            }
        } // end of gen matching

        bRecHitN = bRecHitGx.size();
        fRecHitN = fRecHitGx.size();
        cout << i+1 << "event, " << "N_merged : " << N_merged << endl;

        tt->Fill();
    } // event loop

    //cout << endl << "Before mixing" << endl;
    //for(std::vector<float>::iterator it = egCrysClusterPhi.begin(); it != egCrysClusterPhi.end(); ++it)
    //{
    //    std::cout << ' ' << *it;
    //}
    //cout << endl;

    //for(int k=0; k < 6; k++){
    //    cout << k << " merged event : " << mergedCheck[k] << endl; 
    //}
    fTree1->ResetBranchAddresses();

}
//------------------------------------------------------------------------------
