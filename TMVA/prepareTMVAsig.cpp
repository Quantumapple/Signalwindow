#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "time.h"

#include <vector>
#include <iostream>

void prepareTMVAsig()
{
    TFile *f = TFile::Open("../option1/sig/sig_final.root"); 
    TTreeReader tree_reader("t", f);

    TTreeReaderValue<std::vector<float>> egeta(tree_reader, "ntEgEta");
    TTreeReaderValue<std::vector<float>> egphi(tree_reader, "ntEgPhi");
    TTreeReaderValue<std::vector<float>> eget(tree_reader, "ntEgEt");

    TTreeReaderValue<std::vector<bool>> pixtrkFlag(tree_reader, "ntCl_match");

    TTreeReaderValue<std::vector<float>> isoval(tree_reader, "IsoValue");
    TTreeReaderValue<std::vector<int>> numtrks(tree_reader, "NumOfTrks");

    TTreeReaderValue<std::vector<float>> pixpt1(tree_reader, "track_pT1");

    // Declare branches for output //
    vector<float>   ntEgEt;
    vector<float>   ntEgEta;
    vector<float>   ntEgPhi;
    vector<float>   IsoValue;
    vector<int>     NumOfTrks;
    vector<float>   track_pT1;

    TFile *result = new TFile("tmva_input_r6_sig.root","RECREATE");
    result->mkdir("t");
    result->cd("t");

    TTree *mytree = new TTree("t","t");

    mytree->Branch("ntEgEt", &ntEgEt);
    mytree->Branch("ntEgEta", &ntEgEta);
    mytree->Branch("ntEgPhi", &ntEgPhi);
    mytree->Branch("IsoValue", &IsoValue);
    mytree->Branch("NumOfTrks", &NumOfTrks); 
    mytree->Branch("track_pT1", &track_pT1);

    int eventCount = 0;

    clock_t tStart = clock();
    bool flag = true;
    //// Event loop ////
    while( tree_reader.Next() ) {
        if( eventCount % 500000 == 0 ) cout << "Event: " << eventCount << endl;
        //if( eventCount > 5000 ) break;

        ntEgEt.clear();
        ntEgEta.clear();
        ntEgPhi.clear();
        IsoValue.clear();
        NumOfTrks.clear();
        track_pT1.clear();

        for(int i=0; i < eget->size(); i++) {

            float EgEt  = eget->at(i);
            float EgEta = egeta->at(i);
            float EgPhi = egphi->at(i);

            // Region separation //
            if( fabs(EgEta) > 0.8 ) continue; // region 1
            if( fabs(EgEta) < 2.7 || fabs(EgEta) > 3.0 ) continue; // region 6
            
            if( pixtrkFlag->at(i) ) {
                ntEgEt.push_back(EgEt);
                ntEgEta.push_back(EgEta);
                ntEgPhi.push_back(EgPhi);
                NumOfTrks.push_back(numtrks->at(i));
                track_pT1.push_back(pixpt1->at(i));
                if( numtrks->at(i) != 0 ) {
                    IsoValue.push_back(isoval->at(i));
                }
                else {
                    IsoValue.push_back(0.);
            }


        } // egamma loop
        mytree->Fill();
        eventCount++;

    } // event loop

    cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s" << endl;
    result->Write();
}
