#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include <TMVA/DataLoader.h>
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

void high_granularity_roc_r1() {
    
    TFile *f = TFile::Open("output_R1.root"); 
    TTreeReader tree_reader("Dataset_R1/TestTree", f);
    
    TTreeReaderValue<int>   classID(tree_reader, "classID");
    TTreeReaderValue<float> weight(tree_reader, "weight");
    TTreeReaderValue<float> bdtg(tree_reader, "BDTG");
    TTreeReaderValue<float> dnn(tree_reader, "DNN_CPU");
    TTreeReaderValue<float> ld(tree_reader, "LD");
    TTreeReaderValue<float> mlps(tree_reader, "MLPS");
    
    std::vector<float> in_ld_r1;
    std::vector<float> in_bdtg_r1;
    std::vector<float> in_dnn_r1;
    std::vector<float> in_mlps_r1;
    std::vector<float> weights;
    std::vector<bool>  targets;

    in_ld_r1.clear();
    in_bdtg_r1.clear();
    in_dnn_r1.clear();
    in_mlps_r1.clear();
    weights.clear();
    targets.clear();

    while( tree_reader.Next() ) {
        in_ld_r1.push_back((*ld));
        in_bdtg_r1.push_back((*bdtg));
        in_dnn_r1.push_back((*dnn));
        in_mlps_r1.push_back((*mlps));
        weights.push_back((*weight));
        targets.push_back((*classID));
    }

    targets.flip();

    TMVA::ROCCurve *roc_bdtg = new TMVA::ROCCurve(in_bdtg_r1, targets, weights);
    TMVA::ROCCurve *roc_mlps = new TMVA::ROCCurve(in_mlps_r1, targets, weights);
    TMVA::ROCCurve *roc_dnn = new TMVA::ROCCurve(in_dnn_r1, targets, weights);
    TMVA::ROCCurve *roc_ld = new TMVA::ROCCurve(in_ld_r1, targets, weights);

    auto sensitivity = roc_bdtg->ComputeSensitivity(10);

    /*
    // get a graph given number of points (e.g. 1000)
    
    const UInt_t points = 100;

    
    auto curve1 = roc_bdtg->GetROCCurve(10);
    auto curve2 = roc_mlps->GetROCCurve(10);
    auto curve3 = roc_dnn->GetROCCurve(10);
    auto curve4 = roc_ld->GetROCCurve(10);

    cout << curve1->GetN() << endl;
    
    curve1->SetName("r1_roc_BDTG");
    curve2->SetName("r1_roc_MLPS");
    curve3->SetName("r1_roc_DNN");
    curve4->SetName("r1_roc_Ld");

    TFile *output = new TFile("r1_roc.root", "recreate");
    curve1->Write();
    curve2->Write();
    curve3->Write();
    curve4->Write();

    output->Close();
    */
    
}
