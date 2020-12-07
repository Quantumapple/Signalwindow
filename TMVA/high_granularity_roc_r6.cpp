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

void high_granularity_roc_r6() {
    
    TFile *f = TFile::Open("output_R6.root"); 
    TTreeReader tree_reader("Dataset_R6/TestTree", f);
    
    TTreeReaderValue<int>   classID(tree_reader, "classID");
    TTreeReaderValue<float> weight(tree_reader, "weight");
    TTreeReaderValue<float> bdtg(tree_reader, "BDTG");
    TTreeReaderValue<float> dnn(tree_reader, "DNN_CPU");
    TTreeReaderValue<float> ld(tree_reader, "LD");
    TTreeReaderValue<float> mlps(tree_reader, "MLPS");
    
    std::vector<float> in_ld_r6;
    std::vector<float> in_bdtg_r6;
    std::vector<float> in_dnn_r6;
    std::vector<float> in_mlps_r6;
    std::vector<float> weights;
    std::vector<bool>  targets;

    in_ld_r6.clear();
    in_bdtg_r6.clear();
    in_dnn_r6.clear();
    in_mlps_r6.clear();
    weights.clear();
    targets.clear();

    while( tree_reader.Next() ) {
        in_ld_r6.push_back((*ld));
        in_bdtg_r6.push_back((*bdtg));
        in_dnn_r6.push_back((*dnn));
        in_mlps_r6.push_back((*mlps));
        weights.push_back((*weight));
        targets.push_back((*classID));
    }

    targets.flip();

    TMVA::ROCCurve *roc_bdtg = new TMVA::ROCCurve(in_bdtg_r6, targets, weights);
    TMVA::ROCCurve *roc_mlps = new TMVA::ROCCurve(in_mlps_r6, targets, weights);
    TMVA::ROCCurve *roc_dnn = new TMVA::ROCCurve(in_dnn_r6, targets, weights);
    TMVA::ROCCurve *roc_ld = new TMVA::ROCCurve(in_ld_r6, targets, weights);
    
    // get a graph given number of points (e.g. 1000)
    
    const unsigned int points = 5;

    auto curve1 = roc_bdtg->GetROCCurve(points);
    auto curve2 = roc_mlps->GetROCCurve(points);
    auto curve3 = roc_dnn->GetROCCurve(points);
    auto curve4 = roc_ld->GetROCCurve(points);
    
    curve1->SetName("r6_roc_BDTG");
    curve2->SetName("r6_roc_MLPS");
    curve3->SetName("r6_roc_DNN");
    curve4->SetName("r6_roc_Ld");

    TFile *output = new TFile("r6_roc.root", "recreate");
    curve1->Write();
    curve2->Write();
    curve3->Write();
    curve4->Write();

    output->Close();
    
}
