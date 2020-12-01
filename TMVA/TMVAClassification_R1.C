#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

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

int TMVAClassification_R1(TString myMethodList = "")
{

   // This loads the library
   TMVA::Tools::Instance();

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification_R1" << std::endl;

   // Register the training and test trees
   // Monte Carlo - Signal Sample data input
   TFile* signalfile = new TFile("inputs/tmva_input_r1_sig.root");
   TTree* signaltree = (TTree*)signalfile->Get("t");
   std::cout << "--- TMVAClassification   : Using as signal file: " << signalfile->GetName() << std::endl;
   
   // Monte Carlo - Background Sample data input (Minimum-Bias sample)
   TFile* datafile = new TFile("inputs/tmva_input_r1_bkg.root");
   TTree* datatree = (TTree*)datafile->Get("t");
   std::cout << "--- TMVAClassification   : Using as real data file: " <<datafile->GetName() << std::endl;  

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "output_R1.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   std::cout << "--- TMVAClassification   : Creating as resulting file: " << outputFile->GetName() << std::endl;

   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("Dataset_R1");

   //------------------------------------------------------------------------------------------------------------------------------------
   dataloader->AddVariable( "ntEgEt", "ntEgEt", "units", 'F' );
   dataloader->AddVariable( "ntEgEta", "ntEgEta", "units", 'F' );
   dataloader->AddVariable( "ntEgPhi", "ntEgPhi", "units", 'F' );
   dataloader->AddVariable( "IsoValue", "IsoValue", "units", 'F' );
   dataloader->AddVariable( "NumOfTrks", "NumOfTrks", "units", 'I' );
   dataloader->AddVariable( "track_pT1", "track_pT1", "units", 'F' );

   //------------------------------------------------------------------------------------------------------------------------------
   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;
  
   // You can add an arbitrary number of signal or background trees  //essential part of classification
   dataloader->AddSignalTree    ( signaltree,     signalWeight );
   dataloader->AddBackgroundTree( datatree,   backgroundWeight );
   //------------------------------------------------------------------------------------------------------------------------------

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = "";
   TCut mycutb = "";
   // for example: TCut mycuts = "B_J_mass==3.097 || B_J_chi2>0.1 || mumNHits>5 || mupNHits>5"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";

   // Tell the dataloader how to use the training and testing events
   // To also specify the number of testing events, use:
   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
       "nTrain_Signal=41613:nTrain_Background=118273:SplitMode=Random:NormMode=NumEvents:!V" );

   //------------------------------------------------------------------------------------------------------------------------------
   // Mapping
   std::map<std::string,int> Use;
   
   Use["BDTA"]            = 1; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["MLPS"]            = 1; // Recommended ANN
   Use["MLPR"]            = 1; // Recommended ANN
   Use["LD"]              = 1; // Linear Discriminant identical to Fisher
   Use["DNN_CPU"]         = 1; // Multi-core accelerated DNN.

   //------------------------------------------------------------------------------------------------------------------------------
   // ### Book MVA methods 

   // BDT - Boosted Decision Tree
   if( Use["BDTA"] ) {
       // Adaptive Boost
       factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTA",
               "!H:!V:NTrees=250:MinNodeSize=2.5%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
   }
   
   if( Use["BDTG"] ) {
   // Gradiant Boost
       factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
               "!H:!V:NTrees=250:MinNodeSize=2.5%:MaxDepth=3:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20" );
   }
  
   // MLP - Multi Layer Perceptron
   if( Use["MLPS"] ) {
       // Sigmoid activation function 
       factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPS", "H:!V:NeuronType=sigmoid:VarTransform=N:NCycles=600:HiddenLayers=5:TestRate=1:!UseRegulator");
   }
   if( Use["MLPR"] ) {
       // ReLU activation function
       factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPR", "H:!V:NeuronType=ReLU:VarTransform=N:NCycles=600:HiddenLayers=5:TestRate=1:!UseRegulator" );
   }
   
   // Linear discriminant (same as Fisher discriminant)
   if( Use["LD"] ) {
       factory->BookMethod( dataloader, TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
   }
   if( Use["DNN_CPU"] ) {
      // General layout.
      TString layoutString ("Layout=TANH|128,TANH|128,TANH|128,LINEAR");

      // Training strategies.
      TString training0("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=30,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=None,"
                        "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
      TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
      TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
      TString trainingStrategyString ("TrainingStrategy=");
      trainingStrategyString += training0 + "|" + training1 + "|" + training2;
      
      // General Options.
      TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                          "WeightInitialization=XAVIERUNIFORM");
      dnnOptions.Append (":"); dnnOptions.Append (layoutString);
      dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);
      
      TString cpuOptions = dnnOptions + ":Architecture=CPU";
      factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU", cpuOptions);
   }


//---------------------------------------------------------------------------------------------------------------------------------
   // Train MVAs using the set of training events
   factory->TrainAllMethods();
   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();
   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // Save the output
   outputFile->Close();

   std::cout << "==> TMVAClassification   : Creating as resulting root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification_R1 is done!" << std::endl;

   delete factory;
   delete dataloader;
   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

   return 0;
}

int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   return TMVAClassification_R1(methodList);
}

