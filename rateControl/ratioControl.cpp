#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "time.h"

#include <vector>
#include <iostream>

void ratioControl()
{
    TFile *f = TFile::Open("bkg_final.root"); 
    TTreeReader tree_reader("t", f);

    TTreeReaderValue<int> egn(tree_reader, "ntnEg2");

    TTreeReaderValue<std::vector<float>> egeta(tree_reader, "ntEgEta");
    TTreeReaderValue<std::vector<float>> egphi(tree_reader, "ntEgPhi");
    TTreeReaderValue<std::vector<float>> eget(tree_reader, "ntEgEt");
    TTreeReaderValue<std::vector<float>> egem(tree_reader, "ntEgEem");
    TTreeReaderValue<std::vector<float>> eghad(tree_reader, "ntEgEhad");

    TTreeReaderValue<std::vector<bool>> pixtrkFlag(tree_reader, "ntCl_match");
    TTreeReaderValue<std::vector<bool>> trkisoFlag(tree_reader, "ntCl_iso_match");

    //TTreeReaderValue<std::vector<float>> isoval(tree_reader, "IsoValue");
    //TTreeReaderValue<std::vector<int>> numtrks(tree_reader, "NumOfTrks");

    //TTreeReaderValue<std::vector<float>> pixpt1(tree_reader, "track_pT1");
    //TTreeReaderValue<std::vector<float>> pixpt2(tree_reader, "track_pT2");
    //TTreeReaderValue<std::vector<float>> pixpt3(tree_reader, "track_pT3");
    //TTreeReaderValue<std::vector<float>> pixpt4(tree_reader, "track_pT4");
    //TTreeReaderValue<std::vector<float>> pixpt5(tree_reader, "track_pT5");
    //TTreeReaderValue<std::vector<float>> pixpt6(tree_reader, "track_pT6");

    // Declare branches for output //
    int             ntnEg2;
    vector<float>   ntEgEt;
    vector<float>   ntEgEta;
    vector<float>   ntEgPhi;
    vector<bool>    ntCl_match;
    vector<bool>    ntCl_iso_match;
    //vector<float>   IsoValue;
    //vector<int>     NumOfTrks;
    //vector<float>   track_pT1;
    //vector<float>   track_pT2;
    //vector<float>   track_pT3;
    //vector<float>   track_pT4;
    //vector<float>   track_pT5;
    //vector<float>   track_pT6;

    TFile *result = new TFile("results.root","RECREATE");
    result->mkdir("t");
    result->cd("t");

    TTree *mytree = new TTree("t","t");

    mytree->Branch("ntnEg2", &ntnEg2, "ntnEg2/I");

    mytree->Branch("ntEgEt", &ntEgEt);
    mytree->Branch("ntEgEta", &ntEgEta);
    mytree->Branch("ntEgPhi", &ntEgPhi);

    mytree->Branch("ntCl_match", &ntCl_match);
    mytree->Branch("ntCl_iso_match", &ntCl_iso_match);

    //mytree->Branch("IsoValue", &IsoValue);
    //mytree->Branch("NumOfTrks", &NumOfTrks); 

    //mytree->Branch("track_pT1", &track_pT1);
    //mytree->Branch("track_pT2", &track_pT2);
    //mytree->Branch("track_pT3", &track_pT3);
    //mytree->Branch("track_pT4", &track_pT4);
    //mytree->Branch("track_pT5", &track_pT5);
    //mytree->Branch("track_pT6", &track_pT6);


    int eventCount = 0;

    clock_t tStart = clock();
    bool flag = true;
    //// Event loop ////
    while( tree_reader.Next() ) {
        if( eventCount % 500000 == 0 ) cout << "Event: " << eventCount << endl;
        //if( eventCount > 5000 ) break;

        ntnEg2 = 0;

        ntEgEt.clear();
        ntEgEta.clear();
        ntEgPhi.clear();

        ntCl_match.clear();
        ntCl_iso_match.clear();

        //IsoValue.clear();
        //NumOfTrks.clear();

        //track_pT1.clear();
        //track_pT2.clear();
        //track_pT3.clear();
        //track_pT4.clear();
        //track_pT5.clear();
        //track_pT6.clear();

        for(int i=0; i < eget->size(); i++) {

            float EgEt  = eget->at(i);
            float EgEta = egeta->at(i);
            float EgPhi = egphi->at(i);
            float Eem   = egem->at(i);
            float Ehad  = eghad->at(i);
            float ratio = Ehad/Eem;

            // ratio = HoverE
            // Simplest ratio cut
            //if( (ratio < 0.3 || ratio > 10) || Ehad > 200. ) continue;
            
                
            //if( Ehad > 200. ) continue;

            float rndm = gRandom->Uniform(1.);

            if( EgEt <= 10. ) {
                if( (ratio < 0.01 || ratio > 200) || Ehad > 200. ) continue;
            }
            else if( EgEt < 15. ) {
                //if( (ratio < 0.3 || ratio > 10) || Ehad > 200. ) continue;
                if( (ratio < 0.01 || ratio > 20) || Ehad > 200. ) continue;
                if( rndm < 0.05 ) continue;
            }
            else {
                if( ratio > 400 ) continue;

                else {

                    if( EgEt <= 20. )                    { if( ratio < 0.01 && rndm < 0.88 )  continue; } 
                    else if( EgEt > 20. && EgEt <= 30. ) { if( ratio < 0.01 && rndm < 0.72 )  continue; }
                    else if( EgEt > 30. && EgEt <= 40. ) { if( ratio < 0.01 && rndm < 0.74 )  continue; }
                    else if( EgEt > 40. && EgEt <= 55. ) { if( ratio < 0.01 && rndm < 0.78 )  continue; }
                    else if( EgEt > 55. && EgEt <= 60. ) { if( ratio < 0.01 && rndm < 0.80 )  continue; }
                    else if( EgEt > 60. && EgEt <= 67. ) { if( ratio < 0.01 && rndm < 0.87 )  continue; }
                    else if( EgEt > 67. && EgEt <= 70. ) { if( ratio < 0.01 && rndm < 0.93 )  continue; }
                    else if( EgEt > 70. && EgEt <= 80. ) { if( ratio < 0.01 && rndm < 0.99 )  continue; }
                    else if( EgEt > 80. && EgEt <= 85. ) { if( ratio < 0.01 && rndm < 0.95 )  continue; }
                    else if( EgEt > 85. && EgEt <= 95. ) { if( ratio < 0.01 && rndm < 0.87 )  continue; }
                    else                                 { if( ratio < 0.01 && rndm < 0.90 )  continue; } 
                }
            }
            

            if( flag ) {
                if( EgEt < 18 ) {
                  if( rndm > 0.80 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }

                else if( EgEt < 20. ) {
                  if( rndm > 0.75 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }

                else if( EgEt < 28. ) {
                  if( rndm > 0.76 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }

                else if( EgEt < 30. ) {
                  if( rndm > 0.66 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }

                else if( EgEt < 33. ) {
                  if( rndm > 0.61 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }

                else if( EgEt < 35. ) {
                  if( rndm > 0.75 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }

                else if( EgEt < 40. ) {
                  if( rndm > 0.53 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }

                else if( EgEt < 41. ) {
                  if( rndm > 0.38 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }

                else if( EgEt < 42. ) {
                  if( rndm > 0.52 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }

                else if( EgEt < 43. ) {
                  if( rndm > 0.60 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }

                else if( EgEt < 44. ) {
                  if( rndm > 0.43 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }

                else if( EgEt < 47. ) {
                  if( rndm > 0.47 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }

                else if( EgEt < 52. ) {
                  if( rndm > 0.27 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }

                else if( EgEt < 63. ) {
                  if( rndm > 0.0 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }
                
                else if( EgEt < 65. ) {
                  if( rndm > 0.10 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }
                
                else if( EgEt < 70. ) {
                  if( rndm > 0.40 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }
                
                else {
                  if ( rndm > 0.0 ) {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(pixtrkFlag->at(i));
                    ntCl_iso_match.push_back(trkisoFlag->at(i));
                  }
                  else {
                    ntnEg2++;
                    ntEgEt.push_back(EgEt);
                    ntEgEta.push_back(EgEta);
                    ntEgPhi.push_back(EgPhi);
                    ntCl_match.push_back(0.);
                    ntCl_iso_match.push_back(0.);
                  }
                }
                
            }

            
            // Fill variables
            //ntnEg2++;
            //ntEgEt.push_back(EgEt);
            //ntEgEta.push_back(EgEta);
            //ntEgPhi.push_back(EgPhi);
            //ntCl_match.push_back(pixtrkFlag->at(i));
            //ntCl_iso_match.push_back(trkisoFlag->at(i));

            //if(ntCl_match.at(i) == 1 ){
            //IsoValue.push_back(IsoValue.at(flag));
            //NumOfTrks.push_back(NumOfTrks.at(flag));
            //}

        } // egamma loop
        mytree->Fill();
        eventCount++;

    } // event loop

    cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s" << endl;
    result->Write();
}
