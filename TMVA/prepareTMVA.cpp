#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "time.h"

#include <vector>
#include <iostream>

void prepareTMVA()
{
    TFile *f = TFile::Open("../option1/bkg/bkg_final.root"); 
    TTreeReader tree_reader("t", f);

    TTreeReaderValue<int> egn(tree_reader, "ntnEg2");

    TTreeReaderValue<std::vector<float>> egeta(tree_reader, "ntEgEta");
    TTreeReaderValue<std::vector<float>> egphi(tree_reader, "ntEgPhi");
    TTreeReaderValue<std::vector<float>> eget(tree_reader, "ntEgEt");
    TTreeReaderValue<std::vector<float>> egem(tree_reader, "ntEgEem");
    TTreeReaderValue<std::vector<float>> eghad(tree_reader, "ntEgEhad");

    TTreeReaderValue<std::vector<bool>> pixtrkFlag(tree_reader, "ntCl_match");
    TTreeReaderValue<std::vector<bool>> trkisoFlag(tree_reader, "ntCl_iso_match");

    TTreeReaderValue<std::vector<float>> isoval(tree_reader, "IsoValue");
    TTreeReaderValue<std::vector<int>> numtrks(tree_reader, "NumOfTrks");

    TTreeReaderValue<std::vector<float>> pixpt1(tree_reader, "track_pT1");

    // Declare branches for output //
    //int             ntnEg2;
    vector<float>   ntEgEt;
    vector<float>   ntEgEta;
    vector<float>   ntEgPhi;
    vector<float>   IsoValue;
    vector<int>     NumOfTrks;
    vector<float>   track_pT1;

    TFile *result = new TFile("tmva_input_r6_bkg.root","RECREATE");
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
            float Eem   = egem->at(i);
            float Ehad  = eghad->at(i);
            float ratio = Ehad/Eem;
            float rndm = gRandom->Uniform(1.);

            // Region separation //
            //if( fabs(EgEta) > 0.8 ) continue; // region 1
            if( fabs(EgEta) < 2.7 || fabs(EgEta) > 3.0 ) continue; // region 6


            if( EgEt <= 10. ) {
                if( (ratio < 0.01 || ratio > 200) || Ehad > 200. ) continue;
            }
            else if( EgEt < 15. ) {
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
                        }
                    }
                }
                else if( EgEt < 20. ) {
                    if( rndm > 0.75 ) {
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
                        }
                    }
                }

                else if( EgEt < 28. ) {
                    if( rndm > 0.76 ) {
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
                        }
                    }
                }

                else if( EgEt < 30. ) {
                    if( rndm > 0.66 ) {
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
                        }
                    }
                }

                else if( EgEt < 33. ) {
                    if( rndm > 0.61 ) {
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
                        }
                    }
                }

                else if( EgEt < 35. ) {
                    if( rndm > 0.75 ) {
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
                        }
                    }
                }

                else if( EgEt < 40. ) {
                    if( rndm > 0.53 ) {
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
                        }
                    }
                }

                else if( EgEt < 41. ) {
                    if( rndm > 0.38 ) {
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
                        }
                    }
                }

                else if( EgEt < 42. ) {
                    if( rndm > 0.52 ) {
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
                        }
                    }
                }

                else if( EgEt < 43. ) {
                    if( rndm > 0.60 ) {
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
                        }
                    }
                }

                else if( EgEt < 44. ) {
                    if( rndm > 0.43 ) {
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
                        }
                    }
                }

                else if( EgEt < 47. ) {
                    if( rndm > 0.47 ) {
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
                        }
                    }
                }

                else if( EgEt < 52. ) {
                    if( rndm > 0.27 ) {
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
                        }
                    }
                }

                else if( EgEt < 63. ) {
                    if( rndm > 0.0 ) {
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
                        }
                    }
                }

                else if( EgEt < 65. ) {
                    if( rndm > 0.10 ) {
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
                        }
                    }
                }

                else if( EgEt < 70. ) {
                    if( rndm > 0.40 ) {
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
                        }
                    }
                }

                else {
                    if ( rndm > 0.0 ) {
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
                        }
                    }
                }

            }

        } // egamma loop
        mytree->Fill();
        eventCount++;

    } // event loop

    cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s" << endl;
    result->Write();
}
