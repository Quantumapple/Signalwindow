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
    
    Tower *eCalTower;
    Track *eCalTrack;

    TLorentzVector egVector;

    Long64_t entry;

    bool kFlag = false;

    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {
        if( entry % 10000 == 0 ) cout << "Event: " << entry << endl; 
        
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

        //cout << "Event loop start" << endl;
        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        // Loop over Genparticle
        for(unsigned int k = 0; k < branchParticle->GetEntriesFast(); k++)
        {
            particle = (GenParticle*) branchParticle->At(k);
            propgenElPartPt.push_back(particle->PT);
            propgenElPartEta.push_back(particle->Eta);
            propgenElPartPhi.push_back(particle->Phi);
        }
        
        int TrackN = branchTrack->GetEntriesFast();
        for(unsigned int k = 0; k < TrackN; k++)
        {
            eCalTrack = (Track*) branchTrack->At(k);
            double eCalX, eCalY, eCalZ, eCalEt, eCalEta, eCalPhi;
            eCalX = eCalTrack->XOuter/10.;
            eCalY = eCalTrack->YOuter/10.;

            float R = sqrt(eCalX*eCalX+eCalY*eCalY);

            //ECal track eta
            eCalEta = eCalTrack->EtaOuter;
            eCalEta += -8.084e-5;
            eCalEta = gRandom->Gaus(eCalEta, 0.004276);
            
            if( fabs(eCalEta) < 0.8 )
            {
                //ECal track Et
                eCalEt = eCalTrack->PT;
                float tempEt = eCalEt*0.1012;
                eCalEt -= tempEt;
                eCalEt = gRandom->Gaus(eCalEt, 0.02769);
                
                //ECal track phi
                eCalPhi = eCalTrack->PhiOuter;
                eCalPhi += -5.092e-5;

                if( eCalEt < 30. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.003668-0.0010);
                if( eCalEt >= 30. && eCalEt < 60. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.003668);
                if( eCalEt >= 60. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.003668+0.0010);
                
                eCalX = R*cos(eCalPhi);
                eCalY = R*sin(eCalPhi);
                eCalZ = R*cos(2.*atan(-eCalEta));
            }
            else
            {
                eCalEt = eCalTrack->PT;
                eCalPhi = eCalTrack->PhiOuter;
                eCalX = R*cos(eCalPhi);
                eCalY = R*sin(eCalPhi);
                eCalZ = R*cos(2.*atan(-eCalEta));
            }

            egCrysClusterEta.push_back(eCalEta);
            egCrysClusterPhi.push_back(eCalPhi);
            egCrysClusterEt.push_back(eCalEt);
            egCrysClusterGx.push_back(eCalX);
            egCrysClusterGy.push_back(eCalY);
            egCrysClusterGz.push_back(eCalZ);
        }

        // Loop over 1st pixel layer
        for(unsigned int k = 0; k < branchPixelL1->GetEntriesFast(); k++)
		{
			trackL1 = (Track*) branchPixelL1->At(k);
			if( fabs(trackL1->ZOuter) < 200 )
			{
				if( sqrt(pow(trackL1->XOuter,2)+pow(trackL1->YOuter,2)) < 29.9 ) continue; 

				bRecHitLayer.push_back(1);
                float L1X = trackL1->XOuter/10.;
                float L1Y = trackL1->YOuter/10.;
                float L1Z = trackL1->ZOuter/10.;
				
                bRecHitGx.push_back(gRandom->Gaus(L1X, 5e-6));
				bRecHitGy.push_back(gRandom->Gaus(L1Y, 5e-6));
				bRecHitGz.push_back(gRandom->Gaus(L1Z, 15e-6));
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
                float L2X = trackL2->XOuter/10.;
                float L2Y = trackL2->YOuter/10.;
                float L2Z = trackL2->ZOuter/10.;
				
                bRecHitGx.push_back(gRandom->Gaus(L2X, 5e-6));
				bRecHitGy.push_back(gRandom->Gaus(L2Y, 5e-6));
				bRecHitGz.push_back(gRandom->Gaus(L2Z, 15e-6));
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
                float L3X = trackL3->XOuter/10.;
                float L3Y = trackL3->YOuter/10.;
                float L3Z = trackL3->ZOuter/10.;
				
                bRecHitGx.push_back(gRandom->Gaus(L3X, 5e-6));
				bRecHitGy.push_back(gRandom->Gaus(L3Y, 5e-6));
				bRecHitGz.push_back(gRandom->Gaus(L3Z, 15e-6));
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
                float L4X = trackL4->XOuter/10.;
                float L4Y = trackL4->YOuter/10.;
                float L4Z = trackL4->ZOuter/10.;
				
                bRecHitGx.push_back(gRandom->Gaus(L4X, 5e-6));
				bRecHitGy.push_back(gRandom->Gaus(L4Y, 5e-6));
				bRecHitGz.push_back(gRandom->Gaus(L4Z, 15e-6));
			}
		}
		bRecHitN = bRecHitGx.size();

		// Loop over 1st pixel Disk
		for(unsigned int k = 0; k < branchPixelD1->GetEntriesFast(); k++)
		{
			trackD1 = (Track*) branchPixelD1->At(k);
			float R = sqrt(pow(trackD1->XOuter,2)+pow(trackD1->YOuter,2));
			if( fabs(trackD1->ZOuter) > 248 && R > 28.0 )
                //if( fabs(trackD1->ZOuter) > 248 )
			{
				double j11 = 0, k11, l11;
				double D1X, D1Y, D1Z;
				//trackD1->XOuter
				j11 = trackD1->XOuter;
				k11 = fmod(j11 , 0.012);
				if(fabs(k11) < 0.006) D1X = j11 - k11;
				else{
					if(j11 < 0) D1X = j11 - k11 - 0.012;
					else D1X = j11 - k11 + 0.012;
				}
				//trackD1->YOuter
				j11 = trackD1->YOuter;
				k11 = fmod(j11 , 0.012);
				if(fabs(k11) < 0.006) D1Y = j11 - k11;
				else{
					if(j11 < 0) D1Y = j11 - k11 - 0.012;
					else D1Y = j11 - k11 + 0.012;
				}
				//trackD1->ZOuter
				j11 = trackD1->ZOuter;
				k11 = fmod(j11 , 8.);
				if(fabs(k11) < 4.) D1Z = j11 - k11;
				else{
					if(j11 < 0) D1Z = j11 - k11 - 8.;
					else D1Z = j11 - k11 + 8.;
				}
				fRecHitDisk.push_back(1);
				fRecHitGx.push_back(D1X/10.);
				fRecHitGy.push_back(D1Y/10.);
				fRecHitGz.push_back(D1Z/10.);
			}
		}

		// Loop over 2nd pixel Disk
		for(unsigned int k = 0; k < branchPixelD2->GetEntriesFast(); k++)
		{
			trackD2 = (Track*) branchPixelD2->At(k);
			float R = sqrt(pow(trackD2->XOuter,2)+pow(trackD2->YOuter,2));
			if( fabs(trackD2->ZOuter) > 315 && R > 28.0 )
                //if( fabs(trackD2->ZOuter) > 315 )
			{
				double j11 = 0, k11, l11;
				double D2X, D2Y, D2Z;
				//trackD2->XOuter
				j11 = trackD2->XOuter;
				k11 = fmod(j11 , 0.012);
				if(fabs(k11) < 0.006) D2X = j11 - k11;
				else{
					if(j11 < 0) D2X = j11 - k11 - 0.012;
					else D2X = j11 - k11 + 0.012;
				}
				//trackD2->YOuter
				j11 = trackD2->YOuter;
				k11 = fmod(j11 , 0.012);
				if(fabs(k11) < 0.006) D2Y = j11 - k11;
				else{
					if(j11 < 0) D2Y = j11 - k11 - 0.012;
					else D2Y = j11 - k11 + 0.012;
				}
				//trackD2->ZOuter
				j11 = trackD2->ZOuter;
				k11 = fmod(j11 , 8.);
				if(fabs(k11) < 4.) D2Z = j11 - k11;
				else{
					if(j11 < 0) D2Z = j11 - k11 - 8.;
					else D2Z = j11 - k11 + 4.;
                }
				fRecHitDisk.push_back(2);
				fRecHitGx.push_back(D2X/10.);
				fRecHitGy.push_back(D2Y/10.);
				fRecHitGz.push_back(D2Z/10.);
            }
		}

		// Loop over 3rd pixel Disk
		for(unsigned int k = 0; k < branchPixelD3->GetEntriesFast(); k++)
		{
			trackD3 = (Track*) branchPixelD3->At(k);
			float R = sqrt(pow(trackD3->XOuter,2)+pow(trackD3->YOuter,2));
			if( fabs(trackD3->ZOuter) > 407 && R > 28.0 )
				//if( fabs(trackD3->ZOuter) > 407 )
			{
				double j11 = 0, k11, l11;
				double D3X, D3Y, D3Z;
				//trackD3->XOuter
				j11 = trackD3->XOuter;
				k11 = fmod(j11 , 0.012);
				if(fabs(k11) < 0.006) D3X = j11 - k11;
				else{
					if(j11 < 0) D3X = j11 - k11 - 0.012;
					else D3X = j11 - k11 + 0.012;
				}
				//trackD3->YOuter
				j11 = trackD3->YOuter;
				k11 = fmod(j11 , 0.012);
				if(fabs(k11) < 0.006) D3Y = j11 - k11;
				else{
					if(j11 < 0) D3Y = j11 - k11 - 0.012;
					else D3Y = j11 - k11 + 0.012;
				}
				//trackD3->ZOuter
				j11 = trackD3->ZOuter;
				k11 = fmod(j11 , 8.);
				if(fabs(k11) < 4.) D3Z = j11 - k11;
				else{
					if(j11 < 0) D3Z = j11 - k11 - 8.;
					else D3Z = j11 - k11 + 4.;
                }
				fRecHitDisk.push_back(3);
				fRecHitGx.push_back(D3X/10.);
				fRecHitGy.push_back(D3Y/10.);
				fRecHitGz.push_back(D3Z/10.);
			}
		}

		// Loop over 4th pixel Disk
		for(unsigned int k = 0; k < branchPixelD4->GetEntriesFast(); k++)
		{
			trackD4 = (Track*) branchPixelD4->At(k);
			float R = sqrt(pow(trackD4->XOuter,2)+pow(trackD4->YOuter,2));
			if( fabs(trackD4->ZOuter) > 521 && R > 28.0 )
				//if( fabs(trackD4->ZOuter) > 521 )
			{
				double j11 = 0, k11, l11;
				double D4X, D4Y, D4Z;
				//trackD4->XOuter
				j11 = trackD4->XOuter;
				k11 = fmod(j11 , 0.012);
				if(fabs(k11) < 0.006) D4X = j11 - k11;
				else{
					if(j11 < 0) D4X = j11 - k11 - 0.012;
					else D4X = j11 - k11 + 0.012;
				}
				//trackD4->YOuter
				j11 = trackD4->YOuter;
				k11 = fmod(j11 , 0.012);
				if(fabs(k11) < 0.006) D4Y = j11 - k11;
				else{
					if(j11 < 0) D4Y = j11 - k11 - 0.012;
					else D4Y = j11 - k11 + 0.012;
				}
				//trackD4->ZOuter
				j11 = trackD4->ZOuter;
				k11 = fmod(j11 , 8.);
				if(fabs(k11) < 4.) D4Z = j11 - k11;
				else{
					if(j11 < 0) D4Z = j11 - k11 - 8.;
					else D4Z = j11 - k11 + 8.;
                }
				fRecHitDisk.push_back(4);
				fRecHitGx.push_back(D4X/10.);
				fRecHitGy.push_back(D4Y/10.);
				fRecHitGz.push_back(D4Z/10.);
			}
		}

		// Loop over 5th pixel Disk
		for(unsigned int k = 0; k < branchPixelD5->GetEntriesFast(); k++)
		{
			trackD5 = (Track*) branchPixelD5->At(k);
			float R = sqrt(pow(trackD5->XOuter,2)+pow(trackD5->YOuter,2));
			if( fabs(trackD5->ZOuter) > 667 && R > 28.0 )
				//if( fabs(trackD5->ZOuter) > 667 )
			{
				double j11 = 0, k11, l11;
				double D5X, D5Y, D5Z;
				//trackD5->XOuter
				j11 = trackD5->XOuter;
				k11 = fmod(j11 , 0.012);
				if(fabs(k11) < 0.006) D5X = j11 - k11;
				else{
					if(j11 < 0) D5X = j11 - k11 - 0.012;
					else D5X = j11 - k11 + 0.012;
				}
				//trackD5->YOuter
				j11 = trackD5->YOuter;
				k11 = fmod(j11 , 0.012);
				if(fabs(k11) < 0.006) D5Y = j11 - k11;
				else{
					if(j11 < 0) D5Y = j11 - k11 - 0.012;
					else D5Y = j11 - k11 + 0.012;
				}
				//trackD5->ZOuter
				j11 = trackD5->ZOuter;
				k11 = fmod(j11 , 8.);
				if(fabs(k11) < 4.) D5Z = j11 - k11;
				else{
					if(j11 < 0) D5Z = j11 - k11 - 8.;
					else D5Z = j11 - k11 + 8.;
                }
				fRecHitDisk.push_back(5);
				fRecHitGx.push_back(D5X/10.);
				fRecHitGy.push_back(D5Y/10.);
				fRecHitGz.push_back(D5Z/10.);
			}
		}

		fRecHitN = fRecHitGx.size();

		tt->Fill();
	} // event loop
}

//------------------------------------------------------------------------------

void Example13(const char *inputFile)
{
	gSystem->Load("libDelphes");

	TChain *chain = new TChain("Delphes");
	chain->Add(inputFile);

	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

	TFile *result = new TFile("results.root","RECREATE");
	result->mkdir("l1PiXTRKTree");
	result->cd("l1PiXTRKTree");

	TTree *tree = new TTree("L1PiXTRKTree","L1PiXTRKTree");

	AnalyseEvents(treeReader, result, tree);
	result->Write();

	cout << "** Exiting..." << endl;

	delete treeReader;
	delete chain;

	result->Close();
}

//------------------------------------------------------------------------------
