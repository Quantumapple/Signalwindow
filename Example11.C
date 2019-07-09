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

}

double ReturnCustomGaus(int flag)
{

    int a = gRandom->Integer(2);

    TF1 *EtaGaus = new TF1("EtaGaus","gaus(0)",-0.03,0.03);
    EtaGaus->SetParameters(352.943, -3.9404e-05, 0.00424596);
    
    TF1 *PhiGaus = new TF1("PhiGaus","gaus(0)",-0.03,0.03);
    PhiGaus->SetParameters(421.143, -2.95267e-05, 0.00347567);
    
    TF1 *EtGaus = new TF1("EtGaus","gaus(0)",-0.1,1.0);
    EtGaus->SetParameters(1043.81, 0.0992722, 0.0337839);
    
    TF1 *XGaus = new TF1("XGaus","gaus(0)",-2.0,2.0);
    XGaus->SetParameters(343.316, -0.000918254, 0.347522);
    
    TF1 *YGaus = new TF1("YGaus","gaus(0)",-2.0,2.0);
    YGaus->SetParameters(346.369, -0.00223196, 0.346173);
    
    TF1 *ZGaus = new TF1("ZGaus","gaus(0)",-3.0,3.0);
    ZGaus->SetParameters(222.469, 0.00612539, 0.8371);


    if( flag == 1 && a == 0 ) return EtaGaus->GetRandom();
    if( flag == 1 && a == 1 ) return -1.*EtaGaus->GetRandom();
    
    if( flag == 2 && a == 0 ) return PhiGaus->GetRandom();
    if( flag == 2 && a == 1 ) return -1.*PhiGaus->GetRandom();

    if( flag == 3 ) return fabs(EtGaus->GetRandom());
    
    if( flag == 4 && a == 0 ) return XGaus->GetRandom();
    if( flag == 4 && a == 1 ) return -1.*XGaus->GetRandom();
    
    if( flag == 5 && a == 0 ) return YGaus->GetRandom();
    if( flag == 5 && a == 1 ) return -1.*YGaus->GetRandom();
    
    if( flag == 6 && a == 0 ) return ZGaus->GetRandom();
    if( flag == 6 && a == 1 ) return -1.*ZGaus->GetRandom();
}

//------------------------------------------------------------------------------


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

            //ECal track phi
            eCalPhi = eCalTrack->PhiOuter;
            eCalPhi += ReturnCustomGaus(2);

            //ECal track eta
            eCalEta = eCalTrack->EtaOuter;
            eCalEta += ReturnCustomGaus(1);

            //ECal track Et
            eCalEt = eCalTrack->PT;
            double eCalnum = eCalEt*ReturnCustomGaus(3);
            int a = gRandom->Integer(2);
            if( a == 0 ) eCalEt += eCalnum;
            if( a == 1 ) eCalEt -= eCalnum;

            //ECal track X, Y, Z 
            eCalX = eCalTrack->XOuter/10.;
            eCalX += ReturnCustomGaus(4);
            eCalY = eCalTrack->YOuter/10.;
            eCalY += ReturnCustomGaus(5);
            eCalZ = eCalTrack->ZOuter/10.;
            eCalZ += ReturnCustomGaus(6);

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
				double j11 = 0, k11, l11;
				double L1X, L1Y, L1Z;
				//trackL1->XOuter
				j11 = trackL1->XOuter;
				k11 = fmod(j11 , 0.002);
				if(fabs(k11) < 0.001) L1X = j11 - k11;
				else{
					if(j11 < 0) L1X = j11 - k11 - 0.002;
					else L1X = j11 - k11 + 0.002;
				}
				//trackL1->YOuter
				j11 = trackL1->YOuter;
				k11 = fmod(j11 , 0.002);
				if(fabs(k11) < 0.001) L1Y = j11 - k11;
				else {
					if(j11 < 0) L1Y = j11 - k11 - 0.002;
					else L1Y = j11 - k11 + 0.002;
				}
				//trackL1->ZOuter
				j11 = trackL1->ZOuter;
				k11 = fmod(j11 , 0.008);
				if(fabs(k11) < 0.004) L1Z = j11 - k11;
				else{
					if(j11 < 0) L1Z = j11 - k11 - 0.008;
					else L1Z = j11 - k11 + 0.008;
				}

				bRecHitLayer.push_back(1);
				bRecHitGx.push_back(L1X/10.);
				bRecHitGy.push_back(L1Y/10.);
				bRecHitGz.push_back(L1Z/10.);
			}
		}

		// Loop over 2nd pixel layer
		for(unsigned int k = 0; k < branchPixelL2->GetEntriesFast(); k++)
		{
			trackL2 = (Track*) branchPixelL2->At(k);
			if( fabs(trackL2->ZOuter) < 200 )
			{
				if( sqrt(pow(trackL2->XOuter,2)+pow(trackL2->YOuter,2)) < 69.9 ) continue; 
				double j11 = 0, k11, l11;
				double L1X, L1Y, L1Z;
				//trackL2->XOuter
				j11 = trackL2->XOuter;
				k11 = fmod(j11 , 0.002);
				if(fabs(k11) < 0.001) L1X = j11 - k11;
				else{
					if(j11 < 0) L1X = j11 - k11 - 0.002;
					else L1X = j11 - k11 + 0.002;
				}
				//trackL2->YOuter
				j11 = trackL2->YOuter;
				k11 = fmod(j11 , 0.002);
				if(fabs(k11) < 0.001) L1Y = j11 - k11;
				else{
					if(j11 < 0) L1Y = j11 - k11 - 0.002;
					else L1Y = j11 - k11 + 0.002;
				}
				//trackL2->ZOuter
				j11 = trackL2->ZOuter;
				k11 = fmod(j11 , 0.008);
				if(fabs(k11) < 0.004) L1Z = j11 - k11;
				else{
					if(j11 < 0) L1Z = j11 - k11 - 0.004;
					else L1Z = j11 - k11 + 0.004;
				}

				bRecHitLayer.push_back(2);
				bRecHitGx.push_back(L1X/10.);
				bRecHitGy.push_back(L1Y/10.);
				bRecHitGz.push_back(L1Z/10.);
			}
		}

		// Loop over 3rd pixel layer
		for(unsigned int k = 0; k < branchPixelL3->GetEntriesFast(); k++)
		{
			trackL3 = (Track*) branchPixelL3->At(k);
			if( fabs(trackL3->ZOuter) < 200 )
			{
				if( sqrt(pow(trackL3->XOuter,2)+pow(trackL3->YOuter,2)) < 117.9 ) continue; 
				double j11 = 0, k11, l11;
				double L1X, L1Y, L1Z;
				//trackL3->XOuter
				j11 = trackL3->XOuter;
				k11 = fmod(j11 , 0.002);
				if(fabs(k11) < 0.001) L1X = j11 - k11;
				else{
					if(j11 < 0) L1X = j11 - k11 - 0.002;
					else L1X = j11 - k11 + 0.002;
				}
				//trackL3->YOuter
				j11 = trackL3->YOuter;
				k11 = fmod(j11 , 0.002);
				if(fabs(k11) < 0.001) L1Y = j11 - k11;
				else{
					if(j11 < 0) L1Y = j11 - k11 - 0.002;
					else L1Y = j11 - k11 + 0.002;
				}
				//trackL3->ZOuter
				j11 = trackL3->ZOuter;
				k11 = fmod(j11 , 0.008);
				if(fabs(k11) < 0.004) L1Z = j11 - k11;
				else{
					if(j11 < 0) L1Z = j11 - k11 - 0.008;
					else L1Z = j11 - k11 + 0.008;
				}

				bRecHitLayer.push_back(3);
				bRecHitGx.push_back(L1X/10.);
				bRecHitGy.push_back(L1Y/10.);
				bRecHitGz.push_back(L1Z/10.);
			}
		}

		// Loop over 4th pixel layer
		for(unsigned int k = 0; k < branchPixelL4->GetEntriesFast(); k++)
		{
			trackL4 = (Track*) branchPixelL4->At(k);
			if( fabs(trackL4->ZOuter) < 200 )
			{
				if( sqrt(pow(trackL4->XOuter,2)+pow(trackL4->YOuter,2)) < 157.9 ) continue; 
				double j11 = 0, k11, l11;
				double L1X, L1Y, L1Z;
				//trackL4->XOuter
				j11 = trackL4->XOuter;
				k11 = fmod(j11 , 0.002);
				if(fabs(k11) < 0.001) L1X = j11 - k11;
				else{
					if(j11 < 0) L1X = j11 - k11 - 0.002;
					else L1X = j11 - k11 + 0.002;
				}
				//trackL4->YOuter
				j11 = trackL4->YOuter;
				k11 = fmod(j11 , 0.002);
				if(fabs(k11) < 0.001) L1Y = j11 - k11;
				else{
					if(j11 < 0) L1Y = j11 - k11 - 0.002;
					else L1Y = j11 - k11 + 0.002;
				}
				//trackL4->ZOuter
				j11 = trackL4->ZOuter;
				k11 = fmod(j11 , 0.008);
				if(fabs(k11) < 0.004) L1Z = j11 - k11;
				else{
					if(j11 < 0) L1Z = j11 - k11 - 0.008;
					else L1Z = j11 - k11 + 0.008;
				}

				bRecHitLayer.push_back(4);
				bRecHitGx.push_back(L1X/10.);
				bRecHitGy.push_back(L1Y/10.);
				bRecHitGz.push_back(L1Z/10.);
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

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
	result->Print("png");
}

//------------------------------------------------------------------------------

void Example11(const char *inputFile)
{
	gSystem->Load("libDelphes");

	TChain *chain = new TChain("Delphes");
	chain->Add(inputFile);

	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	//ExRootResult *result = new ExRootResult();

	TFile *result = new TFile("results.root","RECREATE");
	//TFile *result = new TFile("L1X_modified_by_Example9.root","RECREATE");
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
