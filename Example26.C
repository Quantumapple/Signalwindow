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

void AnalyseEvents(ExRootTreeReader *treeReader, TFile *result, TTree *tt);

//------------------------------------------------------------------------------

void Example26(const char *inputFile)
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

void AnalyseEvents(ExRootTreeReader *treeReader, TFile *result, TTree *tt)
{
    int event = 0;

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
    std::vector<int> propgenElPartPdgID;

    std::vector<float> egCrysClusterEt;
    std::vector<float> egCrysClusterEta;
    std::vector<float> egCrysClusterPhi;
    std::vector<float> egCrysClusterGx;
    std::vector<float> egCrysClusterGy;
    std::vector<float> egCrysClusterGz;
    
    std::vector<float> egCrysEt;
    std::vector<float> egCrysEta;
    std::vector<float> egCrysPhi;
    std::vector<float> egCrysGx;
    std::vector<float> egCrysGy;
    std::vector<float> egCrysGz;
   
    tt->Branch("event", &event, "event/I");
    tt->Branch("propgenElPartPhi", &propgenElPartPhi); 
    tt->Branch("propgenElPartEta", &propgenElPartEta); 
    tt->Branch("propgenElPartPt", &propgenElPartPt);
    tt->Branch("propgenElPartPdgID", &propgenElPartPdgID);
    
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
    
    tt->Branch("egCrysEt", &egCrysEt);
    tt->Branch("egCrysEta", &egCrysEta);
    tt->Branch("egCrysPhi", &egCrysPhi);
    tt->Branch("egCrysGx",&egCrysGx);
    tt->Branch("egCrysGy",&egCrysGy);
    tt->Branch("egCrysGz",&egCrysGz);

    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchTower = treeReader->UseBranch("Tower");
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
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");

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
    
    Tower *CalTower;
    Track *CalTrack;

    Long64_t entry;

    bool kFlag = false;

    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {
        event = 1;
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
        propgenElPartPdgID.clear();
    
        egCrysClusterEta.clear();
        egCrysClusterPhi.clear();
        egCrysClusterEt.clear();
        egCrysClusterGx.clear();
        egCrysClusterGy.clear();
        egCrysClusterGz.clear();
        
        egCrysEta.clear();
        egCrysPhi.clear();
        egCrysEt.clear();
        egCrysGx.clear();
        egCrysGy.clear();
        egCrysGz.clear();

        //cout << "Event loop start" << endl;
        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        // Loop over Genparticle
        for(unsigned int k = 0; k < branchParticle->GetEntriesFast(); k++)
        {
            particle = (GenParticle*) branchParticle->At(k);
            propgenElPartPt.push_back(particle->PT);
            propgenElPartPdgID.push_back(particle->PID);
            propgenElPartEta.push_back(particle->Eta);
            propgenElPartPhi.push_back(particle->Phi);
        }

        int TrackN = branchTrack->GetEntriesFast();
        for(unsigned int k = 0; k < TrackN; k++)
        {
            CalTrack = (Track*) branchTrack->At(k);
            double eCalX, eCalY, eCalZ, eCalEt, eCalEta, eCalPhi;
            
            eCalPhi = CalTrack->PhiOuter;
            eCalEta = CalTrack->EtaOuter;
            eCalX = CalTrack->XOuter;
            eCalY = CalTrack->YOuter;
            eCalZ = CalTrack->ZOuter;
            eCalEt = CalTrack->PT;

            egCrysEta.push_back(eCalEta);
            egCrysPhi.push_back(eCalPhi);
            egCrysEt.push_back(eCalEt);
            egCrysGx.push_back(eCalX);
            egCrysGy.push_back(eCalY);
            egCrysGz.push_back(eCalZ);
        }

        int TowerN = branchTower->GetEntriesFast();
        for(unsigned int k = 0; k < TowerN; k++)
        {                                 
            CalTower = (Tower*) branchTower->At(k);

            //if( CalTower->Eem == 0 || CalTower->Ehad != 0 ) continue;
            if( CalTower->Eem/CalTower->Ehad < 0.432 ) continue;  // cut for Minbias
            //if( CalTower->Eem/CalTower->Ehad < 0.810 ) continue;  // cut for old Minbias
            //if( CalTower->Eem/CalTower->Ehad > 0.655 ) continue;  // cut for QCD
            double CalX, CalY, CalZ, CalEt, CalEta, CalPhi;
            CalPhi = CalTower->Phi;     
            CalEta = CalTower->Eta;     
            //CalEt = gRandom->Gaus(CalTower->ET, CalTower->ET*0.04644);
            if( fabs(CalEta) < 0.8 ) CalEt = CalTower->ET;
            if( fabs(CalEta) > 0.8 && fabs(CalEta) < 2.1 ) CalEt = gRandom->Gaus(CalTower->ET, CalTower->ET*0.01439);
            if( fabs(CalEta) > 2.1 ) CalEt = gRandom->Gaus(CalTower->ET, CalTower->ET*0.01445);

            //Cal Phi smearing
            if( fabs(CalEta) < 0.8 ){
                if( CalEt > 15. && CalEt < 20. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0026);
                if( CalEt >= 20. && CalEt < 40. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0021);
                if( CalEt >= 40. && CalEt < 60. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0019);
                if( CalEt >= 60. && CalEt < 80. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0017);
                if( CalEt >= 80. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0016);
            }
            if( fabs(CalEta) > 0.8 && fabs(CalEta) < 1.4 ){
                if( CalEt > 15. && CalEt < 20. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0028);
                if( CalEt >= 20. && CalEt < 40. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0024);
                if( CalEt >= 40. && CalEt < 60. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0022);
                if( CalEt >= 60. && CalEt < 80. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0019);
                if( CalEt >= 80. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0017);
            }
            if( fabs(CalEta) > 1.4 && fabs(CalEta) < 1.7 ){
                if( CalEt > 15. && CalEt < 20. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0031);
                if( CalEt >= 20. && CalEt < 40. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0023);
                if( CalEt >= 40. && CalEt < 60. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0022);
                if( CalEt >= 60. && CalEt < 80. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0020);
                if( CalEt >= 80. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0019);
            }
            if( fabs(CalEta) > 1.7 && fabs(CalEta) < 2.1 ){
                if( CalEt > 15. && CalEt < 20. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0037);
                if( CalEt >= 20. && CalEt < 40. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0035);
                if( CalEt >= 40. && CalEt < 60. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0031);
                if( CalEt >= 60. && CalEt < 80.) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0027);
                if( CalEt >= 80. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0024);
            }
            if( fabs(CalEta) > 2.1 && fabs(CalEta) < 2.7 ){
                if( CalEt > 15. && CalEt < 20. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0033);
                if( CalEt >= 20. && CalEt < 40. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0022);
                if( CalEt >= 40. && CalEt < 60. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0018);
                if( CalEt >= 60. && CalEt < 80. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0017);
                if( CalEt >= 80. && CalEt < 95. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0018);
                if( CalEt >= 95. ) CalPhi = gRandom->Gaus(CalPhi, 0.004000-0.0026);
            }
            if( fabs(CalEta) > 2.7 && fabs(CalEta) < 3.0){
                CalPhi = gRandom->Gaus(CalPhi, 0.00445);
            }

            //Cal Eta smearing
            if( fabs(CalEta) < 0.8 ) CalEta = gRandom->Gaus(CalEta, 0.0034);
            if( fabs(CalEta) > 0.8 && fabs(CalEta) < 1.4) CalEta = gRandom->Gaus(CalEta, 0.0037);
            if( fabs(CalEta) > 1.4 && fabs(CalEta) < 1.7) CalEta = gRandom->Gaus(CalEta, 0.0033);
            if( fabs(CalEta) > 1.7 && fabs(CalEta) < 2.1) CalEta = gRandom->Gaus(CalEta, 0.0059);
            if( fabs(CalEta) > 2.1 && fabs(CalEta) < 2.7) CalEta = gRandom->Gaus(CalEta, 0.0059);
            if( fabs(CalEta) > 2.7 && fabs(CalEta) < 3.0) CalEta = gRandom->Gaus(CalEta, 0.0059);

            float distZ = 330.;
            if( fabs(CalEta) < 1.457 )
            {
                float theta = 2.*atan(exp(-CalEta));
                CalX = 137*cos(CalPhi);   // x = r*cos(phi)
                CalY = 137*sin(CalPhi);   // y = r*sin(phi)
                CalZ = 137*tan(TMath::Pi()/2.-theta);     // z = r*tan(theta)
                //cout << eCalEta << ", " << TMath::RadToDeg*theta << ", " << eCalZ << endl; 
            }   
            if( CalEta > 1.457 )
            {   
                float theta = 2.*atan(exp(-CalEta));
                CalZ = 330.;
                float r = CalZ*tan(theta);
                if( r < 35. || r > 155 ) continue;
                CalX = r*cos(CalPhi);
                CalY = r*sin(CalPhi);
            }   
            if( CalEta < -1.457 )
            {   
                float theta = 2.*atan(exp(-CalEta));
                CalZ = -330.;
                float r = CalZ*tan(theta);
                if( r < 35. || r > 155 ) continue;
                CalX = r*cos(CalPhi);
                CalY = r*sin(CalPhi);
            }

            egCrysClusterEta.push_back(CalEta);
            egCrysClusterPhi.push_back(CalPhi);
            egCrysClusterEt.push_back(CalEt);
            egCrysClusterGx.push_back(CalX);
            egCrysClusterGy.push_back(CalY);
            egCrysClusterGz.push_back(CalZ);
        }

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
				k11 = fmod(j11 , 0.026);
				if(fabs(k11) < 0.013) L1X = j11 - k11;
				else{
					if(j11 < 0) L1X = j11 - k11 - 0.026;
					else L1X = j11 - k11 + 0.026;
				}
				//trackL1->YOuter
				j11 = trackL1->YOuter;
				k11 = fmod(j11 , 0.026);
				if(fabs(k11) < 0.013) L1Y = j11 - k11;
				else {
					if(j11 < 0) L1Y = j11 - k11 - 0.026;
					else L1Y = j11 - k11 + 0.026;
				}
				//trackL1->ZOuter
				j11 = trackL1->ZOuter;
				k11 = fmod(j11 , 0.048);
				if(fabs(k11) < 0.024) L1Z = j11 - k11;
				else{
					if(j11 < 0) L1Z = j11 - k11 - 0.048;
					else L1Z = j11 - k11 + 0.048;
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
				k11 = fmod(j11 , 0.026);
				if(fabs(k11) < 0.013) L1X = j11 - k11;
				else{
					if(j11 < 0) L1X = j11 - k11 - 0.026;
					else L1X = j11 - k11 + 0.026;
				}
				//trackL2->YOuter
				j11 = trackL2->YOuter;
				k11 = fmod(j11 , 0.026);
				if(fabs(k11) < 0.013) L1Y = j11 - k11;
				else{
					if(j11 < 0) L1Y = j11 - k11 - 0.026;
					else L1Y = j11 - k11 + 0.026;
				}
				//trackL2->ZOuter
				j11 = trackL2->ZOuter;
				k11 = fmod(j11 , 0.048);
				if(fabs(k11) < 0.024) L1Z = j11 - k11;
				else{
					if(j11 < 0) L1Z = j11 - k11 - 0.048;
					else L1Z = j11 - k11 + 0.048;
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
				k11 = fmod(j11 , 0.026);
				if(fabs(k11) < 0.013) L1X = j11 - k11;
				else{
					if(j11 < 0) L1X = j11 - k11 - 0.026;
					else L1X = j11 - k11 + 0.026;
				}
				//trackL3->YOuter
				j11 = trackL3->YOuter;
				k11 = fmod(j11 , 0.026);
				if(fabs(k11) < 0.013) L1Y = j11 - k11;
				else{
					if(j11 < 0) L1Y = j11 - k11 - 0.026;
					else L1Y = j11 - k11 + 0.026;
				}
				//trackL3->ZOuter
				j11 = trackL3->ZOuter;
				k11 = fmod(j11 , 0.048);
				if(fabs(k11) < 0.024) L1Z = j11 - k11;
				else{
					if(j11 < 0) L1Z = j11 - k11 - 0.048;
					else L1Z = j11 - k11 + 0.048;
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
				k11 = fmod(j11 , 0.026);
				if(fabs(k11) < 0.013) L1X = j11 - k11;
				else{
					if(j11 < 0) L1X = j11 - k11 - 0.026;
					else L1X = j11 - k11 + 0.026;
				}
				//trackL4->YOuter
				j11 = trackL4->YOuter;
				k11 = fmod(j11 , 0.026);
				if(fabs(k11) < 0.013) L1Y = j11 - k11;
				else{
					if(j11 < 0) L1Y = j11 - k11 - 0.026;
					else L1Y = j11 - k11 + 0.026;
				}
				//trackL4->ZOuter
				j11 = trackL4->ZOuter;
				k11 = fmod(j11 , 0.048);
				if(fabs(k11) < 0.024) L1Z = j11 - k11;
				else{
					if(j11 < 0) L1Z = j11 - k11 - 0.048;
					else L1Z = j11 - k11 + 0.048;
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
				k11 = fmod(j11 , 0.015);
				if(fabs(k11) < 0.0075) D1Z = j11 - k11;
				else{
					if(j11 < 0) D1Z = j11 - k11 - 0.015;
					else D1Z = j11 - k11 + 0.015;
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
				k11 = fmod(j11 , 0.015);
				if(fabs(k11) < 0.0075) D2Z = j11 - k11;
				else{
					if(j11 < 0) D2Z = j11 - k11 - 0.015;
					else D2Z = j11 - k11 + 0.015;
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
				k11 = fmod(j11 , 0.015);
				if(fabs(k11) < 0.0075) D3Z = j11 - k11;
				else{
					if(j11 < 0) D3Z = j11 - k11 - 0.015;
					else D3Z = j11 - k11 + 0.015;
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
				k11 = fmod(j11 , 0.015);
				if(fabs(k11) < 0.0075) D4Z = j11 - k11;
				else{
					if(j11 < 0) D4Z = j11 - k11 - 0.015;
					else D4Z = j11 - k11 + 0.015;
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
				k11 = fmod(j11 , 0.015);
				if(fabs(k11) < 0.0075) D5Z = j11 - k11;
				else{
					if(j11 < 0) D5Z = j11 - k11 - 0.015;
					else D5Z = j11 - k11 + 0.015;
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
