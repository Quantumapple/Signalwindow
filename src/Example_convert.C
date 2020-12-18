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
void calorimeterPosition(double phi, double eta, double e, double &x, double &y, double &z);

//------------------------------------------------------------------------------

void Example_convert(const char *inputFile)
{
	gSystem->Load("libDelphes");

	TChain *chain = new TChain("Delphes");
	chain->Add(inputFile);

	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

    TFile *result = new TFile("result.root","RECREATE");
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
    std::vector<int>   bRecHitLayer;
    std::vector<float> bRecHitGx; 
    std::vector<float> bRecHitGy; 
    std::vector<float> bRecHitGz;
    
    int fRecHitN = 0;
    std::vector<int>   fRecHitDisk;
    std::vector<float> fRecHitGx; 
    std::vector<float> fRecHitGy; 
    std::vector<float> fRecHitGz; 

    std::vector<float> propgenElPartPhi;
    std::vector<float> propgenElPartEta;
    std::vector<float> propgenElPartPt;
    std::vector<float> propgenElPartX;
    std::vector<float> propgenElPartY;
    std::vector<float> propgenElPartZ;

    std::vector<float> egCrysClusterEem;
    std::vector<float> egCrysClusterEhad;
    std::vector<float> egCrysClusterEnergy;
    std::vector<float> egCrysClusterEt;
    std::vector<float> egCrysClusterEta;
    std::vector<float> egCrysClusterPhi;
    std::vector<float> egCrysClusterGx;
    std::vector<float> egCrysClusterGy;
    std::vector<float> egCrysClusterGz;
    
    tt->Branch("event", &event, "event/I");
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

    tt->Branch("egCrysClusterEem", &egCrysClusterEem);
    tt->Branch("egCrysClusterEhad", &egCrysClusterEhad);
    tt->Branch("egCrysClusterEnergy", &egCrysClusterEnergy);
    tt->Branch("egCrysClusterEt", &egCrysClusterEt);
    tt->Branch("egCrysClusterEta", &egCrysClusterEta);
    tt->Branch("egCrysClusterPhi", &egCrysClusterPhi);
    tt->Branch("egCrysClusterGx",&egCrysClusterGx);
    tt->Branch("egCrysClusterGy",&egCrysClusterGy);
    tt->Branch("egCrysClusterGz",&egCrysClusterGz);
    
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchTower = treeReader->UseBranch("Tower");
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
    
    Tower *CalTower;
    //Track *CalTrack;

    Long64_t entry;

    bool kFlag = false;

    int inf_count=0;

    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {
        event = 1;
        if( entry % 100000 == 0 ) cout << "Event: " << entry << endl; 
        
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

        propgenElPartPt.clear();
        propgenElPartPhi.clear();
        propgenElPartEta.clear();
        propgenElPartX.clear();
        propgenElPartY.clear();
        propgenElPartZ.clear();
    
        egCrysClusterEem.clear();
        egCrysClusterEhad.clear();
        egCrysClusterEnergy.clear();
        egCrysClusterEt.clear();
        egCrysClusterEta.clear();
        egCrysClusterPhi.clear();
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
            double j11 = 0;
            double genX, genY, genZ;
            //particle->X
            genX = particle->X;
            //particle->Y
            genY = particle->Y;
            //particle->Z
            genZ = particle->Z;

            propgenElPartPt.push_back(particle->PT);
            propgenElPartEta.push_back(particle->Eta);
            propgenElPartPhi.push_back(particle->Phi);
            propgenElPartX.push_back(genX/10.);
            propgenElPartY.push_back(genY/10.);
            propgenElPartZ.push_back(genZ/10.);
        }

        int TowerN = branchTower->GetEntriesFast();
        for(unsigned int k = 0; k < TowerN; k++)
        {                                 
            CalTower = (Tower*) branchTower->At(k);
 
            double CalX, CalY, CalZ, CalEem, CalEhad, CalEt, CalEta, CalPhi, CalE;
            CalE = CalTower->E;
            CalEem = CalTower->Eem;
            CalEhad = CalTower->Ehad;
            CalPhi = CalTower->Phi;     
            CalEta = CalTower->Eta;     
            
            //calorimeterPosition(CalPhi, CalEta, CalE, CalX, CalY, CalZ);
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

            egCrysClusterEem.push_back(CalEem);
            egCrysClusterEhad.push_back(CalEhad);
            egCrysClusterEnergy.push_back(CalE);
            egCrysClusterEt.push_back(CalEt);
            egCrysClusterEta.push_back(CalEta);
            egCrysClusterPhi.push_back(CalPhi);
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
				double L1X, L1Y, L1Z;
				//trackL1->XOuter
				L1X = trackL1->XOuter;
				//trackL1->YOuter
				L1Y = trackL1->YOuter;
				//trackL1->ZOuter
				L1Z = trackL1->ZOuter;

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
				double L1X, L1Y, L1Z;
				//trackL2->XOuter
				L1X = trackL2->XOuter;
				//trackL2->YOuter
				L1Y = trackL2->YOuter;
				//trackL2->ZOuter
				L1Z = trackL2->ZOuter;

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
				double L1X, L1Y, L1Z;
				//trackL3->XOuter
				L1X = trackL3->XOuter;
				//trackL3->YOuter
				L1Y = trackL3->YOuter;
				//trackL3->ZOuter
				L1Z = trackL3->ZOuter;

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
				double L1X, L1Y, L1Z;
				//trackL4->XOuter
				L1X = trackL4->XOuter;
				//trackL4->YOuter
				L1Y = trackL4->YOuter;
				//trackL4->ZOuter
				L1Z = trackL4->ZOuter;

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
			if( fabs(trackD1->ZOuter) > 249.9 && R > 28.0 )
			{
				double D1X, D1Y, D1Z;
				//trackD1->XOuter
                D1X = trackD1->XOuter;
				//trackD1->YOuter
                D1Y = trackD1->YOuter;
				//trackD1->ZOuter
                D1Z = trackD1->ZOuter;

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
			if( fabs(trackD2->ZOuter) > 319.9 && R > 28.0 )
			{
				double D2X, D2Y, D2Z;
				//trackD2->XOuter
                D2X = trackD2->XOuter;
				//trackD2->YOuter
                D2Y = trackD2->YOuter;
				//trackD2->ZOuter
                D2Z = trackD2->ZOuter;

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
			if( fabs(trackD3->ZOuter) > 408.9 && R > 28.0 )
			{
				double D3X, D3Y, D3Z;
				//trackD3->XOuter
                D3X = trackD3->XOuter;
				//trackD3->YOuter
                D3Y = trackD3->YOuter;
				//trackD3->ZOuter
                D3Z = trackD3->ZOuter;

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
			if( fabs(trackD4->ZOuter) > 522.9 && R > 28.0 )
			{
				double D4X, D4Y, D4Z;
				//trackD4->XOuter
                D4X = trackD4->XOuter;
				//trackD4->YOuter
                D4Y = trackD4->YOuter;
				//trackD4->ZOuter
                D4Z = trackD4->ZOuter;

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
			if( fabs(trackD5->ZOuter) > 668.9 && R > 28.0 )
			{
				double D5X, D5Y, D5Z;
				//trackD5->XOuter
                D5X = trackD5->XOuter;
				//trackD5->YOuter
                D5Y = trackD5->YOuter;
				//trackD5->ZOuter
                D5Z = trackD5->ZOuter;

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

// -------------- get Calorimeter position
void calorimeterPosition(double phi, double eta, double e, double &x, double &y, double &z) {
    double depth = 0.89*(7.7+ log(e) );
    double theta = 2*atan(exp(-1*eta));
    double r = 0;
    if( fabs(eta) > 1.479 )
    {
        double ecalZ = 315.4*fabs(eta)/eta;

        r = ecalZ / cos( 2*atan( exp( -1*eta ) ) ) + depth;
        x = r * cos( phi ) * sin( theta );
        y = r * sin( phi ) * sin( theta );
        z = r * cos( theta );
    }
    else
    {
        double rperp = 129.0;
        double zface =  sqrt( cos( theta ) * cos( theta ) /
                ( 1 - cos( theta ) * cos( theta ) ) *
                rperp * rperp ) * fabs( eta ) / eta;
        r = sqrt( rperp * rperp + zface * zface ) + depth;
        x = r * cos( phi ) * sin( theta );
        y = r * sin( phi ) * sin( theta );
        z = r * cos( theta );
    }
}

//------------------------------------------------------------------------------
