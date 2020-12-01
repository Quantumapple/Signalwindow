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

void Example36(const char *inputFile)
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
    int tmPevent = 0;

    int tmPbRecHitN = 0;
    std::vector<int>   tmPbRecHitLayer;
    std::vector<float> tmPbRecHitGx; 
    std::vector<float> tmPbRecHitGy; 
    std::vector<float> tmPbRecHitGz;
    
    int tmPfRecHitN = 0;
    std::vector<int>   tmPfRecHitDisk;
    std::vector<float> tmPfRecHitGx; 
    std::vector<float> tmPfRecHitGy; 
    std::vector<float> tmPfRecHitGz; 

    std::vector<float> tmPpropgenElPartPhi;
    std::vector<float> tmPpropgenElPartEta;
    std::vector<float> tmPpropgenElPartPt;
    std::vector<float> tmPpropgenElPartX;
    std::vector<float> tmPpropgenElPartY;
    std::vector<float> tmPpropgenElPartZ;

    std::vector<float> tmPegCrysClusterEem;
    std::vector<float> tmPegCrysClusterEhad;
    std::vector<float> tmPegCrysClusterEnergy;
    std::vector<float> tmPegCrysClusterEt;
    std::vector<float> tmPegCrysClusterEta;
    std::vector<float> tmPegCrysClusterPhi;
    std::vector<float> tmPegCrysClusterGx;
    std::vector<float> tmPegCrysClusterGy;
    std::vector<float> tmPegCrysClusterGz;
    
    tt->Branch("tmPevent", &tmPevent, "tmPevent/I");
    tt->Branch("tmPpropgenElPartPhi", &tmPpropgenElPartPhi); 
    tt->Branch("tmPpropgenElPartEta", &tmPpropgenElPartEta); 
    tt->Branch("tmPpropgenElPartPt", &tmPpropgenElPartPt);
    tt->Branch("tmPpropgenElPartX", &tmPpropgenElPartX);
    tt->Branch("tmPpropgenElPartY", &tmPpropgenElPartY);
    tt->Branch("tmPpropgenElPartZ", &tmPpropgenElPartZ);
    
    tt->Branch("tmPbRecHitN", &tmPbRecHitN, "tmPbRecHitN/I");
    tt->Branch("tmPbRecHitLayer", &tmPbRecHitLayer);
    tt->Branch("tmPbRecHitGx", &tmPbRecHitGx);
    tt->Branch("tmPbRecHitGy", &tmPbRecHitGy);
    tt->Branch("tmPbRecHitGz", &tmPbRecHitGz);
    
    tt->Branch("tmPfRecHitN", &tmPfRecHitN, "tmPfRecHitN/I");
    tt->Branch("tmPfRecHitDisk", &tmPfRecHitDisk);
    tt->Branch("tmPfRecHitGx", &tmPfRecHitGx);
    tt->Branch("tmPfRecHitGy", &tmPfRecHitGy);
    tt->Branch("tmPfRecHitGz", &tmPfRecHitGz);

    tt->Branch("tmPegCrysClusterEem", &tmPegCrysClusterEem);
    tt->Branch("tmPegCrysClusterEhad", &tmPegCrysClusterEhad);
    tt->Branch("tmPegCrysClusterEnergy", &tmPegCrysClusterEnergy);
    tt->Branch("tmPegCrysClusterEt", &tmPegCrysClusterEt);
    tt->Branch("tmPegCrysClusterEta", &tmPegCrysClusterEta);
    tt->Branch("tmPegCrysClusterPhi", &tmPegCrysClusterPhi);
    tt->Branch("tmPegCrysClusterGx",&tmPegCrysClusterGx);
    tt->Branch("tmPegCrysClusterGy",&tmPegCrysClusterGy);
    tt->Branch("tmPegCrysClusterGz",&tmPegCrysClusterGz);
    
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
        tmPevent = 1;
        if( entry % 100000 == 0 ) cout << "Event: " << entry << endl; 
        
        tmPbRecHitN = 0;
        tmPbRecHitLayer.clear();
        tmPbRecHitGx.clear();
        tmPbRecHitGy.clear();
        tmPbRecHitGz.clear();
       
        tmPfRecHitN = 0;
        tmPfRecHitDisk.clear();
        tmPfRecHitGx.clear(); 
        tmPfRecHitGy.clear(); 
        tmPfRecHitGz.clear();

        tmPpropgenElPartPt.clear();
        tmPpropgenElPartPhi.clear();
        tmPpropgenElPartEta.clear();
        tmPpropgenElPartX.clear();
        tmPpropgenElPartY.clear();
        tmPpropgenElPartZ.clear();
    
        tmPegCrysClusterEem.clear();
        tmPegCrysClusterEhad.clear();
        tmPegCrysClusterEnergy.clear();
        tmPegCrysClusterEt.clear();
        tmPegCrysClusterEta.clear();
        tmPegCrysClusterPhi.clear();
        tmPegCrysClusterGx.clear();
        tmPegCrysClusterGy.clear();
        tmPegCrysClusterGz.clear();
        
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
            j11 = particle->X;
            genX = gRandom->Gaus(j11,0.0078);
            //particle->Y
            j11 = particle->Y;
            genY = gRandom->Gaus(j11,0.0069);
            //particle->Z
            j11 = particle->Z;
            genZ = gRandom->Gaus(j11,0.0167);

            tmPpropgenElPartPt.push_back(particle->PT);
            tmPpropgenElPartEta.push_back(particle->Eta);
            tmPpropgenElPartPhi.push_back(particle->Phi);
            tmPpropgenElPartX.push_back(genX/10.);
            tmPpropgenElPartY.push_back(genY/10.);
            tmPpropgenElPartZ.push_back(genZ/10.);
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
            
            if( fabs(CalEta) < 0.8 ) CalEt = CalTower->ET;
            if( fabs(CalEta) > 0.8 && fabs(CalEta) < 2.1 ) CalEt = gRandom->Gaus(CalTower->ET, CalTower->ET*0.01439);
            if( fabs(CalEta) > 2.1 ) CalEt = gRandom->Gaus(CalTower->ET, CalTower->ET*0.01445);

            
            double factor_r1 = 2.85;
            double factor_r2 = 3.85;
            //Cal Phi smearing
            if( fabs(CalEta) < 0.8 ){
                //if( CalEt <= 15. ) CalPhi = gRandom->Gaus(CalPhi, (0.004000-0.0026)*factor_r1);
                //if( CalEt > 15. && CalEt < 20. ) CalPhi = gRandom->Gaus(CalPhi, (0.004000-0.0026)*factor_r1);
                if( CalEt < 20. ) CalPhi = gRandom->Gaus(CalPhi, (0.004000-0.0026)*factor_r1);
                if( CalEt >= 20. && CalEt < 40. ) CalPhi = gRandom->Gaus(CalPhi,(0.004000-0.0021)*factor_r1);
                if( CalEt >= 40. && CalEt < 60. ) CalPhi = gRandom->Gaus(CalPhi,(0.004000-0.0019)*factor_r1);
                if( CalEt >= 60. && CalEt < 80. ) CalPhi = gRandom->Gaus(CalPhi,(0.004000-0.0017)*factor_r1);
                if( CalEt >= 80. ) CalPhi = gRandom->Gaus(CalPhi, (0.004000-0.0016)*factor_r1);
            }
            if( fabs(CalEta) > 0.8 && fabs(CalEta) < 1.4 ){
                //if( CalEt > 15. && CalEt < 20. ) CalPhi = gRandom->Gaus(CalPhi, (0.004000-0.0028)*factor_r2);
                if( CalEt < 20. ) CalPhi = gRandom->Gaus(CalPhi, (0.004000-0.0028)*factor_r2);
                if( CalEt >= 20. && CalEt < 40. ) CalPhi = gRandom->Gaus(CalPhi, (0.004000-0.0024)*factor_r2);
                if( CalEt >= 40. && CalEt < 60. ) CalPhi = gRandom->Gaus(CalPhi, (0.004000-0.0022)*factor_r2);
                if( CalEt >= 60. && CalEt < 80. ) CalPhi = gRandom->Gaus(CalPhi, (0.004000-0.0019)*factor_r2);
                if( CalEt >= 80. ) CalPhi = gRandom->Gaus(CalPhi, (0.004000-0.0017)*factor_r2);
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
                if( CalEt >  15. && CalEt < 20. ) CalPhi = gRandom->Gaus(CalPhi, 0.003000-0.0033);
                if( CalEt >= 20. && CalEt < 40. ) CalPhi = gRandom->Gaus(CalPhi, 0.003000-0.0022);
                if( CalEt >= 40. && CalEt < 60. ) CalPhi = gRandom->Gaus(CalPhi, 0.003000-0.0018);
                if( CalEt >= 60. && CalEt < 80. ) CalPhi = gRandom->Gaus(CalPhi, 0.003000-0.0017);
                if( CalEt >= 80. ) CalPhi = gRandom->Gaus(CalPhi, 0.003000-0.0018);
                //if( CalEt >= 80. && CalEt < 95. ) CalPhi = gRandom->Gaus(CalPhi, 0.003000-0.0018);
                //if( CalEt >= 95. ) CalPhi = gRandom->Gaus(CalPhi, 0.003000-0.0026);
            }
            if( fabs(CalEta) > 2.7 && fabs(CalEta) < 3.0){
                CalPhi = gRandom->Gaus(CalPhi, 0.00345);
            }

            //Cal Eta smearing
            //if( fabs(CalEta) < 0.8 ) CalEta = gRandom->Gaus(CalEta, 0.003405); //vX
            if( fabs(CalEta) < 0.8 ) CalEta = gRandom->Gaus(CalEta, 0.001975); //v1
            //if( fabs(CalEta) > 0.8 && fabs(CalEta) < 1.4) CalEta = gRandom->Gaus(CalEta, 0.003392); //v3
            if( fabs(CalEta) > 0.8 && fabs(CalEta) < 1.4) CalEta = gRandom->Gaus(CalEta, 0.002992); //v3
            //if( fabs(CalEta) > 1.4 && fabs(CalEta) < 1.7) CalEta = gRandom->Gaus(CalEta, 0.003551); //v0
            //if( fabs(CalEta) > 1.4 && fabs(CalEta) < 1.7) CalEta = gRandom->Gaus(CalEta, 0.003511); //v3
            if( fabs(CalEta) > 1.4 && fabs(CalEta) < 1.7) CalEta = gRandom->Gaus(CalEta, 0.003491); //v4
            //if( fabs(CalEta) > 1.7 && fabs(CalEta) < 2.1) CalEta = gRandom->Gaus(CalEta, 0.001239); //v0 
            if( fabs(CalEta) > 1.7 && fabs(CalEta) < 2.1) CalEta = gRandom->Gaus(CalEta, 0.002239); 
            //if( fabs(CalEta) > 2.1 && fabs(CalEta) < 2.7) CalEta = gRandom->Gaus(CalEta, 0.002144); // v0
            if( fabs(CalEta) > 2.1 && fabs(CalEta) < 2.7) CalEta = gRandom->Gaus(CalEta, 0.002044); 
            //if( fabs(CalEta) > 2.7 && fabs(CalEta) < 3.0) CalEta = gRandom->Gaus(CalEta, 0.003834); // v0 
            //if( fabs(CalEta) > 2.7 && fabs(CalEta) < 3.0) CalEta = gRandom->Gaus(CalEta, 0.003334); 
            if( fabs(CalEta) > 2.7 && fabs(CalEta) < 3.0) CalEta = gRandom->Gaus(CalEta, 0.003384); 

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

            tmPegCrysClusterEem.push_back(CalEem);
            tmPegCrysClusterEhad.push_back(CalEhad);
            tmPegCrysClusterEnergy.push_back(CalE);
            tmPegCrysClusterEt.push_back(CalEt);
            tmPegCrysClusterEta.push_back(CalEta);
            tmPegCrysClusterPhi.push_back(CalPhi);
            tmPegCrysClusterGx.push_back(CalX);
            tmPegCrysClusterGy.push_back(CalY);
            tmPegCrysClusterGz.push_back(CalZ);
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
                L1X = gRandom->Gaus(j11,0.0073);
				//trackL1->YOuter
				j11 = trackL1->YOuter;
                L1Y = gRandom->Gaus(j11,0.0073);
				//trackL1->ZOuter
				j11 = trackL1->ZOuter;
                L1Z = gRandom->Gaus(j11,0.0073);

				tmPbRecHitLayer.push_back(1);
				tmPbRecHitGx.push_back(L1X/10.);
				tmPbRecHitGy.push_back(L1Y/10.);
				tmPbRecHitGz.push_back(L1Z/10.);
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
                L1X = gRandom->Gaus(j11,0.0073);
				//trackL2->YOuter
				j11 = trackL2->YOuter;
                L1Y = gRandom->Gaus(j11,0.0073);
				//trackL2->ZOuter
				j11 = trackL2->ZOuter;
                L1Z = gRandom->Gaus(j11,0.0083);

				tmPbRecHitLayer.push_back(2);
				tmPbRecHitGx.push_back(L1X/10.);
				tmPbRecHitGy.push_back(L1Y/10.);
				tmPbRecHitGz.push_back(L1Z/10.);
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
                L1X = gRandom->Gaus(j11,0.0073);
				//trackL3->YOuter
				j11 = trackL3->YOuter;
                L1Y = gRandom->Gaus(j11,0.0073);
				//trackL3->ZOuter
				j11 = trackL3->ZOuter;
                L1Z = gRandom->Gaus(j11,0.0093);

				tmPbRecHitLayer.push_back(3);
				tmPbRecHitGx.push_back(L1X/10.);
				tmPbRecHitGy.push_back(L1Y/10.);
				tmPbRecHitGz.push_back(L1Z/10.);
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
                L1X = gRandom->Gaus(j11,0.0073);
				//trackL4->YOuter
				j11 = trackL4->YOuter;
                L1Y = gRandom->Gaus(j11,0.0073);
				//trackL4->ZOuter
				j11 = trackL4->ZOuter;
                L1Z = gRandom->Gaus(j11,0.0103);

				tmPbRecHitLayer.push_back(4);
				tmPbRecHitGx.push_back(L1X/10.);
				tmPbRecHitGy.push_back(L1Y/10.);
				tmPbRecHitGz.push_back(L1Z/10.);
			}
		}
		tmPbRecHitN = tmPbRecHitGx.size();

		// Loop over 1st pixel Disk
		for(unsigned int k = 0; k < branchPixelD1->GetEntriesFast(); k++)
		{
			trackD1 = (Track*) branchPixelD1->At(k);
			float R = sqrt(pow(trackD1->XOuter,2)+pow(trackD1->YOuter,2));
			if( fabs(trackD1->ZOuter) > 249.9 && R > 28.0 )
                //if( fabs(trackD1->ZOuter) > 248 )
			{
				double D1X, D1Y, D1Z;
				//trackD1->XOuter
                D1X = trackD1->XOuter;
				//trackD1->YOuter
                D1Y = trackD1->YOuter;
				//trackD1->ZOuter
                D1Z = trackD1->ZOuter;

				tmPfRecHitDisk.push_back(1);
				tmPfRecHitGx.push_back(D1X/10.);
				tmPfRecHitGy.push_back(D1Y/10.);
				tmPfRecHitGz.push_back(D1Z/10.);
			}
		}

		// Loop over 2nd pixel Disk
		for(unsigned int k = 0; k < branchPixelD2->GetEntriesFast(); k++)
		{
			trackD2 = (Track*) branchPixelD2->At(k);
			float R = sqrt(pow(trackD2->XOuter,2)+pow(trackD2->YOuter,2));
			if( fabs(trackD2->ZOuter) > 319.9 && R > 28.0 )
                //if( fabs(trackD2->ZOuter) > 315 )
			{
				double D2X, D2Y, D2Z;
				//trackD2->XOuter
                D2X = trackD2->XOuter;
				//trackD2->YOuter
                D2Y = trackD2->YOuter;
				//trackD2->ZOuter
                D2Z = trackD2->ZOuter;

				tmPfRecHitDisk.push_back(2);
				tmPfRecHitGx.push_back(D2X/10.);
				tmPfRecHitGy.push_back(D2Y/10.);
				tmPfRecHitGz.push_back(D2Z/10.);
            }
		}

		// Loop over 3rd pixel Disk
		for(unsigned int k = 0; k < branchPixelD3->GetEntriesFast(); k++)
		{
			trackD3 = (Track*) branchPixelD3->At(k);
			float R = sqrt(pow(trackD3->XOuter,2)+pow(trackD3->YOuter,2));
			if( fabs(trackD3->ZOuter) > 408.9 && R > 28.0 )
				//if( fabs(trackD3->ZOuter) > 407 )
			{
				double D3X, D3Y, D3Z;
				//trackD3->XOuter
                D3X = trackD3->XOuter;
				//trackD3->YOuter
                D3Y = trackD3->YOuter;
				//trackD3->ZOuter
                D3Z = trackD3->ZOuter;

				tmPfRecHitDisk.push_back(3);
				tmPfRecHitGx.push_back(D3X/10.);
				tmPfRecHitGy.push_back(D3Y/10.);
				tmPfRecHitGz.push_back(D3Z/10.);
			}
		}

		// Loop over 4th pixel Disk
		for(unsigned int k = 0; k < branchPixelD4->GetEntriesFast(); k++)
		{
			trackD4 = (Track*) branchPixelD4->At(k);
			float R = sqrt(pow(trackD4->XOuter,2)+pow(trackD4->YOuter,2));
			if( fabs(trackD4->ZOuter) > 522.9 && R > 28.0 )
				//if( fabs(trackD4->ZOuter) > 521 )
			{
				double D4X, D4Y, D4Z;
				//trackD4->XOuter
                D4X = trackD4->XOuter;
				//trackD4->YOuter
                D4Y = trackD4->YOuter;
				//trackD4->ZOuter
                D4Z = trackD4->ZOuter;

				tmPfRecHitDisk.push_back(4);
				tmPfRecHitGx.push_back(D4X/10.);
				tmPfRecHitGy.push_back(D4Y/10.);
				tmPfRecHitGz.push_back(D4Z/10.);
			}
		}

		// Loop over 5th pixel Disk
		for(unsigned int k = 0; k < branchPixelD5->GetEntriesFast(); k++)
		{
			trackD5 = (Track*) branchPixelD5->At(k);
			float R = sqrt(pow(trackD5->XOuter,2)+pow(trackD5->YOuter,2));
			if( fabs(trackD5->ZOuter) > 668.9 && R > 28.0 )
				//if( fabs(trackD5->ZOuter) > 667 )
			{
				double D5X, D5Y, D5Z;
				//trackD5->XOuter
                D5X = trackD5->XOuter;
				//trackD5->YOuter
                D5Y = trackD5->YOuter;
				//trackD5->ZOuter
                D5Z = trackD5->ZOuter;

				tmPfRecHitDisk.push_back(5);
				tmPfRecHitGx.push_back(D5X/10.);
				tmPfRecHitGy.push_back(D5Y/10.);
				tmPfRecHitGz.push_back(D5Z/10.);
			}
		}
		tmPfRecHitN = tmPfRecHitGx.size();

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
