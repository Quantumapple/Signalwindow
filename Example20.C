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

double poly(double egEt, double egEta)
{
    double tempEt, over;
    over = 100;
    tempEt = egEt;
    if( fabs(egEta < 0.8) ){
        if( tempEt < 100 ){
            egEt = (0.0106055)*pow(egEt,0)+(-0.007118)*pow(egEt,1)+(0.000153065)*pow(egEt,2)+(5.42911e-08)*pow(egEt,3)+(-1.14254e-08)*pow(egEt,4)+(-9.67285e-11)*pow(egEt,5)+(-1.27576e-13)*pow(egEt,6)+(6.94584e-15)*pow(egEt,7)+(1.06681e-16)*pow(egEt,8)+(8.55015e-19)*pow(egEt,9)+(2.53734e-21)*pow(egEt,10)+(-2.74998e-23)*pow(egEt,11)+(-7.34961e-25)*pow(egEt,12)+(-1.09489e-26)*pow(egEt,13)+(-8.13622e-29)*pow(egEt,14)+(-7.04202e-33)*pow(egEt,15)+(1.3781e-32)*pow(egEt,16);
        }
        if( tempEt >= 100 ){
            egEt = (0.0106055)*pow(over,0)+(-0.007118)*pow(over,1)+(0.000153065)*pow(over,2)+(5.42911e-08)*pow(over,3)+(-1.14254e-08)*pow(over,4)+(-9.67285e-11)*pow(over,5)+(-1.27576e-13)*pow(over,6)+(6.94584e-15)*pow(over,7)+(1.06681e-16)*pow(over,8)+(8.55015e-19)*pow(over,9)+(2.53734e-21)*pow(over,10)+(-2.74998e-23)*pow(over,11)+(-7.34961e-25)*pow(over,12)+(-1.09489e-26)*pow(over,13)+(-8.13622e-29)*pow(over,14)+(-7.04202e-33)*pow(over,15)+(1.3781e-32)*pow(over,16);
        }
    }
    if( fabs(egEta) >= 0.8 && fabs(egEta) < 1.4 ){
        if( tempEt < 100 ){
            egEt = (0.0140984)+(-0.00475802)*egEt+(8.27489e-05)*pow(egEt,2)+(2.27419e-07)*pow(egEt,3)+(-6.84511e-09)*pow(egEt,4)+(-6.28051e-11)*pow(egEt,5)+(6.97593e-13)*pow(egEt,6);
        }
        if( tempEt >= 100 ){
            egEt = (0.0140984)+(-0.00475802)*over+(8.27489e-05)*pow(over,2)+(2.27419e-07)*pow(over,3)+(-6.84511e-09)*pow(over,4)+(-6.28051e-11)*pow(over,5)+(6.97593e-13)*pow(over,6);
        }
    }
    if( fabs(egEta) >= 1.4 && fabs(egEta) < 1.7 ){
        if( tempEt < 100 ){
            egEt = (0.460002)+(-0.0190628)*egEt+(0.00021923)*pow(egEt,2)+(1.0012e-06)*pow(egEt,3)+(-1.61066e-08)*pow(egEt,4)+(-1.80704e-10)*pow(egEt,5)+(1.66094e-12)*pow(egEt,6);
        }
        if( tempEt >= 100 ){
            egEt = (0.460002)+(-0.0190628)*over+(0.00021923)*pow(over,2)+(1.0012e-06)*pow(over,3)+(-1.61066e-08)*pow(over,4)+(-1.80704e-10)*pow(over,5)+(1.66094e-12)*pow(over,6);
        }
    }
    if( fabs(egEta) >= 1.7 && fabs(egEta) < 2.1 ){
        if( tempEt < 100 ){
            egEt = (0.279437)+(-0.00985193)*egEt+(0.000113354)*pow(egEt,2)+(4.95472e-07)*pow(egEt,3)+(-7.90583e-09)*pow(egEt,4)+(-8.71698e-11)*pow(egEt,5)+(7.39395e-13)*pow(egEt,6);
        }
        if( tempEt >= 100 ){
            egEt = (0.279437)+(-0.00985193)*over+(0.000113354)*pow(over,2)+(4.95472e-07)*pow(over,3)+(-7.90583e-09)*pow(over,4)+(-8.71698e-11)*pow(over,5)+(7.39395e-13)*pow(over,6);
        }
    }
    if( fabs(egEta) >= 2.1 && fabs(egEta) < 2.7 ){
        if( tempEt < 100 ){
            egEt = (0.139995)+(-0.00211373)*egEt+(-1.85485e-06)*pow(egEt,2)+(1.73091e-07)*pow(egEt,3)+(1.79371e-09)*pow(egEt,4)+(1.24589e-12)*pow(egEt,5)+(-2.78562e-13)*pow(egEt,6);
        }
        if( tempEt >= 100){
            egEt = (0.139995)+(-0.00211373)*over+(-1.85485e-06)*pow(over,2)+(1.73091e-07)*pow(over,3)+(1.79371e-09)*pow(over,4)+(1.24589e-12)*pow(over,5)+(-2.78562e-13)*pow(over,6);
        }
    }
    if( fabs(egEta) >= 2.7 && fabs(egEta) < 3.0 ){
        if( tempEt < 100 ){
            egEt = (0.12824)+(-0.00186523)*egEt+(-1.66273e-06)*pow(egEt,2)+(1.51476e-07)*pow(egEt,3)+(1.56831e-09)*pow(egEt,4)+(1.07672e-12)*pow(egEt,5)+(-2.44234e-13)*pow(egEt,6);
        }
        if( tempEt >= 100){
            egEt = (0.12824)+(-0.00186523)*over+(-1.66273e-06)*pow(over,2)+(1.51476e-07)*pow(over,3)+(1.56831e-09)*pow(over,4)+(1.07672e-12)*pow(over,5)+(-2.44234e-13)*pow(over,6);
        }
    }
    
    return tempEt+(egEt*tempEt);
}

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
    //TClonesArray *branchTrack = treeReader->UseBranch("EFlowTrack");
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
    //Track *eCalTrack;

    TLorentzVector egVector;

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

        int TrackN = branchTower->GetEntriesFast();
        for(unsigned int k = 0; k < TrackN; k++)
        {                                 
            eCalTower = (Tower*) branchTower->At(k);
            double eCalX, eCalY, eCalZ, eCalEt, eCalEta, eCalPhi;

            //ECal tower phi                
            eCalPhi = eCalTower->Phi;     

            //ECal tower eta                
            eCalEta = eCalTower->Eta;     

            //ECal tower Et                 
            eCalEt = eCalTower->ET;
            eCalEt = poly(eCalEt, eCalEta);

            //ECal Phi smearing
            if( fabs(eCalEta) < 0.8 ){
                if( eCalEt > 15. && eCalEt < 20. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0026);
                if( eCalEt >= 20. && eCalEt < 40. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0021);
                if( eCalEt >= 40. && eCalEt < 60. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0019);
                if( eCalEt >= 60. && eCalEt < 80. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0017);
                if( eCalEt >= 80. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0016);
            }
            if( fabs(eCalEta) > 0.8 && fabs(eCalEta) < 1.4 ){
                if( eCalEt > 15. && eCalEt < 20. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0028);
                if( eCalEt >= 20. && eCalEt < 40. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0024);
                if( eCalEt >= 40. && eCalEt < 60. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0022);
                if( eCalEt >= 60. && eCalEt < 80. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0019);
                if( eCalEt >= 80. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0017);
            }
            if( fabs(eCalEta) > 1.4 && fabs(eCalEta) < 1.7 ){
                if( eCalEt > 15. && eCalEt < 20. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0031);
                if( eCalEt >= 20. && eCalEt < 40. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0023);
                if( eCalEt >= 40. && eCalEt < 60. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0022);
                if( eCalEt >= 60. && eCalEt < 80. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0020);
                if( eCalEt >= 80. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0019);
            }
            if( fabs(eCalEta) > 1.7 && fabs(eCalEta) < 2.1 ){
                if( eCalEt > 15. && eCalEt < 20. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0037);
                if( eCalEt >= 20. && eCalEt < 40. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0035);
                if( eCalEt >= 40. && eCalEt < 60. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0031);
                if( eCalEt >= 60. && eCalEt < 80.) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0027);
                if( eCalEt >= 80. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0024);
            }
            if( fabs(eCalEta) > 2.1 && fabs(eCalEta) < 2.7 ){
                if( eCalEt > 15. && eCalEt < 20. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0033);
                if( eCalEt >= 20. && eCalEt < 40. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0022);
                if( eCalEt >= 40. && eCalEt < 60. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0018);
                if( eCalEt >= 60. && eCalEt < 80. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0017);
                if( eCalEt >= 80. && eCalEt < 95. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0018);
                if( eCalEt >= 95. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.004000-0.0026);
            }
            if( fabs(eCalEta) > 2.7 && fabs(eCalEta) < 3.0){
                eCalPhi = gRandom->Gaus(eCalPhi, 0.00445);
            }

            //ECal Eta smearing
            if( fabs(eCalEta) < 0.8 ) eCalEta = gRandom->Gaus(eCalEta, 0.0034);
            if( fabs(eCalEta) > 0.8 && fabs(eCalEta) < 1.4) eCalEta = gRandom->Gaus(eCalEta, 0.0037);
            if( fabs(eCalEta) > 1.4 && fabs(eCalEta) < 1.7) eCalEta = gRandom->Gaus(eCalEta, 0.0033);
            if( fabs(eCalEta) > 1.7 && fabs(eCalEta) < 2.1) eCalEta = gRandom->Gaus(eCalEta, 0.0059);
            if( fabs(eCalEta) > 2.1 && fabs(eCalEta) < 2.7) eCalEta = gRandom->Gaus(eCalEta, 0.0059);
            if( fabs(eCalEta) > 2.7 && fabs(eCalEta) < 3.0) eCalEta = gRandom->Gaus(eCalEta, 0.0059);

            float distZ = 330.;
            if( fabs(eCalEta) < 1.47 )
            {
                float theta = 2.*atan(exp(-eCalEta));
                eCalX = 137*cos(eCalPhi);   // x = r*cos(phi)
                eCalY = 137*sin(eCalPhi);   // y = r*sin(phi)
                eCalZ = 137*tan(TMath::Pi()/2.-theta);     // z = r*tan(theta)
                if( fabs(eCalZ) > 278. ) continue;
                //cout << eCalEta << ", " << TMath::RadToDeg*theta << ", " << eCalZ << endl; 
            }   
            if( eCalEta > 1.47 )
            {   
                float theta = 2.*atan(exp(-eCalEta));
                eCalZ = 330.;
                float r = eCalZ*tan(theta);
                if( r < 35. || r > 155 ) continue;
                eCalX = r*cos(eCalPhi);
                eCalY = r*sin(eCalPhi);
            }   
            if( eCalEta < -1.47 )
            {   
                float theta = 2.*atan(exp(-eCalEta));
                eCalZ = -330.;
                float r = eCalZ*tan(theta);
                if( r < 35. || r > 155 ) continue;
                eCalX = r*cos(eCalPhi);
                eCalY = r*sin(eCalPhi);
            }

            //1.29m = 129cm                
            //ECal position (unit: cm)      
            //eCalX = 129*cos(eCalPhi);   // x = r*cos(phi)
            //eCalY = 129*sin(eCalPhi);   // y = r*sin(phi)
            //float theta = 2.*atan(exp(-eCalEta));
            //eCalZ = 129./tan(theta);     // z = r*tan(theta)

            egCrysClusterEta.push_back(eCalEta);
            egCrysClusterPhi.push_back(eCalPhi);
            egCrysClusterEt.push_back(eCalEt);
            egCrysClusterGx.push_back(eCalX);
            egCrysClusterGy.push_back(eCalY);
            egCrysClusterGz.push_back(eCalZ);
        }
/*
        for(unsigned int k = 0; k < TrackN; k++)
        {
            eCalTower = (Tower*) branchTower->At(k);
            eCalTrack = (Track*) branchTrack->At(k);
            double eCalX, eCalY, eCalZ, eCalEt, eCalEta, eCalPhi;
            eCalX = eCalTrack->XOuter/10.;
            //eCalX -= 0.001923;
            //eCalX = gRandom->Gaus(eCalX, 0.3729);
            
            eCalY = eCalTrack->YOuter/10.;
            //eCalY -= 0.00132;
            //eCalY = gRandom->Gaus(eCalY, 0.3711);
            
            eCalZ = eCalTrack->ZOuter/10.;
            //eCalZ += 0.004241;
            //eCalZ = gRandom->Gaus(eCalZ, 0.6597-0.15);

            TVector3 v1;
            v1.SetXYZ(eCalX, eCalY, eCalZ);
            eCalEta = v1.Eta();
            eCalPhi = v1.Phi();
            
            //ECal track Et
            //eCalEt = eCalTrack->PT;
            eCalEt = eCalTower->ET;
            cout << "eCalEt :" << eCalEt << endl;
            //float tempEt = eCalEt*0.1012;
            //eCalEt -= tempEt;
            //eCalEt = gRandom->Gaus(eCalEt, eCalEt*0.02769);

            //float R = sqrt(eCalX*eCalX+eCalY*eCalY);
            //float R3 = sqrt(eCalX*eCalX+eCalY*eCalY+eCalZ*eCalZ);
            
            //ECal track eta
            //eCalEta = eCalTrack->EtaOuter;
            //eCalEta += -8.084e-5;
            //if( eCalEt < 30. ) eCalEta = gRandom->Gaus(eCalEta, 0.004276-0.002);
            //if( eCalEt >= 30. && eCalEt < 60. ) eCalEta = gRandom->Gaus(eCalEta, 0.004276-0.001);
            //if( eCalEt >= 60. ) eCalEta = gRandom->Gaus(eCalEta, 0.004276);
            
            //if( fabs(eCalEta) < 0.8 )
            //{
            //    //ECal track Et
            //    eCalEt = eCalTrack->PT;
            //    float tempEt = eCalEt*0.1012;
            //    eCalEt -= tempEt;
            //    eCalEt = gRandom->Gaus(eCalEt, eCalEt*0.02769);
            //    
            //    //ECal track phi
            //    eCalPhi = eCalTrack->PhiOuter;
            //    eCalPhi += -5.092e-5;

            //    if( eCalEt < 30. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.003668-0.0010);
            //    if( eCalEt >= 30. && eCalEt < 60. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.003668);
            //    if( eCalEt >= 60. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.003668+0.0010);
            //    
            //    eCalX = R*cos(eCalPhi);
            //    eCalY = R*sin(eCalPhi);
            //    eCalZ = R3*cos(2.*atan(-eCalEta));
            //}
            //else
            //{
            //    eCalEt = eCalTrack->PT;
            //    eCalPhi = eCalTrack->PhiOuter;
            //    eCalX = R*cos(eCalPhi);
            //    eCalY = R*sin(eCalPhi);
            //    eCalZ = R3*cos(2.*atan(-eCalEta));
            //}
            
            //ECal track Et
            //eCalEt = eCalTrack->PT;
            //float tempEt = eCalEt*0.100361;
            //eCalEt -= tempEt;

            
            //if( eCalEt < 30. ) eCalEt = gRandom->Gaus(eCalEt, 0.0237839);
            //if( eCalEt >= 30. && eCalEt < 60. ) eCalEt = gRandom->Gaus(eCalEt, 0.0337839);
            //if( eCalEt >= 60. ) eCalEt = gRandom->Gaus(eCalEt, 0.0437839);
            
            //ECal track eta
            eCalEta = eCalTrack->EtaOuter;
            eCalEta += -3.9404e-05;
            //eCalEta = gRandom->Gaus(eCalEta, 0.00424596);

            if( eCalEt < 30. ) eCalEta = gRandom->Gaus(eCalEta, 0.0001);
            if( eCalEt >= 30. && eCalEt < 60. ) eCalEta = gRandom->Gaus(eCalEta, 0.0002);
            if( eCalEt >= 60. ) eCalEta = gRandom->Gaus(eCalEta, 0.0003);
            
            //ECal track phi
            eCalPhi = eCalTrack->PhiOuter;
            eCalPhi += -2.95267e-05;
            eCalPhi = gRandom->Gaus(eCalPhi, 0.00347567);

            if( eCalEt < 30. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.00347567-0.0010);
            if( eCalEt >= 30. && eCalEt < 60. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.00347567);
            if( eCalEt >= 60. ) eCalPhi = gRandom->Gaus(eCalPhi, 0.00347567+0.0010);
           

            
            //ECal track X, Y, Z 
            eCalX = eCalTrack->XOuter/10.;
            eCalX += -0.000918254;
            eCalX = gRandom->Gaus(eCalX, 0.347522);
            
            eCalY = eCalTrack->YOuter/10.;
            eCalY += -0.00223196;
            eCalY = gRandom->Gaus(eCalY, 0.346173);
            
            eCalZ = eCalTrack->ZOuter/10.;
            eCalZ += 0.00612539;
            eCalZ = gRandom->Gaus(eCalZ, 0.8371);
            

            egCrysClusterEta.push_back(eCalEta);
            egCrysClusterPhi.push_back(eCalPhi);
            egCrysClusterEt.push_back(eCalEt);
            egCrysClusterGx.push_back(eCalX);
            egCrysClusterGy.push_back(eCalY);
            egCrysClusterGz.push_back(eCalZ);
        }*/

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

void Example20(const char *inputFile)
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
