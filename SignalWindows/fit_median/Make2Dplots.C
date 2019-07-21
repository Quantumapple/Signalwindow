#define Make2Dplots_cxx
#include "Make2Dplots.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGaxis.h>

#include <iostream>
#include <fstream>

#include "./Style/CMS_lumi.C"

class L1eGamma
{
    public:
        float egEt;
        float EgdPhi;

        L1eGamma() { 
            egEt = 0; EgdPhi = 0;
        }

        L1eGamma(float a, float b) {
            egEt = a; EgdPhi = b;
        }
};

double getMedian(const std::vector<float> &vec)
{
    double median=0.;

    //odd number, definate median
    if(vec.size() % 2 !=0) {
        int middleNr = (vec.size()+1)/2;
        median = vec[middleNr];
    }else{ //even number, take median as halfway between the two middle values
        int middleNr = (vec.size()+1)/2;
        median= vec[middleNr];
        if(middleNr+1 <(int) vec.size()) median+= vec[middleNr+1];
        median/=2.;
    }
    return median;
}

int getMedianIndex(const std::vector<float> &vec)
{

    //odd number, definate median
    if(vec.size() % 2 !=0) {
        int middleNr = (vec.size()+1)/2;
        return middleNr;
    }else{ //even number, take median as halfway between the two middle values
        int middleNr = (vec.size()+1)/2;
        return middleNr;
    }
}

void test::Loop(int eta_ = 1)
{
    bool debug = false;

    if (fChain == 0) return;

    const TString eta_region[6] = {"EtaRegion1", "EtaRegion2", "EtaRegion3", "EtaRegion4", "EtaRegion5", "EtaRegion6"};

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    TH2F* pix1egDphi_dist = new TH2F("pix1egDphi_dist","pix1egDphi_dist", 90,10,100,2500,-0.2,0.2);
    TH2F* pix2egDphi_dist = new TH2F("pix2egDphi_dist","pix2egDphi_dist", 90,10,100,2500,-0.2,0.2);
    TH2F* pix3egDphi_dist = new TH2F("pix3egDphi_dist","pix3egDphi_dist", 90,10,100,2500,-0.2,0.2);
    TH2F* pix4egDphi_dist = new TH2F("pix4egDphi_dist","pix4egDphi_dist", 90,10,100,2500,-0.2,0.2);

    TH2F* pix12egDphi_dist = new TH2F("pix12egDphi_dist","pix12egDphi_dist", 90,10,100,2500,-0.2,0.2);
    TH2F* pix13egDphi_dist = new TH2F("pix13egDphi_dist","pix13egDphi_dist", 90,10,100,2500,-0.2,0.2);
    TH2F* pix14egDphi_dist = new TH2F("pix14egDphi_dist","pix14egDphi_dist", 90,10,100,2500,-0.2,0.2);
    TH2F* pix23egDphi_dist = new TH2F("pix23egDphi_dist","pix23egDphi_dist", 90,10,100,2500,-0.2,0.2);
    TH2F* pix24egDphi_dist = new TH2F("pix24egDphi_dist","pix24egDphi_dist", 90,10,100,2500,-0.2,0.2);
    TH2F* pix34egDphi_dist = new TH2F("pix34egDphi_dist","pix34egDphi_dist", 90,10,100,2500,-0.2,0.2);

    TH2F* pix012Dphi_dist = new TH2F("pix012Dphi_dist","pix012Dphi_dist", 90,10,100,5000,-0.02,0.02);
    TH2F* pix013Dphi_dist = new TH2F("pix013Dphi_dist","pix013Dphi_dist", 90,10,100,5000,-0.02,0.02);
    TH2F* pix014Dphi_dist = new TH2F("pix014Dphi_dist","pix014Dphi_dist", 90,10,100,5000,-0.02,0.02);
    TH2F* pix023Dphi_dist = new TH2F("pix023Dphi_dist","pix023Dphi_dist", 90,10,100,5000,-0.02,0.02);
    TH2F* pix024Dphi_dist = new TH2F("pix024Dphi_dist","pix024Dphi_dist", 90,10,100,5000,-0.02,0.02);
    TH2F* pix034Dphi_dist = new TH2F("pix034Dphi_dist","pix034Dphi_dist", 90,10,100,5000,-0.02,0.02);
    TH2F* pix123Dphi_dist = new TH2F("pix123Dphi_dist","pix123Dphi_dist", 90,10,100,5000,-0.02,0.02);
    TH2F* pix124Dphi_dist = new TH2F("pix124Dphi_dist","pix124Dphi_dist", 90,10,100,5000,-0.02,0.02);
    TH2F* pix134Dphi_dist = new TH2F("pix134Dphi_dist","pix134Dphi_dist", 90,10,100,5000,-0.02,0.02);
    TH2F* pix234Dphi_dist = new TH2F("pix234Dphi_dist","pix234Dphi_dist", 90,10,100,5000,-0.02,0.02);

    const int binSize = 90;

    vector<L1eGamma> pix1egDphi_;
    vector<L1eGamma> pix2egDphi_;
    vector<L1eGamma> pix3egDphi_;
    vector<L1eGamma> pix4egDphi_;

    vector<L1eGamma> pix12egDphi_;
    vector<L1eGamma> pix13egDphi_;
    vector<L1eGamma> pix14egDphi_;
    vector<L1eGamma> pix23egDphi_;
    vector<L1eGamma> pix24egDphi_;
    vector<L1eGamma> pix34egDphi_;

    vector<L1eGamma> pix012Dphi_;
    vector<L1eGamma> pix013Dphi_;
    vector<L1eGamma> pix014Dphi_;
    vector<L1eGamma> pix023Dphi_;
    vector<L1eGamma> pix024Dphi_;
    vector<L1eGamma> pix034Dphi_;
    vector<L1eGamma> pix123Dphi_;
    vector<L1eGamma> pix124Dphi_;
    vector<L1eGamma> pix134Dphi_;
    vector<L1eGamma> pix234Dphi_;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if( eta_ == 1 && fabs(ntEgEta->at(0)) > 0.8 ) continue; 
        if( eta_ == 2 && (fabs(ntEgEta->at(0)) < 0.8 || fabs(ntEgEta->at(0)) > 1.4)) continue; 
        if( eta_ == 3 && (fabs(ntEgEta->at(0)) < 1.4 || fabs(ntEgEta->at(0)) > 1.7)) continue; 
        if( eta_ == 4 && (fabs(ntEgEta->at(0)) < 1.7 || fabs(ntEgEta->at(0)) > 2.1)) continue; 
        if( eta_ == 5 && (fabs(ntEgEta->at(0)) < 2.1 || fabs(ntEgEta->at(0)) > 2.7)) continue; 
        if( eta_ == 6 && (fabs(ntEgEta->at(0)) < 2.7 || fabs(ntEgEta->at(0)) > 3.0)) continue; 

        for(unsigned long i = 0; i < ntPix1EGdphi->size(); i++) {
            pix1egDphi_dist->Fill(ntEgEt->at(0), ntPix1EGdphi->at(i));
            pix1egDphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix1EGdphi->at(i) ) ); 
        }          
        for(unsigned long i = 0; i < ntPix2EGdphi->size(); i++) {
            pix2egDphi_dist->Fill(ntEgEt->at(0), ntPix2EGdphi->at(i));
            pix2egDphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix2EGdphi->at(i) ) ); 
        }
        for(unsigned long i = 0; i < ntPix3EGdphi->size(); i++) { 
            pix3egDphi_dist->Fill(ntEgEt->at(0), ntPix3EGdphi->at(i));
            pix3egDphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix3EGdphi->at(i) ) ); 
        }
        for(unsigned long i = 0; i < ntPix4EGdphi->size(); i++) { 
            pix4egDphi_dist->Fill(ntEgEt->at(0), ntPix4EGdphi->at(i));
            pix4egDphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix4EGdphi->at(i) ) ); 
        }

        for(unsigned long i = 0; i < ntPix12EGdphi->size(); i++) { 
            if(fabs(ntPix12EGdphi->at(i)) > 0.2 ) continue;
            pix12egDphi_dist->Fill(ntEgEt->at(0), ntPix12EGdphi->at(i));
            pix12egDphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix12EGdphi->at(i) ) );
        }
        for(unsigned long i = 0; i < ntPix13EGdphi->size(); i++) {
            if(fabs(ntPix13EGdphi->at(i)) > 0.2 ) continue;
            pix13egDphi_dist->Fill(ntEgEt->at(0), ntPix13EGdphi->at(i));
            pix13egDphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix13EGdphi->at(i) ) );
        }
        for(unsigned long i = 0; i < ntPix14EGdphi->size(); i++) {
            if(fabs(ntPix14EGdphi->at(i)) > 0.2 ) continue;
            pix14egDphi_dist->Fill(ntEgEt->at(0), ntPix14EGdphi->at(i));
            pix14egDphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix14EGdphi->at(i) ) );
        }
        for(unsigned long i = 0; i < ntPix23EGdphi->size(); i++) { 
            if(fabs(ntPix23EGdphi->at(i)) > 0.2 ) continue;
            pix23egDphi_dist->Fill(ntEgEt->at(0), ntPix23EGdphi->at(i));
            pix23egDphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix23EGdphi->at(i) ) );
        }
        for(unsigned long i = 0; i < ntPix24EGdphi->size(); i++) { 
            if(fabs(ntPix24EGdphi->at(i)) > 0.2 ) continue;
            pix24egDphi_dist->Fill(ntEgEt->at(0), ntPix24EGdphi->at(i));
            pix24egDphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix24EGdphi->at(i) ) );
        }
        for(unsigned long i = 0; i < ntPix34EGdphi->size(); i++) {
            if(fabs(ntPix34EGdphi->at(i)) > 0.2 ) continue;
            pix34egDphi_dist->Fill(ntEgEt->at(0), ntPix34EGdphi->at(i));
            pix34egDphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix34EGdphi->at(i) ) );
        }

        for(unsigned long i = 0; i < ntPix012dphi->size(); i++) {
            if(fabs(ntPix012dphi->at(i)) > 0.015 ) continue;
            pix012Dphi_dist->Fill(ntEgEt->at(0), ntPix012dphi->at(i));
            pix012Dphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix012dphi->at(i) ) );  
        }
        for(unsigned long i = 0; i < ntPix013dphi->size(); i++) {
            if(fabs(ntPix013dphi->at(i)) > 0.015 ) continue;
            pix013Dphi_dist->Fill(ntEgEt->at(0), ntPix013dphi->at(i));
            pix013Dphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix013dphi->at(i) ) );  
        }
        for(unsigned long i = 0; i < ntPix014dphi->size(); i++) { 
            if(fabs(ntPix014dphi->at(i)) > 0.015 ) continue;
            pix014Dphi_dist->Fill(ntEgEt->at(0), ntPix014dphi->at(i));
            pix014Dphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix014dphi->at(i) ) );  
        }
        for(unsigned long i = 0; i < ntPix023dphi->size(); i++) {
            if(fabs(ntPix023dphi->at(i)) > 0.015 ) continue;
            pix023Dphi_dist->Fill(ntEgEt->at(0), ntPix023dphi->at(i));
            pix023Dphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix023dphi->at(i) ) );  
        }
        for(unsigned long i = 0; i < ntPix024dphi->size(); i++) {
            if(fabs(ntPix024dphi->at(i)) > 0.015 ) continue;
            pix024Dphi_dist->Fill(ntEgEt->at(0), ntPix024dphi->at(i));
            pix024Dphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix024dphi->at(i) ) );  
        }
        for(unsigned long i = 0; i < ntPix034dphi->size(); i++) {
            if(fabs(ntPix034dphi->at(i)) > 0.015 ) continue;
            pix034Dphi_dist->Fill(ntEgEt->at(0), ntPix034dphi->at(i));
            pix034Dphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix034dphi->at(i) ) );  
        }
        for(unsigned long i = 0; i < ntPix123dphi->size(); i++) {
            if(fabs(ntPix123dphi->at(i)) > 0.015 ) continue;
            pix123Dphi_dist->Fill(ntEgEt->at(0), ntPix123dphi->at(i));
            pix123Dphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix123dphi->at(i) ) );  
        }
        for(unsigned long i = 0; i < ntPix124dphi->size(); i++) {
            if(fabs(ntPix124dphi->at(i)) > 0.015 ) continue;
            pix124Dphi_dist->Fill(ntEgEt->at(0), ntPix124dphi->at(i));
            pix124Dphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix124dphi->at(i) ) );  
        }
        for(unsigned long i = 0; i < ntPix134dphi->size(); i++) {
            if(fabs(ntPix134dphi->at(i)) > 0.015 ) continue;
            pix134Dphi_dist->Fill(ntEgEt->at(0), ntPix134dphi->at(i));
            pix134Dphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix134dphi->at(i) ) );  
        }
        for(unsigned long i = 0; i < ntPix234dphi->size(); i++) {
            if(fabs(ntPix234dphi->at(i)) > 0.015 ) continue;
            pix234Dphi_dist->Fill(ntEgEt->at(0), ntPix234dphi->at(i));
            pix234Dphi_.push_back( L1eGamma( ntEgEt->at(0), ntPix234dphi->at(i) ) );  
        }

    } // event loop

    if( debug ) cout << "Pass event loop" << endl;

    float x[binSize];

    // eg-pixel dphi
    float pix1egDphi[binSize], pix1egDphiErr[binSize];
    float pix2egDphi[binSize], pix2egDphiErr[binSize];
    float pix3egDphi[binSize], pix3egDphiErr[binSize];
    float pix4egDphi[binSize], pix4egDphiErr[binSize];

    // eg-pixelpixel dphi
    float pix12egDphi[binSize], pix12egDphiErr[binSize];
    float pix13egDphi[binSize], pix13egDphiErr[binSize];
    float pix14egDphi[binSize], pix14egDphiErr[binSize];
    float pix23egDphi[binSize], pix23egDphiErr[binSize];
    float pix24egDphi[binSize], pix24egDphiErr[binSize];
    float pix34egDphi[binSize], pix34egDphiErr[binSize];

    // pixel-pixel dphi
    float pix012Dphi[binSize], pix012DphiErr[binSize];
    float pix013Dphi[binSize], pix013DphiErr[binSize];
    float pix014Dphi[binSize], pix014DphiErr[binSize];
    float pix023Dphi[binSize], pix023DphiErr[binSize];
    float pix024Dphi[binSize], pix024DphiErr[binSize];
    float pix034Dphi[binSize], pix034DphiErr[binSize];
    float pix123Dphi[binSize], pix123DphiErr[binSize];
    float pix124Dphi[binSize], pix124DphiErr[binSize];
    float pix134Dphi[binSize], pix134DphiErr[binSize];
    float pix234Dphi[binSize], pix234DphiErr[binSize];

    // vectors for TGraphErrors
    vector<float*> pixegDphi, pixegDphiErr;
    vector<float*> pixpixegDphi, pixpixegDphiErr;
    vector<float*> pixpixDphi, pixpixDphiErr;

    pixegDphi.clear();
    pixegDphiErr.clear();
    pixpixegDphi.clear();
    pixpixegDphiErr.clear();
    pixpixDphi.clear(); 
    pixpixDphiErr.clear();

    pixegDphi.push_back(pix1egDphi);
    pixegDphi.push_back(pix2egDphi);
    pixegDphi.push_back(pix3egDphi);
    pixegDphi.push_back(pix4egDphi);

    pixegDphiErr.push_back(pix1egDphiErr);
    pixegDphiErr.push_back(pix2egDphiErr);
    pixegDphiErr.push_back(pix3egDphiErr);
    pixegDphiErr.push_back(pix4egDphiErr);

    pixpixegDphi.push_back(pix12egDphi);
    pixpixegDphi.push_back(pix13egDphi);
    pixpixegDphi.push_back(pix14egDphi);
    pixpixegDphi.push_back(pix23egDphi);
    pixpixegDphi.push_back(pix24egDphi);
    pixpixegDphi.push_back(pix34egDphi);

    pixpixegDphiErr.push_back(pix12egDphiErr);
    pixpixegDphiErr.push_back(pix13egDphiErr);
    pixpixegDphiErr.push_back(pix14egDphiErr);
    pixpixegDphiErr.push_back(pix23egDphiErr);
    pixpixegDphiErr.push_back(pix24egDphiErr);
    pixpixegDphiErr.push_back(pix34egDphiErr);

    pixpixDphi.push_back(pix012Dphi);
    pixpixDphi.push_back(pix013Dphi);
    pixpixDphi.push_back(pix014Dphi);
    pixpixDphi.push_back(pix023Dphi);
    pixpixDphi.push_back(pix024Dphi);
    pixpixDphi.push_back(pix034Dphi);
    pixpixDphi.push_back(pix123Dphi);
    pixpixDphi.push_back(pix124Dphi);
    pixpixDphi.push_back(pix134Dphi);
    pixpixDphi.push_back(pix234Dphi);

    pixpixDphiErr.push_back(pix012DphiErr);
    pixpixDphiErr.push_back(pix013DphiErr);
    pixpixDphiErr.push_back(pix014DphiErr);
    pixpixDphiErr.push_back(pix023DphiErr);
    pixpixDphiErr.push_back(pix024DphiErr);
    pixpixDphiErr.push_back(pix034DphiErr);
    pixpixDphiErr.push_back(pix123DphiErr);
    pixpixDphiErr.push_back(pix124DphiErr);
    pixpixDphiErr.push_back(pix134DphiErr);
    pixpixDphiErr.push_back(pix234DphiErr);

    //vectors for median 
    vector<float> pix1egDphiMed;
    vector<float> pix2egDphiMed;
    vector<float> pix3egDphiMed;
    vector<float> pix4egDphiMed;

    vector<float> pix12egDphiMed;
    vector<float> pix13egDphiMed;
    vector<float> pix14egDphiMed;
    vector<float> pix23egDphiMed;
    vector<float> pix24egDphiMed;
    vector<float> pix34egDphiMed;

    vector<float> pix012DphiMed;
    vector<float> pix013DphiMed;
    vector<float> pix014DphiMed;
    vector<float> pix023DphiMed;
    vector<float> pix024DphiMed;
    vector<float> pix034DphiMed;
    vector<float> pix123DphiMed;
    vector<float> pix124DphiMed;
    vector<float> pix134DphiMed;
    vector<float> pix234DphiMed;

    if( debug ) cout << pix1egDphi_.size() << endl; 
    if( debug ) cout << pix2egDphi_.size() << endl; 
    if( debug ) cout << pix3egDphi_.size() << endl; 
    if( debug ) cout << pix4egDphi_.size() << endl; 
    cout << endl;
    if( debug ) cout << pix12egDphi_.size() << endl; 
    if( debug ) cout << pix13egDphi_.size() << endl; 
    if( debug ) cout << pix14egDphi_.size() << endl; 
    if( debug ) cout << pix23egDphi_.size() << endl; 
    if( debug ) cout << pix24egDphi_.size() << endl; 
    if( debug ) cout << pix34egDphi_.size() << endl; 
    cout << endl;
    if( debug ) cout << pix012Dphi_.size() << endl; 
    if( debug ) cout << pix013Dphi_.size() << endl; 
    if( debug ) cout << pix014Dphi_.size() << endl; 
    if( debug ) cout << pix023Dphi_.size() << endl; 
    if( debug ) cout << pix024Dphi_.size() << endl; 
    if( debug ) cout << pix034Dphi_.size() << endl; 
    if( debug ) cout << pix123Dphi_.size() << endl; 
    if( debug ) cout << pix124Dphi_.size() << endl; 
    if( debug ) cout << pix134Dphi_.size() << endl; 
    if( debug ) cout << pix234Dphi_.size() << endl; 


    //Fill
    for(Int_t nth = 0; nth < binSize; nth++)
    {
        if( debug ) cout << "Which bin: " << nth << " between " << 10.+nth << " and " << 11.+nth << endl;
        cout << "Which bin: " << nth << " between " << 10.+nth << " and " << 11.+nth << endl;
        x[nth] = 10.5 + float(nth);
        pix1egDphiMed.clear();
        pix2egDphiMed.clear();
        pix3egDphiMed.clear();
        pix4egDphiMed.clear();

        pix12egDphiMed.clear();
        pix13egDphiMed.clear();
        pix14egDphiMed.clear();
        pix23egDphiMed.clear();
        pix24egDphiMed.clear();
        pix34egDphiMed.clear();

        pix012DphiMed.clear();
        pix013DphiMed.clear();
        pix014DphiMed.clear();
        pix023DphiMed.clear();
        pix024DphiMed.clear();
        pix034DphiMed.clear();
        pix123DphiMed.clear();
        pix124DphiMed.clear();
        pix134DphiMed.clear();
        pix234DphiMed.clear();

        // Region of Interest
        for(Int_t j = 0; j < pix1egDphi_.size(); j++)
        {
            float tempEt = pix1egDphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix1egDphiMed.push_back(pix1egDphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix2egDphi_.size(); j++)
        {
            float tempEt = pix2egDphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix2egDphiMed.push_back(pix2egDphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix3egDphi_.size(); j++)
        {
            float tempEt = pix3egDphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix3egDphiMed.push_back(pix3egDphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix4egDphi_.size(); j++)
        {
            float tempEt = pix4egDphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix4egDphiMed.push_back(pix4egDphi_[j].EgdPhi);
        }

        // Pixel + egamma matching
        for(Int_t j = 0; j < pix12egDphi_.size(); j++)
        {
            float tempEt = pix12egDphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix12egDphiMed.push_back(pix12egDphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix13egDphi_.size(); j++)
        {
            float tempEt = pix13egDphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix13egDphiMed.push_back(pix13egDphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix14egDphi_.size(); j++)
        {
            float tempEt = pix14egDphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix14egDphiMed.push_back(pix14egDphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix23egDphi_.size(); j++)
        {
            float tempEt = pix23egDphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix23egDphiMed.push_back(pix23egDphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix24egDphi_.size(); j++)
        {
            float tempEt = pix24egDphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix24egDphiMed.push_back(pix24egDphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix34egDphi_.size(); j++)
        {
            float tempEt = pix34egDphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix34egDphiMed.push_back(pix34egDphi_[j].EgdPhi);
        }

        // Pixel-pixel matching
        for(Int_t j = 0; j < pix012Dphi_.size(); j++)
        {
            float tempEt = pix012Dphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix012DphiMed.push_back(pix012Dphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix013Dphi_.size(); j++)
        {
            float tempEt = pix013Dphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix013DphiMed.push_back(pix013Dphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix014Dphi_.size(); j++)
        {
            float tempEt = pix014Dphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix014DphiMed.push_back(pix014Dphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix023Dphi_.size(); j++)
        {
            float tempEt = pix023Dphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix023DphiMed.push_back(pix023Dphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix024Dphi_.size(); j++)
        {
            float tempEt = pix024Dphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix024DphiMed.push_back(pix024Dphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix034Dphi_.size(); j++)
        {
            float tempEt = pix034Dphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix034DphiMed.push_back(pix034Dphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix123Dphi_.size(); j++)
        {
            float tempEt = pix123Dphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix123DphiMed.push_back(pix123Dphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix124Dphi_.size(); j++)
        {
            float tempEt = pix124Dphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix124DphiMed.push_back(pix124Dphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix134Dphi_.size(); j++)
        {
            float tempEt = pix134Dphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix134DphiMed.push_back(pix134Dphi_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix234Dphi_.size(); j++)
        {
            float tempEt = pix234Dphi_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix234DphiMed.push_back(pix234Dphi_[j].EgdPhi);
        }

        if( debug ) cout << "Fill new vecotrs!" << endl;

        std::sort(pix1egDphiMed.begin(), pix1egDphiMed.end()); 
        std::sort(pix2egDphiMed.begin(), pix2egDphiMed.end()); 
        std::sort(pix3egDphiMed.begin(), pix3egDphiMed.end()); 
        std::sort(pix4egDphiMed.begin(), pix4egDphiMed.end()); 

        std::sort(pix12egDphiMed.begin(), pix12egDphiMed.end()); 
        std::sort(pix13egDphiMed.begin(), pix13egDphiMed.end()); 
        std::sort(pix14egDphiMed.begin(), pix14egDphiMed.end()); 
        std::sort(pix23egDphiMed.begin(), pix23egDphiMed.end()); 
        std::sort(pix24egDphiMed.begin(), pix24egDphiMed.end()); 
        std::sort(pix34egDphiMed.begin(), pix34egDphiMed.end()); 

        std::sort(pix012DphiMed.begin(), pix012DphiMed.end()); 
        std::sort(pix013DphiMed.begin(), pix013DphiMed.end()); 
        std::sort(pix014DphiMed.begin(), pix014DphiMed.end()); 
        std::sort(pix023DphiMed.begin(), pix023DphiMed.end()); 
        std::sort(pix024DphiMed.begin(), pix024DphiMed.end()); 
        std::sort(pix034DphiMed.begin(), pix034DphiMed.end()); 
        std::sort(pix123DphiMed.begin(), pix123DphiMed.end()); 
        std::sort(pix124DphiMed.begin(), pix124DphiMed.end()); 
        std::sort(pix134DphiMed.begin(), pix134DphiMed.end()); 
        std::sort(pix234DphiMed.begin(), pix234DphiMed.end()); 

        if( debug ) cout << "vector sorting is done!" << endl;

        int pix1egDphi_size = pix1egDphiMed.size();
        int pix2egDphi_size = pix2egDphiMed.size();
        int pix3egDphi_size = pix3egDphiMed.size();
        int pix4egDphi_size = pix4egDphiMed.size();

        int pix12egDphi_size = pix12egDphiMed.size();
        int pix13egDphi_size = pix13egDphiMed.size();
        int pix14egDphi_size = pix14egDphiMed.size();
        int pix23egDphi_size = pix23egDphiMed.size();
        int pix24egDphi_size = pix24egDphiMed.size();
        int pix34egDphi_size = pix34egDphiMed.size();

        int pix012Dphi_size = pix012DphiMed.size();
        int pix013Dphi_size = pix013DphiMed.size();
        int pix014Dphi_size = pix014DphiMed.size();
        int pix023Dphi_size = pix023DphiMed.size();
        int pix024Dphi_size = pix024DphiMed.size();
        int pix034Dphi_size = pix034DphiMed.size();
        int pix123Dphi_size = pix123DphiMed.size();
        int pix124Dphi_size = pix124DphiMed.size();
        int pix134Dphi_size = pix134DphiMed.size();
        int pix234Dphi_size = pix234DphiMed.size();

        if( debug ) cout << pix1egDphi_size << endl;
        if( debug ) cout << pix2egDphi_size << endl;
        if( debug ) cout << pix3egDphi_size << endl;
        if( debug ) cout << pix4egDphi_size << endl;

        if( debug ) cout << pix12egDphi_size << endl;
        if( debug ) cout << pix13egDphi_size << endl;
        if( debug ) cout << pix14egDphi_size << endl;
        if( debug ) cout << pix23egDphi_size << endl;
        if( debug ) cout << pix24egDphi_size << endl;
        if( debug ) cout << pix34egDphi_size << endl;

        if( debug ) cout << pix012Dphi_size << endl;
        if( debug ) cout << pix013Dphi_size << endl;
        if( debug ) cout << pix014Dphi_size << endl;
        if( debug ) cout << pix023Dphi_size << endl;
        if( debug ) cout << pix024Dphi_size << endl;
        if( debug ) cout << pix034Dphi_size << endl;
        if( debug ) cout << pix123Dphi_size << endl;
        if( debug ) cout << pix124Dphi_size << endl;
        if( debug ) cout << pix134Dphi_size << endl;
        if( debug ) cout << pix234Dphi_size << endl;

        // get index for median
        int pix1egDphi_medianIndex = getMedianIndex(pix1egDphiMed);
        int pix2egDphi_medianIndex = getMedianIndex(pix2egDphiMed);
        int pix3egDphi_medianIndex = getMedianIndex(pix3egDphiMed);
        int pix4egDphi_medianIndex = getMedianIndex(pix4egDphiMed);

        int pix12egDphi_medianIndex = getMedianIndex(pix12egDphiMed);
        int pix13egDphi_medianIndex = getMedianIndex(pix13egDphiMed);
        int pix14egDphi_medianIndex = getMedianIndex(pix14egDphiMed);
        int pix23egDphi_medianIndex = getMedianIndex(pix23egDphiMed);
        int pix24egDphi_medianIndex = getMedianIndex(pix24egDphiMed);
        int pix34egDphi_medianIndex = getMedianIndex(pix34egDphiMed);

        int pix012Dphi_medianIndex = getMedianIndex(pix012DphiMed);
        int pix013Dphi_medianIndex = getMedianIndex(pix013DphiMed);
        int pix014Dphi_medianIndex = getMedianIndex(pix014DphiMed);
        int pix023Dphi_medianIndex = getMedianIndex(pix023DphiMed);
        int pix024Dphi_medianIndex = getMedianIndex(pix024DphiMed);
        int pix034Dphi_medianIndex = getMedianIndex(pix034DphiMed);
        int pix123Dphi_medianIndex = getMedianIndex(pix123DphiMed);
        int pix124Dphi_medianIndex = getMedianIndex(pix124DphiMed);
        int pix134Dphi_medianIndex = getMedianIndex(pix134DphiMed);
        int pix234Dphi_medianIndex = getMedianIndex(pix234DphiMed);

        if( debug ) cout << "Successfully get median index" << endl;

        // get median
        float pix1egDphi_median = getMedian(pix1egDphiMed);
        float pix2egDphi_median = getMedian(pix2egDphiMed);
        float pix3egDphi_median = getMedian(pix3egDphiMed);
        float pix4egDphi_median = getMedian(pix4egDphiMed);

        float pix12egDphi_median = getMedian(pix12egDphiMed);
        float pix13egDphi_median = getMedian(pix13egDphiMed);
        float pix14egDphi_median = getMedian(pix14egDphiMed);
        float pix23egDphi_median = getMedian(pix23egDphiMed);
        float pix24egDphi_median = getMedian(pix24egDphiMed);
        float pix34egDphi_median = getMedian(pix34egDphiMed);

        float pix012Dphi_median = getMedian(pix012DphiMed);
        float pix013Dphi_median = getMedian(pix013DphiMed);
        float pix014Dphi_median = getMedian(pix014DphiMed);
        float pix023Dphi_median = getMedian(pix023DphiMed);
        float pix024Dphi_median = getMedian(pix024DphiMed);
        float pix034Dphi_median = getMedian(pix034DphiMed);
        float pix123Dphi_median = getMedian(pix123DphiMed);
        float pix124Dphi_median = getMedian(pix124DphiMed);
        float pix134Dphi_median = getMedian(pix134DphiMed);
        float pix234Dphi_median = getMedian(pix234DphiMed);

        if( debug ) cout << "Successfully get median" << endl;

        // calculate error
        int pix012Dphi_low = (int)(pix012Dphi_medianIndex - (0.668 * pix012Dphi_medianIndex));
        int pix013Dphi_low = (int)(pix013Dphi_medianIndex - (0.668 * pix013Dphi_medianIndex));
        int pix014Dphi_low = (int)(pix014Dphi_medianIndex - (0.668 * pix014Dphi_medianIndex));
        int pix023Dphi_low = (int)(pix023Dphi_medianIndex - (0.668 * pix023Dphi_medianIndex));
        int pix024Dphi_low = (int)(pix024Dphi_medianIndex - (0.668 * pix024Dphi_medianIndex));
        int pix034Dphi_low = (int)(pix034Dphi_medianIndex - (0.668 * pix034Dphi_medianIndex));
        int pix123Dphi_low = (int)(pix123Dphi_medianIndex - (0.668 * pix123Dphi_medianIndex));
        int pix124Dphi_low = (int)(pix124Dphi_medianIndex - (0.668 * pix124Dphi_medianIndex));
        int pix134Dphi_low = (int)(pix134Dphi_medianIndex - (0.668 * pix134Dphi_medianIndex));
        int pix234Dphi_low = (int)(pix234Dphi_medianIndex - (0.668 * pix234Dphi_medianIndex));

        int pix012Dphi_high = (int)(pix012Dphi_medianIndex + (0.668 * (pix012Dphi_size-1-pix012Dphi_medianIndex)));
        int pix013Dphi_high = (int)(pix013Dphi_medianIndex + (0.668 * (pix013Dphi_size-1-pix013Dphi_medianIndex)));
        int pix014Dphi_high = (int)(pix014Dphi_medianIndex + (0.668 * (pix014Dphi_size-1-pix014Dphi_medianIndex)));
        int pix023Dphi_high = (int)(pix023Dphi_medianIndex + (0.668 * (pix023Dphi_size-1-pix023Dphi_medianIndex)));
        int pix024Dphi_high = (int)(pix024Dphi_medianIndex + (0.668 * (pix024Dphi_size-1-pix024Dphi_medianIndex)));
        int pix034Dphi_high = (int)(pix034Dphi_medianIndex + (0.668 * (pix034Dphi_size-1-pix034Dphi_medianIndex)));
        int pix123Dphi_high = (int)(pix123Dphi_medianIndex + (0.668 * (pix123Dphi_size-1-pix123Dphi_medianIndex)));
        int pix124Dphi_high = (int)(pix124Dphi_medianIndex + (0.668 * (pix124Dphi_size-1-pix124Dphi_medianIndex)));
        int pix134Dphi_high = (int)(pix134Dphi_medianIndex + (0.668 * (pix134Dphi_size-1-pix134Dphi_medianIndex)));
        int pix234Dphi_high = (int)(pix234Dphi_medianIndex + (0.668 * (pix234Dphi_size-1-pix234Dphi_medianIndex)));

        float pix012Dphi_medianErr = pix012DphiMed.at(pix012Dphi_high) - pix012DphiMed.at(pix012Dphi_low);
        float pix013Dphi_medianErr = pix013DphiMed.at(pix013Dphi_high) - pix013DphiMed.at(pix013Dphi_low);
        float pix014Dphi_medianErr = pix014DphiMed.at(pix014Dphi_high) - pix014DphiMed.at(pix014Dphi_low);
        float pix023Dphi_medianErr = pix023DphiMed.at(pix023Dphi_high) - pix023DphiMed.at(pix023Dphi_low);
        float pix024Dphi_medianErr = pix024DphiMed.at(pix024Dphi_high) - pix024DphiMed.at(pix024Dphi_low);
        float pix034Dphi_medianErr = pix034DphiMed.at(pix034Dphi_high) - pix034DphiMed.at(pix034Dphi_low);
        float pix123Dphi_medianErr = pix123DphiMed.at(pix123Dphi_high) - pix123DphiMed.at(pix123Dphi_low);
        float pix124Dphi_medianErr = pix124DphiMed.at(pix124Dphi_high) - pix124DphiMed.at(pix124Dphi_low);
        float pix134Dphi_medianErr = pix134DphiMed.at(pix134Dphi_high) - pix134DphiMed.at(pix134Dphi_low);
        float pix234Dphi_medianErr = pix234DphiMed.at(pix234Dphi_high) - pix234DphiMed.at(pix234Dphi_low);

        float pix1egDphi_medianErr = 0.;
        float pix2egDphi_medianErr = 0.;
        float pix3egDphi_medianErr = 0.;
        float pix4egDphi_medianErr = 0.;

        float pix12egDphi_medianErr = 0.;
        float pix13egDphi_medianErr = 0.;
        float pix14egDphi_medianErr = 0.;
        float pix23egDphi_medianErr = 0.;
        float pix24egDphi_medianErr = 0.;
        float pix34egDphi_medianErr = 0.;

        pix1egDphi[nth] = pix1egDphi_median;
        pix2egDphi[nth] = pix2egDphi_median;
        pix3egDphi[nth] = pix3egDphi_median;
        pix4egDphi[nth] = pix4egDphi_median;

        pix1egDphiErr[nth] = pix1egDphi_medianErr/2.;
        pix2egDphiErr[nth] = pix2egDphi_medianErr/2.;
        pix3egDphiErr[nth] = pix3egDphi_medianErr/2.;
        pix4egDphiErr[nth] = pix4egDphi_medianErr/2.;

        pix12egDphi[nth] = pix12egDphi_median;
        pix13egDphi[nth] = pix13egDphi_median;
        pix14egDphi[nth] = pix14egDphi_median;
        pix23egDphi[nth] = pix23egDphi_median;
        pix24egDphi[nth] = pix24egDphi_median;
        pix34egDphi[nth] = pix34egDphi_median;

        pix12egDphiErr[nth] = pix12egDphi_medianErr/2.;
        pix13egDphiErr[nth] = pix13egDphi_medianErr/2.;
        pix14egDphiErr[nth] = pix14egDphi_medianErr/2.;
        pix23egDphiErr[nth] = pix23egDphi_medianErr/2.;
        pix24egDphiErr[nth] = pix24egDphi_medianErr/2.;
        pix34egDphiErr[nth] = pix34egDphi_medianErr/2.;

        pix012Dphi[nth] = pix012Dphi_median;
        pix013Dphi[nth] = pix013Dphi_median;
        pix014Dphi[nth] = pix014Dphi_median;
        pix023Dphi[nth] = pix023Dphi_median;
        pix024Dphi[nth] = pix024Dphi_median;
        pix034Dphi[nth] = pix034Dphi_median;
        pix123Dphi[nth] = pix123Dphi_median;
        pix124Dphi[nth] = pix124Dphi_median;
        pix134Dphi[nth] = pix134Dphi_median;
        pix234Dphi[nth] = pix234Dphi_median;

        pix012DphiErr[nth] = pix012Dphi_medianErr/2.;
        pix013DphiErr[nth] = pix013Dphi_medianErr/2.;
        pix014DphiErr[nth] = pix014Dphi_medianErr/2.;
        pix023DphiErr[nth] = pix023Dphi_medianErr/2.;
        pix024DphiErr[nth] = pix024Dphi_medianErr/2.;
        pix034DphiErr[nth] = pix034Dphi_medianErr/2.;
        pix123DphiErr[nth] = pix123Dphi_medianErr/2.;
        pix124DphiErr[nth] = pix124Dphi_medianErr/2.;
        pix134DphiErr[nth] = pix134Dphi_medianErr/2.;
        pix234DphiErr[nth] = pix234Dphi_medianErr/2.;

    } // Et scanning loop

    // save the fit parameters
    ofstream fit_result1;

    TString prefix1 = "ROI_"; 
    TString postfix1 = ".txt"; 

    //setTDRStyle();
    fit_result1.open(prefix1+eta_region[eta_-1]+postfix1);

    // draw eg-pix dphi 
    for(int i = 0; i < 4; i++){

        writeExtraText = true;
        extraText  = "Phase-2 simulation";
        gROOT->SetBatch();
        TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
        gStyle->SetOptStat(0);
        gStyle->SetPalette(1);
        c1->SetTopMargin(0.07);
        c1->SetRightMargin(0.05);
        c1->SetGridy();
        c1->SetGridx();
        TGaxis::SetMaxDigits(3);
        c1->cd();

        //

        TH2F* scatter_plot = NULL;
        if(i==0) scatter_plot = pix1egDphi_dist;
        if(i==1) scatter_plot = pix2egDphi_dist;
        if(i==2) scatter_plot = pix3egDphi_dist;
        if(i==3) scatter_plot = pix4egDphi_dist;

        if(scatter_plot != NULL){
            scatter_plot->GetYaxis()->SetRangeUser(-0.1,0.1);
            scatter_plot->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
            scatter_plot->GetYaxis()->SetTitle("#Delta#phi [rad.]");
            scatter_plot->SetMarkerColor(1);
            //scatter_plot->SetMarkerSize(.25);
            scatter_plot->SetMarkerStyle(6);
            scatter_plot->Draw("scat=1.");
        }
        if(scatter_plot == NULL) std::cout << "NULL" << std::endl;

        TGraphErrors* gr1 = new TGraphErrors(binSize, x, pixegDphi.at(i), 0, pixegDphiErr.at(i));

        TString nthsw;
        nthsw.Form("%d", i+1);
        gr1->SetMarkerColor(kBlue);
        gr1->SetLineColor(kBlue);
        gr1->SetMarkerSize(.8);
        gr1->SetMarkerStyle(20);
        //gr1->SetFillColorAlpha(kBlue, 0.4); // for error
        //gr1->Draw("samepE3");

        TF1 *median_fitFunc = new TF1("func","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])) )", 10., 100.);
        median_fitFunc->SetLineColor(kRed);
        median_fitFunc->SetLineStyle(7);
        median_fitFunc->SetLineWidth(4);

        gr1->Fit(median_fitFunc,"0");
        double ab[4]={0};
        median_fitFunc->GetParameters(ab);


        double x1[200],y1[200],x2[200],y2[200];
        for ( int j = 0;j<110;j++){
            x1[j]=j+10;
            x2[j]=0.5;
            y1[j]=ab[0]*pow(j+10,0) + ab[1]*pow(j+10,ab[2])*exp(-pow(j+10,ab[3]));
            y2[j]=0.017;  // eCal-Pix matching shadow width, cmssw 0.017
        }
        TGraphErrors* gr2 = new TGraphErrors(110,x1,y1,x2,y2); //0.0017 0.003

        gr2->SetLineColor(kYellow);
        gr2->SetFillColorAlpha(kYellow, 0.40);
        //    
        gr2->Draw("sameE3");
        median_fitFunc->Draw("lsame");

        TLegend* leg = new TLegend(0.3,0.7,0.5,0.9);
        //leg->AddEntry(scatter_plot, "Scattered #Delta#phi distribution as a function of L1 EG E_{T}","pl");
        //leg->AddEntry(gr1, "Median point of #Delta#phi","fp");
        leg->AddEntry(median_fitFunc,"Fit of median point", "l");
        leg->AddEntry(gr2,"Signal window", "f");
        leg->AddEntry((TObject*)0, eta_region[eta_-1], "");
        leg->SetTextFont(42);
        leg->SetTextSize(0.035);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();

        CMS_lumi( c1, 4, 0 );
        c1->Update();

        c1->SaveAs("PixEG"+nthsw+"_"+eta_region[eta_-1]+".png");

        fit_result1 << endl;
        fit_result1 << "if( region == " << eta_ << " && i == " << i << " ){" <<endl;
        for( int j=0; j < 4; j++){
            fit_result1 << "p[" << j << "] = " << median_fitFunc->GetParameter(j) << ";" << endl;
        }
        fit_result1 << "}" << endl;
        fit_result1 << endl;

        delete gr1;
        delete c1;
    }

    // save the fit parameters
    ofstream fit_result2;

    TString prefix2 = "EGmatching_";
    TString postfix2 = ".txt";

    fit_result2.open(prefix2+eta_region[eta_-1]+postfix2);

    const TString eg_pixpixSW[6] = {"12","13","14","23","24","34"};  
    // draw eg-pixpix dphi 
    for(int i = 0; i < 6; i++){

        writeExtraText = true;
        extraText  = "Phase-2 simulation";
        gROOT->SetBatch();
        TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
        gStyle->SetOptStat(0);
        gStyle->SetPalette(1);
        c1->SetTopMargin(0.07);
        c1->SetRightMargin(0.05);
        c1->SetGridy();
        c1->SetGridx();
        TGaxis::SetMaxDigits(3);
        c1->cd();

        //
        TGraphErrors* gr1 = new TGraphErrors(binSize, x, pixpixegDphi.at(i), 0, pixpixegDphiErr.at(i));

        TH2F* scatter_plot = NULL;
        if(i==0) scatter_plot = pix12egDphi_dist;
        if(i==1) scatter_plot = pix13egDphi_dist;
        if(i==2) scatter_plot = pix14egDphi_dist;
        if(i==3) scatter_plot = pix23egDphi_dist;
        if(i==4) scatter_plot = pix24egDphi_dist;
        if(i==5) scatter_plot = pix34egDphi_dist;

        if(scatter_plot != NULL){
            scatter_plot->GetYaxis()->SetRangeUser(-0.15,0.15);
            scatter_plot->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
            scatter_plot->GetYaxis()->SetTitle("#Delta#phi [rad.]");
            scatter_plot->SetMarkerColor(1);
            //scatter_plot->SetMarkerSize(.25);
            scatter_plot->SetMarkerStyle(6);
            scatter_plot->Draw("scat=1.");

        }

        TString nthsw;
        nthsw.Form("%d", i+1);
        gr1->SetLineColor(kBlue);
        gr1->SetMarkerColor(kBlue);
        gr1->SetMarkerSize(.8);
        gr1->SetMarkerStyle(20);
        gr1->SetFillColorAlpha(kBlue, 0.4);
        gr1->SetLineWidth(1);
        //gr1->Draw("samep");

        TF1 *median_fitFunc = new TF1("func","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])) )", 10., 100.);
        median_fitFunc->SetLineColor(kRed);
        median_fitFunc->SetLineStyle(7);
        median_fitFunc->SetLineWidth(4);

        gr1->Fit(median_fitFunc,"0");
        median_fitFunc->Draw("lsame");

        double ab[4]={0};
        median_fitFunc->GetParameters(ab);

        //MARK
        double x1[200],y1[200],x2[200],y2[200];
        for ( int j = 0;j<110;j++){
            x1[j]=j;
            x2[j]=0.5;
            y1[j]=ab[0]*pow(j,0) + ab[1]*pow(j,ab[2])*exp(-pow(j,ab[3]));
            y2[j]=0.017;  // eCal-PixPix matching shadow width  cmssw 0.017
        }
        TGraphErrors* gr2 = new TGraphErrors(110,x1,y1,x2,y2); //0.0017 0.003

        gr2->SetLineColor(kYellow);
        gr2->SetFillColorAlpha(kYellow, 0.40);
        //    

        gr2->Draw("sameE3");
        median_fitFunc->Draw("lsame");
        TLegend* leg = new TLegend(0.3,0.2,0.5,0.4);
        //leg->AddEntry(scatter_plot, "Scattered #Delta#phi distribution as a function of L1 EG E_{T}","pl");
        //leg->AddEntry(gr1, "Median point of #Delta#phi","lp");
        leg->AddEntry(median_fitFunc,"Fit of median point", "l");
        leg->AddEntry(gr2,"Signal window", "f");
        leg->AddEntry((TObject*)0, eta_region[eta_-1], "");
        leg->SetTextFont(42);
        leg->SetTextSize(0.035);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();

        CMS_lumi( c1, 4, 0 );
        c1->Update();

        //    c1->SaveAs("PixPixEG"+nthsw+"_"+eta_region[eta_-1]+".pdf");
        c1->SaveAs("PixPixEG"+nthsw+"_"+eta_region[eta_-1]+".png");

        fit_result2 << endl;
        fit_result2 << "if( region == " << eta_ << " && i == " << i << " ){" <<endl;
        for( int j=0; j < 4; j++){
            fit_result2 << "p[" << j << "] = " << median_fitFunc->GetParameter(j) << ";" << endl;
        }
        fit_result2 << "}" << endl;
        fit_result2 << endl;
        delete gr2;
        delete gr1;
        delete c1;
    }

    // save the fit parameters
    ofstream fit_result3;

    TString prefix3 = "Pixelmatching_";
    TString postfix3 = ".txt";

    fit_result3.open(prefix3+eta_region[eta_-1]+postfix3);

    const TString pixpix1_SW[10] = {"01","01","01","02","02","03","12","12","13","23"};  
    const TString pixpix2_SW[10] = {"12","13","14","23","24","34","23","24","34","34"};  

    // draw pix-pixpix dphi 
    for(int i = 0; i < 10; i++){

        writeExtraText = true;
        extraText  = "Phase-2 simulation";
        gROOT->SetBatch();
        TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
        gStyle->SetOptStat(0);
        gStyle->SetPalette(1);
        c1->SetTopMargin(0.07);
        c1->SetRightMargin(0.05);
        c1->SetGridy();
        c1->SetGridx();
        TGaxis::SetMaxDigits(4);

        //  
        TGraphErrors* gr1 = new TGraphErrors(binSize, x, pixpixDphi.at(i), 0, pixpixDphiErr.at(i));

        TH2F* scatter_plot = NULL;
        if(i==0) scatter_plot = pix012Dphi_dist;
        if(i==1) scatter_plot = pix013Dphi_dist;
        if(i==2) scatter_plot = pix014Dphi_dist;
        if(i==3) scatter_plot = pix023Dphi_dist;
        if(i==4) scatter_plot = pix024Dphi_dist;
        if(i==5) scatter_plot = pix034Dphi_dist;
        if(i==6) scatter_plot = pix123Dphi_dist;
        if(i==7) scatter_plot = pix124Dphi_dist;
        if(i==8) scatter_plot = pix134Dphi_dist;
        if(i==9) scatter_plot = pix234Dphi_dist;

        if(scatter_plot != NULL){

            scatter_plot->GetYaxis()->SetRangeUser(-0.015,0.015);
            scatter_plot->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
            scatter_plot->GetYaxis()->SetTitle("#Delta#phi [rad.]");
            scatter_plot->SetMarkerColor(1);
            //scatter_plot->SetMarkerSize(.25);
            scatter_plot->SetMarkerStyle(6);
            scatter_plot->Draw("scat=1.");
        }

        TString nthsw;
        nthsw.Form("%d", i+1);
        gr1->SetLineColor(kBlue);
        gr1->SetMarkerColor(kBlue);
        gr1->SetMarkerSize(.8);
        gr1->SetMarkerStyle(20);
        gr1->SetFillColorAlpha(kYellow, 0.4);
        gr1->SetLineWidth(1);
        //gr1->Draw("samep");

        TF1 *median_fitFunc = new TF1("func","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])) )", 10., 100.);

        median_fitFunc->SetParLimits(0, -0.05, 0.001);
        median_fitFunc->SetParLimits(1, -0.001, 5);
        median_fitFunc->SetParLimits(2, -2., 0.2);
        median_fitFunc->SetParLimits(3, -0.7, 0.5);

        median_fitFunc->SetLineColor(kRed);
        median_fitFunc->SetLineStyle(7);
        median_fitFunc->SetLineWidth(4);

        gr1->Fit(median_fitFunc,"0");

        double ab[4]={0};
        median_fitFunc->GetParameters(ab);


        double x1[200],y1[200],x2[200],y2[200];
        for ( int j = 0;j<110;j++){
            x1[j]=j+10;
            x2[j]=0.5;
            y1[j]=ab[0]*pow(j+10,0) + ab[1]*pow(j+10,ab[2])*exp(-pow(j+10,ab[3]));
            y2[j]=0.003;  // Pix-Pix matching shadow width  cmssw 0.003
        }

        TGraphErrors* gr2 = new TGraphErrors(110,x1,y1,x2,y2); //0.0017 0.003

        gr2->SetLineColor(kYellow);
        gr2->SetFillColorAlpha(kYellow, 0.40);
        gr2->Draw("sameE3");  //pix-pix
        median_fitFunc->Draw("lsame");

        TLegend* leg = new TLegend(0.2,0.2,0.4,0.4);
        //leg->AddEntry(scatter_plot, "Scattered #Delta#phi distribution as a function of L1 EG E_{T}","peZ");
        //leg->AddEntry(scatter_plot, "Scattered #Delta#phi distribution as a function of L1 EG E_{T}","pl");
        //leg->AddEntry(gr1, "Median point of #Delta#phi", "lp");
        //leg->AddEntry(gr1, "Median point of #Delta#phi", "pe3");
        leg->AddEntry(median_fitFunc,"Fit of median point", "l");
        leg->AddEntry(gr2,"Signal window", "f");
        leg->AddEntry((TObject*)0, eta_region[eta_-1], "");
        leg->SetTextFont(42);
        leg->SetTextSize(0.035);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();

        CMS_lumi( c1, 4, 0 );
        c1->Update();

        c1->SaveAs("PixPix"+nthsw+"_"+eta_region[eta_-1]+".png");

        fit_result3 << endl;
        fit_result3 << "if( region == " << eta_ << " && i == " << i << " ){" <<endl;
        for( int j=0; j < 4; j++){
            fit_result3 << "p[" << j << "] = " << median_fitFunc->GetParameter(j) << ";" << endl;
        }
        fit_result3 << "}" << endl;
        fit_result3 << endl;
        delete gr2;
        delete gr1;
        delete c1; 
    }

    delete pix1egDphi_dist;
    delete pix2egDphi_dist;
    delete pix3egDphi_dist;
    delete pix4egDphi_dist;

    delete pix12egDphi_dist;
    delete pix13egDphi_dist;
    delete pix14egDphi_dist;
    delete pix23egDphi_dist;
    delete pix24egDphi_dist;
    delete pix34egDphi_dist;

    delete pix012Dphi_dist;
    delete pix013Dphi_dist;
    delete pix014Dphi_dist;
    delete pix023Dphi_dist;
    delete pix024Dphi_dist;
    delete pix034Dphi_dist;
    delete pix123Dphi_dist;
    delete pix124Dphi_dist;
    delete pix134Dphi_dist;
    delete pix234Dphi_dist;
    std::cout << "DONE" << std::endl;

} 

