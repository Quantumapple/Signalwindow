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

        //static bool EtComp(const L1eGamma &t1, const L1eGamma &t2)
        //    return t1.EgEt < t2.EgEt;
};

double getMedian(const std::vector<float> &vec)
{
    double median=0.;

    //odd number, definate median
    if(vec.size() % 2 !=0) {
        int middleNr = (vec.size()+1)/2;
        median = vec[middleNr];
    }
    else{ //even number, take median as halfway between the two middle values
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
    }
    else{ //even number, take median as halfway between the two middle values
        int middleNr = (vec.size()+1)/2;
        return middleNr;
    }
}

void Make2Dplots::Loop(int eta_ = 1)
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    bool debug = false;
    const TString eta_region[6] = {"EtaRegion1", "EtaRegion2", "EtaRegion3", "EtaRegion4", "EtaRegion5", "EtaRegion6"};

    TH2F* pix12egDeta_dist = new TH2F("pix12egDeta_dist","pix12egDeta_dist", 90,10,100,100,-0.2,0.2);
    TH2F* pix13egDeta_dist = new TH2F("pix13egDeta_dist","pix13egDeta_dist", 90,10,100,100,-0.2,0.2);
    TH2F* pix14egDeta_dist = new TH2F("pix14egDeta_dist","pix14egDeta_dist", 90,10,100,100,-0.2,0.2);
    TH2F* pix23egDeta_dist = new TH2F("pix23egDeta_dist","pix23egDeta_dist", 90,10,100,100,-0.2,0.2);
    TH2F* pix24egDeta_dist = new TH2F("pix24egDeta_dist","pix24egDeta_dist", 90,10,100,100,-0.2,0.2);
    TH2F* pix34egDeta_dist = new TH2F("pix34egDeta_dist","pix34egDeta_dist", 90,10,100,100,-0.2,0.2);
    TH2F* pix123Deta_dist = new TH2F("pix123Deta_dist","pix123Deta_dist", 90,10,100,100,-0.02,0.02);
    TH2F* pix124Deta_dist = new TH2F("pix124Deta_dist","pix124Deta_dist", 90,10,100,100,-0.02,0.02);
    TH2F* pix134Deta_dist = new TH2F("pix134Deta_dist","pix134Deta_dist", 90,10,100,100,-0.02,0.02);
    TH2F* pix234Deta_dist = new TH2F("pix234Deta_dist","pix234Deta_dist", 90,10,100,100,-0.02,0.02);

    const int binSize = 90;
    vector<L1eGamma> pix12egDeta_;
    vector<L1eGamma> pix13egDeta_;
    vector<L1eGamma> pix14egDeta_;
    vector<L1eGamma> pix23egDeta_;
    vector<L1eGamma> pix24egDeta_;
    vector<L1eGamma> pix34egDeta_;

    vector<L1eGamma> pix123Deta_;
    vector<L1eGamma> pix124Deta_;
    vector<L1eGamma> pix134Deta_;
    vector<L1eGamma> pix234Deta_;

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if( eta_ == 1 && fabs(ntEgEta->at(0)) > 0.8 ) continue; // test with barrel only
        if( eta_ == 2 && (fabs(ntEgEta->at(0)) < 0.8 || fabs(ntEgEta->at(0)) > 1.4)) continue; 
        if( eta_ == 3 && (fabs(ntEgEta->at(0)) < 1.4 || fabs(ntEgEta->at(0)) > 1.7)) continue; 
        if( eta_ == 4 && (fabs(ntEgEta->at(0)) < 1.7 || fabs(ntEgEta->at(0)) > 2.1)) continue; 
        if( eta_ == 5 && (fabs(ntEgEta->at(0)) < 2.1 || fabs(ntEgEta->at(0)) > 2.7)) continue; 
        if( eta_ == 6 && (fabs(ntEgEta->at(0)) < 2.7 || fabs(ntEgEta->at(0)) > 3.)) continue; 

        for(unsigned long i = 0; i < ntPix12EGdeta->size(); i++) { 
            if(fabs(ntPix12EGdeta->at(i)) > 0.2 ) continue;
            pix12egDeta_.push_back( L1eGamma( ntEgEt->at(0), ntPix12EGdeta->at(i) ) ); 
            pix12egDeta_dist->Fill(ntEgEt->at(0), ntPix12EGdeta->at(i));
        }
        for(unsigned long i = 0; i < ntPix13EGdeta->size(); i++) {
            if(fabs(ntPix13EGdeta->at(i)) > 0.2 ) continue;
            pix13egDeta_.push_back( L1eGamma( ntEgEt->at(0), ntPix13EGdeta->at(i) ) ); 
            pix13egDeta_dist->Fill(ntEgEt->at(0), ntPix13EGdeta->at(i));
        }
        for(unsigned long i = 0; i < ntPix14EGdeta->size(); i++) {
            if(fabs(ntPix14EGdeta->at(i)) > 0.2 ) continue;
            pix14egDeta_.push_back( L1eGamma( ntEgEt->at(0), ntPix14EGdeta->at(i) ) ); 
            pix14egDeta_dist->Fill(ntEgEt->at(0), ntPix14EGdeta->at(i));
        }
        for(unsigned long i = 0; i < ntPix23EGdeta->size(); i++) { 
            if(fabs(ntPix23EGdeta->at(i)) > 0.2 ) continue;
            pix23egDeta_.push_back( L1eGamma( ntEgEt->at(0), ntPix23EGdeta->at(i) ) ); 
            pix23egDeta_dist->Fill(ntEgEt->at(0), ntPix23EGdeta->at(i));
        }
        for(unsigned long i = 0; i < ntPix24EGdeta->size(); i++) { 
            if(fabs(ntPix24EGdeta->at(i)) > 0.2 ) continue;
            pix24egDeta_.push_back( L1eGamma( ntEgEt->at(0), ntPix24EGdeta->at(i) ) ); 
            pix24egDeta_dist->Fill(ntEgEt->at(0), ntPix24EGdeta->at(i));
        }
        for(unsigned long i = 0; i < ntPix34EGdeta->size(); i++) {
            if(fabs(ntPix34EGdeta->at(i)) > 0.2 ) continue;
            pix34egDeta_.push_back( L1eGamma( ntEgEt->at(0), ntPix34EGdeta->at(i) ) ); 
            pix34egDeta_dist->Fill(ntEgEt->at(0), ntPix34EGdeta->at(i));
        }
        for(unsigned long i = 0; i < ntPix123deta->size(); i++) {
            if(fabs(ntPix123deta->at(i)) > 0.015 ) continue;
            pix123Deta_.push_back( L1eGamma( ntEgEt->at(0), ntPix123deta->at(i) ) ); 
            pix123Deta_dist->Fill(ntEgEt->at(0), ntPix123deta->at(i));
        }
        for(unsigned long i = 0; i < ntPix124deta->size(); i++) {
            if(fabs(ntPix124deta->at(i)) > 0.015 ) continue;
            pix124Deta_.push_back( L1eGamma( ntEgEt->at(0), ntPix124deta->at(i) ) ); 
            pix124Deta_dist->Fill(ntEgEt->at(0), ntPix124deta->at(i));
        }
        for(unsigned long i = 0; i < ntPix134deta->size(); i++) {
            if(fabs(ntPix134deta->at(i)) > 0.015 ) continue;
            pix134Deta_.push_back( L1eGamma( ntEgEt->at(0), ntPix134deta->at(i) ) ); 
            pix134Deta_dist->Fill(ntEgEt->at(0), ntPix134deta->at(i));
        }
        for(unsigned long i = 0; i < ntPix234deta->size(); i++) {
            if(fabs(ntPix234deta->at(i)) > 0.015 ) continue;
            pix234Deta_.push_back( L1eGamma( ntEgEt->at(0), ntPix234deta->at(i) ) ); 
            pix234Deta_dist->Fill(ntEgEt->at(0), ntPix234deta->at(i));
        }

    }// event loop

    if( debug ) cout << "Pass event loop" << endl;

    float x[binSize];

    // eg-pixelpixel deta
    float pix12egDeta[binSize], pix12egDetaErr[binSize];
    float pix13egDeta[binSize], pix13egDetaErr[binSize];
    float pix14egDeta[binSize], pix14egDetaErr[binSize];
    float pix23egDeta[binSize], pix23egDetaErr[binSize];
    float pix24egDeta[binSize], pix24egDetaErr[binSize];
    float pix34egDeta[binSize], pix34egDetaErr[binSize];

    // pixel-pixel deta

    float pix123Deta[binSize], pix123DetaErr[binSize];
    float pix124Deta[binSize], pix124DetaErr[binSize];
    float pix134Deta[binSize], pix134DetaErr[binSize];
    float pix234Deta[binSize], pix234DetaErr[binSize];

    vector<float*> pixpixegDeta, pixpixegDetaErr;
    vector<float*> pixpixDeta, pixpixDetaErr;

    pixpixegDeta.push_back(pix12egDeta);
    pixpixegDeta.push_back(pix13egDeta);
    pixpixegDeta.push_back(pix14egDeta);
    pixpixegDeta.push_back(pix23egDeta);
    pixpixegDeta.push_back(pix24egDeta);
    pixpixegDeta.push_back(pix34egDeta);

    pixpixegDetaErr.push_back(pix12egDetaErr);
    pixpixegDetaErr.push_back(pix13egDetaErr);
    pixpixegDetaErr.push_back(pix14egDetaErr);
    pixpixegDetaErr.push_back(pix23egDetaErr);
    pixpixegDetaErr.push_back(pix24egDetaErr);
    pixpixegDetaErr.push_back(pix34egDetaErr);

    pixpixDeta.push_back(pix123Deta);
    pixpixDeta.push_back(pix124Deta);
    pixpixDeta.push_back(pix134Deta);
    pixpixDeta.push_back(pix234Deta);

    pixpixDetaErr.push_back(pix123DetaErr);
    pixpixDetaErr.push_back(pix124DetaErr);
    pixpixDetaErr.push_back(pix134DetaErr);
    pixpixDetaErr.push_back(pix234DetaErr);

    vector<float> pix12egDetaMed; 
    vector<float> pix13egDetaMed; 
    vector<float> pix14egDetaMed; 
    vector<float> pix23egDetaMed; 
    vector<float> pix24egDetaMed; 
    vector<float> pix34egDetaMed; 

    vector<float> pix123DetaMed; 
    vector<float> pix124DetaMed; 
    vector<float> pix134DetaMed; 
    vector<float> pix234DetaMed; 

    //Fill
    for(Int_t nth = 0; nth < binSize; nth++)
    {
        if( debug ) cout << "Which bin: " << nth << " between " << 10.+nth << " and " << 11.+nth << endl;

        x[nth] = 10.5 + float(nth);
        pix12egDetaMed.clear(); 
        pix13egDetaMed.clear(); 
        pix14egDetaMed.clear(); 
        pix23egDetaMed.clear(); 
        pix24egDetaMed.clear(); 
        pix34egDetaMed.clear(); 

        pix123DetaMed.clear(); 
        pix124DetaMed.clear(); 
        pix134DetaMed.clear(); 
        pix234DetaMed.clear(); 

        // Pixel + egamma matching
        for(Int_t j = 0; j < pix12egDeta_.size(); j++)
        {
            float tempEt = pix12egDeta_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix12egDetaMed.push_back(pix12egDeta_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix13egDeta_.size(); j++)
        {
            float tempEt = pix13egDeta_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix13egDetaMed.push_back(pix13egDeta_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix14egDeta_.size(); j++)
        {
            float tempEt = pix14egDeta_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix14egDetaMed.push_back(pix14egDeta_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix23egDeta_.size(); j++)
        {
            float tempEt = pix23egDeta_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix23egDetaMed.push_back(pix23egDeta_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix24egDeta_.size(); j++)
        {
            float tempEt = pix24egDeta_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix24egDetaMed.push_back(pix24egDeta_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix34egDeta_.size(); j++)
        {
            float tempEt = pix34egDeta_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix34egDetaMed.push_back(pix34egDeta_[j].EgdPhi);
        }

        // Pixel + pixel matching
        for(Int_t j = 0; j < pix123Deta_.size(); j++)
        {
            float tempEt = pix123Deta_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix123DetaMed.push_back(pix123Deta_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix124Deta_.size(); j++)
        {
            float tempEt = pix124Deta_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix124DetaMed.push_back(pix124Deta_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix134Deta_.size(); j++)
        {
            float tempEt = pix134Deta_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix134DetaMed.push_back(pix134Deta_[j].EgdPhi);
        }
        for(Int_t j = 0; j < pix234Deta_.size(); j++)
        {
            float tempEt = pix234Deta_[j].egEt;
            if( tempEt < 10.+nth || tempEt > 11.+nth ) continue;
            pix234DetaMed.push_back(pix234Deta_[j].EgdPhi);
        }

        if( debug ) cout << "Fill new vecotrs!" << endl;

        std::sort(pix12egDetaMed.begin(), pix12egDetaMed.end()); 
        std::sort(pix13egDetaMed.begin(), pix13egDetaMed.end()); 
        std::sort(pix14egDetaMed.begin(), pix14egDetaMed.end()); 
        std::sort(pix23egDetaMed.begin(), pix23egDetaMed.end()); 
        std::sort(pix24egDetaMed.begin(), pix24egDetaMed.end()); 
        std::sort(pix34egDetaMed.begin(), pix34egDetaMed.end()); 

        std::sort(pix123DetaMed.begin(), pix123DetaMed.end()); 
        std::sort(pix124DetaMed.begin(), pix124DetaMed.end()); 
        std::sort(pix134DetaMed.begin(), pix134DetaMed.end()); 
        std::sort(pix234DetaMed.begin(), pix234DetaMed.end()); 

        if( debug ) cout << "vector sorting is done!" << endl;

        int pix12egDeta_size = pix12egDetaMed.size();
        int pix13egDeta_size = pix13egDetaMed.size();
        int pix14egDeta_size = pix14egDetaMed.size();
        int pix23egDeta_size = pix23egDetaMed.size();
        int pix24egDeta_size = pix24egDetaMed.size();
        int pix34egDeta_size = pix34egDetaMed.size();

        int pix123Deta_size = pix123DetaMed.size();
        int pix124Deta_size = pix124DetaMed.size();
        int pix134Deta_size = pix134DetaMed.size();
        int pix234Deta_size = pix234DetaMed.size();

        // get index for median

        int pix12egDeta_medianIndex = getMedianIndex(pix12egDetaMed);
        int pix13egDeta_medianIndex = getMedianIndex(pix13egDetaMed);
        int pix14egDeta_medianIndex = getMedianIndex(pix14egDetaMed);
        int pix23egDeta_medianIndex = getMedianIndex(pix23egDetaMed);
        int pix24egDeta_medianIndex = getMedianIndex(pix24egDetaMed);
        int pix34egDeta_medianIndex = getMedianIndex(pix34egDetaMed);

        int pix123Deta_medianIndex = getMedianIndex(pix123DetaMed);
        int pix124Deta_medianIndex = getMedianIndex(pix124DetaMed);
        int pix134Deta_medianIndex = getMedianIndex(pix134DetaMed);
        int pix234Deta_medianIndex = getMedianIndex(pix234DetaMed);

        if( debug ) cout << "Successfully get median index" << endl;

        // get median

        float pix12egDeta_median = getMedian(pix12egDetaMed);
        float pix13egDeta_median = getMedian(pix13egDetaMed);
        float pix14egDeta_median = getMedian(pix14egDetaMed);
        float pix23egDeta_median = getMedian(pix23egDetaMed);
        float pix24egDeta_median = getMedian(pix24egDetaMed);
        float pix34egDeta_median = getMedian(pix34egDetaMed);

        float pix123Deta_median = getMedian(pix123DetaMed);
        float pix124Deta_median = getMedian(pix124DetaMed);
        float pix134Deta_median = getMedian(pix134DetaMed);
        float pix234Deta_median = getMedian(pix234DetaMed);

        if( debug ) cout << "Successfully get median" << endl;

        pix12egDeta[nth] = pix12egDeta_median;
        pix13egDeta[nth] = pix13egDeta_median;
        pix14egDeta[nth] = pix14egDeta_median;
        pix23egDeta[nth] = pix23egDeta_median;
        pix24egDeta[nth] = pix24egDeta_median;
        pix34egDeta[nth] = pix34egDeta_median;

        pix12egDetaErr[nth] = 0.;
        pix13egDetaErr[nth] = 0.;
        pix14egDetaErr[nth] = 0.;
        pix23egDetaErr[nth] = 0.;
        pix24egDetaErr[nth] = 0.;
        pix34egDetaErr[nth] = 0.;

        pix123Deta[nth] = pix123Deta_median;
        pix124Deta[nth] = pix124Deta_median;
        pix134Deta[nth] = pix134Deta_median;
        pix234Deta[nth] = pix234Deta_median;

        pix123DetaErr[nth] = 0.;
        pix124DetaErr[nth] = 0.;
        pix134DetaErr[nth] = 0.;
        pix234DetaErr[nth] = 0.;

    } // Et scanning loop

    const TString eg_pixpixSW[6] = {"12","13","14","23","24","34"};
    // draw eg-pixpix deta 
    for(int i = 0; i < 6; i++){

        setTDRStyle();
        writeExtraText = true;
        extraText  = "Phase-2 simulation";
        gROOT->SetBatch();
        TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
        gStyle->SetOptStat(0);
        gStyle->SetPalette(1);
        c1->SetTopMargin(0.07);
        c1->SetGridy();
        c1->SetGridx();
        TGaxis::SetMaxDigits(3);

        
        TGraphErrors* gr1 = new TGraphErrors(binSize, x, pixpixegDeta.at(i), 0, pixpixegDetaErr.at(i)); //0.015 0.01

        TH2F* scatter_plot = NULL;
        if(i==0) scatter_plot = pix12egDeta_dist;
        if(i==1) scatter_plot = pix13egDeta_dist;
        if(i==2) scatter_plot = pix14egDeta_dist;
        if(i==3) scatter_plot = pix23egDeta_dist;
        if(i==4) scatter_plot = pix24egDeta_dist;
        if(i==5) scatter_plot = pix34egDeta_dist;

        if(scatter_plot != NULL){
            scatter_plot->GetYaxis()->SetRangeUser(-0.05,0.05);//0.15
            scatter_plot->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
            scatter_plot->GetYaxis()->SetTitle("#Delta#eta [rad.]");
            scatter_plot->GetYaxis()->SetNoExponent();
            scatter_plot->SetMarkerColor(1);
            scatter_plot->SetMarkerStyle(6);
            scatter_plot->Draw("scat=1.");

        }

        TString nthsw;
        nthsw.Form("%d", i+1);
        gr1->SetTitle("#Delta#eta(L1 EG, pixel-pixel "+nthsw+")");
        gr1->GetYaxis()->SetRangeUser(-0.05,0.05);//0.15
        gr1->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
        gr1->GetYaxis()->SetTitle("#Delta#eta");
        gr1->GetYaxis()->SetNoExponent();
        gr1->SetLineColor(kYellow);
        gr1->SetMarkerColor(kYellow);
        gr1->SetMarkerSize(1.);
        gr1->SetMarkerStyle(20);
        gr1->SetFillColorAlpha(kYellow, 0.40);
        gr1->SetLineWidth(1);
        TF1 *median_fitFunc = new TF1("func","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]) )", 10., 100.);
        median_fitFunc->SetLineColor(kRed);
        median_fitFunc->SetLineStyle(7);
        median_fitFunc->SetLineWidth(4);

        median_fitFunc->GetYaxis()->SetNoExponent();
        gr1->Fit(median_fitFunc,"0");

        double ab[4]={0};
        median_fitFunc->GetParameters(ab);
        cout<<ab[0]*pow(2,0) + ab[1]*pow(2,ab[2])*exp(-pow(2,ab[3]))<<endl;


        double x1[200],y1[200],x2[200],y2[200];
        for ( int j = 10;j<120;j++){
            x1[j-10]=j;
            x2[j-10]=0.5;
            y1[j-10]=ab[0]*pow(j,0) + ab[1]*pow(j,ab[2])*exp(-pow(j,ab[3]));
            y2[j-10]=0.01;
        }
        TGraphErrors* gr2 = new TGraphErrors(110,x1,y1,x2,y2); //0.0017 0.003

        gr2->SetLineColor(kYellow);
        gr2->SetFillColorAlpha(kYellow, 0.40);
        gr2->GetYaxis()->SetNoExponent();
        gr2->Draw("sameE3");
        median_fitFunc->Draw("lsame");

        TLegend* leg = new TLegend(0.3,0.2,0.5,0.4);
        leg->AddEntry(median_fitFunc,"Fit of Median point", "l");
        leg->AddEntry(gr2,"Signal window", "f");
        leg->AddEntry((TObject*)0, eta_region[eta_-1], "");
        leg->SetTextFont(42);
        leg->SetTextSize(0.035);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();

        CMS_lumi( c1, 4, 0 );
        c1->Update();

        c1->SaveAs("PixPixEG"+nthsw+"_"+eta_region[eta_-1]+"_eta.png");
        c1->SaveAs("testEG"+nthsw+"_"+eta_region[eta_-1]+"_eta.root");


        delete gr2;
        delete gr1;
        delete c1;
    }

    const TString pixpix1_SW[4] = {"12","12","13","23"};
    const TString pixpix2_SW[4] = {"23","24","34","34"};
    
    // draw pix-pixpix deta 
    for(int i = 0; i < 4; i++){

        setTDRStyle();
        writeExtraText = true;
        extraText  = "Phase-2 simulation";
        gROOT->SetBatch();
        TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
        gStyle->SetOptStat(0);
        gStyle->SetPalette(1);
        c1->SetTopMargin(0.07);
        c1->SetGridy();
        c1->SetGridx();
        TGaxis::SetMaxDigits(4);

        //  
        TGraphErrors* gr1 = new TGraphErrors(binSize, x, pixpixDeta.at(i), 0, pixpixDetaErr.at(i)); //0.0017 0.003
        TH2F* scatter_plot = NULL;
        if(i==0) scatter_plot = pix123Deta_dist;
        if(i==1) scatter_plot = pix124Deta_dist;
        if(i==2) scatter_plot = pix134Deta_dist;
        if(i==3) scatter_plot = pix234Deta_dist;

        if(scatter_plot != NULL){

            scatter_plot->GetYaxis()->SetRangeUser(-0.015,0.015);
            scatter_plot->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
            scatter_plot->GetYaxis()->SetTitle("#Delta#eta [rad.]");
            scatter_plot->SetMarkerColor(1);
            scatter_plot->SetMarkerStyle(6);
            scatter_plot->Draw("scat=1.");
        }


        TString nthsw;
        nthsw.Form("%d", i+7);
        gr1->SetTitle("#Delta#eta(pixel-pixel, pixel-pixel "+nthsw+")");
        gr1->GetYaxis()->SetRangeUser(-0.015,0.015);
        gr1->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
        gr1->GetYaxis()->SetTitle("#Delta#eta");
        gr1->SetLineColor(kYellow);
        gr1->SetMarkerColor(kYellow);
        gr1->SetMarkerSize(1.);
        gr1->SetMarkerStyle(20);
        gr1->SetLineWidth(1);

        TF1 *median_fitFunc = new TF1("func","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])) )", 10., 100.);

        median_fitFunc->SetParLimits(0, -0.05, 0.001);
        median_fitFunc->SetParLimits(1, -0.001, 5);
        median_fitFunc->SetParLimits(2, -2., 0.2);
        median_fitFunc->SetParLimits(3, -0.7, 0.5);

        median_fitFunc->SetLineColor(kRed);
        median_fitFunc->SetLineStyle(7);
        median_fitFunc->SetLineWidth(4);

        gr1->Fit(median_fitFunc,"0");
        median_fitFunc->Draw("lsame");

        double ab[4]={0};
        median_fitFunc->GetParameters(ab);


        double x1[200],y1[200],x2[200],y2[200];
        for ( int j = 10;j<120;j++){
            x1[j]=j;
            x2[j]=0.5;
            y1[j]=ab[0]*pow(j,0) + ab[1]*pow(j,ab[2])*exp(-pow(j,ab[3]));
            y2[j]=0.0017;
        }
        TGraphErrors* gr2 = new TGraphErrors(110,x1,y1,x2,y2); //0.0017 0.003

        gr2->SetFillColorAlpha(kYellow, 0.40);
            
        gr2->Draw("sameE3");
        TLegend* leg = new TLegend(0.2,0.2,0.4,0.4);
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

        c1->SaveAs("PixPix"+nthsw+"_"+eta_region[eta_-1]+"_eta.png");


        delete gr2;
        delete gr1;
        delete c1; 

    }
    delete pix12egDeta_dist;
    delete pix13egDeta_dist;
    delete pix14egDeta_dist;
    delete pix23egDeta_dist;
    delete pix24egDeta_dist;
    delete pix34egDeta_dist;

    delete pix123Deta_dist;
    delete pix124Deta_dist;
    delete pix134Deta_dist;
    delete pix234Deta_dist;
    std::cout << "DONE" << std::endl;
}
