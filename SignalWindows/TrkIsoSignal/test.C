#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <string>
#include <TLorentzVector.h>
#include <TGraph.h>
#include "../../../trkIsoPtFitFunction.h"

using namespace std;

void test::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    float dr_cut = 0.1;

    bit1 = 0x1;
    bit2 = 0x1;

    debug = false;

    const double EM_PiX_dphi_width_[9] = {0.02, 0.03, 0.03, 0.033, 0.035, 0.037, 0.039, 0.04, 0.05}; // 
    const double EM_PiX_deta_width_[9] = {0.01, 0.015, 0.015, 0.033, 0.035, 0.037, 0.039, 0.04, 0.05}; // 

    const double PiX_PiX_dphi_width_[9] = {0.0017, 0.003, 0.003, 0.0033, 0.0035, 0.0037, 0.0039, 0.004, 0.005};
    const double PiX_PiX_deta_width_[9] = {0.0017, 0.003, 0.005, 0.0033, 0.0035, 0.0037, 0.0039, 0.004, 0.005};

    TH1F *aa1 = new TH1F("aa1","z_{vtx}-z_{gen}",200,-0.04,0.04);
    TH1F *aa2 = new TH1F("aa2","z_{vtx}-z_{gen}",200,-0.04,0.04);
    TH1F *aa3 = new TH1F("aa3","z_{vtx}-z_{gen}",200,-0.06,0.06);
    TH1F *aa4 = new TH1F("aa4","z_{vtx}-z_{gen}",200,-0.1,0.1);
    TH1F *aa5 = new TH1F("aa5","z_{vtx}-z_{gen}",200,-0.2,0.2);
    TH1F *aa6 = new TH1F("aa6","z_{vtx}-z_{gen}",200,-0.2,0.2);
    
    TH1F *h1 = new TH1F("h1","; Isolation value; ",100,0,1);
    TH1F *h2 = new TH1F("h2","; Isolation value; ",100,0,1);
    TH1F *h3 = new TH1F("h3","; Isolation value; ",100,0,1);
    TH1F *h4 = new TH1F("h4","; Isolation value; ",100,0,1);
    TH1F *h5 = new TH1F("h5","; Isolation value; ",100,0,1);
    TH1F *h6 = new TH1F("h6","; Isolation value; ",100,0,1);
    
    TH1F *nh1 = new TH1F("nh1","; Number of pixel tracks; ",10,0,10);
    TH1F *nh2 = new TH1F("nh2","; Number of pixel tracks; ",10,0,10);
    TH1F *nh3 = new TH1F("nh3","; Number of pixel tracks; ",10,0,10);
    TH1F *nh4 = new TH1F("nh4","; Number of pixel tracks; ",10,0,10);
    TH1F *nh5 = new TH1F("nh5","; Number of pixel tracks; ",10,0,10);
    TH1F *nh6 = new TH1F("nh6","; Number of pixel tracks; ",10,0,10);

    TH1F *h_tracketa_all = new TH1F("h_tracketa_all","All track eta",160,-4.,4.);
    TH1F *h_tracketa_1 = new TH1F("h_tracketa_1","leading track eta",160,-4.,4.);
    TH1F *h_tracketa_2 = new TH1F("h_tracketa_2","second track eta",160,-4.,4.);
    TH1F *h_tracketa_3 = new TH1F("h_tracketa_3","third track eta",160,-4.,4.);
    TH1F *h_tracketa_4 = new TH1F("h_tracketa_4","fourth track eta",160,-4.,4.);

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if (!(jentry%1) ) cout << "Processing entry " << jentry << "/" << nentries << endl;
        FillCutFlow("NoCut", 1.);

        matchedEgEt  = -999.;
        matchedEgEta = -999.;
        matchedEgPhi = -999.;

        ntnEg2 = 0;
        ntEgEt.clear();
        ntEgEta.clear();
        ntEgPhi.clear();
        
        ntnEg2Iso = 0;
        ntEgEtIso.clear();
        ntEgEtaIso.clear();
        ntEgPhiIso.clear();

        PiXTRKbit.clear();
        trigger_bit_width.clear();
        trigger_bit_width_iso.clear();
        pix_comb.clear();

        ntCl_match.clear();
        withoutEM_match.clear();
        withEM_match.clear();

        NumOfTrks.clear();
        IsoValue.clear();
        track_pT1.clear();
        track_pT2.clear();
        track_pT3.clear();
        track_pT4.clear();
        track_pT5.clear();
        track_pT6.clear();
        track_pT7.clear();
        track_pT8.clear();
        track_pT9.clear();
        track_pT10.clear();
        track_pT11.clear();
        track_pT12.clear();
        track_pT13.clear();
        track_pT14.clear();
        track_pT15.clear();
        track_eta1.clear();
        track_eta2.clear();
        track_eta3.clear();

        nt_genPhi = propgenElPartPhi->at(0);
        nt_genEta = propgenElPartEta->at(0);
        nt_genPt = propgenElPartPt->at(0);
        float nt_genZ = propgenElPartZ->at(0);

        float tempDR = 999.;
        int   indx = -1;
        int   indx_closestEg = -1;

        //find closest egamma object to the gen electron
        EgN=egCrysClusterEt->size();

        float closest_dr = 9999.;
        int closest_eg = 0;
        int egCount = 0;

        // loop over barrel egamma objects
        for(int i=0; i < EgN;i++){

            float dPhi = deltaPhi(propgenElPartPhi->at(0), egCrysClusterPhi->at(i));

            float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-egCrysClusterEta->at(i),2));
            if(egCrysClusterEt->at(i) < 10) continue;
            egCount++;
            if(current_dr < closest_dr){
                closest_dr = current_dr;
                closest_eg = i;
            }
        }// end of loop to find the closest egamma to gen electron 
        pix_comb_ = 0x0;

        nPix123_segments = 0;
        nPix124_segments = 0;
        nPix134_segments = 0;
        nPix234_segments = 0;

        // find egamma objects passing pixtrk signal windows
        if(closest_dr < dr_cut && closest_dr != 9999.){

            debug = false;

            EgEt =egCrysClusterEt ->at(closest_eg);
            EgEta=egCrysClusterEta->at(closest_eg);
            EgPhi=egCrysClusterPhi->at(closest_eg);

            cout << "Gen matched" << endl;
            cout << "et: " << EgEt << endl;
            cout << "eta: " << EgEta << endl;
            cout << "phi: " << EgPhi << endl;
            cout << endl;


            float EgGx = egCrysClusterGx->at(closest_eg);
            float EgGy = egCrysClusterGy->at(closest_eg);
            float EgGz = egCrysClusterGz->at(closest_eg);
            emvector.SetXYZ(EgGx,EgGy,EgGz);

            tempDR = closest_dr; 
            matchedEgEta = EgEta;
            matchedEgPhi = EgPhi;
            matchedEgEt  = EgEt;
            indx_closestEg = indx;

            eta_region = 0; // intialize
            if( fabs(EgEta) <= 0.8 ) eta_region =1;
            if( fabs(EgEta) <= 1.4 && fabs(EgEta) > 0.8 ) eta_region =2;
            if( fabs(EgEta) <= 1.7 && fabs(EgEta) > 1.4 ) eta_region =3;
            if( fabs(EgEta) <= 2.1 && fabs(EgEta) > 1.7 ) eta_region =4;
            if( fabs(EgEta) <= 2.7 && fabs(EgEta) > 2.1 ) eta_region =5;
            if( fabs(EgEta) <= 3.0 && fabs(EgEta) > 2.7 ) eta_region =6;

            if( fabs(EgEta) > 3.0 ) continue;
            //if( eta_region != 6 ) continue;

            Bool_t flag123 = false;
            Bool_t flag124 = false;
            Bool_t flag134 = false;
            Bool_t flag234 = false;

            Float_t recoPV = 0.;

            Float_t zp1 = -99.;
            Float_t zp2 = -99.;
            Float_t zp3 = -99.;
            Float_t zp4 = -99.;

            ntnEg2++;
            ntEgEt.push_back(EgEt);
            ntEgEta.push_back(EgEta);
            ntEgPhi.push_back(EgPhi);

            // set regin of interest
            SetROI(eta_region);

            // initialize pixel hit variables
            first_layer_hits.clear();
            second_layer_hits.clear();
            third_layer_hits.clear();
            fourth_layer_hits.clear();

            first_layer_hits_Ele_or_Pos.clear();
            second_layer_hits_Ele_or_Pos.clear();
            third_layer_hits_Ele_or_Pos.clear();
            fourth_layer_hits_Ele_or_Pos.clear();
            hitted_layers.clear();

            layers[0] = 1; // beam spot
            layers[1] = 0; layers[2] = 0; layers[3] = 0; layers[4] = 0;
            r = 0;

            StorePixelHit(eta_region); // save pixel hits in Region of Interest for the given eta region`

            // check which pixel has hits
            for( int i=1; i < 5; i++){ 
                if( layers[i] != 0 ){ 
                    hitted_layers.push_back(i); 
                }
            }

            trigger_bit_width_ = 0x0;
            trigger_bit_width_iso_ = 0x0;

            // set pixtrk signal boundary
            for(int nth_eg_pix_deta = 0; nth_eg_pix_deta < 9; nth_eg_pix_deta++){
                if( nth_eg_pix_deta != 0) continue;
                if(eta_region == 1) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta], EM_PiX_deta_width_[nth_eg_pix_deta], PiX_PiX_dphi_width_[nth_eg_pix_deta], PiX_PiX_deta_width_[nth_eg_pix_deta]);
                else if(eta_region == 2) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta+1], EM_PiX_deta_width_[nth_eg_pix_deta+1], PiX_PiX_dphi_width_[nth_eg_pix_deta+1], PiX_PiX_deta_width_[nth_eg_pix_deta+1]);
                else if(eta_region == 3) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta+1], EM_PiX_deta_width_[nth_eg_pix_deta+1], PiX_PiX_dphi_width_[nth_eg_pix_deta+1], PiX_PiX_deta_width_[nth_eg_pix_deta+1]);
                else if(eta_region == 4) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta+1], EM_PiX_deta_width_[nth_eg_pix_deta+1], PiX_PiX_dphi_width_[nth_eg_pix_deta+1], PiX_PiX_deta_width_[nth_eg_pix_deta+1]);
                else if(eta_region == 5) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta+1], EM_PiX_deta_width_[nth_eg_pix_deta+1], PiX_PiX_dphi_width_[nth_eg_pix_deta+1], PiX_PiX_deta_width_[nth_eg_pix_deta+1]);
                else if(eta_region == 6) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta+1], EM_PiX_deta_width_[nth_eg_pix_deta+1], PiX_PiX_dphi_width_[nth_eg_pix_deta+1], PiX_PiX_deta_width_[nth_eg_pix_deta+2]);
                else SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta], EM_PiX_deta_width_[nth_eg_pix_deta], PiX_PiX_dphi_width_[nth_eg_pix_deta], PiX_PiX_deta_width_[nth_eg_pix_deta]);

                // PixTRK algorithm 
                PixTrkPassed = false;
                withoutEM_count_Ele = 0, withEM_count_Ele = 0;

                fourth_layer_missing = 0;
                third_layer_missing = 0;
                second_layer_missing = 0;
                first_layer_missing = 0;

                // loop over every 3 out of 4 pixel combination 
                for( std::vector<int>::iterator first_hit = hitted_layers.begin(); first_hit != hitted_layers.end(); first_hit++){
                    for ( std::vector<int>::iterator second_hit = first_hit+1; second_hit != hitted_layers.end(); second_hit++){
                        for ( std::vector<int>::iterator third_hit = second_hit+1; third_hit != hitted_layers.end(); third_hit++){


                            // loop over every pixel hits in the given pixel combination
                            for( int k=0; k < layers[*first_hit]; k++){
                                for( int i=0; i < layers[*second_hit]; i++){
                                    _pass_Ele = 0, _pass_Pos = 0;
                                    L012_pass_Ele = 0, L012_pass_Pos = 0;
                                    L013_pass_Ele = 0, L013_pass_Pos = 0;
                                    L023_pass_Ele = 0, L023_pass_Pos = 0;

                                    if( *first_hit == 1 && *second_hit == 2 )
                                        TriggeringWith_1st2ndPixel(k,i);

                                    if( *first_hit == 1 && *second_hit == 3 )
                                        TriggeringWith_1st3rdPixel(k,i);

                                    if( *first_hit == 2 && *second_hit == 3 )
                                        TriggeringWith_2nd3rdPixel(k,i);

                                    // skip only if both _pass_Ele and _pass_Pos are 0 i.e., both electron and positron signal window are not satisfied
                                    if( !_pass_Ele && !_pass_Pos ) continue;

                                    for( int j=0; j < layers[*third_hit]; j++){
                                        all_cut_pass_Ele = 0, all_cut_pass_Pos = 0;
                                        withoutEM_pass_Ele = 0, withEM_pass_Ele = 0;

                                        L012_pass_Ele = 0, L012_pass_Pos = 0;
                                        L013_pass_Ele = 0, L013_pass_Pos = 0;
                                        L014_pass_Ele = 0, L014_pass_Pos = 0;
                                        L023_pass_Ele = 0, L023_pass_Pos = 0;
                                        L024_pass_Ele = 0, L024_pass_Pos = 0;
                                        L034_pass_Ele = 0, L034_pass_Pos = 0;
                                        L123_pass_Ele = 0, L123_pass_Pos = 0;
                                        L124_pass_Ele = 0, L124_pass_Pos = 0;
                                        L134_pass_Ele = 0, L134_pass_Pos = 0;
                                        L234_pass_Ele = 0, L234_pass_Pos = 0;

                                        L12_EM_Ele = 0, L12_EM_Pos = 0;
                                        L13_EM_Ele = 0, L13_EM_Pos = 0;
                                        L14_EM_Ele = 0, L14_EM_Pos = 0;
                                        L23_EM_Ele = 0, L23_EM_Pos = 0;
                                        L24_EM_Ele = 0, L24_EM_Pos = 0;
                                        L34_EM_Ele = 0, L34_EM_Pos = 0;

                                        dPhi = StandaloneDPhi( *first_hit, *second_hit, *third_hit, k, i, j );
                                        dEta = StandaloneDEta( *first_hit, *second_hit, *third_hit, k, i, j );

                                        if( *first_hit == 1 && *second_hit == 2 && *third_hit == 3 ){ // for efficiency counting  !!caution of the position of this codition
                                            // This is for the case that the first hit is in the first pixel layer and the second hit is in the second pixel layer and the third hit is in the third layer. 
                                            TriggeringWithout_4thPixel(k, i, j);

                                            if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
                                                    (second_layer_hits_Ele_or_Pos[i] == 1 || second_layer_hits_Ele_or_Pos[i] ==3) &&
                                                    (third_layer_hits_Ele_or_Pos[j] == 1 || third_layer_hits_Ele_or_Pos[j] ==3)){
                                                if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele && L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
                                                    all_cut_pass_Ele = 1; 
                                                if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele)
                                                    withoutEM_pass_Ele = 1;
                                                if(L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
                                                    withEM_pass_Ele = 1;
                                            }

                                            if( L012_pass_Pos && L013_pass_Pos && L023_pass_Pos && L123_pass_Pos && L12_EM_Pos && L13_EM_Pos && L23_EM_Pos &&
                                                    (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
                                                    (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
                                                    (third_layer_hits_Ele_or_Pos[j] == 2 || third_layer_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

                                            if( all_cut_pass_Ele == 1){
                                                pix_comb_ = pix_comb_ | (bit1 << 1);
                                                nPix123_segments++;
                                            }
                                            
                                            // Save vertex coordinate and check whether pass or not
                                            if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ) {
                                                flag123 = true;
                                                Float_t R1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
                                                Float_t R2 = sqrt(pow(second_layer_hits[i].X(),2)+pow(second_layer_hits[i].Y(),2));
                                                Float_t R3 = sqrt(pow(third_layer_hits[j].X(),2)+pow(third_layer_hits[j].Y(),2));
                                                Float_t Z1 = first_layer_hits[k].Z();
                                                Float_t Z2 = second_layer_hits[i].Z();
                                                Float_t Z3 = third_layer_hits[j].Z();
                                                if( eta_region <= 2 || eta_region >= 4 ) zp1 = (R3*Z1 - R1*Z3)/(R3-R1);
                                                if( eta_region == 3 ) zp1 = (R3*Z2 - R2*Z3)/(R3-R2);
                                            }

                                            if(skip){
                                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit loop
                                                    k = layers[*first_hit];
                                                    i = layers[*second_hit];
                                                    j = layers[*third_hit]; 
                                                    fourth_layer_missing = 1;
                                                }
                                            }
                                        }
                                        if( *first_hit == 1 && *second_hit == 2 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                                            // This is for the case that the first hit is in the first pixel layer and the second is in the second layer and the third hit is in the fourth layer.
                                            TriggeringWithout_3rdPixel(k, i, j);

                                            if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
                                                    (second_layer_hits_Ele_or_Pos[i] == 1 || second_layer_hits_Ele_or_Pos[i] ==3) &&
                                                    (fourth_layer_hits_Ele_or_Pos[j] == 1 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                                if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele && L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) all_cut_pass_Ele = 1;
                                                if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele) withoutEM_pass_Ele = 1;
                                                if(L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) withEM_pass_Ele = 1;
                                            } 

                                            if( L012_pass_Pos && L014_pass_Pos && L024_pass_Pos && L124_pass_Pos && L12_EM_Pos && L14_EM_Pos && L24_EM_Pos &&
                                                    (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
                                                    (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
                                                    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

                                            if( all_cut_pass_Ele == 1){
                                                pix_comb_ = pix_comb_ | (bit1 << 2);
                                                nPix124_segments++;
                                            }
                                            
                                            // Save vertex coordinate and check whether pass or not
                                            if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ) {
                                                flag124 = true;
                                                Float_t R1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
                                                Float_t R2 = sqrt(pow(second_layer_hits[i].X(),2)+pow(second_layer_hits[i].Y(),2));
                                                Float_t R4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
                                                Float_t Z1 = first_layer_hits[k].Z();
                                                Float_t Z2 = second_layer_hits[i].Z();
                                                Float_t Z4 = fourth_layer_hits[j].Z();
                                                if( eta_region == 1 || eta_region >= 4 ) zp2 = (R4*Z1 - R1*Z4)/(R4-R1);
                                                if( eta_region == 2 || eta_region == 3 ) zp2 = (R4*Z2 - R2*Z4)/(R4-R2);
                                            }

                                            if(skip){
                                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
                                                    k = layers[*first_hit];
                                                    i = layers[*second_hit];
                                                    j = layers[*third_hit]; 
                                                    third_layer_missing = 1;
                                                }
                                            }
                                        }
                                        if( *first_hit == 1 && *second_hit == 3 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                                            TriggeringWithout_2ndPixel(k, i, j);

                                            if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
                                                    (third_layer_hits_Ele_or_Pos[i] == 1 || third_layer_hits_Ele_or_Pos[i] ==3) &&
                                                    (fourth_layer_hits_Ele_or_Pos[j] == 1 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                                if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele && L13_EM_Ele && L14_EM_Ele && L34_EM_Ele)all_cut_pass_Ele = 1; 
                                                if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele) withoutEM_pass_Ele = 1;
                                                if(L13_EM_Ele && L14_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
                                            }

                                            if( L013_pass_Pos && L014_pass_Pos && L034_pass_Pos && L134_pass_Pos && L13_EM_Pos && L14_EM_Pos && L34_EM_Pos &&
                                                    (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
                                                    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
                                                    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

                                            if( all_cut_pass_Ele == 1){
                                                pix_comb_ = pix_comb_ | (bit1 << 3);
                                                nPix134_segments++;
                                            }
                                
                                            // Save vertex coordinate and check whether pass or not
                                            if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ) {
                                                flag134 = true;
                                                Float_t R1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
                                                Float_t R3 = sqrt(pow(third_layer_hits[i].X(),2)+pow(third_layer_hits[i].Y(),2));
                                                Float_t R4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
                                                Float_t Z1 = first_layer_hits[k].Z();
                                                Float_t Z3 = third_layer_hits[i].Z();
                                                Float_t Z4 = fourth_layer_hits[j].Z();
                                                if( eta_region == 1 || eta_region >= 3 ) zp3 = (R4*Z1 - R1*Z4)/(R4-R1);
                                                if( eta_region == 2 ) zp3 = (R4*Z3 - R3*Z4)/(R4-R3);
                                            }

                                            if(skip){
                                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
                                                    k = layers[*first_hit];
                                                    i = layers[*second_hit];
                                                    j = layers[*third_hit]; 
                                                    second_layer_missing = 1; 
                                                }
                                            }
                                        }
                                        if( *first_hit == 2 && *second_hit == 3 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                                            TriggeringWithout_1stPixel(k, i, j);

                                            if( (second_layer_hits_Ele_or_Pos[k] == 1 || second_layer_hits_Ele_or_Pos[k] ==3) &&
                                                    (third_layer_hits_Ele_or_Pos[i] == 1 || third_layer_hits_Ele_or_Pos[i] ==3) &&
                                                    (fourth_layer_hits_Ele_or_Pos[j] == 1 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                                if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele && L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) all_cut_pass_Ele = 1;
                                                if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele) withoutEM_pass_Ele = 1;
                                                if(L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
                                            } 

                                            if( L023_pass_Pos && L024_pass_Pos && L034_pass_Pos && L234_pass_Pos && L23_EM_Pos && L24_EM_Pos && L34_EM_Pos &&
                                                    (second_layer_hits_Ele_or_Pos[k] == 2 || second_layer_hits_Ele_or_Pos[k] ==3) &&
                                                    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
                                                    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

                                            if( all_cut_pass_Ele == 1){
                                                pix_comb_ = pix_comb_ | (bit1 << 4);
                                                nPix234_segments++;
                                            } 
                                
                                            // Save vertex coordinate and check whether pass or not
                                            if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ) {
                                                flag234 = true;
                                                Float_t R2 = sqrt(pow(second_layer_hits[k].X(),2)+pow(second_layer_hits[k].Y(),2));
                                                Float_t R3 = sqrt(pow(third_layer_hits[i].X(),2)+pow(third_layer_hits[i].Y(),2));
                                                Float_t R4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
                                                Float_t Z2 = second_layer_hits[k].Z();
                                                Float_t Z3 = third_layer_hits[i].Z();
                                                Float_t Z4 = fourth_layer_hits[j].Z();
                                                if( eta_region == 1 || eta_region >= 3 ) zp4 = (R4*Z2 - R2*Z4)/(R4-R2);
                                                if( eta_region == 2 ) zp4 = (R4*Z3 - R3*Z4)/(R4-R3);
                                            }

                                            if(skip){
                                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
                                                    k = layers[*first_hit];
                                                    i = layers[*second_hit];
                                                    j = layers[*third_hit]; 
                                                    first_layer_missing = 1;
                                                } 
                                            }
                                        }

                                        //if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ) { PixTrkPassed = true;}
                                        if( all_cut_pass_Ele == 1 ) { PixTrkPassed = true;}
                                        if( withoutEM_pass_Ele == 1 ) withoutEM_count_Ele = 1;
                                        if( withEM_pass_Ele == 1 ) withEM_count_Ele = 1;
                                    } // loop for third layer hits
                                } // loop for second layer hits       
                            } // loop for first layer hits
                        }          
                    }
                }

                if( PixTrkPassed ) { 
                    trigger_bit_width_ = trigger_bit_width_| (bit1 << nth_eg_pix_deta);
                }
            
                // %%%%%%%%%%%%%%%%%%%%% Track isolation algorithm %%%%%%%%%%%%%%%%%%%%%

                TrkIsoPassed = false;
                if( !PixTrkPassed ) {
                    continue; // Skip when L1 Egamma doesn't pass PixTRK algorithm
                }

                // Select which combination will be used to calculate reconstructed vertex
                if( eta_region == 1 || eta_region >= 4 ) {
                    if( flag124 || flag134 ) { 
                        if( zp2 != -99.) recoPV = zp2;
                        if( zp3 != -99.) recoPV = zp3;
                    }
                    if( !flag124 && !flag134 && flag123 ) recoPV = zp1;
                    if( !flag124 && !flag134 && !flag123 && flag234 ) recoPV = zp4;
                }
                if( eta_region == 2 ) {
                    if( flag234 || flag134 ) { 
                        if( zp4 != -99. ) recoPV = zp4;
                        if( zp3 != -99. ) recoPV = zp3;
                    }
                    if( !flag234 && !flag134 && flag123 ) recoPV = zp1;
                    if( !flag234 && !flag134 && !flag123 && flag124 ) recoPV = zp2;
                }
                if( eta_region == 3 ) {
                    if( flag124 || flag234 ) {
                        if( zp4 != -99. ) recoPV = zp4;
                        if( zp2 != -99. ) recoPV = zp2;
                    }
                    if( !flag124 && !flag234 && flag123 ) recoPV = zp1;
                    if( !flag124 && !flag234 && !flag123 && flag134 ) recoPV = zp3;
                }

                if( eta_region == 1 ) aa1->Fill(recoPV-nt_genZ);
                if( eta_region == 2 ) aa2->Fill(recoPV-nt_genZ);
                if( eta_region == 3 ){
                    //cout << "flag123 : " << flag123 << ", flag124 : " << flag124 << ", flag 134 : " << flag134 << ", flag234 : " << flag234 << endl;
                    aa3->Fill(recoPV-nt_genZ);
                }
                if( eta_region == 4 ) aa4->Fill(recoPV-nt_genZ);
                if( eta_region == 5 ) aa5->Fill(recoPV-nt_genZ);
                if( eta_region == 6 ) aa6->Fill(recoPV-nt_genZ);
     
                // initialize pixel hit variables to use in track isolation algorithm
                first_hits.clear();
                second_hits.clear();
                third_hits.clear();
                fourth_hits.clear();

                // Store pixel clusters in vectors
                StorePixelHitForIso(eta_region, recoPV); 

                L123.clear();
                L124.clear();
                L134.clear();
                L234.clear();

                // 1st step of the track isolation - filter out with kinematic parameters
                IsoWith_1st2nd3rd(eta_region, recoPV);
                IsoWith_1st2nd4th(eta_region, recoPV);
                IsoWith_1st3rd4th(eta_region, recoPV);
                IsoWith_2nd3rd4th(eta_region, recoPV);

                // Erase duplication in each combination
                if( L123.size() >= 2 ) 
                {
                    sort(L123.begin(), L123.end(), track::comp3);
                    L123.erase(unique(L123.begin(), L123.end(), track::uni3),L123.end());
                    sort(L123.begin(), L123.end(), track::comp2);
                    L123.erase(unique(L123.begin(), L123.end(), track::uni2),L123.end());
                    sort(L123.begin(), L123.end(), track::comp1);
                    L123.erase(unique(L123.begin(), L123.end(), track::uni1),L123.end());
                }
                if( L124.size() >= 2 ) 
                {
                    sort(L124.begin(), L124.end(), track::comp3);
                    L124.erase(unique(L124.begin(), L124.end(), track::uni3),L124.end());
                    sort(L124.begin(), L124.end(), track::comp2);
                    L124.erase(unique(L124.begin(), L124.end(), track::uni2),L124.end());
                    sort(L124.begin(), L124.end(), track::comp1);
                    L124.erase(unique(L124.begin(), L124.end(), track::uni1),L124.end());
                }
                if( L134.size() >= 2 ) 
                {
                    sort(L134.begin(), L134.end(), track::comp3);
                    L134.erase(unique(L134.begin(), L134.end(), track::uni3),L134.end());
                    sort(L134.begin(), L134.end(), track::comp2);
                    L134.erase(unique(L134.begin(), L134.end(), track::uni2),L134.end());
                    sort(L134.begin(), L134.end(), track::comp1);
                    L134.erase(unique(L134.begin(), L134.end(), track::uni1),L134.end());
                }
                if( L234.size() >= 2 ) 
                {
                    sort(L234.begin(), L234.end(), track::comp3);
                    L234.erase(unique(L234.begin(), L234.end(), track::uni3),L234.end());
                    sort(L234.begin(), L234.end(), track::comp2);
                    L234.erase(unique(L234.begin(), L234.end(), track::uni2),L234.end());
                    sort(L234.begin(), L234.end(), track::comp1);
                    L234.erase(unique(L234.begin(), L234.end(), track::uni1),L234.end());
                }
                // Make vector to contain all pixel clusters combinations from different layer combinations and erase duplication
                vector<track> all;
                all.clear();

                for(vector<track>::iterator a1 = L123.begin(); a1 != L123.end(); ++a1) 
                    all.push_back(track((*a1).pos_x3, (*a1).pos_x2, (*a1).pos_x1, (*a1).pos_y3, (*a1).pos_y2, (*a1).pos_y1, (*a1).pos_z3, (*a1).pos_z2, (*a1).pos_z1, 1 ));
                for(vector<track>::iterator a2 = L124.begin(); a2 != L124.end(); ++a2) 
                    all.push_back(track((*a2).pos_x3, (*a2).pos_x2, (*a2).pos_x1, (*a2).pos_y3, (*a2).pos_y2, (*a2).pos_y1, (*a2).pos_z3, (*a2).pos_z2, (*a2).pos_z1, 2 ));
                for(vector<track>::iterator a3 = L134.begin(); a3 != L134.end(); ++a3) 
                    all.push_back(track((*a3).pos_x3, (*a3).pos_x2, (*a3).pos_x1, (*a3).pos_y3, (*a3).pos_y2, (*a3).pos_y1, (*a3).pos_z3, (*a3).pos_z2, (*a3).pos_z1, 3 ));
                for(vector<track>::iterator a4 = L234.begin(); a4 != L234.end(); ++a4) 
                    all.push_back(track((*a4).pos_x3, (*a4).pos_x2, (*a4).pos_x1, (*a4).pos_y3, (*a4).pos_y2, (*a4).pos_y1, (*a4).pos_z3, (*a4).pos_z2, (*a4).pos_z1, 4 ));

                // Erase vector compare in the same queue
                sort(all.begin(), all.end(), track::comp3);
                all.erase(unique(all.begin(), all.end(), track::uni3),all.end());
                sort(all.begin(), all.end(), track::comp2);
                all.erase(unique(all.begin(), all.end(), track::uni2),all.end());
                sort(all.begin(), all.end(), track::comp1);
                all.erase(unique(all.begin(), all.end(), track::uni1),all.end());

                // Erase vector compare in different queues
                all.erase(unique(all.begin(), all.end(), track::uni12),all.end());
                all.erase(unique(all.begin(), all.end(), track::uni13),all.end());
                all.erase(unique(all.begin(), all.end(), track::uni23),all.end());
                all.erase(unique(all.begin(), all.end(), track::uni21),all.end());
                all.erase(unique(all.begin(), all.end(), track::uni31),all.end());
                all.erase(unique(all.begin(), all.end(), track::uni32),all.end());

                // For distribution, we consider L1 e/gamma larger than 20 GeV
                //if( EgEt < 20. ) continue;

                ntnEg2Iso++;
                ntEgEtIso.push_back(EgEt);
                ntEgEtaIso.push_back(EgEta);
                ntEgPhiIso.push_back(EgPhi);

                vector<leading> pT_vector;
                Int_t all_size = all.size();

                //cout << "Number of tracks: " << all_size << endl;
                //cout << endl;

                if( all.size() == 0 || all.size() >= 16) { 
                    TrkIsoPassed = true;
                    
                    if( EgEt > 20. ) {
                        h1->Fill(0);
                        h2->Fill(0);
                        h3->Fill(0);
                        h4->Fill(0);
                        h5->Fill(0);
                        h6->Fill(0);
                    }
                    NumOfTrks.push_back(0);
                    IsoValue.push_back(0);
                    track_pT1.push_back(-999.);
                    track_pT2.push_back(-999.);
                    track_pT3.push_back(-999.);
                    track_pT4.push_back(-999.);
                    track_pT5.push_back(-999.);
                    track_pT6.push_back(-999.);
                    track_pT7.push_back(-999.);
                    track_pT8.push_back(-999.);
                    track_pT9.push_back(-999.);
                    track_pT10.push_back(-999.);
                    track_pT11.push_back(-999.);
                    track_pT12.push_back(-999.);
                    track_pT13.push_back(-999.);
                    track_pT14.push_back(-999.);
                    track_pT15.push_back(-999.);
                    track_eta1.push_back(-999.);
                    track_eta2.push_back(-999.);
                    track_eta3.push_back(-999.);
                } // When the last vector size <= 1

                else {

                    pT_vector.clear();
                    for(Int_t cur = 0; cur < all_size; cur++)
                    {
                        if( all[cur].index == 1 )
                        {
                            TVector3 pixel1; pixel1.SetXYZ( all[cur].pos_x1, all[cur].pos_y1, all[cur].pos_z1 - recoPV );
                            TVector3 pixel2; pixel2.SetXYZ( all[cur].pos_x3 - all[cur].pos_x1, all[cur].pos_y3 - all[cur].pos_y1, all[cur].pos_z3 - all[cur].pos_z1 );

                            Float_t phi1 = pixel1.Phi(); Float_t phi2 = pixel2.Phi();
                            Float_t dPhi = deltaPhi(phi1,phi2);
                            Float_t recopT = pT_fit(eta_region, all[cur].index, dPhi);
                            Float_t trackEta = pixel2.Eta();
                            if( recopT > 0. ) pT_vector.push_back(leading(recopT, trackEta));
                        }
                        if( all[cur].index == 2 )
                        {
                            TVector3 pixel1; pixel1.SetXYZ( all[cur].pos_x2, all[cur].pos_y2, all[cur].pos_z2 - recoPV );
                            TVector3 pixel2; pixel2.SetXYZ( all[cur].pos_x3 - all[cur].pos_x2, all[cur].pos_y3 - all[cur].pos_y2, all[cur].pos_z3 - all[cur].pos_z2 );

                            Float_t phi1 = pixel1.Phi(); Float_t phi2 = pixel2.Phi();
                            Float_t dPhi = deltaPhi(phi1,phi2);
                            Float_t recopT = pT_fit(eta_region, all[cur].index, dPhi);
                            Float_t trackEta = pixel2.Eta();
                            if( recopT > 0. ) pT_vector.push_back(leading(recopT, trackEta));
                        }
                        if( all[cur].index == 3 )
                        {
                            TVector3 pixel1; pixel1.SetXYZ( all[cur].pos_x2, all[cur].pos_y2, all[cur].pos_z2 - recoPV );
                            TVector3 pixel2; pixel2.SetXYZ( all[cur].pos_x3 - all[cur].pos_x2, all[cur].pos_y3 - all[cur].pos_y2, all[cur].pos_z3 - all[cur].pos_z2 );

                            Float_t phi1 = pixel1.Phi(); Float_t phi2 = pixel2.Phi();
                            Float_t dPhi = deltaPhi(phi1,phi2);
                            Float_t recopT = pT_fit(eta_region, all[cur].index, dPhi);
                            Float_t trackEta = pixel2.Eta();
                            if( recopT > 0. ) pT_vector.push_back(leading(recopT, trackEta));
                        }
                        if( all[cur].index == 4 )
                        {
                            TVector3 pixel1; pixel1.SetXYZ( all[cur].pos_x1, all[cur].pos_y1, all[cur].pos_z1 - recoPV );
                            TVector3 pixel2; pixel2.SetXYZ( all[cur].pos_x3 - all[cur].pos_x2, all[cur].pos_y3 - all[cur].pos_y2, all[cur].pos_z3 - all[cur].pos_z2 );

                            Float_t phi1 = pixel1.Phi(); Float_t phi2 = pixel2.Phi();
                            Float_t dPhi = deltaPhi(phi1,phi2);
                            Float_t recopT = pT_fit(eta_region, all[cur].index, dPhi);
                            Float_t trackEta = pixel2.Eta();
                            if( recopT > 0. ) pT_vector.push_back(leading(recopT, trackEta));
                        }
                    }

                    sort(pT_vector.begin(), pT_vector.end(), leading::ptsort );
                    
                    Float_t denomi = 0.; Float_t nomi = 0.;
                    Int_t vec_size = pT_vector.size();
                    for(Int_t k = 0; k < vec_size; k++) {
                        denomi += pT_vector[k].pt;
                        if( k < 4 ) h_tracketa_all->Fill(pT_vector[k].eta);
                    }
                    for(Int_t k = 1; k < vec_size; k++) nomi += pT_vector[k].pt;

                    Float_t ratio = nomi/denomi;
                    NumOfTrks.push_back(vec_size);

                    switch( vec_size ) {
                        case 1:
                            IsoValue.push_back(0);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(-999.);
                            track_pT3.push_back(-999.);
                            track_pT4.push_back(-999.);
                            track_pT5.push_back(-999.);
                            track_pT6.push_back(-999.);
                            track_pT7.push_back(-999.);
                            track_pT8.push_back(-999.);
                            track_pT9.push_back(-999.);
                            track_pT10.push_back(-999.);
                            track_pT11.push_back(-999.);
                            track_pT12.push_back(-999.);
                            track_pT13.push_back(-999.);
                            track_pT14.push_back(-999.);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            track_eta1.push_back(pT_vector[0].eta);
                            track_eta2.push_back(-999.);
                            track_eta3.push_back(-999.);
                            break;
                        case 2:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(-999.);
                            track_pT4.push_back(-999.);
                            track_pT5.push_back(-999.);
                            track_pT6.push_back(-999.);
                            track_pT7.push_back(-999.);
                            track_pT8.push_back(-999.);
                            track_pT9.push_back(-999.);
                            track_pT10.push_back(-999.);
                            track_pT11.push_back(-999.);
                            track_pT12.push_back(-999.);
                            track_pT13.push_back(-999.);
                            track_pT14.push_back(-999.);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            track_eta1.push_back(pT_vector[0].eta);
                            track_eta2.push_back(pT_vector[1].eta);
                            track_eta3.push_back(-999.);
                            break;
                        case 3:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(pT_vector[2].pt);
                            track_pT4.push_back(-999.);
                            track_pT5.push_back(-999.);
                            track_pT6.push_back(-999.);
                            track_pT7.push_back(-999.);
                            track_pT8.push_back(-999.);
                            track_pT9.push_back(-999.);
                            track_pT10.push_back(-999.);
                            track_pT11.push_back(-999.);
                            track_pT12.push_back(-999.);
                            track_pT13.push_back(-999.);
                            track_pT14.push_back(-999.);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            h_tracketa_3->Fill(pT_vector[2].eta);
                            track_eta1.push_back(pT_vector[0].eta);
                            track_eta2.push_back(pT_vector[1].eta);
                            track_eta3.push_back(pT_vector[2].eta);
                            break;
                        case 4:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(pT_vector[2].pt);
                            track_pT4.push_back(pT_vector[3].pt);
                            track_pT5.push_back(-999.);
                            track_pT6.push_back(-999.);
                            track_pT7.push_back(-999.);
                            track_pT8.push_back(-999.);
                            track_pT9.push_back(-999.);
                            track_pT10.push_back(-999.);
                            track_pT11.push_back(-999.);
                            track_pT12.push_back(-999.);
                            track_pT13.push_back(-999.);
                            track_pT14.push_back(-999.);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            h_tracketa_3->Fill(pT_vector[2].eta);
                            h_tracketa_4->Fill(pT_vector[3].eta);
                            track_eta1.push_back(pT_vector[0].eta);
                            track_eta2.push_back(pT_vector[1].eta);
                            track_eta3.push_back(pT_vector[2].eta);
                            break;
                        case 5:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(pT_vector[2].pt);
                            track_pT4.push_back(pT_vector[3].pt);
                            track_pT5.push_back(pT_vector[4].pt);
                            track_pT6.push_back(-999.);
                            track_pT7.push_back(-999.);
                            track_pT8.push_back(-999.);
                            track_pT9.push_back(-999.);
                            track_pT10.push_back(-999.);
                            track_pT11.push_back(-999.);
                            track_pT12.push_back(-999.);
                            track_pT13.push_back(-999.);
                            track_pT14.push_back(-999.);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            h_tracketa_3->Fill(pT_vector[2].eta);
                            h_tracketa_4->Fill(pT_vector[3].eta);
                            track_eta1.push_back(pT_vector[0].eta);
                            track_eta2.push_back(pT_vector[1].eta);
                            track_eta3.push_back(pT_vector[2].eta);
                            break;
                        case 6:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(pT_vector[2].pt);
                            track_pT4.push_back(pT_vector[3].pt);
                            track_pT5.push_back(pT_vector[4].pt);
                            track_pT6.push_back(pT_vector[5].pt);
                            track_pT7.push_back(-999.);
                            track_pT8.push_back(-999.);
                            track_pT9.push_back(-999.);
                            track_pT10.push_back(-999.);
                            track_pT11.push_back(-999.);
                            track_pT12.push_back(-999.);
                            track_pT13.push_back(-999.);
                            track_pT14.push_back(-999.);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            h_tracketa_3->Fill(pT_vector[2].eta);
                            h_tracketa_4->Fill(pT_vector[3].eta);
                            track_eta1.push_back(pT_vector[0].eta);
                            track_eta2.push_back(pT_vector[1].eta);
                            track_eta3.push_back(pT_vector[2].eta);
                            break;
                        case 7:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(pT_vector[2].pt);
                            track_pT4.push_back(pT_vector[3].pt);
                            track_pT5.push_back(pT_vector[4].pt);
                            track_pT6.push_back(pT_vector[5].pt);
                            track_pT7.push_back(pT_vector[6].pt);
                            track_pT8.push_back(-999.);
                            track_pT9.push_back(-999.);
                            track_pT10.push_back(-999.);
                            track_pT11.push_back(-999.);
                            track_pT12.push_back(-999.);
                            track_pT13.push_back(-999.);
                            track_pT14.push_back(-999.);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            h_tracketa_3->Fill(pT_vector[2].eta);
                            h_tracketa_4->Fill(pT_vector[3].eta);
                            track_eta1.push_back(pT_vector[0].eta);
                            track_eta2.push_back(pT_vector[1].eta);
                            track_eta3.push_back(pT_vector[2].eta);
                            break;
                        case 8:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(pT_vector[2].pt);
                            track_pT4.push_back(pT_vector[3].pt);
                            track_pT5.push_back(pT_vector[4].pt);
                            track_pT6.push_back(pT_vector[5].pt);
                            track_pT7.push_back(pT_vector[6].pt);
                            track_pT8.push_back(pT_vector[7].pt);
                            track_pT9.push_back(-999.);
                            track_pT10.push_back(-999.);
                            track_pT11.push_back(-999.);
                            track_pT12.push_back(-999.);
                            track_pT13.push_back(-999.);
                            track_pT14.push_back(-999.);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            h_tracketa_3->Fill(pT_vector[2].eta);
                            h_tracketa_4->Fill(pT_vector[3].eta);
                            track_eta1.push_back(pT_vector[0].eta);
                            track_eta2.push_back(pT_vector[1].eta);
                            track_eta3.push_back(pT_vector[2].eta);
                            break;
                        case 9:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(pT_vector[2].pt);
                            track_pT4.push_back(pT_vector[3].pt);
                            track_pT5.push_back(pT_vector[4].pt);
                            track_pT6.push_back(pT_vector[5].pt);
                            track_pT7.push_back(pT_vector[6].pt);
                            track_pT8.push_back(pT_vector[7].pt);
                            track_pT9.push_back(pT_vector[8].pt);
                            track_pT10.push_back(-999.);
                            track_pT11.push_back(-999.);
                            track_pT12.push_back(-999.);
                            track_pT13.push_back(-999.);
                            track_pT14.push_back(-999.);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            h_tracketa_3->Fill(pT_vector[2].eta);
                            h_tracketa_4->Fill(pT_vector[3].eta);
                            track_eta1.push_back(pT_vector[0].eta);
                            track_eta2.push_back(pT_vector[1].eta);
                            track_eta3.push_back(pT_vector[2].eta);
                            break;
                        case 10:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(pT_vector[2].pt);
                            track_pT4.push_back(pT_vector[3].pt);
                            track_pT5.push_back(pT_vector[4].pt);
                            track_pT6.push_back(pT_vector[5].pt);
                            track_pT7.push_back(pT_vector[6].pt);
                            track_pT8.push_back(pT_vector[7].pt);
                            track_pT9.push_back(pT_vector[8].pt);
                            track_pT10.push_back(pT_vector[9].pt);
                            track_pT11.push_back(-999.);
                            track_pT12.push_back(-999.);
                            track_pT13.push_back(-999.);
                            track_pT14.push_back(-999.);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            h_tracketa_3->Fill(pT_vector[2].eta);
                            h_tracketa_4->Fill(pT_vector[3].eta);
                            break;
                        case 11:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(pT_vector[2].pt);
                            track_pT4.push_back(pT_vector[3].pt);
                            track_pT5.push_back(pT_vector[4].pt);
                            track_pT6.push_back(pT_vector[5].pt);
                            track_pT7.push_back(pT_vector[6].pt);
                            track_pT8.push_back(pT_vector[7].pt);
                            track_pT9.push_back(pT_vector[8].pt);
                            track_pT10.push_back(pT_vector[9].pt);
                            track_pT11.push_back(pT_vector[10].pt);
                            track_pT12.push_back(-999.);
                            track_pT13.push_back(-999.);
                            track_pT14.push_back(-999.);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            h_tracketa_3->Fill(pT_vector[2].eta);
                            h_tracketa_4->Fill(pT_vector[3].eta);
                            break;
                        case 12:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(pT_vector[2].pt);
                            track_pT4.push_back(pT_vector[3].pt);
                            track_pT5.push_back(pT_vector[4].pt);
                            track_pT6.push_back(pT_vector[5].pt);
                            track_pT7.push_back(pT_vector[6].pt);
                            track_pT8.push_back(pT_vector[7].pt);
                            track_pT9.push_back(pT_vector[8].pt);
                            track_pT10.push_back(pT_vector[9].pt);
                            track_pT11.push_back(pT_vector[10].pt);
                            track_pT12.push_back(pT_vector[11].pt);
                            track_pT13.push_back(-999.);
                            track_pT14.push_back(-999.);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            h_tracketa_3->Fill(pT_vector[2].eta);
                            h_tracketa_4->Fill(pT_vector[3].eta);
                            break;
                        case 13:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(pT_vector[2].pt);
                            track_pT4.push_back(pT_vector[3].pt);
                            track_pT5.push_back(pT_vector[4].pt);
                            track_pT6.push_back(pT_vector[5].pt);
                            track_pT7.push_back(pT_vector[6].pt);
                            track_pT8.push_back(pT_vector[7].pt);
                            track_pT9.push_back(pT_vector[8].pt);
                            track_pT10.push_back(pT_vector[9].pt);
                            track_pT11.push_back(pT_vector[10].pt);
                            track_pT12.push_back(pT_vector[11].pt);
                            track_pT13.push_back(pT_vector[12].pt);
                            track_pT14.push_back(-999.);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            h_tracketa_3->Fill(pT_vector[2].eta);
                            h_tracketa_4->Fill(pT_vector[3].eta);
                            break;
                        case 14:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(pT_vector[2].pt);
                            track_pT4.push_back(pT_vector[3].pt);
                            track_pT5.push_back(pT_vector[4].pt);
                            track_pT6.push_back(pT_vector[5].pt);
                            track_pT7.push_back(pT_vector[6].pt);
                            track_pT8.push_back(pT_vector[7].pt);
                            track_pT9.push_back(pT_vector[8].pt);
                            track_pT10.push_back(pT_vector[9].pt);
                            track_pT11.push_back(pT_vector[10].pt);
                            track_pT12.push_back(pT_vector[11].pt);
                            track_pT13.push_back(pT_vector[12].pt);
                            track_pT14.push_back(pT_vector[13].pt);
                            track_pT15.push_back(-999.);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            h_tracketa_3->Fill(pT_vector[2].eta);
                            h_tracketa_4->Fill(pT_vector[3].eta);
                            break;
                        case 15:
                            IsoValue.push_back(ratio);
                            track_pT1.push_back(pT_vector[0].pt);
                            track_pT2.push_back(pT_vector[1].pt);
                            track_pT3.push_back(pT_vector[2].pt);
                            track_pT4.push_back(pT_vector[3].pt);
                            track_pT5.push_back(pT_vector[4].pt);
                            track_pT6.push_back(pT_vector[5].pt);
                            track_pT7.push_back(pT_vector[6].pt);
                            track_pT8.push_back(pT_vector[7].pt);
                            track_pT9.push_back(pT_vector[8].pt);
                            track_pT10.push_back(pT_vector[9].pt);
                            track_pT11.push_back(pT_vector[10].pt);
                            track_pT12.push_back(pT_vector[11].pt);
                            track_pT13.push_back(pT_vector[12].pt);
                            track_pT14.push_back(pT_vector[13].pt);
                            track_pT15.push_back(pT_vector[14].pt);
                            h_tracketa_1->Fill(pT_vector[0].eta);
                            h_tracketa_2->Fill(pT_vector[1].eta);
                            h_tracketa_3->Fill(pT_vector[2].eta);
                            h_tracketa_4->Fill(pT_vector[3].eta);
                            break;
                    } // switch function
                    
                    if( eta_region == 1 ) nh1->Fill(vec_size);
                    if( eta_region == 2 ) nh2->Fill(vec_size);
                    if( eta_region == 3 ) nh3->Fill(vec_size);
                    if( eta_region == 4 ) nh4->Fill(vec_size);
                    if( eta_region == 5 ) nh5->Fill(vec_size);
                    if( eta_region == 6 ) nh6->Fill(vec_size);

                    if( EgEt > 20 ) {
                        if( eta_region == 1 ) h1->Fill(ratio);
                        if( eta_region == 2 ) h2->Fill(ratio);
                        if( eta_region == 3 ) h3->Fill(ratio);
                        if( eta_region == 4 ) h4->Fill(ratio);
                        if( eta_region == 5 ) h5->Fill(ratio);
                        if( eta_region == 6 ) h6->Fill(ratio);
                    }

                    if( eta_region == 1 || eta_region == 4 ) {
                        if( ratio < 0.18 ) TrkIsoPassed = true;
                    }
                    if( eta_region == 2 ) {
                        if( ratio < 0.14 ) TrkIsoPassed = true;
                    }
                    if( eta_region == 3 ) {
                        if( ratio < 0.08 ) TrkIsoPassed = true;
                    }
                    if( eta_region == 5 ) {
                        if( ratio < 0.42 ) TrkIsoPassed = true;
                    }
                    if( eta_region == 6 ) {
                        if( ratio < 0.23 ) TrkIsoPassed = true;
                    }


                }

                if( PixTrkPassed && TrkIsoPassed ) {
                    trigger_bit_width_iso_ = trigger_bit_width_iso_| (bit2 << nth_eg_pix_deta);
                }


            } // nth_eg_pix_deta loop
        
        } // end of egamma (if function)

        trigger_bit_width.push_back(trigger_bit_width_);
        trigger_bit_width_iso.push_back(trigger_bit_width_iso_);
        pix_comb.push_back(pix_comb_);

        pixtrk_tree->Fill();
    } // event loop

    outfile->Write();
}
