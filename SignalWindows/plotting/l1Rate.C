#define l1Rate_cxx
#include "l1Rate.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void l1Rate::Loop()
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   
   float genN = 488500.00; 
   int eta_r = 1;

   int nbins = 151; float x1 = 9. ; float x2 = 160. ;
   TH1F* hEG = new TH1F("hEG",";E_{T} threshold (GeV); Rate (kHz)",nbins,x1,x2);
   TH1F* hTKEG = new TH1F("hTKEG",";E_{T} threshold (GeV); Rate (kHz)",nbins,x1,x2);
   TH1F* hPXEG = new TH1F("hPXEG",";E_{T} threshold (GeV); Rate (kHz)",nbins,x1,x2);
   TH1F* hPXIsoEG = new TH1F("hPXIsoEG",";E_{T} threshold (GeV); Rate (kHz)",nbins,x1,x2);

   TH1F* r_PXEG = new TH1F("L1PXEG_r_factor",";E_{T} threshold (GeV); Rate reduction factor",nbins,x1,x2);
   TH1F* r_PXIsoEG = new TH1F("L1PXIsoEG_r_factor",";E_{T} threshold (GeV); Rate reduction factor",nbins,x1,x2);

   int width_bit[27] = {0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80, 0x100, 0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x8000, 0x10000, 0x20000, 0x40000, 0x80000,
                        0x100000, 0x200000, 0x400000, 0x800000, 0x1000000, 0x2000000, 0x4000000};
   
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      if(!ntnEg2) continue;

      float max_eget = -999.;
      for(int i=0; i < ntnEg2; i++){

         //if(eta_r == 1 && ntEgEt->at(i) > max_eget && (fabs(ntEgEta->at(i)) < 1.5) ) max_eget = ntEgEt->at(i);
         //if(eta_r == 1 && ntEgEt->at(i) > max_eget && (fabs(ntEgEta->at(i)) > 1.5 && fabs(ntEgEta->at(i)) < 2.5) ) max_eget = ntEgEt->at(i);
         //if(eta_r == 1 && ntEgEt->at(i) > max_eget && (fabs(ntEgEta->at(i)) > 2.5 && fabs(ntEgEta->at(i)) < 3.0) ) max_eget = ntEgEt->at(i);
         //if(eta_r == 1 && ntEgEt->at(i) > max_eget && (fabs(ntEgEta->at(i)) < 2.5) ) max_eget = ntEgEt->at(i);
         if(eta_r == 1 && ntEgEt->at(i) > max_eget && (fabs(ntEgEta->at(i)) < 3.0) ) max_eget = ntEgEt->at(i);
           
      }// eg loop

      float max_pxeget = -999.;
      for(int i=0; i < ntnEg2; i++){
         //if(eta_r == 1 && ntEgEt->at(i) > max_pxeget && ntCl_match->at(i) && (fabs(ntEgEta->at(i)) < 1.5)) max_pxeget = ntEgEt->at(i);
         //if(eta_r == 1 && ntEgEt->at(i) > max_pxeget && ntCl_match->at(i) && (fabs(ntEgEta->at(i)) > 1.5 && fabs(ntEgEta->at(i)) < 2.5) ) max_pxeget = ntEgEt->at(i);
         //if(eta_r == 1 && ntEgEt->at(i) > max_pxeget && ntCl_match->at(i) && (fabs(ntEgEta->at(i)) > 2.5 && fabs(ntEgEta->at(i)) < 3.0) ) max_pxeget = ntEgEt->at(i);
         //if(eta_r == 1 && ntEgEt->at(i) > max_pxeget && ntCl_match->at(i) && (fabs(ntEgEta->at(i)) < 2.5) ) max_pxeget = ntEgEt->at(i);
         if(eta_r == 1 && ntEgEt->at(i) > max_pxeget && ntCl_match->at(i) && (fabs(ntEgEta->at(i)) < 3.0) ) max_pxeget = ntEgEt->at(i);
           
      }// pxeg loop
      
      float max_pxeget_iso = -999.;
      for(int i=0; i < ntnEg2; i++){
         //if(eta_r == 1 && ntEgEt->at(i) > max_pxeget_iso && ntCl_iso_match->at(i) && (fabs(ntEgEta->at(i)) < 1.5)) max_pxeget_iso = ntEgEt->at(i);
         //if(eta_r == 1 && ntEgEt->at(i) > max_pxeget_iso && ntCl_iso_match->at(i) && (fabs(ntEgEta->at(i)) > 1.5 && fabs(ntEgEta->at(i)) < 2.5) ) max_pxeget_iso = ntEgEt->at(i);
         //if(eta_r == 1 && ntEgEt->at(i) > max_pxeget_iso && ntCl_iso_match->at(i) && (fabs(ntEgEta->at(i)) > 2.5 && fabs(ntEgEta->at(i)) < 3.0) ) max_pxeget_iso = ntEgEt->at(i);
         //if(eta_r == 1 && ntEgEt->at(i) > max_pxeget_iso && ntCl_iso_match->at(i) && (fabs(ntEgEta->at(i)) < 2.5) ) max_pxeget_iso = ntEgEt->at(i);
         if(eta_r == 1 && ntEgEt->at(i) > max_pxeget_iso && ntCl_iso_match->at(i) && (fabs(ntEgEta->at(i)) < 3.0) ) max_pxeget_iso = ntEgEt->at(i);
           
      }// pxegIso loop
   
      if(max_eget > 0 ) hEG->Fill(max_eget);
      
      if(max_pxeget > 0 ){ 
        bool check_l1eg = false;
        for(int i=0; i < ntnEg2; i++){
           if(ntEgEt->at(i) == max_eget){
             check_l1eg = true;    
             continue;
           }
        }
        
        if(check_l1eg) hPXEG->Fill(max_pxeget);
      }
      
      if(max_pxeget_iso > 0 ){ 
        bool check_l1eg_iso = false;
        for(int i=0; i < ntnEg2; i++){
           if(ntEgEt->at(i) == max_eget){
             check_l1eg_iso = true;    
             continue;
           }
        }
        
        if(check_l1eg_iso) hPXIsoEG->Fill(max_pxeget_iso);
      }

   
   } // event loop

   for (int i=0; i<= nbins; i++) {
       hEG->SetBinContent(i, hEG -> Integral(i, nbins+1) );
       hEG->SetBinError(i, sqrt( hEG -> GetBinContent(i) ) ); 

       hPXEG->SetBinContent(i, hPXEG -> Integral(i, nbins+1) );
       hPXEG->SetBinError(i, sqrt( hPXEG -> GetBinContent(i) ) ); 
       
       hPXIsoEG->SetBinContent(i, hPXIsoEG -> Integral(i, nbins+1) );
       hPXIsoEG->SetBinError(i, sqrt( hPXIsoEG -> GetBinContent(i) ) ); 

   }

   cout << "0 bin: " << hEG->GetBinContent(0) << endl;
   cout << "1 bin: " << hEG->GetBinContent(1) << endl;
   cout << "11 bin: " << hEG->GetBinContent(16) * 30000./genN<< endl;
   cout << "11px bin: " << hPXEG->GetBinContent(16) * 30000./genN<< endl;
   cout << "11px iso bin: " << hPXIsoEG->GetBinContent(16) * 30000./genN<< endl;
   cout << "nbin bin: " << hEG->GetBinContent(nbins) << endl;
   cout << "nbin+1 bin: " << hEG->GetBinContent(nbins+1) << endl;

   hEG -> Scale(30000./genN);
   hPXEG -> Scale(30000./genN);
   hPXIsoEG -> Scale(30000./genN);

   TCanvas *c1 = new TCanvas("c1","c1",800,700);
   gStyle->SetOptStat(0);
   gStyle->SetLineWidth(1); // axis width, default is 1
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.12);
   c1->SetRightMargin(0.03);
   c1->SetLeftMargin(0.2);
   c1->SetLogy();
   c1->SetGrid();
   c1->SetTicky(1);
   c1->SetTickx(1);

   hEG->SetTitle("");
   hEG->GetXaxis()->CenterTitle(true);
   hEG->GetXaxis()->SetTitleOffset(1.);
   hEG->GetXaxis()->SetTitleSize(0.055);
   hEG->GetXaxis()->SetNdivisions(515);
   hEG->GetYaxis()->SetNdivisions(546);
   hEG->GetXaxis()->SetLabelSize(0.05);
   hEG->GetYaxis()->SetLabelSize(0.05);
   hEG->GetXaxis()->SetRangeUser(9.5,100.5);
   hEG->GetYaxis()->SetRangeUser(1,35000);
   hEG->GetYaxis()->SetTitleOffset(1.2);
   hEG->GetYaxis()->SetTitleSize(0.055);
   hEG->SetMinimum(1.);

   hEG->SetMarkerColor(2);
   hEG->SetLineColor(2);
   hEG->SetLineWidth(2);
   hEG->SetMarkerStyle(20);
   hEG->SetMarkerSize(1.);
   hEG->Draw("HIST E");

   hPXEG->SetMarkerColor(4);
   hPXEG->SetLineColor(4);
   hPXEG->SetLineWidth(2);
   hPXEG->SetMarkerStyle(20);
   hPXEG->SetMarkerSize(1.);
   hPXEG->Draw("HIST same E");
   
   hPXIsoEG->SetMarkerColor(8);
   hPXIsoEG->SetLineColor(8);
   hPXIsoEG->SetLineWidth(2);
   hPXIsoEG->SetMarkerStyle(20);
   hPXIsoEG->SetMarkerSize(1.);
   hPXIsoEG->Draw("HIST same E");

   TLegend *Lgd = new TLegend(0.35, 0.55, 0.65, 0.9);
   Lgd->SetFillColor(0);
   Lgd->SetTextFont(42);
   Lgd->SetTextSize(0.03);
   Lgd->SetBorderSize(0);
   Lgd->SetFillStyle(0);
   Lgd->AddEntry(hEG,"Phase-2 L1 EG","lp");
   Lgd->AddEntry(hPXEG,"Phase-2 L1 EG + Pixel","lp");
   Lgd->AddEntry(hPXIsoEG,"Phase-2 L1 EG + Pixel + Isolation","lp");
   Lgd->Draw();

   float r_ = c1->GetRightMargin();
   float t_ = c1->GetTopMargin();

   TLatex t(1-r_,1-t_+0.2*t_,"Delphes Simulation, <PU>=200");
   t.SetNDC();
   t.SetTextFont(42);
   t.SetTextAlign(31);
   t.SetTextSize(0.6*t_);
   t.Draw();

   TString eta_str;

   //if(eta_r == 1) eta_str = "#it{#lbar#bf{#eta}#lbar < 1.5}";
   //if(eta_r == 1) eta_str = "#it{1.5 < #lbar#bf{#eta}#lbar < 2.5}";
   //if(eta_r == 1) eta_str = "#it{2.5 < #lbar#bf{#eta}#lbar < 3.0}";
   //if(eta_r == 1) eta_str = "#it{#lbar#bf{#eta}#lbar < 2.5}";
   if(eta_r == 1) eta_str = "#it{#lbar#bf{#eta}#lbar < 3.0}";

   //TLatex eta_range(70,4500,eta_str);
   TLatex eta_range(66.5,4500,eta_str);
   eta_range.SetTextSize(0.05);
   eta_range.Draw();

   TString l1rate = "l1rate";

   //if(eta_r == 0) l1rate = l1rate + "_all.png";
   //if(eta_r == 1) l1rate = l1rate + "_r1.png";
   if(eta_r == 1) l1rate = l1rate + "_all.png";
   //if(eta_r == 2) l1rate = l1rate + "_r2.png";
   //if(eta_r == 3) l1rate = l1rate + "_r3.png";
   //if(eta_r == 4) l1rate = l1rate + "_r4.png";

   c1->Print(l1rate);
   c1->SaveAs("l1rate_all.pdf");

   ////
   TCanvas *c2 = new TCanvas("c2","c2",700,700);
   gStyle->SetOptStat(0);
   gStyle->SetLineWidth(2); // axis width, default is 1
   c2->SetTopMargin(0.07);
   c2->SetBottomMargin(0.12);
   c2->SetRightMargin(0.03);
   c2->SetLeftMargin(0.2);
   c2->SetLogy();
   c2->SetGrid();
   c2->SetTicky(1);
   c2->SetTickx(1);

   c2->cd();


   r_PXEG->Add(hEG);
   r_PXEG->Divide(hPXEG);
   
   r_PXIsoEG->Add(hEG);
   r_PXIsoEG->Divide(hPXIsoEG);

   r_PXEG->SetTitle("");
   r_PXEG->GetXaxis()->SetTitleOffset(1.3);
   r_PXEG->GetXaxis()->SetTitleSize(0.045);
   r_PXEG->GetXaxis()->SetNdivisions(515);
   r_PXEG->GetYaxis()->SetNdivisions(546);
   r_PXEG->GetXaxis()->SetLabelSize(0.05);
   r_PXEG->GetYaxis()->SetLabelSize(0.05);
   r_PXEG->GetXaxis()->SetRangeUser(10.5,60.5);
   r_PXEG->GetYaxis()->SetRangeUser(1,60);
   r_PXEG->GetYaxis()->SetTitleOffset(1.5);
   r_PXEG->GetYaxis()->SetTitleSize(0.050);

   r_PXEG->SetMarkerColor(1);
   r_PXEG->SetLineColor(1);
   r_PXEG->SetLineWidth(2);
   r_PXEG->SetMarkerStyle(20);
   r_PXEG->SetMarkerSize(1.);

   r_PXIsoEG->SetMarkerColor(4);
   r_PXIsoEG->SetLineColor(4);
   r_PXIsoEG->SetLineWidth(2);
   r_PXIsoEG->SetMarkerStyle(20);
   r_PXIsoEG->SetMarkerSize(1.);

   r_PXEG->Draw("HIST e");
   r_PXIsoEG->Draw("HIST same e");

   cout << endl;
   cout << "11px bin: " << r_PXEG->GetBinContent(16) << " error: " << r_PXEG->GetBinError(11) << endl;
   cout << "11px iso bin: " << r_PXIsoEG->GetBinContent(16) << " error: " << r_PXIsoEG->GetBinError(11) << endl;

   //TLatex t1(12,65,"CMS Preliminary Simulation, Phase 2, <PU>=200");
   TLatex t1(1-r_,1-t_+0.2*t_,"Delphes Simulation, <PU>=200");
   t1.SetNDC();
   t1.SetTextFont(42);
   t1.SetTextAlign(31);
   t1.SetTextSize(0.6*t_);
   //t1.SetTextSize(0.035);
   t1.Draw();

   TLatex eta_range1(12,45,eta_str);
   eta_range1.SetTextSize(0.035);
   eta_range1.Draw();

   TLegend *Lgd1 = new TLegend(0.25, 0.15, 0.75, 0.3); 
   //TLegend *Lgd1 = new TLegend(0.25, 0.65, 0.85, 0.8);
   //TLegend *Lgd1 = new TLegend(0.25, 0.71, 0.85, 0.86);
   Lgd1->SetFillColor(0);
   Lgd1->SetTextFont(42);
   Lgd1->SetTextSize(0.03);
   Lgd1->SetBorderSize(1);
   Lgd1->AddEntry(r_PXEG,"L1 Pixel Detector","lp");
   Lgd1->AddEntry(r_PXIsoEG,"L1 Pixel Detector + isolation","lp");

   Lgd1->Draw();

   TString factor = "reduction_factor";

   //if(eta_r == 0) factor = factor + "_all.png";
   if(eta_r == 1) factor = factor + "_all.png";
   //if(eta_r == 1) factor = factor + "_r1.png";
   //if(eta_r == 2) factor = factor + "_r2.png";
   //if(eta_r == 3) factor = factor + "_r3.png";  
   //if(eta_r == 4) factor = factor + "_r4.png";

   c2->Print(factor);   
   c2->Print("reduction_factor_all.pdf");


}
