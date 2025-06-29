#include "TROOT.h" 
#include "TStyle.h" 
#include "TCanvas.h" 
#include "TMarker.h" 
#include "TH2F.h" 
#include "TLatex.h" 
#include "TPave.h" 
#include "TGraphAsymmErrors.h" 

#include <array>
#include <sstream>
#include <iomanip>

void blue_summary_plot_v2(TString pdffile="B_CMS_ISR_Run2_Electron_m50to64_5_DisRes_Obs_0.pdf",
                                        TString mass_window="[55, 64]", TString AXIS = "p_{T}^{ee} [GeV]", 
                                        TString period1 = "2016a", TString period2 = "2016b",
                                        const std::vector<Float_t>& valu  = {14.76, 14.56, 15.25},
                                        const std::vector<Float_t>& stat  = { 0.06,  0.09,  0.08},
                                        const std::vector<Float_t>& syst  = { 0.14,  0.14,  0.26},
                                        const std::vector<Float_t>& full  = { 0.15,  0.17,  0.27}
                                        ){
	 const Int_t N =  3; 
	 Float_t Indx[N] = { 0, 1, 2,};
	 Int_t   Colo[N] = {632, 1, 1,};
	 TString Name[N] = {mass_window, period1, period2};
	 Float_t Valu[N];
	 Float_t Stat[N];
	 Float_t Syst[N];
	 Float_t Full[N];

     for (int i = 0; i < N; i++){
        Valu[i] = valu[i];
        Stat[i] = stat[i];
        Syst[i] = syst[i];
        Full[i] = full[i];

     }

	 // Text Strings 
	 TString CanTit = "The combination You always wanted to do";
     TString UNCB = "                       (stat)        (syst)";

	 // Logo text
	 const Int_t    EmbDim = 2;
	 const TString  EmbNam[EmbDim] = {"BLUE", "2.4.0"};
	 const Int_t    EmbFnt[EmbDim] = {72, 42};
	 const Double_t EmbXva =  0.08;
	 const Double_t EmbYva =  0.88;
	 const Int_t    EmbCol = 600;
	 TLatex*  EmbTxt[EmbDim];
	 for(Int_t k = 0; k < EmbDim; k++){
		 EmbTxt[k] = new TLatex();
		 EmbTxt[k]->SetNDC();
		 EmbTxt[k]->SetTextColor(EmbCol);
		 EmbTxt[k]->SetTextFont(EmbFnt[k]);
	}

	 // Canvas parameters
	 Int_t Canx = 800, Cany = 650;

	 // Histogram x-axis is booked depending on Nxmi and Nxma with
	 // x=(Valu[0]-Fxmi*largest extension of unc bar to tge left,
	 //    Valu[0]+Fxma*largest extension of unc bar to the right) and y(-2, N+2)
	 // The offset of values (unc text) from combined valu can be steered with Fxva(Fxun)
	 //---------------------------------------------------------------------- 
	 // MOST CASES SHOULD BE ADJUSTABLE BY CHANGING THE FOLLOWING FOUR VALUES 
	 //---------------------------------------------------------------------- 
	 Float_t Fxmi = 2.0, Fxma = 3.0;
	 Float_t Fxva = 1.1, Fxun = 2.0;


	 //----------------------------------------------------------------------------
	 // No number or text needed below this line
	 //----------------------------------------------------------------------------

	 // Get min and max in x
	 Float_t xmi = 0, xma = 0, xfu = 0;
	 Float_t xText = 0, xValu = 0, xUnco = 0;
	 Float_t xmin = Valu[0], xmax = xmin;
	 for(Int_t j = 0; j < N; j++) {
		 if(Valu[j]-Full[j] < xmin)xmin = Valu[j]-Full[j];
		 if(Valu[j]+Full[j] > xmax)xmax = Valu[j]+Full[j];
	 }
	 xfu = (xmax - xmin);
	 xmi = Valu[0] - Fxmi*xfu;
	 xma = Valu[0] + Fxma*xfu;
	 xText = xmi + xfu/10.;
	 xValu = Valu[0] + Fxva*xfu;
	 xUnco = Valu[0] + Fxun*xfu;

	 // Set up graphs with stat and full uncertainties
	 TGraphAsymmErrors *grStat = new TGraphAsymmErrors(N, Valu, Indx, Stat, Stat, 0, 0);
	 TGraphAsymmErrors *grFull = new TGraphAsymmErrors(N, Valu, Indx, Full, Full, 0, 0);

	 // Now draw everything
	 gROOT->Reset(); 
	 gROOT->SetStyle("Plain");
	 gStyle->SetOptStat(0);
	 gStyle->SetPadTickX(1);
	 gStyle->SetPadTickY(1);
	 gStyle->SetEndErrorSize(6);

	 // Make canvas
	 TCanvas *CanDel;
	 CanDel = (TCanvas*) gROOT->FindObject("myCanvas");
	 if(CanDel) delete CanDel;
	 TCanvas *myCanvas = new TCanvas("myCanvas",CanTit, Canx, Cany);
	 myCanvas->SetTopMargin(0.025);
	 myCanvas->SetBottomMargin(0.18);
	 myCanvas->SetLeftMargin(0.04);
	 myCanvas->SetRightMargin(0.04);

	 // Book the TH2F plot
	 TH2F *H2FDel = (TH2F*) gROOT->FindObject("myPlot");
	 if(H2FDel)delete H2FDel;
	 TH2F *myPlot = new TH2F("myPlot", "", 20, xmi, xma, 20, -2.0, N + 1.0);
	 myPlot->GetXaxis()->CenterTitle();
	 myPlot->GetXaxis()->SetNdivisions(6);
	 myPlot->GetXaxis()->SetTitle(AXIS);
	 myPlot->GetXaxis()->SetTitleFont(42);
	 myPlot->GetXaxis()->SetTitleSize(0.05);
	 myPlot->GetXaxis()->SetTitleOffset(1.4);
	 myPlot->GetXaxis()->SetLabelFont(42);
	 myPlot->GetXaxis()->SetLabelSize(0.05);
	 myPlot->GetYaxis()->SetLabelColor(0);
	 myPlot->GetYaxis()->SetNdivisions(1);
	 myPlot->Draw();  

	 // Put the vertical line around the combined value
	 TPave *Comb = new TPave(Valu[0]-Full[0], Indx[0]-0.5, Valu[0]+Full[0], Indx[N-1]+0.5, 0, "br");
	 Comb->SetFillColor(17);
	 Comb->Draw();

	 // Draw results as points
	 grFull->SetLineColor(4);
	 grFull->SetLineWidth(3);
	 grFull->SetMarkerColor(2);
	 grFull->SetMarkerStyle(20);
	 grFull->SetMarkerSize(1.3);
	 grFull->Draw("P");
	 grStat->SetLineColor(2);
	 grStat->SetLineWidth(3);
	 grStat->SetMarkerColor(2);
	 grStat->SetMarkerStyle(20);
	 grStat->SetMarkerSize(1.3);
	 grStat->Draw("Psame");

	 // Store numerical values in a TString
	 // Print points and write values
	 TLatex *Text = new TLatex();
	 char cx[200];
	 for(Int_t j = 0; j < N; j++) {
		 //snprintf(cx,"%5.2f #pm %5.2f #pm %5.2f",Valu[j],Stat[j],Syst[j]);

         std::snprintf(cx, sizeof(cx),
              "%5.2f #pm %5.2f #pm %5.2f",
              Valu[j], Stat[j], Syst[j]);
		 TString Resu = cx;
		 Text->SetTextSize(0.015);
		 Text->SetTextColor(Colo[j]);
		 Text->SetTextFont(42);
		 Text->SetTextSize(0.035);
		 Text->DrawLatex(xText, Indx[j]-0.12, Name[j]);
		 Text->SetTextSize(0.045);
		 Text->DrawLatex(xValu, Indx[j]-0.12, Resu);
	}

	 // Print the header and the uncertainty breakdown
	 Text->SetTextSize(0.04);
	 Text->SetTextFont(42);
	 Text->SetTextColor(kBlack);
		 // The Logo 
		 Double_t EmbXdi = EmbNam[0].Sizeof()*0.019*696* 
			 gPad->GetWh()/(472*gPad->GetWw()); 
		 EmbTxt[0]->DrawLatex(       EmbXva,EmbYva,EmbNam[0].Data());
		 EmbTxt[1]->DrawLatex(EmbXva+EmbXdi,EmbYva,EmbNam[1].Data()); 

	 Text->SetTextSize(0.03);
	 Text->SetTextColor(kBlack);
     Text->DrawLatex(xValu,Indx[0]-1.0, UNCB);

	 // Write the file
	 myCanvas->Print(pdffile);
}
