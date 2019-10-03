#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include <TFile.h>
#include "TFile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TString.h"
#include <fstream>
#include <string>
#include <cmath>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TRandom.h>
#include <TFitResultPtr.h>
#include <RConfig.h>
#include <TStopwatch.h>
#include <TLegend.h>
#include <TH2D.h>


double    P_max = 2;        //In GeV
int       nbins = 100;
int       nbinsp = 10000;
double    sigma_P_CM = 0.143;   //In GeV
double    kf = 0.221;
double    m=0.939;



void SF_p1_fast(){
    

    TStopwatch * Clock1 = new TStopwatch();
    
    // I open the file with the Wave Function Squared and save it int arrays
    ifstream myfile;
    myfile.open("all.txt");
    double k[100], phi2_nn[100], phi2_pn0[100], phi2_pn1[100];
    int i=0, j=0;
    
    for(i=0;i<100;i++){
        myfile >> k[i];
       k[i] = k[i] * (TMath::H()) * TMath::C() * 1e6/ (TMath::Qe()) / 2/ TMath::Pi(); // Conversion
        myfile >> phi2_nn[i];
        myfile >> phi2_pn0[i];
        myfile >> phi2_pn1[i];
    }
    myfile.close();
    
    //Open the histogram and start working for the wavefunction
    
    TH1D* histo_p_nn = new TH1D("histo_p_nn","histo_p_nn",nbins,k[0]/2,k[nbins-1]+k[0]/2);
    TH1D* histo_p_pn0 = new TH1D("histo_p_pn0","histo_p_pn0",nbins,k[0]/2,k[nbins-1]+k[0]/2);
    TH1D* histo_p_pn1 = new TH1D("histo_p_pn1","histo_p_pn1",nbins,k[0]/2,k[nbins-1]+k[0]/2);
    
    TCanvas* canvas1 = new TCanvas("phi2","phi2",0);
    canvas1->SetLogy();
    canvas1 -> cd();

    

    
    
    for(i=0;i<100;i++){
            histo_p_nn->Fill(k[i],phi2_nn[i]);
            histo_p_pn0->Fill(k[i],phi2_pn0[i]);
            histo_p_pn1->Fill(k[i],phi2_pn1[i]);
    }
    
   //  Normalization and drawing
    histo_p_nn->Scale(1/histo_p_nn->Integral(),"width");
    histo_p_nn->Draw("L HIST");
    
    histo_p_pn0->Scale(1/histo_p_pn0->Integral(),"width");
    histo_p_pn0->SetLineColor(2);
    histo_p_pn0->Draw("same L HIST");
    
    histo_p_pn1->Scale(1/histo_p_pn1->Integral(),"width");
    histo_p_pn1->SetLineColor(8);
    histo_p_pn1->Draw("same L HIST");
    
    TLegend* legend = new TLegend();
    legend->SetHeader("Wavefunction squared for relative momentum","C"); // option "C" allows to center the header
    legend->AddEntry(histo_p_nn,"NN");
    legend->AddEntry(histo_p_pn0,"PN Singlet");
    legend->AddEntry(histo_p_pn1,"PN Triplet");
    legend->Draw();
    
    //Here comes the hard work
    
    
    TH2D*   Histo_SF_nn = new TH2D("Spectral Function nn","Spectral Function nn",nbinsp,0,P_max,nbinsp,kf,k[99]);
    //TH2D*   Histo_SF_pn0 = new TH2D("Spectral Function pn0","Spectral Function pn0",nbinsp,kf,P_max,nbins,k[0]/2,k[nbins-1]+k[0]/2);
   // TH2D*   Histo_SF_pn1 = new TH2D("Spectral Function pn1","Spectral Function pn1",nbinsp,kf,P_max,nbins,k[0]/2,k[nbins-1]+k[0]/2);
    
    
    
    //Definitions of variable and functions
    double dp =(P_max)/(double)nbinsp; //* (TMath::H()) * TMath::C() * 1e6/ (TMath::Qe()) / 2/ TMath::Pi();
    double dPr = (k[99]-kf)/nbinsp, Pr=0;
    double p1=0;
    int    A = 12;     //C12
    double DB = 0.033;      //Binding energy
    double Weight_nn = 0, Weight_pn0 = 0, Weight_pn1 = 0;
    double Norm_SF_nn =0, Norm_SF_pn0=0, Norm_SF_pn1 =0;
    double Cost_min = -1, Cost_max = 1;
    
    int p1_i=0, i_stop=0 ;
    
    for(p1_i=0; p1_i<nbinsp; p1_i++){
    p1 = dp * ((double)(p1_i + 0.5));
        Cost_min = ((2*m - DB)*(2*m - DB) - m*m - p1*p1 - Pr*Pr)/(2*p1*Pr);
        if (Cost_min<-1){Cost_min=-1;}
        cout << p1_i<<endl;
        
           for(i=0; i<nbinsp; i++){
               Pr= dPr * ((double)(i + 0.5)) + kf;
                        //Condizioni di non validitá
                    if (Pr<kf||Pr>k[99]||Cost_min>Cost_max){
                        Weight_nn = 0;        //Not in domain
                        Weight_pn0 = 0;
                        Weight_pn1 = 0;
                            
                    }
                    else {
                        
  
                        
                        
                        Weight_nn = Pr/p1  * (histo_p_nn -> Interpolate(Pr) )/ sigma_P_CM *(TMath::Exp((4*p1*Pr*Cost_max-2*(p1*p1 + Pr*Pr))/(sigma_P_CM*sigma_P_CM))-TMath::Exp((4*p1*Pr*Cost_min-2*(p1*p1 + Pr*Pr))/(sigma_P_CM*sigma_P_CM)));
                            /*
                        Weight_pn0 = k[i]/p1  * (histo_p_pn0 -> GetBinContent(i) )/ sigma_P_CM * (TMath::Exp((-2)*(p1-k[i])*(p1-k[i])/(sigma_P_CM*sigma_P_CM)) - TMath::Exp((-2)*(p1*p1 - k[i]*k[i])/(sigma_P_CM*sigma_P_CM)));
                            
                        Weight_pn1 = k[i]/p1  * (histo_p_pn1 -> GetBinContent(i) )/ sigma_P_CM * (TMath::Exp((-2)*(p1-k[i])*(p1-k[i])/(sigma_P_CM*sigma_P_CM)) - TMath::Exp((-2)*(p1*p1 - k[i]*k[i])/(sigma_P_CM*sigma_P_CM)));
*/
                    }
                    Histo_SF_nn ->Fill (p1, Pr, Weight_nn);
              // Histo_SF_pn0 ->Fill (p1, k[i], Weight_pn0);
                   // Histo_SF_pn1 ->Fill (p1, k[i], Weight_pn1);
                       
                    Norm_SF_nn += Weight_nn;
                   // Norm_SF_pn0+= Weight_pn0;
                   // Norm_SF_pn1+=Weight_pn1;
            }
    }
    
    //Normalizations

    Histo_SF_nn -> Scale (1/Norm_SF_nn,"width");
   // Histo_SF_pn0 -> Scale (1/Norm_SF_pn0,"width");
   // Histo_SF_pn1 -> Scale (1/Norm_SF_pn1,"width");
    
    

    
    TCanvas* canvas3 = new TCanvas("Momentum","Momentum",0);
    canvas3 -> Divide (2,1);
    canvas3 -> cd (1);
        gPad ->SetLogy();
       gStyle ->SetOptStat(""); //For Drawing options

    
    //For the momentum
    TH1D* Histo_p1_nn = Histo_SF_nn->ProjectionX("Histo_p1_nn",0,-1);
    Histo_p1_nn->SetTitle("P_{1}");
    Histo_p1_nn->GetXaxis()->SetTitle("P_{1} (GeV/c)");
    Histo_p1_nn->GetYaxis()->SetTitle("S_{P_{1}} (a.u)");
    
 //   TH1D* Histo_p1_pn0 = Histo_SF_pn0->ProjectionX("Histo_p1_pn0",0,-1);
   // TH1D* Histo_p1_pn1 = Histo_SF_pn1->ProjectionX("Histo_p1_pn1",0,-1);
    
    //Histo_p1_pn0 ->SetLineColor(2);
    //Histo_p1_pn1 ->SetLineColor(8);
    
    Histo_p1_nn -> Draw(" HIST");
   // Histo_p1_pn0 -> Draw("same  HIST");
   // Histo_p1_pn1 -> Draw("same  HIST");
    
    /*
    TLegend* legendp1 = new TLegend(0.6,0.7,0.9,0.9);
    legendp1->SetTextSize(.02);
    legendp1->SetHeader("Distribution p1","C"); // option "C" allows to center the header
    legendp1->AddEntry(Histo_p1_nn,"NN");
    legendp1->AddEntry(Histo_p1_pn0,"PN Singlet");
    legendp1->AddEntry(Histo_p1_pn1,"PN Triplet");
    legendp1->Draw();*/
    
    canvas3 -> cd (2);
     gPad    ->SetLogy();
    
    TH1D* Histo_PR_nn = Histo_SF_nn->ProjectionY("Histo_PR_nn",0,-1);
    Histo_PR_nn->GetXaxis()->SetTitle("P_{R} (GeV/c)");
    Histo_PR_nn->GetYaxis()->SetTitle("S_{P_{R}} (a.u)");
    
  //  TH1D* Histo_PR_pn0 = Histo_SF_pn0->ProjectionY("Histo_PR_pn0",0,-1);
   // TH1D* Histo_PR_pn1 = Histo_SF_pn1->ProjectionY("Histo_PR_pn1",0,-1);
    
   // Histo_PR_pn0 ->SetLineColor(2);
   // Histo_PR_pn1 ->SetLineColor(8);
    
    Histo_PR_nn -> Draw("HIST");
   // Histo_PR_pn0 -> Draw("same  HIST");
   // Histo_PR_pn1 -> Draw("same  HIST");
    
    /*
    TLegend* legendPR = new TLegend(0.6,0.7,0.9,0.9);
    legendPR->SetTextSize(.02);
    legendPR->SetHeader("Distribution for PR","C"); // option "C" allows to center the header
    legendPR->AddEntry(Histo_PR_nn,"NN");
    legendPR->AddEntry(Histo_PR_pn0,"PN Singlet");
    legendPR->AddEntry(Histo_PR_pn1,"PN Triplet");
    legendPR->Draw();*/
    
       TFile f1("Contact_formalism.root", "UPDATE");
    
    //Saving into file
    Histo_p1_nn -> Write("P1_nn_1");
    Histo_PR_nn -> Write("PR_nn_1");
    
  //  Histo_p1_pn0 -> Write("P1_pn0_1");
  //  Histo_PR_pn0 -> Write("PR_pn0_1");
    
   // Histo_p1_pn1 -> Write("P1_pn1_1");
   // Histo_PR_pn1 -> Write("PR_pn1_1");
    
    
    
   f1.Close();
    
    cout << "eseguito in: " << Clock1->RealTime()<< " s" << endl;

    
}
