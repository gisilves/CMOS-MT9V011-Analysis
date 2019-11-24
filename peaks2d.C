#include "TSpectrum2.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TH2.h"
#include "TF2.h"
#include "TMath.h"
#include "TROOT.h"
#include "TMarker.h"

TSpectrum2 *s;
TH2F *h2 = 0;

Int_t npeaks = 30;
Double_t fpeaks2(Double_t *x, Double_t *par) {
   Double_t result = 0.1;
   for (Int_t p=0;p<npeaks;p++) {
      Double_t norm   = par[5*p+0];
      Double_t mean1  = par[5*p+1];
      Double_t sigma1 = par[5*p+2];
      Double_t mean2  = par[5*p+3];
      Double_t sigma2 = par[5*p+4];
      result += norm*TMath::Gaus(x[0],mean1,sigma1)*TMath::Gaus(x[1],mean2,sigma2);
   }
   return result;
}


double findPeak2() {
   printf("Generating histogram with %d peaks\n",npeaks);
   Int_t nbinsx = 128;
   Int_t nbinsy = 128;
   Double_t xmin   = 0;
   Double_t xmax   = (Double_t)nbinsx;
   Double_t ymin   = 0;
   Double_t ymax   = (Double_t)nbinsy;
   Double_t dx = (xmax-xmin)/nbinsx;
   Double_t dy = (ymax-ymin)/nbinsy;
   Int_t i, j;

   Double_t** source = new Double_t*[nbinsx+10];
   for (i=0;i<nbinsx+10;i++)
     source[i]=new Double_t[nbinsy+10];
   Double_t** dest = new Double_t*[nbinsx+10];
   for (i=0;i<nbinsx+10;i++)
     dest[i]=new Double_t[nbinsy+10];
   
   // Double_t dx = 10;                                                              
   // Double_t dy = 10; 
   
   delete h2;
   h2 = new TH2F("h2","test",nbinsx,xmin,xmax,nbinsy,ymin,ymax);
   h2->SetStats(0);

   /*
   //generate n peaks at random
   gRandom->SetSeed(100);
   Double_t par[3000];
   Int_t p;
   for (p=0;p<npeaks;p++) {
      par[5*p+0] = gRandom->Uniform(0.2,1);
      par[5*p+1] = gRandom->Uniform(xmin,xmax);
      par[5*p+2] = gRandom->Uniform(dx,5*dx);
      par[5*p+3] = gRandom->Uniform(ymin,ymax);
      par[5*p+4] = gRandom->Uniform(dy,5*dy);
   }
   
   TF2 *f2 = new TF2("f2",fpeaks2,xmin,xmax,ymin,ymax,5*npeaks);
   f2->SetNpx(100);
   f2->SetNpy(100);
   f2->SetParameters(par);
   TCanvas *c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c1");
   if (!c1) c1 = new TCanvas("c1","c1",10,10,1000,700);
   h2->FillRandom("f2",500000);
   */
   for (i = 0; i < nbinsx+10; i++){
     for (j = 0; j < nbinsy+10; j++){
       if(i >= nbinsx || j >= nbinsy){
	 source[i][j] = 0;
       } else {
	 source[i][j] = h2->GetBinContent(i + 1,j + 1);
       }
     }
   }

   //now the real stuff: Finding the peaks
   Int_t nfound = s->SearchHighRes(source, dest, nbinsx, nbinsy, 2, 4, kFALSE, 25, kTRUE, 5);
   
   //searching good and ghost peaks (approximation)
   Int_t pf,ngood = 0;
   Double_t *xpeaks = s->GetPositionX();
   Double_t *ypeaks = s->GetPositionY();

   for (p=0;p<npeaks;p++) {
      for (pf=0;pf<nfound;pf++) {
         Double_t diffx = TMath::Abs(xpeaks[pf] - par[5*p+1]);
         Double_t diffy = TMath::Abs(ypeaks[pf] - par[5*p+3]);
         if (diffx < 2*dx && diffy < 2*dy) ngood++;
      }
   }
   if (ngood > nfound) ngood = nfound;

   //Search ghost peaks (approximation)
   Int_t nghost = 0;
   for (pf=0;pf<nfound;pf++) {
      Int_t nf=0;
      for (p=0;p<npeaks;p++) {
         Double_t diffx = TMath::Abs(xpeaks[pf] - par[5*p+1]);
         Double_t diffy = TMath::Abs(ypeaks[pf] - par[5*p+3]);
         if (diffx < 3*dx && diffy < 3*dy) nf++;
      }
      if (nf == 0){
	nghost++;
      }else{
	cout<< "Peak value [" << pf << "] " <<(Int_t)source[(Int_t)(xpeaks[pf])][(Int_t)(ypeaks[pf])] << endl;

	double cluster = 0;
	int max_widthX = 10;
	int max_widthY = 10;
	int threshold = 100;

	for(int i=-max_widthX/2; i<=max_widthX/2; i++){
	  for(int j=-max_widthY/2; j<=max_widthY/2; j++){
	    int signal = (Int_t)source[(Int_t)(xpeaks[0])+i][(Int_t)(ypeaks[0])+j];
	    if (signal < threshold) continue;
	    cluster+=signal;
	  }
	}
	cout<< "Cluster value [" << pf << "] "  << cluster << endl;
      }
   }
   
   h2->Draw("colz");
   s->Print();
   c1->Update();

   printf("Gener=%d, Found=%d, Good=%d, Ghost=%d\n",npeaks,nfound,ngood,nghost);
   
   //if (!gROOT->IsBatch()) {
   // printf("\nDouble click in the bottom right corner of the pad to continue\n");
   // c1->WaitPrimitive();
   //}
   double perc = 1.0 - double((npeaks-ngood))/(double)npeaks;
   return perc;
}

void peaks2d(Int_t maxpeaks=5) {
   s = new TSpectrum2(2*maxpeaks);
   double mean = 0;
   
   int frames = 1;

   for (int i=0; i<frames; ++i) {
     npeaks = (Int_t)gRandom->Uniform(1,maxpeaks);
     double perc = findPeak2();
     mean += perc;
   }
   
   cout << endl;
   cout << "percentage found " << mean/frames*100  << "%"<< endl;
}


