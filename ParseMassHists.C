


float OptimizeMassWindow(float sigma){
  //ratio=sig/sqrt(bg)
  //we assume bg is flat, so ratio =K*sig_in_window/sqrt(window width);
  //this is a partial integral of a gaussian = an erf function
  //and has no analytical optimal solution, so we scan:
  static TF1 r("analyticalratio","TMath::Erf(x/[0])/sqrt(x)",0.05,10);
  r.SetParameter(0,sigma);
  float bestw=r.GetMaximumX(0.05,10);
  new TCanvas("junk","junk canvas",400,600);
  r.Draw();
  TLine j;
  j.DrawLine(bestw,0,bestw,100);
  return bestw;
}  


void ParseMassHists(){

  const float targetlumi=100e9;//100fb-1, expressed in ub-1
  const float alpha_d_ref=1e-8;//assumed alpha_dark in the code.
  const float alpha_qed=1.0/137.0;
  const float minimumevents=10;//minimnum number of A' events in window for perfect case.
  float targetsig=2;//significance units of sigma=sqrt(Bg)
  //oldstyle:  char * inputpattern="localMassSpectrumM*GeV.hist.root";
  char * inputpattern="otreePlotsM*GeV.hist.root";//caution that the _bg.hist.root are the bg files.
  char * bgname="uncutMassSpectrumBg.hist.root";
   
  //find all files that match the input string
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern);
  filelist->Print();
  printf("found: %s\n",((TFileInfo*)(filelist->GetList()->At(0)))->GetCurrentUrl()->GetFile());//Title());//Print();
  
  TFile *bgfile=TFile::Open(bgname,"READ");
  TH1F *hBgMass=(TH1F*)bgfile->Get("hRecoMassWide");
  TH1F *hBgMassYield=new TH1F(*hBgMass);
  hBgMassYield->SetTitle(Form("Smeared QED e+e- invMass yield @ L=%1.2E/fb;Mass (GeV); Yield per bin",targetlumi/1.0e9));
  hBgMassYield->Scale(targetlumi);
  TF1 *fBgModel=new TF1("fBgModel","[0]*(x^(-2))",2.0,100.0);
  fBgModel->SetParameter(0,1e3);
  fBgModel->SetParameter(1,-2);
  
  //hBgMass->Draw();
  //doesn't work hBgMass->Fit(fBgModel);

  int nMasses=filelist->GetNFiles();
  const int nBinsToAverage=2;//+/- 2 bins from center.
  float mass[nMasses];
  float masswidth[nMasses];//gaus parameter
  float masswindow[nMasses];
  float totBg[nMasses];//integrated smeared bg cross section in window
  float totSig[nMasses];//integrated smeared signal cross section in window
  float fullSigYield[nMasses];//signal yield before smearing
  float smearSigYield[nMasses];//signal yield in hist range (+/-5% from nominal mass, particles in detector range) with smearing
  float windowSigYield[nMasses];//signal yield with smearing and window cut.
  float smearOverFull[nMasses];
  float windowOverSmear[nMasses];
  float ratio[nMasses];
  float epsilon_reach[nMasses];
  float epsilon2_reach[nMasses];
  float epsilon2_max[nMasses];//maximum, zero-background reach where we hit the defined number of events.
  TH1F *hTrueMass[nMasses];
  TH1F *hMass[nMasses];
    TF1 *fGaus=new TF1("fGaus","gaus(0)",0,200);

  for (int i=0;i<nMasses;i++){
    TString filename=((TFileInfo*)(filelist->GetList()->At(i)))->GetCurrentUrl()->GetUrl();
    TFile *infile=TFile::Open(filename,"READ");

    hTrueMass[i]=(TH1F*)infile->Get("hAmass");
    fullSigYield[i]=hTrueMass[i]->Integral()*targetlumi;
    
    hMass[i]=(TH1F*)infile->Get("hRecoMass");
    smearSigYield[i]=hMass[i]->Integral()*targetlumi;
    smearOverFull[i]=smearSigYield[i]/fullSigYield[i];

    hMass[i]->Fit(fGaus);
    masswidth[i]=fGaus->GetParameter(2);
    masswindow[i]=OptimizeMassWindow(fGaus->GetParameter(2));
    mass[i]=0.5*(hMass[i]->GetXaxis()->GetXmax()+hMass[i]->GetXaxis()->GetXmin());
    //masswindow[i]=(hMass[i]->GetXaxis()->GetXmax()-hMass[i]->GetXaxis()->GetXmin());
    int bgbin=hBgMass->GetXaxis()->FindBin(mass[i]);
    int sigbinmin=hMass[i]->GetXaxis()->FindBin(mass[i]-masswindow[i]);
    int sigbinmax=hMass[i]->GetXaxis()->FindBin(mass[i]+masswindow[i]);
    totBg[i]=hBgMass->Integral(bgbin-nBinsToAverage,bgbin+nBinsToAverage)*(1.0/(1.0+2*nBinsToAverage))*(2*masswindow[i]);

    
    totSig[i]=hMass[i]->Integral(sigbinmin,sigbinmax);

    
    windowSigYield[i]=totSig[i]*targetlumi;
    windowOverSmear[i]=windowSigYield[i]/smearSigYield[i];
    ratio[i]=totSig[i]/sqrt(totBg[i])*sqrt(targetlumi);

    epsilo2n_reach[i]=targetsig*alpha_d_ref/alpha_qed*sqrt(totBg[i])/sqrt(targetlumi)/totSig[i];
    epsilon_reach[i]=sqrt(epsilon_reach[i]);
    //mineve=epsmax*alpha*sigtot*L
    //mineve=sigtot*L*alpha_d/alpha_0
    //alpha_d=mineve*alpha_0/(sigtot*L)
    //epsmax=alpha_d/alpha=mineve*alpha_0/(sigtot*L)/alpha
    float epsmax=minimumevents*alpha_d_ref/(totSig[i]*targetlumi*alpha_qed);
    epsilon2_max[i]=epsmax*epsmax;
    //hMass[i]->Draw("same");
    printf("M=%E\tRaw=%E\tSmear=%E\tWindow=%E\n",mass[i],fullSigYield[i],smearSigYield[i],windowSigYield[i]);

  }

  TCanvas *c; int nc=0;
  TVirtualPad *p;
  TGraph *g;
  TLegend *leg;
  TLine *line=new TLine();
  line->SetLineColor(kBlue);
  TLatex *tex=new TLatex();
  tex->SetTextAlign(12);

  if (1){ //plot the sample invM and fits to them
    c=new TCanvas(Form("c%d",nc),"canvas",1200,800);
    nc++;
    c->Divide(3,4);
    for (int i=0;i<nMasses;i++){
      p=c->cd(i+1);
      p->SetLogy();
      gStyle->SetTitleFontSize(0.1);
      hMass[i]->SetTitle(Form("Smeared e+e- invariant mass for Mtrue=%1.1fGeV;Mass (GeV);#sigma (ub/bin)",mass[i]));
      hMass[i]->Draw();      
      p->Modified();
      //hMass[i]->Draw();
      //((TPaveText*)(p->GetPrimitive("title")))->SetTextSize(0.1);
      //per https://root.cern.ch/root/roottalk/roottalk03/1099.html
      // y axes on the histogram aren't resolved until it draws.
      float xmin=hMass[i]->GetXaxis()->GetXmin();
      float xmax=hMass[i]->GetXaxis()->GetXmax();
      float ymin=0;hMass[i]->GetMinimum()/5;//p->GetY1();//hMass[i]->GetYaxis()->GetXmin();
      float ymax=hMass[i]->GetMaximum()*2;//p->GetY2();//hMass[i]->GetYaxis()->GetXmax();
      line->DrawLine(mass[i]-masswindow[i],ymin,mass[i]-masswindow[i],ymax);
      line->DrawLine(mass[i]+masswindow[i],ymin,mass[i]+masswindow[i],ymax);
      //tex->DrawLatex(mass[i]-masswindow[i],pow(10.0,0.5*(log10(ymax)+log10(ymin))),"  hello");
      // tex->DrawLatex(mass[i]-masswindow[i],pow(10.0,0.5*(log10(ymax)+log10(ymin))),Form("#Delta M=%E GeV",2*masswindow[i]));
      tex->DrawLatex(mass[i]-masswindow[i],ymax/1.0e2,Form("#Delta M=%1.2f GeV",2*masswindow[i]));
      //tex->DrawLatex(mass[i]-masswindow[i],pow(10.0,0.55*(log10(ymax)+log10(ymin))),Form("max=%E GeV",ymax));
    }
  }

  if (1){ //plot the SM yield, s/sqrt(B), and the resulting reach vs mass
    c=new TCanvas(Form("c%d",nc),"canvas",1200,400);
    nc++;
    c->Divide(2);
    p=c->cd(1);
    g=new TGraph(nMasses,mass,masswindow);
    g->SetTitle("Window width vs Mass;Mass (GeV);width (GeV)");
    g->Draw("A*");
     g=new TGraph(nMasses,mass,masswidth);
    g->SetTitle("Mass width");
    g->SetMarkerColor(kRed);
    g->Draw("*");
    leg=p->BuildLegend();
    ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel("window width");
   p=c->cd(2);
    p->SetLogy();
    g=new TGraph(nMasses,mass,smearOverFull);
    g->SetTitle("A' Yield Ratios");
    g->SetMarkerColor(kRed);
    g->GetHistogram()->SetMaximum(1.1);
    g->GetHistogram()->SetMinimum(1e-5);
    g->Draw("A*");
    p->SetTitle("test");
    g=new TGraph(nMasses,mass,windowOverSmear);
    g->SetTitle("window/smear");
    g->SetMarkerColor(kBlack);
    g->Draw("*");
    leg=p->BuildLegend();
    ((TLegendEntry*)leg->GetListOfPrimitives()->At(0))->SetLabel("smear/full");
  }
  if (1){ //plot the SM yield, s/sqrt(B), and the resulting reach vs mass
    c=new TCanvas(Form("c%d",nc),"canvas",1200,400);
    nc++;
    c->Divide(2,1);
    p=c->cd(1);
    p->SetLogy();
    //p->SetLogx();
    hBgMassYield->Draw();
    tex->DrawLatex(80,1.0e7,"Contains both pairs");
    tex->DrawLatex(80,0.15e7,"if 2 e- detected");
    for (int i=0;i<nMasses;i++){
      line->SetLineColor(kRed);
      //line->DrawLine(mass[i]-masswindow[i],1e2,mass[i]-masswindow[i],1e8);
      //line->DrawLine(mass[i]+masswindow[i],1e2,mass[i]+masswindow[i],1e8);
    }
    g=new TGraph(nMasses,mass,fullSigYield);
    g->SetTitle("A' truth yield");
    g->SetMarkerColor(kRed);
    g->Draw("*");
    g=new TGraph(nMasses,mass,smearSigYield);
    g->SetTitle("A' reco yield");
    g->SetMarkerColor(kBlue);
    g->Draw("*");
    g=new TGraph(nMasses,mass,windowSigYield);
    g->SetTitle("A' window yield");
    g->SetMarkerColor(kBlack);
    g->Draw("*");
    p=c->cd(2);//->cd(1);
    p->SetLogy();
    //p->SetLogx();
    //g=new TGraph(nMasses,mass,ratio);
    //g->SetTitle(Form("S/sqrt(B) vs Mass (L=%Eub, alpha_d=1e-8);mass(GeV);S/sqrt(B)",targetlumi));
    //g->Draw("A*");
    //p=c->cd(2)->cd(2);
    p->SetLogy();
    p->SetLogx();
    g=new TGraph(nMasses,mass,epsilon2_reach);
    g->SetTitle(Form("#epsilon^2 reach vs Mass (L=%1.2E/ub, S=%1.0f);mass(GeV);#epsilon^2",targetlumi,targetsig));
    //g->Draw("A*");
    g->GetXaxis()->SetLimits(0.1,1000);
    g->GetHistogram()->SetMaximum(1);
    g->GetHistogram()->SetMinimum(1.0e-8);
    g->Draw("A*");
    c->Update();
    g=new TGraph(nMasses,mass,epsilon2_max);
    g->SetTitle(Form("%1.0f eve",minimumevents));
    g->SetMarkerColor(kGreen);
    g->Draw("*");
    p->BuildLegend();
  }
  
}
