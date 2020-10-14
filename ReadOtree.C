
#include "TLorentzVector.h"

struct Dataset{
  TString filename;
  TFile *file;
  TTree *tree;
  float mass;
  TString name;
  int n;
};
struct FloatErr{
  float v;
  float err;
};

Dataset MakeData(const char* f, float m, const char *tname="oTree",const char *shortname="unset"){
  Dataset temp;
  temp.filename=f;
  temp.mass=m;
  temp.file=NULL;
  temp.file=TFile::Open(f,"R");
  if (temp.file==NULL) {
    printf("couldn't find file: %s\n",f);
    assert(false);
  }
  temp.tree=(TTree*)temp.file->Get(tname);
  temp.name=shortname;
  temp.n=temp.tree->GetEntries();

  //get the mass from the first event:
  //TLorentzVector *Atemp=new TLorentzVector(0,0,0,0);
  //temp.tree->SetBranchAddress("A4",&Atemp);
  //temp.tree->GetEntry(0);
  //temp.mass=Atemp->M();
  
  
  return temp;
}


float Epsilon(float tau, float theta){
  float tanth=tan(theta/2.0);
  float epsinv=1+2*(1+tau)*tanth*tanth;
  return 1/epsinv;
}
float Tau(float Q2, float M2){
  return Q2/(4*M2);
}
float GD2(float Q2){
  return 1/pow((1+Q2/0.71),4);
}
float GeGmFF(float tau, float theta, float Q2){
  static const float musquared=2.79*2.79;
  //return GD2(Q2)+tau/eps*musquared*GD2(Q2); //slightly inefficient
  return (1+2*tau/Epsilon(tau,theta)*musquared*GD2(Q2));
}
float PointFF(float tau, float theta){
  float tanth=tan(theta/2.0);
  //return GD2(Q2)+tau/eps*musquared*GD2(Q2); //slightly inefficient
  return (1+2*tau*tanth*tanth);
}
 

 float MottWithRecoil(float E, float M2, float Q2, float theta){
  static const float alpha2=1.0/137.0/137.0;
  static const float GeVtoUbarn=389.4;
  float costh=cos(theta/2.0);
  float costh2=costh*costh;
  float sinth=sin(theta/2.0);
  float sinth2=sinth*sinth;
  float A=(alpha2*costh2)/(4*E*E*sinth2*sinth2);
  float B=1/(1+Tau(Q2,M2));
  return GeVtoUbarn*A*B;
 }

float MottFF(float E, float M2, float Q2, float theta){
  return MottWithRecoil(E,M2,Q2,theta)*GeGmFF(Tau(Q2,M2), theta, Q2);
}
float MottPoint(float E, float M2, float Q2, float theta){
  return MottWithRecoil(E,M2,Q2,theta)*PointFF(Tau(Q2,M2), theta);
}


void ReadOtree(char *treefile){
  Dataset d=MakeData(treefile,0,"oTree","data");

  //initial state:
  TLorentzVector *e04=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("e04",&e04);
  TLorentzVector *P04=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("P04",&P04);

  //intermediate state:  
  TLorentzVector *A4=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("A4",&A4);
  TLorentzVector *e14=new TLorentzVector(0,0,0,0); //beam electron after ISR. not from the tree.  We have to build this ourselves
  TLorentzVector *ef4=new TLorentzVector(0,0,0,0); //final state electron before FSR.  not from the tree.  We have to build this ourselves
  float Q2; d.tree->SetBranchAddress("Q2",&Q2);
  
  //final state:
  TLorentzVector *P4=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("P4",&P4);
  TLorentzVector *es4=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("es4",&es4);
  TLorentzVector *p4=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("p4",&p4);
  TLorentzVector *e4=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("e4",&e4);

  //event weight
  float w; d.tree->SetBranchAddress("w",&w);
  
  d.tree->GetEntry(0);
  d.mass=A4->M();
  float wscale=15.0/d.n; //temporary to check the scale per discussion with Jan
  //d.tree->Draw("es.Z()");


  //outputs:
  vector<double> ve,vq2,vq2alt,vq2best,vq2check,vq2checkISR,vq2checkFSR,vQ2,vxs,vw,vth,vcorr;
  TH1F *hXsComparison=new TH1F("hXsComparison","total weight of events vs xs from ep scatter;xs;weight",1e5,1e-12,1e12);
  TH1F *hXsLogComparison=new TH1F("hXsLogComparison","total weight of events vs Log10(xs) from ep scatter;log10(xs);weight/barn-width",100,-12,12);
  TH1F *hAngle=new TH1F("hANgle","total weight of events vs electron theta from ep scatter;theta (rad);weight",50,0,1e-5);

  

  for (int i=0;i<d.n;i++){
    d.tree->GetEntry(i);
    //w=15*w; //hack to check a weighting issue.
    //construct the intermediate electron assuming the A' radiates off the initial state:
    *e14=*e04-*A4;
    //construct the interfinal electron assuming the A' radiates off the pre-final state:
    *ef4=*es4+*A4;

    
    //compute the angles in the fixed target frame:
    //first, boost into the fixed target frame:
    TVector3 boostToFixed=P04->BoostVector();
    TLorentzVector *P04f,*e04f,*P4f,*es4f,*e14f, *ef4f,*A4f;
    //initial state
    P04f=new TLorentzVector(*P04); P04f->Boost(-1*boostToFixed);
    e04f=new TLorentzVector(*e04); e04f->Boost(-1*boostToFixed);
    //final states:
    P4f=new TLorentzVector(*P4); P4f->Boost(-1*boostToFixed);
    es4f=new TLorentzVector(*es4); es4f->Boost(-1*boostToFixed);
    //tentative intermediate states
    ef4f=new TLorentzVector(*ef4); ef4f->Boost(-1*boostToFixed);//before FSR
    e14f=new TLorentzVector(*e14); e14f->Boost(-1*boostToFixed);//after ISR
    A4f=new TLorentzVector(*A4); A4f->Boost(-1*boostToFixed);
    //now these are in the fixed target frame.
    double q2check=(*P4f-*P04f)*(*P4f-*P04f);
    double q2checkISR=(*es4f-*e14f)*(*es4f-*e14f);
    double q2checkFSR=(*ef4f-*e04f)*(*ef4f-*e04f);

    
    //the scattering angle is the angle between the incoming and outgoing electron:
    double theta=es4f->Angle(e14f->Vect());
    double sinth=TMath::Sin(theta/2.0);
    double costh=TMath::Cos(theta/2.0);

    //the incoming energy is just the E term off the incoming electron:
    double energy=e14f->E();


    //the 'eprime' recoil-corrected energy is from the eqn:
    double eprime=energy/(1+2*energy/P4f->M()*sinth*sinth);
    double q2=-4*energy*eprime*sinth*sinth;
    if (i<10){
      printf( "i=%d e04f.E=%1.2E  e14f.E=%1.2E  theta=%1.2E  E'=%1.2E\n",i,e04f->E(),energy,theta,eprime);
    }

    //assuming FSR:
    //the scattering angle is the angle between the incoming and outgoing electron:
    double alt_theta=ef4f->Angle(e04f->Vect());
    double alt_sinth=TMath::Sin(alt_theta/2.0);
    double alt_costh=TMath::Cos(alt_theta/2.0);

    //the incoming energy is just the E term off the incoming electron.  in FSR, this is the beam energy:
    double alt_energy=e04f->E();

    //the 'eprime' recoil-corrected energy is from the eqn:
    double alt_eprime=alt_energy/(1+2*alt_energy/P4f->M()*alt_sinth*alt_sinth);
    double alt_q2=-4*alt_energy*alt_eprime*alt_sinth*alt_sinth;

    double bestq2=q2;
    if (abs(Q2+alt_q2)<abs(Q2+q2)) bestq2=alt_q2;
    

    double alpha=1.0/137.0;


    //double xs_0=alpha
    double xs=alpha*alpha*costh*costh/(4*energy*energy*pow(sinth,4));
    xs*=eprime/energy;
    xs*=1 - q2/(2*P4f->M2())*(sinth*sinth)/(costh*costh);

    float tau=Tau(Q2,P4f->M2());
    
    double ffcorrection=GeGmFF(tau,theta,Q2)/PointFF(tau,theta);

    
    //printf("cross section: xs=%E\n",xs);
    ve.push_back(energy);
    vq2.push_back(-q2);
    vq2alt.push_back(-alt_q2);
    vq2check.push_back(-q2check);
    vq2checkISR.push_back(-q2checkISR);
    vq2checkFSR.push_back(-q2checkFSR);
    vq2best.push_back(-bestq2);
    vQ2.push_back(Q2);
    vxs.push_back(xs);
    vw.push_back(w*wscale);
    vth.push_back(theta);
    vcorr.push_back(ffcorrection);

    hXsComparison->Fill(xs,w*wscale);
    hAngle->Fill(theta,w*wscale);

    int bin=hXsLogComparison->FindBin(log10(xs));
    double width=pow(10,hXsLogComparison->GetBinLowEdge(bin+1))-pow(10,hXsLogComparison->GetBinLowEdge(bin));
    //printf("bin width=%E, w/width=%E\n",width,w/width);
    hXsLogComparison->Fill(log10(xs),w*wscale/width);
    
    //P04f->Print();
    //boostToFixed.Print();
    

    
  }

  TCanvas *c;
  int nc=0;//number of canvases
  TGraph *g;
  TGraph2D *g2;

  
  if (0){ //plot the Q^2 comparison
    c=new TCanvas(Form("c%d",nc),"canvas",1200,800);
    nc++;
    c->Divide(3,2);
    c->cd(1);
    g=new TGraph(vq2.size(),&(vq2[0]),&(vQ2[0]));
    g->SetTitle("Q^2 in fixed target (assuming ISR A') vs Q^2 from proton vectors;Q^2 from EE' (GeV^2);Q^2 from proton (GeV^2)");
    g->Draw("A*");
    c->cd(2);
    g=new TGraph(vq2.size(),&(vq2best[0]),&(vQ2[0]));
    g->SetTitle("Q^2 in fixed target (picking best option) vs Q^2 from proton vectors;Q^2 from electrons (GeV^2);Q^2 from proton (GeV^2)");
    g->Draw("A*");
    
    c->cd(3);
    g=new TGraph(vq2.size(),&(vq2alt[0]),&(vQ2[0]));
    g->SetTitle("Q^2 in fixed target (assuming FSR A') vs Q^2 from proton vectors;Q^2 from EE' (GeV^2);Q^2 from proton (GeV^2)");
    g->Draw("A*");
    c->cd(4);
    g=new TGraph(vq2.size(),&(vq2checkISR[0]),&(vQ2[0]));
    g->SetTitle("Q^2 in fixed target ISR electrons vs Q^2 from otree;Q^2 from fixed (GeV^2);Q^2 from otree (GeV^2)");
    g->Draw("A*");
    c->cd(5);
    g=new TGraph(vq2.size(),&(vq2check[0]),&(vQ2[0]));
    g->SetTitle("Q^2 in fixed target protons vs Q^2 from otree;Q^2 from fixed (GeV^2);Q^2 from otree (GeV^2)");
    g->Draw("A*");
    c->cd(6);
    g=new TGraph(vq2.size(),&(vq2checkFSR[0]),&(vQ2[0]));
    g->SetTitle("Q^2 in fixed target FSR electrons vs Q^2 from otree;Q^2 from fixed (GeV^2);Q^2 from otree (GeV^2)");
    g->Draw("A*");
  }
  if (0){ //plot calc'd xs vs theta, compare that to the 
    c=new TCanvas(Form("c%d",nc),"canvas",800,800);
    nc++;

    c->Divide(2,2);
    c->cd(1)->SetLogx();
    c->cd(1)->SetLogy();
    c->cd(1)->SetLogz();
    g=new TGraph(vq2.size(),&(vth[0]),&(vxs[0]));
    g->SetTitle("cross section vs theta (fixed target frame, using ISR assumption);theta (rad);xs arb*barns");
    g->Draw("A*");
    c->cd(2)->SetLogx();
    c->cd(2)->SetLogy();
    c->cd(2)->SetLogz();
    g=new TGraph(vq2.size(),&(vth[0]),&(ve[0]));
    g->SetTitle("energy vs theta (fixed target frame, using ISR assumption);theta (rad);electron energy (GeV)");
    g->Draw("A*");
    c->cd(3)->SetLogx();
    c->cd(3)->SetLogy();
    c->cd(3)->SetLogz();    g2=new TGraph2D(vq2.size(),&(ve[0]),&(vth[0]),&(vxs[0]));
    g2->SetTitle("cross section vs electron energy and scattering angle (fixed target frame);electron E (GeV);theta (rad); xs(arb*barns)");
    g2->Draw("P");
    c->cd(4)->SetLogx();
    c->cd(4)->SetLogy();
    c->cd(4)->SetLogz();
    g=new TGraph(vq2.size(),&(vw[0]),&(vxs[0]));
    g->SetTitle("cross section vs mg weight; mg weight (ub);xs arb*barns");
    g->Draw("A*");
  }

  if (0){ //compare calc'd theta to the madgraph weight. 
    c=new TCanvas(Form("c%d",nc),"canvas",1200,400);
    nc++;

    c->Divide(3,1);
     c->cd(1)->SetLogx();
    c->cd(1)->SetLogy();
    c->cd(1)->SetLogz();
    g=new TGraph(vq2.size(),&(vw[0]),&(vxs[0]));
    g->SetTitle("cross section vs mg weight; mg weight (ub);xs arb*barns");
    g->Draw("A*");
    c->cd(2)->SetLogy();
    hXsLogComparison->Draw();
    c->cd(3)->SetLogx();
    c->cd(3)->SetLogy();
    hXsComparison->Draw("hist");

  }


  
  if (1){ //compare weight vs theta to the cross section calculated for same 
    c=new TCanvas(Form("c%d",nc),"canvas",1200,400);
    nc++;


        c->Divide(3,1);

    c->cd(1)->SetLogx();
    //c->cd(1)->SetLogy();
    //c->cd(1)->SetLogz();
    g=new TGraph(vq2.size(),&(vth[0]),&(vw[0]));
    g->SetTitle("event weight vs fixed-target angle; angle (rad); mg weight (ub)");
    //g->Draw("A*");
    //c->cd(2)->SetLogx();
    c->cd(2)->SetLogy();
    //c->cd(2)->SetLogz();
   //generate the cross section for each bin:
    double centers[hAngle->GetNbinsX()];
    double xscalc[hAngle->GetNbinsX()];
    double xsratio[hAngle->GetNbinsX()];
    double xsrecalc[hAngle->GetNbinsX()];
    double xsffcalc[hAngle->GetNbinsX()];
    hAngle->GetXaxis()->GetCenter(centers);
  

    
    //hAngle->Scale(1/(centers[1]-centers[0]));
    hAngle->Draw();
  

   for (int i=0;i<hAngle->GetNbinsX();i++){
     if (i<10){
       printf("computing center i=%d:  th=%E rads  ve[i]=%E\n",i,centers[i],ve[i]);
     }
      //the 'eprime' recoil-corrected energy is from the eqn:
      double sinth2=pow(sin(centers[i]/2),2);
      double costh2=pow(cos(centers[i]/2),2);
      double eprime=ve[0]/(1+2*ve[0]/0.938*sinth2);
      double q2=-4*ve[0]*eprime*sinth2;
      double alpha=1.0/137.0;
      xscalc[i]=2*TMath::Pi()*sin(centers[i])*(centers[1]-centers[0])*389.4*( (alpha*alpha*costh2/(4*ve[0]*ve[0]*sinth2*sinth2))
									      *eprime/ve[0]
									      *(1 - q2/(2*0.938*0.938)*sinth2/costh2));
      xsrecalc[i]=2*TMath::Pi()*sin(centers[i])*(centers[1]-centers[0])* MottPoint(ve[0],0.938*0.938,-q2,centers[i]);
      xsffcalc[i]=2*TMath::Pi()*sin(centers[i])*(centers[1]-centers[0])* MottFF(ve[0],0.938*0.938,-q2,centers[i]);
      xsratio[i]=xsffcalc[i]/xsrecalc[i];

    }
   g=new TGraph(hAngle->GetNbinsX(),centers,xscalc);
    g->SetTitle("cross section;angle (rad);xs (ub)");
    g->Draw("C");
	g=new TGraph(hAngle->GetNbinsX(),centers,xsrecalc);
    g->SetTitle("Pointlike");
    g->SetLineColor(kRed);
    g->Draw("C");
	g=new TGraph(hAngle->GetNbinsX(),centers,xsffcalc);
    g->SetTitle("Dipole");
    g->SetLineColor(kBlue);
    g->Draw("C");
    c->cd(1);//->SetLogy();
    g=new TGraph(hAngle->GetNbinsX(),xscalc,xsrecalc);
    g->SetTitle("calc vs recalc");
    //g->SetLineColor(kBlue);
    //g->Draw("AC*");
    g=new TGraph(hAngle->GetNbinsX(),centers,xsratio);
    g->SetTitle("ratio of ff vs point;scattering angle (rad);xsff/xspoint");
    g->SetLineColor(kBlue);
    g->SetMarkerColor(kBlue);
    g->Draw("AC*");
      //c->cd(2)->SetLogx();
    c->cd(3);//->SetLogy();
    TH1F *hXsRatio=new TH1F(*hAngle);
    hXsRatio->Reset();
    hXsRatio->SetTitle("Ratio of CalcXs/MadWeight;theta(rad);xs/w");
    for (int i=0;i<hAngle->GetNbinsX();i++){
      hXsRatio->Fill(centers[i],xscalc[i]);
    }
    hXsRatio->Divide(hAngle);
    hXsRatio->Draw("hist");

    //c->cd(2)->SetLogz();
  }




  
  if (0){ //plot the size of the correction vs Q2 for the collision:
    c=new TCanvas(Form("c%d",nc),"canvas",1200,500);
    nc++;
    c->Divide(3,1);
    c->cd(1)->SetLogx();
    g=new TGraph(vq2.size(),&(vQ2[0]),&(vcorr[0]));
    g->SetTitle("Size of FF correction vs Q^2 (assuming ISR A') ;Q^2 (GeV^2);cross section correction");
    g->Draw("A*");
    c->cd(2)->SetLogy();
    TH1F *hCorrection=new TH1F("hCorrection","Form Factor Correction;ratio of FF/pointlike;event weight",200,0,3);
    for (int i=0;i<vq2.size();i++){
      hCorrection->Fill(vcorr[i],vw[i]);
    }
    hCorrection->Draw();
  }
  return;
}
