
#include "TLorentzVector.h"

TRandom3 rng;

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

struct EventFourVectors{
  TLorentzVector *P0;//beam proton
  TLorentzVector *e0;//beam electron
  TLorentzVector *Pf;//final state proton
  TLorentzVector *es;//final state spectator e-
  TLorentzVector *A;//intermediate A' or gammastar
  TLorentzVector *e;//decay electron
  TLorentzVector *p;//decay positron
};
EventFourVectors det;//=MakeEventFourVectors(); //detector and fixed-target frame four vectors.
EventFourVectors fix;//=MakeEventFourVectors(); //these are manually associated with the pointers below:
EventFourVectors smear;//=MakeEventFourVectors(); //smeared four vectors


void SetFixedTargetVectors(){
  TVector3 boostToFixed=det.P0->BoostVector();//get the boost from lab to proton beam frame
  *(fix.P0)=*(det.P0); fix.P0->Boost(-1*boostToFixed);
  *(fix.e0)=*(det.e0); fix.e0->Boost(-1*boostToFixed);
  *(fix.Pf)=*(det.Pf); fix.Pf->Boost(-1*boostToFixed);
  *(fix.es)=*(det.es); fix.es->Boost(-1*boostToFixed);
  *(fix.A)=*(det.A); fix.A->Boost(-1*boostToFixed);
  *(fix.e)=*(det.e); fix.e->Boost(-1*boostToFixed);
  *(fix.p)=*(det.p); fix.p->Boost(-1*boostToFixed);
};



TLorentzVector SmearHadron(TLorentzVector x){
  return x; //I don't care about the hadron for now.
}



TLorentzVector SmearLepton(TLorentzVector x){
  if (x.Vect().Mag()==0) return x; //can't smear a zero.
  //use the handbook detector to smear the momentum vector.  WE assume this always washes out the mass of the particle, so we make it massless.
  TLorentzVector smeared(0,0,0,0);
  float eta=x.Eta();
  float sigma=0;
  if (eta<-4.5 || eta>4.5){
    //no tracking or cal
    return smeared;
  }
  smeared=x;
  TVector3 mom=x.Vect();

  if (eta<-3.5){   //EMcal Zone: -- should this override the tracking as well?
    smeared.SetPhi(rng.Gaus(smeared.Phi(),0.001));
    smeared.SetTheta(rng.Gaus(smeared.Theta(),0.001));
    //this applies to all angles, but tracking will do better, so only apply if outside of tracking.
  } else if (eta <-2.5){//backward tracking
    mom.SetMag(rng.Gaus(mom.Mag(),sqrt(pow(0.001*mom.Mag2(),2)+pow(0.02*mom.Mag(),2))));
    smeared.SetVectM(mom,0);
  } else if (eta <-1.0){//backward tracking
    mom.SetMag(rng.Gaus(mom.Mag(),sqrt(pow(0.0005*mom.Mag2(),2)+pow(0.01*mom.Mag(),2))));
    smeared.SetVectM(mom,0);
  } else if (eta <1.0){//mid tracking
    mom.SetMag(rng.Gaus(mom.Mag(),sqrt(pow(0.0005*mom.Mag2(),2)+pow(0.005*mom.Mag(),2))));
    smeared.SetVectM(mom,0);
  } else if (eta <2.5){//forward tracking
    mom.SetMag(rng.Gaus(mom.Mag(),sqrt(pow(0.0005*mom.Mag2(),2)+pow(0.01*mom.Mag(),2))));
    smeared.SetVectM(mom,0);
  } else if (eta <3.5){
    mom.SetMag(rng.Gaus(mom.Mag(),sqrt(pow(0.001*mom.Mag2(),2)+pow(0.02*mom.Mag(),2))));
    smeared.SetVectM(mom,0);
  } if (eta<4.5){    //EMcal Zone:  (see question above.)
    smeared.SetPhi(rng.Gaus(smeared.Phi(),0.001));
    smeared.SetTheta(rng.Gaus(smeared.Theta(),0.001));
    //this applies to all angles, but tracking will do better, so only apply if outside of tracking.
  }

  //we ought to also smear the energy, but the tracking precision will beat that as long as we have tracking, so we only need to smear in the range outside that.
  if (eta<-3.5){
    mom.SetMag(rng.Gaus(mom.Mag(),sqrt(pow(0.02*mom.Mag(),2)+pow(0.01,2)*mom.Mag())));
    smeared.SetVectM(mom,0);
  } else if (eta>3.5){
    mom.SetMag(rng.Gaus(mom.Mag(),sqrt(pow(0.02*mom.Mag(),2)+pow(0.12,2)*mom.Mag())));
    smeared.SetVectM(mom,0);
  } 
  return smeared;
}

void SmearVectors(){
  *(smear.P0)=*(det.P0);//input beam
  *(smear.e0)=*(det.e0);//input beam
  *(smear.Pf)=SmearHadron(*(det.Pf));
  *(smear.es)=SmearLepton(*(det.es));
  *(smear.e)=SmearLepton(*(det.e));
  *(smear.p)=SmearLepton(*(det.p));
  *(smear.A)=*(smear.e)+*(smear.p);
  return;
}

EventFourVectors MakeEventFourVectors(){
  //this is so obviously a class...
  EventFourVectors x;
  x.P0=new TLorentzVector(0,0,0,0);
  x.e0=new TLorentzVector(0,0,0,0);
  x.Pf=new TLorentzVector(0,0,0,0);
  x.es=new TLorentzVector(0,0,0,0);
  x.A=new TLorentzVector(0,0,0,0);
  x.e=new TLorentzVector(0,0,0,0);
  x.p=new TLorentzVector(0,0,0,0);
  return x;
}
  
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
  smear=MakeEventFourVectors();
 
  //initial state:
  TLorentzVector *e04=det.e0=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("e04",&e04);
  TLorentzVector *P04=det.P0=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("P04",&P04);

  //intermediate state:  
  TLorentzVector *A4=det.A=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("A4",&A4);
  TLorentzVector *e14=new TLorentzVector(0,0,0,0); //beam electron after ISR. not from the tree.  We have to build this ourselves
  TLorentzVector *ef4=new TLorentzVector(0,0,0,0); //final state electron before FSR.  not from the tree.  We have to build this ourselves
  float Q2; d.tree->SetBranchAddress("Q2",&Q2);
  
  //final state:
  TLorentzVector *P4=det.Pf=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("P4",&P4);
  TLorentzVector *es4=det.es=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("es4",&es4);
  TLorentzVector *p4=det.p=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("p4",&p4);
  TLorentzVector *e4=det.e=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("e4",&e4);
  //for backward compatibility, we manually define the boosted as well:
  TLorentzVector *P04f,*e04f,*P4f,*p4f,*e4f,*es4f,*e14f, *ef4f,*A4f;
  //initial state
  P04f=fix.P0=new TLorentzVector(*P04);
  e04f=fix.e0=new TLorentzVector(*e04);
  //final states:
  P4f=fix.Pf=new TLorentzVector(*P4);
  es4f=fix.es=new TLorentzVector(*es4);
  p4f=fix.p=new TLorentzVector(0,0,0,0);
  e4f=fix.e=new TLorentzVector(0,0,0,0);

  //tentative intermediate states
  ef4f=new TLorentzVector(0,0,0,0);
  e14f=new TLorentzVector(0,0,0,0);
  A4f=fix.A=new TLorentzVector(0,0,0,0);

  //event weight
  float w; d.tree->SetBranchAddress("w",&w);
  
  d.tree->GetEntry(0);
  d.mass=A4->M();
  float wscale=1.0/d.n; //temporary to check the scale per discussion with Jan
  //d.tree->Draw("es.Z()");






  
  //outputs:
  vector<double> ve,vq2,vq2alt,vq2best,vq2check,vq2checkISR,vq2checkFSR,vQ2,vxs,vw,vth,vcorr;
  TH1F *hXsComparison=new TH1F("hXsComparison","total weight of events vs xs from ep scatter;xs;weight",1e5,1e-12,1e12);
  TH1F *hXsLogComparison=new TH1F("hXsLogComparison","total weight of events vs Log10(xs) from ep scatter;log10(xs);weight/barn-width",100,-12,12);
  TH1F *hAngle=new TH1F("hAngle","total weight of events vs electron theta from ep scatter;theta (rad);weight",50,0,1e-3);
  TH2F *hPositronPolar=new TH2F("hPositronPolar","positron energy and direction;theta (rad);energy (GeV)",60,-TMath::Pi(),TMath::Pi(),50,0,30);
  TH2F *hParticleMom[5];
  TH2F *hPositronMom=hParticleMom[0]=new TH2F("hPositronMom","positron long vs transverse momentum;pZ;pT",100,-30,30,50,0,30);
  TH2F *hElectronMom=hParticleMom[1]=new TH2F("hElectronMom","decay electron long vs transverse momentum;pZ;pT",100,-30,30,50,0,30);
  TH2F *hSpectatorMom=hParticleMom[2]=new TH2F("hSpectatorMom","spectator electron long vs transverse momentum;pZ;pT",100,-30,30,50,0,30);
  TH2F *hAPrimeMom=hParticleMom[3]=new TH2F("hAPrimeMom","heavy photon long vs transverse momentum;pZ;pT",100,-30,30,50,0,30);
  TH2F *hProtonMom=hParticleMom[4]=new TH2F("hProtonMom","proton long vs transverse momentum;pZ;pT",100,-250,250,50,0,30);
  TH2F *hFixParticleMom[5];
  TH2F *hFixPositronMom=hFixParticleMom[0]=new TH2F("hFixPositronMom","positron pT vs pZ in Fixed Target;pZ;pT",100,-1000,12000,50,0,30);
  TH2F *hFixElectronMom=hFixParticleMom[1]=new TH2F("hFixElectronMom","decay electron pT vs pZ in Fixed Target;pZ;pT",100,-1000,12000,50,0,30);
  TH2F *hFixSpectatorMom=hFixParticleMom[2]=new TH2F("hFixSpectatorMom","spectator electron pT vs pZ in Fixed Target;pZ;pT",100,-1000,12000,50,0,30);
  TH2F *hFixAPrimeMom=hFixParticleMom[3]=new TH2F("hFixAPrimeMom","heavy photon pT vs pZ in Fixed Target;pZ;pT",100,-1000,12000,50,0,30);
  TH2F *hFixProtonMom=hFixParticleMom[4]=new TH2F("hFixProtonMom","proton pT vs pZ in Fixed Target;pZ;pT",100,-10,100,50,0,30);

  TH2F* hQ2Angle=new TH2F("hQ2Angle","Q2 and fixed-target scattering angle of elastic e-;theta (rad);Q2 (GeV^2)",100,0,3.14,100,-10,100);
  TH2F* hLogsQ2Angle=new TH2F("hLogsQ2Angle","log10(Q2) and fixed-target scattering angle of elastic e-;log10(theta) (rad);log10(Q2) (GeV^2)",100,-7,1,100,-4,5);
  TH2F* hISRLogsQ2Angle=new TH2F("hISRLogsQ2Angle","log10(Q2) and fixed-target scattering angle assuming ISR A';log10(theta) (rad);log10(Q2) (GeV^2)",100,-7,1,100,-4,5);
  TH2F* hFSRLogsQ2Angle=new TH2F("hFSRLogsQ2Angle","log10(Q2) and fixed-target scattering angle assuming FSR A';log10(theta) (rad);log10(Q2) (GeV^2)",100,-7,1,100,-4,5);
  

  for (int i=0;i<d.n;i++){
    d.tree->GetEntry(i);
    //w=15*w; //hack to check a weighting issue.
    w*=wscale; //scale weight tot he number of entries.  Careful if adding disjoint sets!
    
    //generate the correct fixed target fourvectors:
    SetFixedTargetVectors();//boosts the tree vectors to the proton state.
    //smear the detector vectors into the 'smear' counterparts
    SmearVectors();
    
    //construct the pre-elastic scatter electron assuming the A' radiates off the initial state:
    *e14=*e04-*A4;
    *e14f=*e04f-*A4f;
    
    //construct the post-elastic scatter electron assuming the A' radiates off that pre-final state:
    *ef4=*es4+*A4;
    *ef4f=*es4f+*A4f;

    //plot momentum and direction for the particles:
    hPositronPolar->Fill(det.p->Theta(),det.p->E(),w);
    hPositronMom->Fill(det.p->Vect().Z(),det.p->Vect().Perp(),w);
    hElectronMom->Fill(det.e->Vect().Z(),det.e->Vect().Perp(),w);
    hAPrimeMom->Fill(det.A->Vect().Z(),det.A->Vect().Perp(),w);
    hProtonMom->Fill(det.Pf->Vect().Z(),det.Pf->Vect().Perp(),w);
    hSpectatorMom->Fill(det.es->Vect().Z(),det.es->Vect().Perp(),w);
    hFixPositronMom->Fill(fix.p->Vect().Z(),fix.p->Vect().Perp(),w);
    hFixElectronMom->Fill(fix.e->Vect().Z(),fix.e->Vect().Perp(),w);
    hFixAPrimeMom->Fill(fix.A->Vect().Z(),fix.A->Vect().Perp(),w);
    hFixProtonMom->Fill(fix.Pf->Vect().Z(),fix.Pf->Vect().Perp(),w);
    hFixSpectatorMom->Fill(fix.es->Vect().Z(),fix.es->Vect().Perp(),w);

    //plot Q2 and final-state electron angle for the spectator:
    //q=e_before-e_after=(e_beam-A)-e_final if ISR = e_beam-(e_final+A) if FSR, so q same either way:)
    // TLorentzVector elasticq=*(det.e0)-*(det.A)-*(det.es);
    TLorentzVector protonq=*(det.Pf)-*(det.P0);
    float elasticq2=-protonq*protonq;//elasticq*elasticq;
    hQ2Angle->Fill(fix.es->Theta(),elasticq2,w);
    hLogsQ2Angle->Fill(log10(fix.es->Theta()),log10(elasticq2),w);
    hISRLogsQ2Angle->Fill(log10(fix.es->Angle(e14f->Vect())),log10(elasticq2),w);
    hFSRLogsQ2Angle->Fill(log10(ef4f->Angle(fix.e0->Vect())),log10(elasticq2),w);
    /* outdated:
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
    */
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

    hXsComparison->Fill(xs,w);
    hAngle->Fill(theta,w);

    int bin=hXsLogComparison->FindBin(log10(xs));
    double width=pow(10,hXsLogComparison->GetBinLowEdge(bin+1))-pow(10,hXsLogComparison->GetBinLowEdge(bin));
    //printf("bin width=%E, w/width=%E\n",width,w/width);
    hXsLogComparison->Fill(log10(xs),w/width);
    
    //P04f->Print();
    //boostToFixed.Print();
    

    
  }

  TCanvas *c;
  int nc=0;//number of canvases
  TGraph *g;
  TGraph2D *g2;

 if (0){ //plot the particle directions in the lab frame
    c=new TCanvas(Form("c%d",nc),"canvas",1200,800);
    nc++;
    c->Divide(1,2);
    c->cd(1)->Divide(3,1);
    //hPositronPolar->Draw("lego2 pol");
    //hPositronPolar->Draw("colz");
    //c->cd(2);
    //hPositronMom->Draw("colz");

    for (int i=0;i<3;i++){
      c->cd(1)->cd(i+1)->SetLogz();
      hParticleMom[i]->Draw("colz");
    }
    c->cd(2)->Divide(2,1);
    for (int i=0;i<2;i++){
      c->cd(2)->cd(i+1)->SetLogz();
      hParticleMom[i+3]->Draw("colz");
    }
      
 }

 if (0){ //plot the elastic scatter components in the fixed target frame
    c=new TCanvas(Form("c%d",nc),"canvas",1200,400);
    nc++;
    c->Divide(5,1);
    c->cd(1)->SetLogz();//->Divide(3,1);
    hQ2Angle->Draw("colz");
    c->cd(2)->SetLogz();//->Divide(3,1);
    hLogsQ2Angle->Draw("colz");
    c->cd(3)->SetLogz();//->Divide(3,1);
    hISRLogsQ2Angle->Draw("colz");
    c->cd(4)->SetLogz();//->Divide(3,1);
    hFSRLogsQ2Angle->Draw("colz");
    c->cd(5)->SetLogz();
    int nx=hLogsQ2Angle->GetNbinsX();
    int ny=hLogsQ2Angle->GetNbinsY();
    double logthetas[nx],logq2[ny];
    hLogsQ2Angle->GetXaxis()->GetCenter(logthetas);
    hLogsQ2Angle->GetYaxis()->GetCenter(logq2);
    TH2F *hCorrTemp=new TH2F(*hLogsQ2Angle);
    hCorrTemp->SetTitle("XS Form Factor Correction vs Q2 and scatter angle;log10(theta) (rad);log10(Q2) (GeV^2)");
    for (int i=0;i<nx;i++){
      float th=pow(10,logthetas[i]);
      for (int j=0;j<ny;j++){
	float q2=pow(10,logq2[j]);
	float tau=Tau(Q2,0.938*0.938);
	hCorrTemp->Fill(logthetas[i],logq2[j],GeGmFF(tau,th,q2)/PointFF(tau,th));
      }
    }
    hCorrTemp->Draw("colz");
	  
 }

 if (0){ //plot the particle directions in the fixed target frame
    c=new TCanvas(Form("c%d",nc),"canvas",1200,800);
    nc++;
    c->Divide(1,2);
    c->cd(1)->Divide(3,1);
    //hPositronPolar->Draw("lego2 pol");
    //hPositronPolar->Draw("colz");
    //c->cd(2);
    //hPositronMom->Draw("colz");

    for (int i=0;i<3;i++){
      c->cd(1)->cd(i+1)->SetLogz();
      hFixParticleMom[i]->Draw("colz");
    }
    c->cd(2)->Divide(2,1);
    for (int i=0;i<2;i++){
      c->cd(2)->cd(i+1)->SetLogz();
      hFixParticleMom[i+3]->Draw("colz");
    }
      
 }
  
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


  
  if (0){ //compare weight vs theta to the cross section calculated for same 
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
