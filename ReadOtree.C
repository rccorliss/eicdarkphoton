
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

bool PartOk(TLorentzVector four){
  if (four.Vect().Mag()>0) return true;
  if (four.E()>0) return true;
  return false;
}

bool PartOk(TLorentzVector *four){
  return PartOk(*four);
}

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
    //degraded by 10 to prove it works:
    //mom.SetMag(rng.Gaus(mom.Mag(),sqrt(pow(0.005*mom.Mag2(),2)+pow(0.05*mom.Mag(),2))));
    mom.SetMag(rng.Gaus(mom.Mag(),sqrt(pow(0.0005*mom.Mag2(),2)+pow(0.005*mom.Mag(),2))));
    smeared.SetVectM(mom,0);
  } else if (eta <2.5){//forward tracking
    mom.SetMag(rng.Gaus(mom.Mag(),sqrt(pow(0.0005*mom.Mag2(),2)+pow(0.01*mom.Mag(),2))));
    smeared.SetVectM(mom,0);
  } else if (eta <3.5){
    mom.SetMag(rng.Gaus(mom.Mag(),sqrt(pow(0.001*mom.Mag2(),2)+pow(0.02*mom.Mag(),2))));
    smeared.SetVectM(mom,0);
  } else if (eta<4.5){    //EMcal Zone:  (see question above.)
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


void ReadOtree(char *treefile, bool isBg=false){
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
  float wscale=1.0;//Now handled in ReadMGsimple.C.  for old signal, I need to divide here by: /d.n;// per discussion with Jan.
  if (isBg){
    wscale=1.0; //for bg sets where I have summed together many disjoint sets, I do the per-set norm outside of this code
    d.mass=25;
  }
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
  TH2F *hRecoAPrimeMom=new TH2F("hRecoAPrimeMom","Reco A' pT vs pZ;pZ;pT",100,-30,30,50,0,30);

  TH2F* hQ2Angle=new TH2F("hQ2Angle","Q2 and fixed-target scattering angle of elastic e-;theta (rad);Q2 (GeV^2)",100,0,3.14,100,-10,100);
  TH2F* hLogsQ2Angle=new TH2F("hLogsQ2Angle","log10(Q2) and fixed-target scattering angle of elastic e-;log10(theta) (rad);log10(Q2) (GeV^2)",100,-7,1,100,-4,5);
  TH2F* hISRLogsQ2Angle=new TH2F("hISRLogsQ2Angle","log10(Q2) and fixed-target scattering angle assuming ISR A';log10(theta) (rad);log10(Q2) (GeV^2)",100,-7,1,100,-4,5);
  TH2F* hFSRLogsQ2Angle=new TH2F("hFSRLogsQ2Angle","log10(Q2) and fixed-target scattering angle assuming FSR A';log10(theta) (rad);log10(Q2) (GeV^2)",100,-7,1,100,-4,5);

  TH1F* hAmass=new TH1F("hAmass","True intermediate particle mass from MG record;Mass (GeV)",100,d.mass*0.95,d.mass*1.05);
  TH1F* hPairMass=new TH1F("hPairMass","Unsmeared e+e- invariant mass;Mass (GeV)",100,d.mass*0.95,d.mass*1.05);
  TH1F* hRecoMass=new TH1F("hRecoMass","Smeared e+e- invariant mass;Mass (GeV)",100,d.mass*0.95,d.mass*1.05);
  TH1F* hRecoMassWide=new TH1F("hRecoMass","Smeared e+e- invariant mass;Mass (GeV)",200,0,200);
  TH2F* hCompareMass=new TH2F("hCompareMass","Unsmeared e+e- invariant mass vs MG truth;True Mass (GeV);Pair Mass (GeV)",100,d.mass*0.95,d.mass*1.05,100,d.mass*0.95,d.mass*1.05);
  TH2F* hCompareRecoMass=new TH2F("hCompareRecoMass","Smeared e+e- invariant mass vs MG truth;True Mass (GeV);Pair Mass (GeV)",100,d.mass*0.95,d.mass*1.05,100,d.mass*0.95,d.mass*1.05);
  TH2F* hRecoMassPosEta=new TH2F("hRecoMassPosEta","Smeared e+e- invariant mass vs positron eta;eta;Mass (GeV)",20,-5,5,100,d.mass*0.95,d.mass*1.05);
  TH2F* hRecoMassSumEta=new TH2F("hRecoMassSumEta","Smeared e+e- invariant mass vs sum of etas;eta e+ + eta e-;Mass (GeV)",20,-5,5,100,d.mass*0.95,d.mass*1.05);
  TH2F* hTrueMassPosEta=new TH2F("hTrueMassPosEta","true e+e- invariant mass vs positron eta;eta;Mass (GeV)",20,-5,5,100,d.mass*0.95,d.mass*1.05);
  TH2F* hTrueMassSumEta=new TH2F("hTrueMassSumEta","true e+e- invariant mass vs sum of etas;eta e+ + eta e-;Mass (GeV)",20,-5,5,100,d.mass*0.95,d.mass*1.05);
  
  TH2F * hLogDeltaEtaVsEta=new TH2F("hLogDeltaEtaVsEta","Log Delta Eta vs Eta;eta;abs(log10(eta-etatrue))",50,-6,6,50,-6,2);
  TH2F * hLogDeltaPtVsEta=new TH2F("hLogDeltaPtVsEta","Log (#Delta Pt)/Pt vs Eta;eta;abs(log10(pt-pttrue/pttrue))",50,-6,6,50,-6,2);
  TH2F * hDeltaEtaVsEta=new TH2F("hDeltaEtaVsEta","Delta Eta vs Eta;eta;eta-etatrue",50,-6,6,50,-0.04,0.04);
  TH2F * hDeltaPtVsEta=new TH2F("hDeltaPtVsEta"," (#Delta Pt)/Pt vs Eta;eta;pt-pttrue/pttrue",50,-6,6,50,-0.1,0.1);


  TH2F * hTruePairEtas=new TH2F("hTruePairEtas","true e- eta vs e+ eta;eta e+;eta e-",50,-6,6,50,-6,6);
  TH2F * hRecoPairEtas=new TH2F("hRecoPairEtas","smeared e- eta vs e+ eta;eta e+;eta e-",50,-6,6,50,-6,6);
  TH2F * hTruePairPts=new TH2F("hTruePairPts","true e- pT vs e+ pT;pT e+  (GeV);pT e-  (GeV)",50,0,40,50,0,40);
  TH2F * hRecoPairPts=new TH2F("hRecoPairPts","smeared e- pT vs e+ pT;pT e+  (GeV);pT e-  (GeV)",50,0,40,50,0,40);
  TH1F* hTrueDeltaPhi=new TH1F("hTrueDeltaPhi","true e- phi relative to e+;#Delta #phi (rad)", 100,0,2*TMath::Pi());
  TH1F* hRecoDeltaPhi=new TH1F("hRecoDeltaPhi","smeared e- phi relative to e+;#Delta #phi (rad)", 100,0,2*TMath::Pi());

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
    hRecoAPrimeMom->Fill(smear.A->Vect().Z(),smear.A->Vect().Perp(),w);

    hLogDeltaEtaVsEta->Fill(det.p->Eta(),log10(abs(smear.p->Eta()-det.p->Eta())),w);
    hLogDeltaPtVsEta->Fill(det.p->Eta(),log10(abs((smear.p->Pt()-det.p->Pt())/det.p->Pt())),w);
    hDeltaEtaVsEta->Fill(det.p->Eta(),smear.p->Eta()-det.p->Eta(),w);
    hDeltaPtVsEta->Fill(det.p->Eta(),(smear.p->Pt()-det.p->Pt())/det.p->Pt(),w);

    //plot Q2 and final-state electron angle for the spectator:
    //q=e_before-e_after=(e_beam-A)-e_final if ISR = e_beam-(e_final+A) if FSR, so q same either way:)
    // TLorentzVector elasticq=*(det.e0)-*(det.A)-*(det.es);
    TLorentzVector protonq=*(det.Pf)-*(det.P0);
    float elasticq2=-protonq*protonq;//elasticq*elasticq;
    hQ2Angle->Fill(fix.es->Theta(),elasticq2,w);
    hLogsQ2Angle->Fill(log10(fix.es->Theta()),log10(elasticq2),w);
    hISRLogsQ2Angle->Fill(log10(fix.es->Angle(e14f->Vect())),log10(elasticq2),w);
    hFSRLogsQ2Angle->Fill(log10(ef4f->Angle(fix.e0->Vect())),log10(elasticq2),w);


    //compare true A' mass to reco mass:
    TLorentzVector tempA=*(det.p)+*(det.e);
    hAmass->Fill(det.A->M(),w);
    if (PartOk(det.p) && PartOk(det.e)){
      hPairMass->Fill(tempA.M(),w);
      hCompareMass->Fill(det.A->M(),tempA.M(),w);
    }
    if (PartOk(smear.p) && PartOk(smear.e)) {//only plot these if the smear particle exists (is in detector range)
      hRecoMass->Fill(smear.A->M(),w);
      hRecoMassWide->Fill(smear.A->M(),w);
      hCompareRecoMass->Fill(det.A->M(),smear.A->M(),w);
      hRecoMassPosEta->Fill(det.p->Eta(),smear.A->M(),w);
    }
    if (isBg){
      tempA=*(det.p)+*(det.es);
      TLorentzVector tempSmearA=*(smear.p)+*(smear.es);
      if (PartOk(det.p) && PartOk(det.es)) hPairMass->Fill(tempA.M(),w);
      if (PartOk(smear.p) && PartOk(smear.es)){
	  hRecoMass->Fill(tempSmearA.M(),w);
	  hRecoMassWide->Fill(tempSmearA.M(),w);
	  hCompareRecoMass->Fill(tempA.M(),tempSmearA.M(),w);
	  hRecoMassPosEta->Fill(det.p->Eta(),tempSmearA.M(),w);
	}
    }
 

    float detdeltaphi=det.e->Phi()-det.p->Phi();
    if (detdeltaphi<0) detdeltaphi+=2*TMath::Pi();
    float smeardeltaphi=smear.e->Phi()-smear.p->Phi();
    if (smeardeltaphi<0) smeardeltaphi+=2*TMath::Pi();
    if (PartOk(det.p) && PartOk(det.e)){
      hTrueDeltaPhi->Fill(detdeltaphi,w);
      hTruePairEtas->Fill(det.p->Eta(),det.e->Eta(),w);
      hTruePairPts->Fill(det.p->Pt(),det.e->Pt(),w);
    }
    if (PartOk(smear.p) && PartOk(smear.e)){
      hRecoDeltaPhi->Fill(smeardeltaphi,w);
      hRecoPairEtas->Fill(smear.p->Eta(),smear.e->Eta(),w);
      hRecoPairPts->Fill(smear.p->Pt(),smear.e->Pt(),w);
    }


    if (isBg){
      detdeltaphi=det.es->Phi()-det.p->Phi();
      if (detdeltaphi<0) detdeltaphi+=2*TMath::Pi();
      smeardeltaphi=smear.es->Phi()-smear.p->Phi();
      if (smeardeltaphi<0) smeardeltaphi+=2*TMath::Pi();
      if (det.p->Vect().Mag()>0 && det.es->Vect().Mag()>0){
	hTrueDeltaPhi->Fill(detdeltaphi,w);
	hTruePairEtas->Fill(det.p->Eta(),det.es->Eta(),w);
	hTruePairPts->Fill(det.p->Pt(),det.es->Pt(),w);
      }
      if (smear.p->Vect().Mag()>0 && smear.es->Vect().Mag()>0){
	hRecoDeltaPhi->Fill(smeardeltaphi,w);
	hRecoPairEtas->Fill(smear.p->Eta(),smear.es->Eta(),w);
	hRecoPairPts->Fill(smear.p->Pt(),smear.es->Pt(),w);
      }
    }

    
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


 if (1){ //plot the basic summary of the kinematics for the A'
    c=new TCanvas(Form("c%d",nc),"canvas",1200,600);
    nc++;
    c->Divide(1,2);
    c->cd(1)->Divide(5,1);
    c->cd(1)->cd(1)->SetLogy();
    hTrueDeltaPhi->Draw();
    c->cd(1)->cd(2)->SetLogz();
    hTruePairEtas->Draw("colz");
    c->cd(1)->cd(3)->SetLogz();
    hTruePairPts->Draw("colz");
    c->cd(1)->cd(4)->SetLogz();
    hAPrimeMom->Draw("colz");
    c->cd(1)->cd(5)->SetLogy();
    hPairMass->Draw("hist");
    c->cd(2)->Divide(5,1);
    c->cd(2)->cd(1)->SetLogy();
    hRecoDeltaPhi->Draw();
    c->cd(2)->cd(2)->SetLogz();
    hRecoPairEtas->Draw("colz");
    c->cd(2)->cd(3)->SetLogz();
    hRecoPairPts->Draw("colz");
    c->cd(2)->cd(4)->SetLogz();
    hRecoAPrimeMom->Draw("colz");
    c->cd(2)->cd(5)->SetLogy();
    hRecoMass->Draw();
    if (isBg){
      c->SaveAs("overviewBg.pdf");
    } else {
      c->SaveAs(Form("overviewM%2.2fGeV.pdf",d.mass));
    }
    c=new TCanvas(Form("c%d",nc),"canvas",600,400);
    nc++;
    c->cd(0)->SetLogy();
    hRecoMassWide->Draw();
    if (isBg){
      hRecoMassWide->SaveAs("uncutMassSpectrumBg.hist.root");
    } else {
      hRecoMassWide->SaveAs(Form("uncutMassSpectrumM%2.2fGeV.hist.root",d.mass));
      hRecoMass->SaveAs(Form("localMassSpectrumM%2.2fGeV.hist.root",d.mass));
      
    }

 }


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

 if (0){ //plot the elastic scatter quantities in the fixed target frame
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

if (0){ //plot the A' mass with and without detector smearing:
    c=new TCanvas(Form("c%d",nc),"canvas",1200,400);
    nc++;
    c->Divide(4,1);
    c->cd(1)->SetLogy();
    //hAmass->Fill(det.A->M(),w);
    hPairMass->Draw("hist");
    hRecoMass->SetLineColor(kRed);
    hRecoMass->SetMarkerColor(kRed);
    hRecoMass->Draw("same hist");
    c->cd(2)->SetLogz();
    hCompareMass->Draw("colz");
    c->cd(3)->SetLogz();
    hCompareRecoMass->Draw("colz");
    c->cd(4)->SetLogz();
    hRecoMassPosEta->Draw("colz");
 }

 if (0){ //plot the eta and pt differences due to smearing:
    c=new TCanvas(Form("c%d",nc),"canvas",600,600);
    nc++;
    c->Divide(2,2);
    c->cd(1)->SetLogz();
    hLogDeltaEtaVsEta->Draw("colz");
    c->cd(2)->SetLogz();
    hLogDeltaPtVsEta->Draw("colz");
    c->cd(3)->SetLogz();
    hDeltaEtaVsEta->Draw("colz");
    c->cd(4)->SetLogz();
    hDeltaPtVsEta->Draw("colz");

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
