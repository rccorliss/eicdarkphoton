
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


void ReadOtree(char *treefile){
  Dataset d=MakeData(treefile,0,"oTree","data");

  //initial state:
  TLorentzVector *e04=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("e04",&e04);
  TLorentzVector *P04=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("P04",&P04);

  //intermediate state:  
  TLorentzVector *A4=new TLorentzVector(0,0,0,0); d.tree->SetBranchAddress("A4",&A4);
  TLorentzVector *e14=new TLorentzVector(0,0,0,0); //not from the tree.  We have to build this ourselves
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
  d.tree->Draw("es.Z()");

  vector<double> vq2,vQ2,vxs,vw,vth;

  for (int i=0;i<d.n;i++){
    d.tree->GetEntry(i);

    //construct the intermediate electron assuming the A' radiates off the initial state:
    *e14=*e04-*A4;

    //compute the angles in the fixed target frame:
    //first, boost into the fixed target frame:
    TVector3 boostToFixed=P04->BoostVector();
    TLorentzVector *P04f,*P4f,*es4f,*e14f, *A4f;
    P04f=new TLorentzVector(*P04); P04f->Boost(-1*boostToFixed);
    P4f=new TLorentzVector(*P4); P4f->Boost(-1*boostToFixed);
    es4f=new TLorentzVector(*es4); es4f->Boost(-1*boostToFixed);
    e14f=new TLorentzVector(*e14); e14f->Boost(-1*boostToFixed);
    A4f=new TLorentzVector(*A4); A4f->Boost(-1*boostToFixed);
    //now these are in the fixed target frame.
    //the scattering angle is the angle between the incoming and outgoing electron:
    double theta=es4f->Angle(e14f->Vect());
    double sinth=sin(theta/2);
    double costh=cos(theta/2);

    //the incoming energy is just the E term off the incoming electron:
    double energy=e14f->E();

    //the 'eprime' recoil-corrected energy is from the eqn:
    double eprime=energy/(1+2*energy/P4f->M()*sinth*sinth);
    double q2=-4*energy*eprime*sinth*sinth;

    double alpha=1/137;

    double xs=alpha*alpha*costh*costh/(4*energy*energy*pow(sinth,4));
    xs*=eprime/energy;
    xs*=(1-q2/(2*P4f->M()*P4f->M())*(sinth*sinth)/(costh*costh));

    vq2.push_back(-q2);
    vQ2.push_back(Q2);
    vxs.push_back(xs);
    vw.push_back(w);
    vth.push_back(theta);

    
    //P04f->Print();
    //boostToFixed.Print();
    

    
  }
  TGraph *g;
  g=new TGraph(vq2.size(),&(vq2[0]),&(vQ2[0]));
  g->SetTitle("Q^2 in fixed target (assuming ISR A') vs Q^2 from proton vectors;Q^2 from EE' (GeV^2);Q^2 from proton (GeV^2)");
  g->Draw("A*");
  return;
}
