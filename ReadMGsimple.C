
/* Read in my simple Madgraph output and parse it into a ttree for my use, and into a DJANGOH-like EIC file for EICsmear

 */

float GuessMass(TTree * t,int *npart, float *px,float *py,float *pz,int *pid);


void ReadMGsimple(const char* filename="eic20x250_ep_epee.1.ttree.root"){

  const float mA=500;//MeV;
  const float mElec=0.511;//MeV;
  const float elBeamEnergy=20;//GeV
  const float pBeamEnergy=250;//GeV
  const float pBeamMass=0.938;//GeV
  float ensum=elBeamEnergy+pBeamEnergy;
  float endiff=elBeamEnergy-pBeamEnergy;
  const float beamS=ensum*ensum-endiff*endiff;
  

  //todo:
  /*
 - get beam energy from filename

   */

  int nruns=1;//number of runs over the same kinematic range that were used to produce the source file
  //kludgy, but I store this information in the filename:
  TString filenameT=filename;
  if (filenameT.Contains("sum")){
    int startat=filenameT.Index("sum");
    std::sscanf(filename+startat,"sum%d_%*s",&nruns);
  }
  float weightscale=1/(1.0*nruns);
  
  
  TFile *file=TFile::Open(filename,"READ");
  TTree* mTree=(TTree*)(file->Get("madTree"));
  TString dj=filename;
  dj.ReplaceAll("ttree.root","djangoh.txt");
  FILE *djangoh=fopen(dj.Data(),"w");

  const int MaxPart=10;
  Int_t npart;
  Float_t weight_ub; //weight in ub.
  Float_t px[MaxPart];
  Float_t py[MaxPart];
  Float_t pz[MaxPart];
  Float_t ene[MaxPart];
  Int_t pid[MaxPart];
  Char_t pidc[MaxPart];//single character to make it easier to recognize a particle. 
  Int_t charge[MaxPart];//particle charge +1,-1, or zero in some cases.

  mTree->SetBranchAddress("npart",&npart);
  mTree->SetBranchAddress("wub",&weight_ub);
  mTree->SetBranchAddress("pid",pid);
  mTree->SetBranchAddress("pidc",pidc);
  mTree->SetBranchAddress("charge",charge);
  mTree->SetBranchAddress("px",px);
  mTree->SetBranchAddress("py",py);
  mTree->SetBranchAddress("pz",pz);
  mTree->SetBranchAddress("ene",ene);


  int neve=mTree->GetEntries();

  //detect the A' mass by looking for commonalities in the first few events:
  float bestGuessMass;
  bestGuessMass=GuessMass(mTree,&npart,px,py,pz,pid);

  printf("Building DJANGOH and oTree outputs from %s, which is a sum of %d MG runs with a guessed mass of %fMeV\n",filename,nruns,bestGuessMass);
  


  //open the file we'll be writing the parsed tree to:
  TFile *ofile=TFile::Open(Form("%s.otree.root",filename),"RECREATE");

  TTree *oTree=new TTree("oTree","parsed tree for specific event structure");
  TVector3 p,e,es,P;//positron,electron,spec. electron, Proton
  TVector3 e0[2];//unsorted electrons;
  float m[2];//unsorted candidate masses.
  float mA0,mA1,mA2;//inv. mass of three candidates (two good ones first, then same-charge one last)
  float weight_scaled; //scaled for the number of files we combined.  still in ub units.
  int ne=0;//

  oTree->Branch("p",&p);
  oTree->Branch("P",&P);
  oTree->Branch("e",&e);
  oTree->Branch("es",&es);
  oTree->Branch("mA0",&mA0);
  oTree->Branch("mA1",&mA1);
  oTree->Branch("mA2",&mA2);
  oTree->Branch("w",&weight_scaled); //still in ub!



  //now that we have all the branches we need, prime the djangoh output with the proper header:
  fprintf(djangoh,"DJANGOH EVENT FILE -- ACTUALLY MADGRAPH PRETENDING TO BE DJANGOH\n");
  fprintf(djangoh,"============================================\n");
  fprintf(djangoh,"I,ievent,IChannel,process,subprocess,nucleon,struckparton,partontrack,y,Q2,x,W2,Nu,truey,trueQ2,truex,trueW2,trueNu,SIGtot (fb),errSIGtot (fb),Depol,9 structures, nrTracks\n");
  fprintf(djangoh,"============================================\n");
  fprintf(djangoh,"I, 1=stable, CHEP PID, parent I, first daughter, last daughter, px, py, pz, E, m, x0,y0,z0\n");
  fprintf(djangoh,"============================================\n");



  TVector3 zero(0,0,0);
  for (int i=0;i<neve;i++){
    mTree->GetEntry(i);
    p=zero;
    e=zero;
    es=zero;
    P=zero;
    mA2=mA1=mA0=0;
    e0[0]=zero;
    e0[1]=zero;
   
    int ne=0;//number of unsorted electrons in this event 0.
    for (int j=0;j<npart;j++){
      if (pid[j]==2212)	P.SetXYZ(px[j],py[j],pz[j]); //only one proton per event in this assumption;
      if (pid[j]==-11)	p.SetXYZ(px[j],py[j],pz[j]); //only one positron per event in this assumption;
      if (pid[j]==11){
	e0[ne].SetXYZ(px[j],py[j],pz[j]);
	ne++;
      }
    }

    //sort which electron is closer to the correct mass:
    for (int j=0;j<2;j++){
      m[j]=sqrt(-2*p.Dot(e0[j])+2*p.Mag()*e0[j].Mag());
    }
    if (abs(m[0]-bestGuessMass)<abs(m[1]-bestGuessMass)){
      mA0=m[0];
      mA1=m[1];
      e=e0[0];
      es=e0[1];
    } else {
      mA0=m[1];
      mA1=m[0];
      e=e0[1];
      es=e0[0];
    }
    
    weight_scaled=weight_ub*weightscale;
    mA2=sqrt(-2*e.Dot(es)+2*e.Mag()*es.Mag());
    oTree->Fill();

    //now that we have all the variables sorted, calc djangoh variabels as best we can, and print to the djangoh file:

    int dj_I=0;
    int dj_chan=0; //pick a valid entry.
    int dj_process=0;
    int dj_subproc=0;
    int dj_nucleon=0;
    int dj_struck=0;
    int dj_parton=1;//just in case it panics... this is the particle line that is the struck quark
    float costh=cos(es.Theta());
    float coshalf=cos(es.Theta()/2);
    float dj_Q2=2*elBeamEnergy*es.Mag()*(1-costh);
    float dj_y=1-es.Mag()/elBeamEnergy*(coshalf*coshalf);
    float dj_x=dj_Q2/(beamS*dj_y);
    float dj_nu=elBeamEnergy-es.Mag();
    float dj_W2=pBeamMass*pBeamMass;
    
 
    // format:   fprintf(djangoh,"I,ievent,IChannel,process,subprocess,nucleon,struckparton,partontrack,y,Q2,x,W2,Nu,truey,trueQ2,truex,trueW2,trueNu,SIGtot (fb),errSIGtot,Depol,9 structures, nrTracks\n");
    fprintf(djangoh,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",dj_I,i+1,dj_chan,dj_process,dj_subproc,dj_nucleon,dj_struck,dj_parton);
    fprintf(djangoh,"\t%f\t%f\t%f\t%f\t%f",dj_y,dj_Q2,dj_x,dj_W2,dj_nu);//reco'd with fs effects?
    fprintf(djangoh,"\t%f\t%f\t%f\t%f\t%f",dj_y,dj_Q2,dj_x,dj_W2,dj_nu);//true hard collision without radiative process?
    fprintf(djangoh,"\t%f\t%f\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t%d\n",weight_scaled*1e9,0.0,npart);
      fprintf(djangoh,"============================================\n");


    for (int j=0;j<npart;j++){
      // format:   fprintf(djangoh,"I, 1=stable, CHEP PID, parent I, 0=stable, 0=stable, px, py, pz, E, m, x0,y0,z0\n");
      fprintf(djangoh,"%d\t1\t%d\t0\t0\t0\t%f\t%f\t%f\t%f\t%f\t0\t0\t0\n",j+1,pid[j],px[j]/1e3,py[j]/1e3,pz[j]/1e3,ene[j]/1e3,pid[j]==2212?pBeamMass:(mElec/1e3));
    }
  fprintf(djangoh,"=============== Event finished ===============\n");




    
  }

  fclose (djangoh);
  oTree->Write();
  file->Close();
  ofile->Close();

return;
}


float GuessMass(TTree * t,int *npart, float *px,float *py,float *pz,int *pid){
  //look at the first few events to guess the most likely mass, if this is a signal set.

  int neve=t->GetEntries();
  float massguess[2];//both options from the first set
  const float guessrange=1;//MeV we're allowed to be away and still be a good guess
  const int winthreshold=5;//want >5 events within range to establish our 'right mass'
  
  int guesswins[2];//number of events that are within range of this guess.
  TVector3 p,e[2];
  
  for (int i=0;i<neve;i++){
    t->GetEntry(i);
    p.SetXYZ(0,0,0);
    e[0].SetXYZ(0,0,0);
    e[1].SetXYZ(0,0,0);
   
    int ne=0;//number of unsorted electrons in this event 0.
    for (int j=0;j<*npart;j++){
      if (pid[j]==2212)	continue; //don't need the proton.
      if (pid[j]==-11)	p.SetXYZ(px[j],py[j],pz[j]); //only one positron per event in this assumption;
      if (pid[j]==11){
	e[ne].SetXYZ(px[j],py[j],pz[j]);
	ne++;
      }
    }

    for (int j=0;j<2;j++){
      //float mtest=sqrt(mElec*mElec*2+2*p.Dot(e0[j]));//
      float mtest=sqrt(-2*p.Dot(e[j])+2*p.Mag()*e[j].Mag());//neglecting rest mass of electrons

      if (i==0){ //fill mass guesses on first event:
	massguess[j]=mtest;
	guesswins[j]=0;
      } else{
	
	for (int k=0;k<2;k++){
	  if (abs(mtest-massguess[k])<guessrange){
	    guesswins[k]++;	  
	  }
	}
      }
    }
    for (int k=0;k<2;k++){
      if (guesswins[k]>winthreshold) return massguess[k];
    }
  }

  return 0; //if we didn't ever find a mass, assume it's zero.
}
  
