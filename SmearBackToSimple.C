
float GuessMass(TTree * t,int *npart, float *px,float *py,float *pz,int *pid);


void SmearBackToSimple(const char* filename="sum100_eic20x250_ep_epee_m5GeV_th_1deglab.djangoh.root"){

  //const float mA=500;//MeV;
  const float mElec=0.511;//MeV;
  const float elBeamEnergy=20;//GeV
  const float pBeamEnergy=250;//GeV
  const float pBeamMass=0.938;//GeV
  float ensum=elBeamEnergy+pBeamEnergy;
  float endiff=elBeamEnergy-pBeamEnergy;
  const float beamS=ensum*ensum-endiff*endiff;
     
  // Load the shared library, if not done automaticlly:
   gSystem->Load("libeicsmear.so" );

   gROOT->ProcessLine(".L smearHandBook.cxx");
   TString outputname=filename;
   outputname.ReplaceAll("djangoh.root","smeared.root");
   SmearTree(BuildHandBookDetector(),filename,outputname.Data());
   
   //guess our mass from the raw file:
   TChain unsmeared("EICTree");
   unsmeared.Add(filename.Data());
   erhic::EventDjangoh* unEve(NULL);
   float bestGuessMass=GuessMassFromUnsmeared(unEve,unsmeared);

     printf("Building smeared oTree output from %s, which has a guessed mass of %fMeV\n",filename,bestGuessMass);

   
   // The TTrees are named EICTree. -- nope, smeared trees are called 'Smeared'
   // Create a TChain for trees with this name.
   TChain tree("Smeared");
   
   // Add the file(s) we want to analyse to the chain.
   // We could add multiple files if we wanted.
   tree.Add(outputname.Data()); // Wild cards are allowed e.g. tree.Add("*.root" );
   
   // Create an object to store the current event from the tree.
   // This is how we access the values in the tree.
   // If you want to use generator-specific values, then
   // the event type should match the type in the TTree. Valid types are
   // EventPythia, EventPepsi, EventRapgap, EventDjangoh, EventMilou.
   // If you only need shared quantities like x, Q2 and the particle list
   // you can use EventBase and the macro will be general for any Monte Carlo.
   Smear::Event* event(NULL);// = new EventPythia;
// EventBase* event(NULL);
   
   // Now associate the contents of the branch with the buffer.
   // The events are stored in a branch named event:
   tree.SetBranchAddress("event", &event ); // Note &event, not event.
   
   int nEvents=tree.GetEntries();

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



  
  //open the file we'll be writing the parsed tree to:
  outputname.ReplaceAll("smeared.root","smeared.otree.root");

  TFile *ofile=TFile::Open(outputname.Data(),"RECREATE");

  TTree *oTree=new TTree("oTree","parsed tree for specific event structure");
  TVector3 p,e,es,P;//positron,electron,spec. electron, Proton
  TVector3 e0[2];//unsorted electrons;
  float mA0,mA1,mA2;//inv. mass of three candidates (two good ones first, then same-charge one last)
  float weight_scaled; //scaled for the number of files we combined.  still in ub units.
  int ne=0;//

  oTree->Branch("p",&p);
  oTree->Branch("P",&P);
  oTree->Branch("e",&e);
  oTree->Branch("es",&es);
  oTree->Branch("mA0",&mA0);
  oTree->Branch("mA1",&mA1);
  oTree->Branch("mA2",&mA1);
  oTree->Branch("w",&weight_scaled); //still in ub!


  // Loop over events:
   for(int i=0; i<nEvents; i++){  
      // Read the next entry from the tree.
      tree.GetEntry(i);
      unsmeared.GetEntry(i);//follow along with our unsmeared tree.  If we didn't need it separately, we could have made it a friend.
      
      int npart = event->GetNTracks();
      weight_scaled=unEve->sigTot*1e-9;//convert fb back into ub for oTree convention
      
      //identify and sort the particles 
      int ne=0;//no electrons found to start
      for(int j=0; j < npart; ++j ) {
         const Particle* particle = event->GetTrack(j);
         // Let's just select charged pions for this example:
         int pid = particle->GetPdgCode();
	 TVector3 mom(particle->GetPx(),particle->GetPy(),particle->GetPz());
	 if (pid==2212) P=mom; //proton;
	 if (pid==-11) p=mom; //only one positron per event in this assumption;
	 if (pid==11){
	   e[ne]=mom;
	   ne++;
	 }
      }


      //sort electrons by which electron+positron pair is closer to the correct mass:
      float mdiff=9999; 
      for (int j=0;j<2;j++){
	//float mtest=sqrt(mElec*mElec*2+2*p.Dot(e0[j]));//
	float mtest=sqrt(-2*p.Dot(e0[j])+2*p.Mag()*e0[j].Mag());//neglecting rest mass of electrons
	//printf("found rest mass=%f\n",mtest);
	if (abs(mtest-bestGuessMass)<mdiff){//if we're closer in this pair than the other pair...
	  mdiff=abs(mtest-bestGuessMass);
	  mA1=mA0;
	  mA0=mtest;
	  e=e0[j];
	  es=e0[(j+1)%2];
	} else {
	  mA1=mtest;
	  es=e0[j];
	}
      }
      weight_scaled=weight_ub*weightscale;
      mA2=sqrt(-2*e.Dot(es)+2*e.Mag()*es.Mag());
    oTree->Fill();

      }
   }





  oTree->Write();
  ofile->Close();

return;
}


float GuessMassFromUnsmeared(TChain *t, erhic::EventMC* eve){
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

    //identify our particles:
    int npart=t->GetNTracks();
    int ne=0;//number of unsorted electrons in this event 0.
    for (int j=0;j<npart;j++){
      const Particle* part = eve->GetTrack(j);
      int pid = particle->GetPdgCode();
      TVector3 mom(particle->GetPx(),particle->GetPy(),particle->GetPz());
      if (pid==2212) continue; //don't need the proton.
      if (pid==-11) p=mom; //only one positron per event in this assumption;
      if (pid==11){
	e[ne]=mom;
	ne++;
      }
   
      //calculate the masses of each combination:
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
  
