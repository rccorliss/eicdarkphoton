/*
To avoid dependencies on anything complicated, we extract the TRootLHEF* stuff into 
a TTree composed only of standard root classes:

*/

R__LOAD_LIBRARY(../ExRootAnalysis/lib/libExRootAnalysis.dylib)


void translateMGsimple(TString core="ZZZ", TString evePath="EventsX/", int MDver=11, int nGoal=0) {
  
  //the particle conventions in LHEF are:
  int n_codes=6;
  int lhef_code[]={-11,11,-13,13,22,100};//where the proton is hacked in, here, as '100'.
  int hep_code[]={-11,11,-13,13,22,2212};//HEP conventions, I hope.
  char lhef_letter[]={'e','e','u','u','g','p'};
  int lhef_charge[]={+1,-1,+1,-1,0,+1};
   
  //----This block loads the ExRootAnalysis library and reads in the tree----
  // Load shared library to read MadGraph-TTree
  assert(gSystem->Load("../ExRootAnalysis/lib/libExRootAnalysis.so")==0);
   
  // Create chain of root trees
  TChain chain("LHEF"); // name of TTree trunk
  chain.Add(evePath+core+"_weighted_events.root");
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  int maxEve = treeReader->GetEntries();
  printf("see %d events in the tree, MDver=%d, nGoal=%d\n",maxEve,MDver,nGoal);
  assert(maxEve>10);


  // Get pointers to branches used in this analysis
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  TFile *mHfile=new TFile(evePath+core+".ttree.root","RECREATE","file with madgraph event tree");
  assert(mHfile->IsOpen());

  //some accumulators to make sure our math's right:
  float totXsec=0;
  float totXsec2=0;

  
  //construct our simple output tree branches:
  TTree *mTree = new TTree("madTree","plain tree with madgraph events");
  const Int_t MaxPart=10;
  Int_t npart;
  Float_t weight_ub; //weight in ub.
  Float_t px[MaxPart];
  Float_t py[MaxPart];
  Float_t pz[MaxPart];
  Float_t ene[MaxPart];
  Int_t pid[MaxPart];
  Char_t pidc[MaxPart];//single character to make it easier to recognize a particle. 
  Int_t charge[MaxPart];//particle charge +1,-1, or zero in some cases.

  mTree->Branch("npart",&npart,"npart/I");
  mTree->Branch("wub",&weight_ub,"wub/F");
  mTree->Branch("pid",pid,"pid[npart]/I");
  mTree->Branch("pidc",pidc,"pidc[npart]/B");
  mTree->Branch("charge",charge,"charge[npart]/I");
  //I wish I could use TVector3, but we can only do variable-length arrays of primitives, not class objects, apparently.
  mTree->Branch("px",px,"px[npart]/F");
  mTree->Branch("py",py,"py[npart]/F");
  mTree->Branch("pz",pz,"pz[npart]/F");
  mTree->Branch("ene",ene,"ene[npart]/F");
  


  //calculate the proper lorentz transformations:
  double prmass=938;
  double prene=250e3;//note this needs to be synched up to the run_card.dat energy!
  double dummymass=1;

  //note that for XYZM the fourth coordinate is rest mass, not energy!
  ROOT::Math::PxPyPzMVector pr; //proton beam in lab frame
  pr.SetCoordinates(0,0,-prene,prmass);
  ROOT::Math::PxPyPzMVector dummy; //proton beam in lab frame
  dummy.SetCoordinates(0,0,0,dummymass);
  
  //find boost that gets us from lab frame to proton rest frame:
  ROOT::Math::PxPyPzMVector::BetaVector betaToFixedTarget=pr.BoostToCM();
  ROOT::Math::Boost boostLabToFixed(betaToFixedTarget);

  //compute the boost back to the lab frame (these should be the  
  ROOT::Math::PxPyPzMVector dummy_fixed=boostLabToFixed(dummy);
  ROOT::Math::PxPyPzMVector::BetaVector betaFixedToLab=dummy_fixed.BoostToCM();
  ROOT::Math::Boost boostFixedToLab(betaFixedToLab);

  ROOT::Math::PxPyPzEVector p_fixed,p_lab; //particle 4vectors for use in the loop.  Note this is an 'E' not an 'M' assignment
  

  for(int ieve=0;ieve<maxEve;ieve++) {
    // unpack MadGraph  TTree ......
    treeReader->ReadEntry(ieve); // This reads in the first event

    TRootLHEFEvent *event = (TRootLHEFEvent*) branchEvent->At(0);

    
    // new MadGraph: MeV & W in ub
    /* Jesse Thaler suggested weight renormalization:
       W= w_madGr*n_madGr
       with the use case:  h->Fill(pT,W)/ sumW
      - it allows concatenating multiple files
      - weight of each event is approximately equal the x-section of given type of events
    */   
    assert(event->Weight>0);
    if (nGoal<=0)
      weight_ub=event->Weight*maxEve;  
    else // use it to rescale the weights
      weight_ub=event->Weight*nGoal;  



    npart=branchParticle->GetEntries()-2; //the first two particles are inbound, so we skip them.
    if(ieve<5) printf("eve=%d  weight=%e  nPart=%d\n",ieve,event->Weight,npart);


    //loop over particles STARTING AT THE OUTBOUND ONES
    for(int ipart=0;ipart<npart;ipart++) {
      TRootLHEFParticle *particle = (TRootLHEFParticle*) branchParticle->At(ipart+2); // creates a pointer to the 1st particle 
      int myPID=particle->PID;
      if(myPID==100) myPID=2212; // change proton to Geant HEP convention.  Is this still needed?  If proton was specially defined here, it wouldn't have gotten the proper ID, so maybe.
      // http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
      // 11=e-,  gam=22, 2212=p
      // commented out to allow other final states: assert(myPID==11 || myPID==11 || myPID==-11 || myPID==2212 || myPID==101 || myPID==22);

      
      pid[ipart]=myPID;
      for (int i=0;i<n_codes;i++){
	if (myPID==hep_code[i]){
	  pidc[ipart]=lhef_letter[i];
	  charge[ipart]=lhef_charge[i];
	  break; //if we found it, no need to keep looking
	}
      }

      //in GeV, presumably
      p_fixed.SetCoordinates(particle->Px,particle->Py,particle->Pz,particle->E);
      p_lab=boostFixedToLab(p_fixed);
      
      px[ipart]=p_lab.X();
      py[ipart]=p_lab.Y();
      pz[ipart]=p_lab.Z();
      ene[ipart]=p_lab.E();
      

    }// one event is unpacked & saved
    
    //if(ieve<5)daliEve.print();
    mTree->Fill();
    totXsec+=event->Weight;
    totXsec2+=weight_ub;
    //break;
  } // end of loop over events
  //mTree->Show(1);// event id startting from 0
  //mTree->Print();

  mTree->Write();

  double eps=totXsec2/totXsec-1.;
  printf("\n test of weights acc Xsect2(pb)=%.4g  eps=%.4g\n",totXsec2,eps);
  assert( fabs(eps)<0.01);

}


