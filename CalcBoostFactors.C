

void CalcBoostFactors(double elE=20e3, double prE=250e3){
  double elmass=0.511;
  double prmass=938;
  double aprimemass=500; //note that if you set this too high, math will break down.  must be <<~ CM energy.
  double dummymass=1;
  double thetamindeg=1;//angle wrt beam axis in lab frame, in degrees
  double thetamaxdeg=179;//ditto, max angle.
  double radperdeg=TMath::Pi()/180;

  //note that for XYZM the fourth coordinate is rest mass, not energy!
  ROOT::Math::PxPyPzMVector el; //electron beam in lab frame
  el.SetCoordinates(0,0,elE,elmass);
  ROOT::Math::PxPyPzMVector pr; //proton beam in lab frame
  pr.SetCoordinates(0,0,-prE,prmass);
  ROOT::Math::PxPyPzMVector cm; // center of mass in lab frame
  cm=el+pr;
  ROOT::Math::PxPyPzMVector dummy; // particle at rest in the lab frame.
  dummy.SetCoordinates(0,0,0,dummymass);


  //calculate some boosts we'll need to convert from frame to frame:

  //find boost that gets us from lab frame to proton rest frame:
  ROOT::Math::PxPyPzMVector::BetaVector betaToFixedTarget=pr.BoostToCM();
  ROOT::Math::Boost boostLabToFixed(betaToFixedTarget);
  ROOT::Math::PxPyPzMVector dummy_fixed=boostLabToFixed(dummy);
  ROOT::Math::PxPyPzMVector::BetaVector betaFixedToLab=dummy_fixed.BoostToCM();
  ROOT::Math::Boost boostFixedToLab(betaFixedToLab);

  
  //find boost that gets us from lab frame to cm frame:
  ROOT::Math::PxPyPzMVector::BetaVector betaToCM=cm.BoostToCM();
  ROOT::Math::Boost boostLabToCM(betaToCM);
  ROOT::Math::PxPyPzMVector cm_in_cm=boostLabToCM(cm);
  ROOT::Math::PxPyPzMVector dummy_in_cm=boostLabToCM(dummy);
  ROOT::Math::PxPyPzMVector::BetaVector betaFromCM=dummy_in_cm.BoostToCM();
  ROOT::Math::Boost boostCmToLab(betaFromCM);

  printf("cm E~%3.0fGeV cm4vec=(%f,%f,%f; %f)\n", cm_in_cm.T()/1000,cm_in_cm.X(), cm_in_cm.Y(), cm_in_cm.Z(), cm_in_cm.T());


 //find boost that gets us from cm frame to aprime frame in the simplest beam-->A' assumption:
    ROOT::Math::PxPyPzMVector Aprime_cm;//Aprime vector in cm frame, for kinematic calculations. 
  //assuming E>>mass, easy soln:
  //pA=pP+pE (mag of momentum)
  //pA+pP+pE=pCM (energy terms)
  //hence pA=0.5*pCM 
    Aprime_cm.SetCoordinates(0,0,cm_in_cm.T()/2,aprimemass);
  ROOT::Math::PxPyPzMVector::BetaVector betaToAp=Aprime_cm.BoostToCM();
  ROOT::Math::Boost boostCmToAp(betaToAp);
  ROOT::Math::PxPyPzMVector cm_in_ap=boostCmToAp(cm_in_cm);
  ROOT::Math::PxPyPzMVector::BetaVector betaFromAp=cm_in_ap.BoostToCM();
  ROOT::Math::Boost boostApToCM(betaFromAp);



  //compute the symmetric decay angle in the lab frame:
  ROOT::Math::PxPyPzMVector decay_ap;//symmetric decay electron in aprime frame.
  decay_ap.SetCoordinates(0,250,0,elmass);//at right angles to the boost direction
  //now boost that decay particle all the way back to the lab frame:
  ROOT::Math::PxPyPzMVector decay_lab=boostCmToLab(boostApToCM(decay_ap));


  //compute the fixed-target angles corresponding to practical lab limits:
  //minimum and maximum electron angles in lab frame
  float minrad=thetamindeg*radperdeg;
  float maxrad=thetamaxdeg*radperdeg;
  ROOT::Math::PxPyPzMVector th_min; //electron angular bounds in lab frame
  th_min.SetCoordinates(100*sin(minrad),0,100*cos(minrad),elmass);
  ROOT::Math::PxPyPzMVector th_max; //electron angular bounds in lab frame
  th_max.SetCoordinates(100*sin(maxrad),0,100*cos(maxrad),elmass);
  ROOT::Math::PxPyPzMVector th_min_fixed, th_max_fixed; //electron angular bounds in proton frame
  th_min_fixed=boostLabToFixed(th_min);
  th_max_fixed=boostLabToFixed(th_max);
  double thetamin_fixed=atan2(th_min_fixed.X(),th_min_fixed.Z())/radperdeg;
  double thetamax_fixed=atan2(th_max_fixed.X(),th_max_fixed.Z())/radperdeg;

  
  //compute the fixed-target energies equivalent to the lab frame beams
  ROOT::Math::PxPyPzMVector el_fixed=boostLabToFixed(el);
  printf("electron pz for fixed proton=%f\n",el_fixed.Z());
  printf("electron lab4vec=(%f,%f,%f; %f)\n",el.X(), el.Y(), el.Z(), el.T());
  printf("electron fix4vec=(%f,%f,%f; %f)\n",el_fixed.X(), el_fixed.Y(), el_fixed.Z(), el_fixed.T());
  ROOT::Math::PxPyPzMVector pr_fixed=boostLabToFixed(pr);
  printf("proton pz for fixed proton=%f\n",pr_fixed.Z());
  printf("proton lab4vec=(%f,%f,%f; %f)\n",pr.X(), pr.Y(), pr.Z(), pr.T());
  printf("proton fix4vec=(%f,%f,%f; %f)\n",pr_fixed.X(), pr_fixed.Y(), pr_fixed.Z(), pr_fixed.T());

  
  printf("lab rest pz in fixed frame=%f\n",dummy_fixed.Z());
  printf("dummy lab4vec=(%f,%f,%f; %f)\n",pr.X(), pr.Y(), pr.Z(), pr.T());
  printf("dummy fix4vec=(%f,%f,%f; %f)\n",pr_fixed.X(), pr_fixed.Y(), pr_fixed.Z(), pr_fixed.T());

  el=boostFixedToLab(el_fixed);
  pr=boostFixedToLab(pr_fixed);

  printf("boosting back:\n");
  printf("beta=(%f,%f,%f)\n",betaFixedToLab.X(),betaFixedToLab.Y(),betaFixedToLab.Z());
  printf("proton lab4vec=(%f,%f,%f; %f)\n",pr.X(), pr.Y(), pr.Z(), pr.T());
  printf("electron lab4vec=(%f,%f,%f; %f)\n",el.X(), el.Y(), el.Z(), el.T());

  printf("\nsummary:\n");
  printf("ebeam_e=%f\tebeam_p=%f\n",el_fixed.Z(),pr_fixed.Z());
  printf("min angle %2.2fdeg (lab) = %2.4fdeg (fixed)\n",thetamindeg,thetamin_fixed);
  printf("max angle %2.2fdeg (lab) = %2.4fdeg (fixed)\n",thetamaxdeg,thetamax_fixed);
  printf("symmetric decay in lab frame for mA=%f is theta(deg)=%f\n",
	 aprimemass,atan2(decay_lab.Y(),decay_lab.Z())/radperdeg);


}
