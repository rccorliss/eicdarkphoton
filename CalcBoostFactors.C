

void CalcBoostFactors(double elE=20e3, double prE=250e3){
  double elmass=0.511;
  double prmass=938;
  double dummymass=1;

  //note that for XYZM the fourth coordinate is rest mass, not energy!
  ROOT::Math::PxPyPzMVector el; //electron beam in lab frame
  el.SetCoordinates(0,0,elE,elmass);
  ROOT::Math::PxPyPzMVector pr; //proton beam in lab frame
  pr.SetCoordinates(0,0,-prE,prmass);
  ROOT::Math::PxPyPzMVector cm; // center of mass in lab frame
  cm=el+pr;
  ROOT::Math::PxPyPzMVector dummy; // particle at rest in the lab frame.
  dummy.SetCoordinates(0,0,0,dummymass);

  //find boost that gets us from lab frame to proton rest frame:
  ROOT::Math::PxPyPzMVector::BetaVector betaToFixedTarget=pr.BoostToCM();
  ROOT::Math::Boost boostLabToFixed(betaToFixedTarget);

    
  ROOT::Math::PxPyPzMVector el_fixed=boostLabToFixed(el);
  printf("electron pz for fixed proton=%f\n",el_fixed.Z());
  printf("electron lab4vec=(%f,%f,%f; %f)\n",el.X(), el.Y(), el.Z(), el.T());
  printf("electron fix4vec=(%f,%f,%f; %f)\n",el_fixed.X(), el_fixed.Y(), el_fixed.Z(), el_fixed.T());
  ROOT::Math::PxPyPzMVector pr_fixed=boostLabToFixed(pr);
  printf("proton pz for fixed proton=%f\n",pr_fixed.Z());
  printf("proton lab4vec=(%f,%f,%f; %f)\n",pr.X(), pr.Y(), pr.Z(), pr.T());
  printf("proton fix4vec=(%f,%f,%f; %f)\n",pr_fixed.X(), pr_fixed.Y(), pr_fixed.Z(), pr_fixed.T());

  //compute the boost back to the lab frame:
  ROOT::Math::PxPyPzMVector dummy_fixed=boostLabToFixed(dummy);
  ROOT::Math::PxPyPzMVector::BetaVector betaFixedToLab=dummy_fixed.BoostToCM();
  ROOT::Math::Boost boostFixedToLab(betaFixedToLab);

  printf("lab rest pz in fixed frame=%f\n",dummy_fixed.Z());
  printf("dummy lab4vec=(%f,%f,%f; %f)\n",pr.X(), pr.Y(), pr.Z(), pr.T());
  printf("dummy fix4vec=(%f,%f,%f; %f)\n",pr_fixed.X(), pr_fixed.Y(), pr_fixed.Z(), pr_fixed.T());

  el=boostFixedToLab(el_fixed);
  pr=boostFixedToLab(pr_fixed);

  printf("boosting back:\n");
  printf("beta=(%f,%f,%f)\n",betaFixedToLab.X(),betaFixedToLab.Y(),betaFixedToLab.Z());
  printf("proton lab4vec=(%f,%f,%f; %f)\n",pr.X(), pr.Y(), pr.Z(), pr.T());
  printf("electron lab4vec=(%f,%f,%f; %f)\n",el.X(), el.Y(), el.Z(), el.T());




}
