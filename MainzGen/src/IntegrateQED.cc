#include <stdlib.h>
#include <fstream>
#include <cstring>
#include <random>
#include <pthread.h>
#include <unistd.h>
#include "FourVector.h"
#include "QuasiRandom.h"
#include "Histograms.h"
#include "QEDCalculation.h"

using namespace std;

//generator used for smearing only:
std::default_random_engine randgen(1.); //prime it with a fixed seed of 1. for reproducibility.
std::normal_distribution<float> gaus(0,1); //mean zero, width of 1.

long double m_heavytarget = 168.514728;      //Ta181
const long double m_electron    = 0.000510998928;
const long double deg           = M_PI/180;

char   setupname[255];
long double E0, Eion,
  spec1, spec1mom, spec1thacc, spec1phiacc, spec1momacc,
  spec2, spec2mom, spec2thacc, spec2phiacc, spec2momacc,
  minct, maxct, minctD, maxctD, minphiD, maxphiD, minE, maxE, minm, maxm;

long double maxctScatter,minctScatter,
  minScatterAngle,maxScatterAngle; //to subdivide a rapidly falling weight.

long double events, allevents, accepted = 0;
long double sum = 0;

const int nHists=27;
Hist *id[nHists];

FourVector smear(const FourVector &x)
{//take in a four vector and smear it with random values.
  //NOTE:  here we don't want quasi-random numbers from sobol, since those are designed to be ~evenly spaced.
  FourVector smeared=FourVector(0,0,0,0);

  if (x.abs()==0) return x; //can't smear a zero.
  //use the handbook detector to smear the momentum vector.  WE assume this always washes out the mass of the particle, so we make it massless.

  //components as useful:
  float theta=(float)x.theta();
  float phi=(float)x.phi();
  // don't need:  float ene=(float)x.energy();
  float mom=(float)x.momentum();
  float momSqr=(float)x.momentumSqr();
  float eta=-log(tan(theta/2.)); // we don't need long double precision here, so snuff it early.

  if (eta<-4.5 || eta>4.5){
    //no tracking or cal.  return a zero vector to indicate the particle was lost.
    
    printf("Particle lost!  This should have triggered the earlier cuts!\n");
    printf("eta=%f, theta=%f\n",eta,theta);
    exit(-2);
    return smeared;
  }
  smeared=x;
 
  float momshift=0,phishift=0,thetashift=0;
  //smear all angles, but only smear the energy/momentum if we don't have tracking:
  if (eta<-3.5){   //EMcal Zone: -- should this override the tracking as well?
    momshift=gaus(randgen)*sqrt( 4e-4*momSqr+1e-4*mom ); //sqrt( pow(0.02*mom.Mag(),2)+pow(0.01,2)*mom.Mag()
    phishift=gaus(randgen)*0.001;
    thetashift=gaus(randgen)*0.001;
    //this applies to all angles, but tracking will do better, so only apply if outside of tracking.
  } else if (eta <-2.5){//backward tracking
    momshift=gaus(randgen)*sqrt( 1e-6*momSqr*momSqr+4e-4*momSqr);//sqrt(pow(0.001*mom.Mag2(),2)+pow(0.02*mom.Mag(),2));
  } else if (eta <-1.0){//backward tracking
    momshift=gaus(randgen)*sqrt( 2.5e-7*momSqr*momSqr+1e-4*momSqr);//sqrt(pow(0.0005*mom.Mag2(),2)+pow(0.01*mom.Mag(),2))
  } else if (eta <1.0){//mid tracking
    momshift=gaus(randgen)*sqrt( 2.5e-7*momSqr*momSqr+2.5e-5*momSqr);//sqrt(pow(0.0005*mom.Mag2(),2)+pow(0.005*mom.Mag(),2))
  } else if (eta <2.5){//forward tracking
    momshift=gaus(randgen)*sqrt( 2.5e-7*momSqr*momSqr+1e-4*momSqr);//sqrt(pow(0.0005*mom.Mag2(),2)+pow(0.01*mom.Mag(),2))));
  } else if (eta <3.5){
    momshift=gaus(randgen)*sqrt( 1e-6*momSqr*momSqr+4e-4*momSqr);//sqrt(pow(0.001*mom.Mag2(),2)+pow(0.02*mom.Mag(),2));
  } else if (eta<4.5){    //EMcal Zone:  (see question above.)
    momshift=gaus(randgen)*sqrt( 4e-4*momSqr+1e-4*mom ); //sqrt( pow(0.02*mom.Mag(),2)+pow(0.01,2)*mom.Mag()
    phishift=gaus(randgen)*0.001;
    thetashift=gaus(randgen)*0.001;
     //this applies to all angles, but tracking will do better, so only apply if outside of tracking.
  }

  smeared=Polar(energy(mom+momshift,m_electron),mom+momshift, theta+thetashift, phi+phishift);

  return smeared;
}

void *integrationpart(void *seed)
{
  SobolSequence sobol(8); 
  sobol.init(*(int *) seed * events);

  char outfilename[20];
  sprintf(outfilename,"part_%d.cvs",*(int *) seed);
  std::ofstream outfile(outfilename,std::ofstream::out);
  bool talkingThread=( (*(int *) seed)==0); //only one thread has to do sanity checks.

const FourVector
  e_in_coll=FourVector(E0,0,0,momentum(E0,m_electron)),
  p_in_coll=FourVector(Eion,0,0,-momentum(Eion,m_heavytarget));

//_coll denotes wrt the collider frame.  The rest proceeds assuming wrt the heavy target rest frame unless otherwise noted.
 
 const FourVector
   e_in=e_in_coll.Lorentz(-p_in_coll),
   p_in=p_in_coll.Lorentz(-p_in_coll);

 if (talkingThread){
   printf("sanity check: p_in_coll=(%2.2LE,%2.2LE,%2.2LE,%2.2LE)\n",p_in_coll[0],p_in_coll[1],p_in_coll[2],p_in_coll[3]);
   printf("sanity check: e_in_coll=(%2.2LE,%2.2LE,%2.2LE,%2.2LE)\n",e_in_coll[0],e_in_coll[1],e_in_coll[2],e_in_coll[3]);
   printf("sanity check: p_in=(%2.2LE,%2.2LE,%2.2LE,%2.2LE)\n",p_in[0],p_in[1],p_in[2],p_in[3]);
   printf("sanity check: e_in=(%2.2LE,%2.2LE,%2.2LE,%2.2LE)\n",e_in[0],e_in[1],e_in[2],e_in[3]);

   FourVector specBounds[4];
   specBounds[0]=FourVector(0,sin(spec1-spec1thacc),0,cos(spec1-spec1thacc));
   specBounds[1]=FourVector(0,sin(spec1+spec1thacc),0,cos(spec1+spec1thacc));
   specBounds[2]=FourVector(0,sin(spec2-spec2thacc),0,cos(spec2-spec2thacc));
   specBounds[3]=FourVector(0,sin(spec2+spec2thacc),0,cos(spec2+spec2thacc));

   Momentum m;
   printf("spectrometer bounds:\n");
   printf(" (collider) e-:  %.3LE<TH<%.3LE\t e+:  %.3LE<TH<%.3LE  (degrees)\n",
	  specBounds[0].theta()/deg,specBounds[1].theta()/deg,
	  specBounds[2].theta()/deg,specBounds[3].theta()/deg);
   long double beta=p_in_coll.beta()[2];
   long double gamma=p_in_coll.gamma();
   //angles for beta~1 particles transform as tan(th')=sin(th)/(gamma(cos(th)-beta):
   // thetaprime=atan2(sin(specBounds[0].theta()),gamma*(cos(specBounds[0].theta())-beta));
  printf(" (fixedtar) e-:  %.3LE<TH<%.3LE\t e+:  %.3LE<TH<%.3LE  (degrees)\n",
	 atan2(sin(specBounds[0].theta()),gamma*(cos(specBounds[0].theta())-beta))/deg,
	 atan2(sin(specBounds[1].theta()),gamma*(cos(specBounds[1].theta())-beta))/deg,
	 atan2(sin(specBounds[2].theta()),gamma*(cos(specBounds[2].theta())-beta))/deg,
	 atan2(sin(specBounds[3].theta()),gamma*(cos(specBounds[3].theta())-beta))/deg);
 }

 

 bool spectrometer_mode=true;
 if (Eion>=m_heavytarget) spectrometer_mode=false;

 if (spectrometer_mode){
   if (talkingThread) printf("in spectrometer_mode  (mass=%.1LE>Eion=%.1LE)\n",m_heavytarget,Eion);

 }
 if (!spectrometer_mode){
   if (talkingThread) printf("in !spectrometer_mode  (mass=%.1LE<Eion=%.1LE)\n",m_heavytarget,Eion);
    //because we are in a boosted frame, our energy bounds can't really be imported from the setup file without more consideration.  They need to be the CM energies.
   maxE=e_in[0];
   minE=m_electron;
   if (talkingThread) printf("scattered electron bounds in fixed target frame are %.2LE < E < %.2LE\n",minE,maxE);
   if (talkingThread) {
     printf("scattered electron theta bounds in fixed target frame are %LE < cos(th) < %LE\n",minctScatter,maxctScatter);
     printf("scattered electron theta bounds in fixed target frame are %LE > 1-cos(th) > %LE\n",1.-minctScatter,1.-maxctScatter);
   }
 }
//before my meddling:
// const FourVector 
//    e_in = FourVector(E0, 0, 0, momentum(E0,m_electron)),
//    p_in(m_heavytarget,0,0,0);//energy(ionMom,m_heavytarget), 0, 0, ionMom);

  
  FourVector e_out, q_out, p_out, e1out, e2out, e1spec, e2spec;
  FourVector e_out_coll, q_out_coll, p_out_coll, e1out_coll, e2out_coll;
  long double lepton = m_electron; 


  //previously this was repeatedly computed inside the loop, but it's a constant, so let's just do it once (c++ standards say this will fall out of scope each time through the loop, so we were repeating this fixed calculation unnecessarily.)
  const long double Solidangle = 1/allevents 
    * (maxE-minE)                // dE
    * (maxm-minm)                // dm
    * (maxct-minct)*2*M_PI           // d\Omega
    * 4 * M_PI                   // d\Omega_e
    * (maxctD-minctD) * (maxphiD-minphiD);// d\Omega_Decay

  //originally 'maxct' above was 1

  bool scatteringAngleRangeSmall=(maxctScatter-minctScatter)<1e-14;

  long double minScatterAngleSq=pow(minScatterAngle,2)/2.;
  long double maxScatterAngleSq=pow(maxScatterAngle,2)/2.;
  if (scatteringAngleRangeSmall & talkingThread){
    printf("Scattering Range is very small! \n If the scattering bounds %LE(rad) and %LE(rad) are not both small, consider selecting a different range.\n >>>>  Using cos(x)=1-x^2/2 expansion...\n",minScatterAngle,maxScatterAngle);
  }
 
  
  long double nfail[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
  
  for (long double i=1;i<=events;i++) {
    
    double rndm[8];
    sobol(rndm);

    long double m          = minm + rndm[1]*(maxm-minm); //select a mass for the virtual (or dark) photon
    if (m < 2*lepton) {nfail[0]++;continue;} //if the mass is below the decay lepton mass, veto

    long double E          = minE + rndm[0]*(maxE-minE); //select a total energy for the electron in the fixed target frame
    if (E<m_electron)  {nfail[1]++;continue;} //skip if the electron total energy is below the electron mass.  (daughter mass, or is this the spectator?)

    //determine the direction of the scattered spectator electron in the fixed target frame
    long double cosThetaScatter=0;
    long double thetae=0;
    if (scatteringAngleRangeSmall){
      //if small, then costh terms will be close to 1.0.  Use cos(x)~1-x^2/2:
      //     cosThetaScatter =minctScatter + rndm[4] *(maxctScatter-minctScatter);
      //     cosThetaScatter-1 =minctScatter-1 + rndm[4] *(maxctScatter-minctScatter);
      //     cosThetaScatter-1 =minctScatter-1 + rndm[4] *(maxctScatter+(1-1)-minctScatter);
      //     cosThetaScatter-1 =minctScatter-1 + rndm[4] *( (maxctScatter-1)-(minctScatter-1) );
      //     -th^2/2=-thmax^2/2+rndm*(-thmin^2/2-(-thmax^2/2) );
      //     th=sqrt(thmax^2+rndm*(thmin^2-thmax^2);
      thetae=sqrt(maxScatterAngleSq+rndm[4] *(minScatterAngleSq-maxScatterAngleSq));
      cosThetaScatter=cos(thetae);
      if (0 && log10(thetae/deg)<-10.) {
	printf("log theta <-10:  log=%LE, theta=%LE maxSQ=%LE, minSQ=%LE, (minSQ-maxSQ)=%LE, rnd=%E)\n",
					 log10(thetae/deg),thetae,maxScatterAngleSq,minScatterAngleSq,
					 (minScatterAngleSq-maxScatterAngleSq),rndm[4]);
	exit(-1);
      }
    }else{
      //generate flat in cosTheta and compute the acos to get the scatter angle:
     cosThetaScatter =minctScatter + rndm[4] *(maxctScatter-minctScatter);
      thetae     = acos(cosThetaScatter); //theta between forward and max allowed
    }
    
    long double phie       = rndm[5]*2*M_PI; //phi between zero and 2pi
    e_out = Polar(E,momentum(E,m_electron), thetae, phie);
    if (!spectrometer_mode) e_out_coll=e_out.Lorentz(p_in_coll);

    //determine the direction of the dark/virtual photon in the center of mass frame
    long double theta      = acos(minct+rndm[2]*(maxct-minct)); //previously max was hardcoded to 1
    long double phi        = rndm[3]*2*M_PI-M_PI;

    //determine the direction of the decay positron in the dark/virtual photon frame
    long double thetadecay = acos(minctD + rndm[6]*(maxctD-minctD));
    long double phidecay   = minphiD+rndm[7]*(maxphiD-minphiD);
   
    FourVector cms = e_in + p_in - e_out;
    long double s = cms.square();
    if (s < pow(m + m_heavytarget, 2))  {
      // printf("s=%.2LE but wanted to make mass %.2LE.  electron scatter (fixed frame) E = %.2LE (max=%.2LE, frac=%.2LE)\n",s,m,e_out[0],maxE,e_out[0]/maxE);
      //printf("sanity check: e_out=(%2.2LE,%2.2LE,%2.2LE,%2.2LE)\n",e_out[0],e_out[1],e_out[2],e_out[3]);
      //printf("sanity check: cms=(%2.2LE,%2.2LE,%2.2LE,%2.2LE)\n",cms[0],cms[1],cms[2],cms[3]);
      // break;
      nfail[2]++;continue;}//skip if center of mass energy is not enough to make the dark photon + target

    //calculate the momentum of the virtual/dark photon in the enter of mass frame, using the Kallen triangle function.
    long double kallenTriangle=(s - pow(m + m_heavytarget,2)) * (s - pow(m - m_heavytarget,2));
    long double pAcms = sqrt(kallenTriangle/(4*s) );

    q_out = Polar(energy(m,pAcms),pAcms, theta, phi).Lorentz(cms);
    if (!spectrometer_mode) q_out_coll=q_out.Lorentz(p_in_coll);

    long double peq = sqrt(m*m/4 - lepton*lepton);

    // e1out is positron e2out is electron 
    e1out = Polar(energy(lepton,peq), peq, 
		  thetadecay, phidecay).Lorentz(q_out);


    //generate detector frame quantities and check acceptance of e1:
    //if spectrometer, check spectrometer acceptance.  If not, require minimum momentum and eta cut
    if (spectrometer_mode){
      if ( fabs((e1out.momentum()-spec1mom)/spec1mom)>spec1momacc)  {nfail[3]++;continue;}
       e1spec = FourVector(e1out[0], e1out[1] * cos(spec1) - e1out[3]*sin(spec1),  
			e1out[2], e1out[1] * sin(spec1) + e1out[3]*cos(spec1));
       
      if (fabs(atan2(e1spec[1],e1spec[3])) > spec2thacc)  {nfail[4]++;continue;}
      if (fabs(atan2(e1spec[2],e1spec[3])) > spec2phiacc )  {nfail[5]++;continue;}
    }
    if (!spectrometer_mode) {
      e1out_coll=e1out.Lorentz(p_in_coll);
      if (e1out_coll.momentum()<spec1mom){nfail[3]++;continue;}
      if (fabs(e1out_coll.theta()-spec1)>spec1thacc) {nfail[4]++;continue;}
    }
    
    e2out = q_out-e1out;

    //generate detector frame quantities and check acceptance of e2:
   //if spectrometer, check spectrometer acceptance.  If not, require minimum momentum and eta cut
    if (spectrometer_mode){
      if (fabs((e2out.momentum()-spec2mom)/spec2mom)>spec2momacc) {nfail[6]++;continue;}
      e2spec = FourVector(e2out[0], e2out[1] * cos(spec2) - e2out[3]*sin(spec2),
			  e2out[2], e2out[1] * sin(spec2) + e2out[3]*cos(spec2));
      
      if (fabs(atan2(e2spec[1],e2spec[3])) > spec2thacc)  {nfail[7]++;continue;}
      if (fabs(atan2(e2spec[2],e2spec[3])) > spec2phiacc )  {nfail[8]++;continue;}
    }
    if (!spectrometer_mode) {
      e2out_coll=e2out.Lorentz(p_in_coll);
      if (e2out_coll.momentum()<spec2mom){ nfail[6]++;continue;}
            if (fabs(e2out_coll.theta()-spec2)>spec2thacc) {nfail[7]++;continue;}
    }
    
    FourVector e1_coll_smear=smear(e1out_coll);
    FourVector e2_coll_smear=smear(e2out_coll);
    FourVector q_coll_smear=e1_coll_smear+e2_coll_smear;
    long double smeared_mass=q_coll_smear.mass();
    
    long double thetaD = e2out.Lorentz(-q_out).rotate(q_out).theta();
    long double phiD   = e2out.Lorentz(-q_out).rotate(q_out).phi();
    long double weight = QEDBackground(e_in,e_out,q_out,m,thetaD,phiD)*Solidangle;

    if (isnan(weight)) 
      cout << "WARNING: "<<setprecision(10)<<weight<<" "<<thetae<<" "<<phie<<endl;
    {
      if (0){
       outfile<<m<<","<<weight<<","<<e1out_coll[0]<<","<<e1out_coll[1]<<","<<e1out_coll[2]<<","<<e1out_coll[3]
	     <<","<<e2out_coll[0]<<","<<e2out_coll[1]<<","<<e2out_coll[2]<<","<<e2out_coll[3]<<"\n";//<<","<<
      }
      static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
      pthread_mutex_lock(&mutex); // filling shared histograms must be locked!

      id[ 0]->fill(m, weight);
      id[ 1]->fill(q_out.energy()/E0,                   weight);
      id[ 2]->fill(thetadecay/deg,                      weight);
      id[ 3]->fill(phidecay/deg,                        weight);
      id[ 4]->fill(thetae/deg,                          weight);
      id[ 5]->fill(phie/deg,                            weight);
      id[ 6]->fill(theta/deg,                           weight);
      id[ 7]->fill(phi/deg,                             weight);
      id[ 8]->fill2d(e1out.theta()/deg, e1out.phi()/deg,    weight);
      id[ 9]->fill2d(e2out.theta()/deg, e2out.varPhi()/deg, weight);
      id[10]->fill(100*(e2out.momentum()-spec2mom)/spec2mom, weight);
      id[11]->fill(100*(e1out.momentum()-spec1mom)/spec1mom, weight);
      id[12]->fill(e_out_coll.momentum(),weight);
      id[13]->fill2d(log10(e1out.theta()/deg),log10(weight),1.);
      id[14]->fill2d(log10(e2out.theta()/deg),log10(weight),1.);
      id[15]->fill(log10(weight),1.);
      id[16]->fill2d(log10(e_out.theta()/deg),log10(weight),1.);
      id[17]->fill2d(log10(m),log10(weight),1.);
      id[18]->fill2d(log10(q_out.theta()/deg),log10(weight),1.);
      id[19]->fill2d(log10(angle(q_out,e_out)/deg),log10(weight),1.);
      id[20]->fill2d(log10(angle(e1out,e_out)/deg),log10(weight),1.);
      id[21]->fill2d(log10(angle(e2out,e_out)/deg),log10(weight),1.);
      id[22]->fill2d(log10(thetae/deg),rndm[4],1.);
      id[23]->fill2d(log10(e_out[0]),log10(weight),1.);
      id[24]->fill2d(log10(q_out[0]),log10(weight),1.);
      id[25]->fill2d(log10(q_out_coll[0]),log10(q_out_coll.theta()/deg),weight);
      id[26]->fill2d(m,smeared_mass,weight);
     sum += weight;

      if (!fmod(++accepted,1000)) 
      	cout << fixed<< setprecision(0)<<"\r"<<accepted<<"/" << i << flush;
      pthread_mutex_unlock(&mutex);
    }
  }
  cout << "\r" << fixed<< setprecision(0)<<"\r"<<accepted<<"/" << allevents<< 
    " Thread " <<*(int *)seed<<" done.\n";
  cout <<"Fail matrix: \n";
  for (int i=0;i<9;i++){
    cout << i<< " " << "\r" << fixed<< setprecision(0)<<"\r"<<nfail[i]<<"  \n";
  }
  cout <<"\n";
  return 0;
  outfile.close();
    }

int main(int argc, char * argv[])
{
  
  
  int jobs=1;
  char opt;
  char *filename = NULL;
  char *prefix = (char *) "Integration_";
  events = 1e6;

  // process command line arguments

  while ((opt = getopt(argc,argv, "j:he:f:p:")) != EOF) {
    switch (opt) {
    case 'j' : jobs = atoi(optarg); break;
    case 'h' : cout << "Usage: " << argv[0] << "[options]\n\nOptions:\n"	
		    << "\t-j N        : nr of jobs (threads), used for several cores\n"
		    << "\t-e N        : total number of events\n"	
		    << "\t-f filename : setup file\n"
		    << "\t-p prefix   : prefix for output files\n"
		    << "\t-h          : this help\n"; exit(-1);
    case 'e' : events = atof(optarg); break;
    case 'f' : filename = new char[strlen(optarg)+1];
               strcpy(filename,optarg); break;
    case 'p' : prefix = new char[strlen(optarg)+1];
               strcpy(prefix, optarg); break;
    default: cerr << "Usage: "<<argv[0]<< " events setup\n";exit(-1);
    }
  }


  maxct=0;//default the to angle 0 degrees.  We'll turn it into a costheta later, if it's not overwritten.
  // read parameter file
  if (!filename) {cerr << "No input file given!" << endl; exit(-1);}
  ifstream in(filename);
  while (!in.eof()) {
    char token[255], dummy[255];
    in >> token >> dummy;
    if (!strcmp(token, "setup"))                  in >> setupname;
    if (!strcmp(token, "beam_energy"))            in >> E0;
    if (!strcmp(token, "ion_energy"))             in >> Eion;
    if (!strcmp(token, "ion_mass"))               in >> m_heavytarget;
    if (!strcmp(token, "e-_angle"))               in >> spec2;
    if (!strcmp(token, "e-_momentum"))            in >> spec2mom;
    if (!strcmp(token, "e-_acceptance_phi"))      in >> spec2phiacc;
    if (!strcmp(token, "e-_acceptance_theta"))    in >> spec2thacc;
    if (!strcmp(token, "e-_acceptance_momentum")) in >> spec2momacc;
    if (!strcmp(token, "e+_angle"))               in >> spec1;
    if (!strcmp(token, "e+_momentum"))            in >> spec1mom;
    if (!strcmp(token, "e+_acceptance_phi"))      in >> spec1phiacc;
    if (!strcmp(token, "e+_acceptance_theta"))    in >> spec1thacc;
    if (!strcmp(token, "e+_acceptance_momentum")) in >> spec1momacc;
    if (!strcmp(token, "min_theta"))              in >> maxct;
    if (!strcmp(token, "max_theta"))              in >> minct;
    if (!strcmp(token, "range_mass"))             in >> minm >> maxm;
    if (!strcmp(token, "range_energy"))           in >> minE >> maxE;
    if (!strcmp(token, "range_decay_phi"))        in >> minphiD >> maxphiD;
    if (!strcmp(token, "range_decay_theta"))      in >> maxctD >> minctD;
  }
  in.close();

  spec1   *= deg;              // convert to radian
  spec1phiacc   *= deg;              // convert to radian
  spec1thacc   *= deg;              // convert to radian
  spec2   *= deg;              // convert to radian
  spec2phiacc   *= deg;              // convert to radian
  spec2thacc   *= deg;              // convert to radian
   minphiD *= deg;              // convert to radian
  maxphiD *= deg;              // convert to radian
  minctD  = cos(minctD * deg); // we only need the cosine for generator
  maxctD  = cos(maxctD * deg);
  minct   = cos(minct * deg);
  maxct   = cos(maxct * deg);

  events = floor(events/jobs); // events per thread
  allevents = events*jobs;     // and the sum of all events
  
  // initialize histograms  


  const double minLogAngle=-3, maxLogAngle=1,
    minLogWeight=-20, maxLogWeight=10;
  //  id[ 0]= new Hist("Dark Photon Mass", "m_{{/Symbol g}''}", "", 
  //		   "GeV", "{/Symbol m}b", (int) ((maxm-minm)/0.0005), minm, maxm);
id[ 0]= new Hist("Dark Photon Mass", "$m_{\\gamma}$", "", 
		   "GeV", "$\\mu b$", (int) ((maxm-minm)/0.000025), minm, maxm);

  
  id[ 1]= new Hist("x","", "","","{/Symbol m}b",100,(E0-maxE)/E0,1);
  id[ 2]= new Hist("Decay Angle {/Symbol q}_D","{/Symbol q}_D","","^o","",
		   200,acos(maxctD)/deg, acos(minctD)/deg);
  id[ 3]= new Hist("Decay Angle {/Symbol f}_D","{/Symbol q}_D","","^o","",
	       180,minphiD/deg,maxphiD/deg);
  id[ 4]= new Hist("{/Symbol q}_{e+}", "", "", "^o", "{/Symbol m}b", 180, 0, 180);
  id[ 5]= new Hist("{/Symbol f}_{e+}", "", "", "^o", "{/Symbol m}b", 180, 0, 360);
  id[ 6]= new Hist("{/Symbol q}", "", "", "^o", "{/Symbol m}b", 
		   200, 0, acos(minct)/deg);
  id[ 7]= new Hist("{/Symbol f}", "", "", "^o", "{/Symbol m}b", 180, -180, 180);
  id[ 8]= new Hist("Positron Acceptance", 
		   "{/Symbol q}_{e+}", "{/Symbol f}_e+", "", "^o", "^o", "", 
		   100, fabs(spec1)/deg-180, fabs(spec1)/deg+180, 100, -180, 180);
  id[ 9]= new Hist("Electron Acceptance","{/Symbol q}","{/Symbol f}",
	       "","^o","^o","",
		   100, fabs(spec2)/deg-180, fabs(spec2)/deg+180, 100, -180, 180);
  id[10]= new Hist("Positron Momentum","{/Symbol D}p_e", "","%","",100,-25,25);
  id[11]= new Hist("Electron Momentum","{/Symbol D}p_e", "","%","",100,-25,25);
  id[12]= new Hist("Spectator Momentum (collider frame)","E", "","GeV","",100,0,E0);
  id[13]= new Hist("e1 Angle (fixed target frame) vs Event Weight","log(theta)","log10(weight)",
	       "","log(deg.)","log(mb)","",
		  100, minLogAngle,maxLogAngle, 100, minLogWeight, maxLogWeight);  
  id[14]= new Hist("e2 Angle (fixed target frame) vs Event Weight","log(theta)","log10(weight)",
	       "","log(deg.)","log(mb)","",
		  100, minLogAngle,maxLogAngle, 100, minLogWeight, maxLogWeight);  
  id[15]= new Hist("Event weight","log10(w)", "","","",100,minLogWeight,5);
  id[16]= new Hist("Spectator Angle (fixed target frame) vs Event Weight","log(theta)","log10(weight)",
	       "","log(deg.)","log(mb)","",
		  100, -8,3, 100, minLogWeight, maxLogWeight);
  id[17]= new Hist("Aprime mass vs Event Weight","log(mass)","log10(weight)",
	       "","log(GeV)","log(mb)","",
		   100, -5,2, 100, minLogWeight, maxLogWeight);
  id[18]= new Hist("qout Angle (fixed target frame) vs Event Weight","log(theta)","log10(weight)",
	       "","log(deg.)","log(mb)","",
		  100, minLogAngle,maxLogAngle, 100, minLogWeight, maxLogWeight);  
  id[19]= new Hist("qout-eout Angle (fixed target frame) vs Event Weight","log(theta)","log10(weight)",
	       "","log(deg.)","log(mb)","",
		  100, minLogAngle,maxLogAngle, 100, minLogWeight, maxLogWeight);  
  id[20]= new Hist("e1-eout Angle (fixed target frame) vs Event Weight","log(theta)","log10(weight)",
	       "","log(deg.)","log(mb)","",
		  100, minLogAngle,maxLogAngle, 100, minLogWeight, maxLogWeight);
  id[21]= new Hist("e2-eout Angle (fixed target frame) vs Event Weight","log(theta)","log10(weight)",
	       "","log(deg.)","log(mb)","",
		  100, minLogAngle,maxLogAngle, 100, minLogWeight, maxLogWeight);  
  id[22]= new Hist("log scattering angle vs rndm[4] (are we able to generate a dist?)","log10(th)","rndm[4]",
	       "","","","",
		   100, -10,1, 100, -1, 2.);  
  id[23]= new Hist("Spectator Energy (fixed target frame) vs Event Weight","log10(E)","log10(weight)",
	       "","log(GeV)","log(mb)","",
		  100, -4,4, 100, minLogWeight, maxLogWeight);
  id[24]= new Hist("Aprime Energy (fixed target frame) vs Event Weight","log10(E)","log10(weight)",
	       "","log(GeV)","log(mb)","",
		  100, -4,4, 100, minLogWeight, maxLogWeight);
  id[25]= new Hist("Aprime Energy and angle (collider frame)","log10(E)","log10(theta)",
	       "weight","log(GeV)","log(deg)","mb",
		  100, -4,4, 100, -3, 3);
  id[26]= new Hist("Gen and smeared Aprime mass","mGen","mReco",
	       "weight","GeV","GeV","mb",
		  100, 0,10, 100, 0, 10);

  // start threads
  pthread_t thread[jobs];
  int seed[jobs];


  int nScatterBins=1;
  long double scatterBin[]={1e-9,180,1e-7,1e-6,1e-4,1e-2,1.,10.,180.};
  for (int j=0;j<nScatterBins;j++){
    minScatterAngle=scatterBin[j]*deg;
    maxScatterAngle=scatterBin[j+1]*deg;
    maxctScatter=cos(scatterBin[j]*deg);
    minctScatter=cos(scatterBin[j+1]*deg);
    if (maxctScatter==minctScatter){
      printf("range indistinguishable.  exiting.\n"); exit(-1);
    }
      for (int i=0;i<jobs;i++) {
	seed[i]=i;
	pthread_create(&thread[i], NULL, integrationpart, (void*) &seed[i]);
      }
    for (int i=0;i<jobs;i++) pthread_join(thread[i],NULL);
    cout << "Sum: " << fixed<< setprecision(30)<<sum<<endl;
  }
  // output histograms

  char fn1[255],fn2[255];
  sprintf(fn1,"%s%s.dat",prefix,setupname);
  sprintf(fn2,"%s%s.gnuplot",prefix,setupname);
  ofstream out(fn1), gnuplot(fn2);
  gnuplot << "set terminal postscript enhanced color\n"
	  << "set output '"<<prefix<<setupname<<".ps'\n";
  for (int i=0; i<nHists; i++) {
    out << *id[i];
    id[i]->writeGnuplotCommands(gnuplot, i, fn1);
  }
}
