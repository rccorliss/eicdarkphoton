//                                                                    -*-c++-*-
// $Id: Generate.h 2469 2012-07-04 06:50:19Z merkel $

#ifndef __Generate_h__
#define __Generate_h__





#include "FourVector/FourVector.h"   // relativistic four vector arithmetics
#include "Cola/QuasiRandom.h"             // Random Generator
#include "Cola/Reactions.h"                // Definition of Reactions
#ifdef __ColaMIT__
#include "ColaMIT/rundatabase.h"	// only for struct SIMUL
#else
#include "Cola/rundatabase.h"	// only for struct SIMUL
#endif
#include "Model/Model.h"             // Model library
#include "simDetectorBase.h"
#include "Cola/Masses.h"
#include "QED.h"
#include "PWIA.h"
#include "GridInterpolation.h"
#include "DM_QED.h"
#include <cstdlib>


#include <fstream>
#include <cstring>


namespace He3HACK
{
  extern ifstream csfile;
  extern ofstream listfile;
  extern int mode;
  extern int num;
  extern float buf[16];

}


const double rad     = M_PI/180;
//test
void plotPolarization();

double ElasticCrossSection(const FourVector &e_in, 
			   const FourVector &e_out,
			   const FourVector &q_out,
			   const double helicity = 0,
			   const unsigned short a = 1,
			   const unsigned short b = 0,
			   class Momentum *P_Spin = NULL,
			   class Momentum *P_Spin_CM = NULL,
			   double *Msqr = NULL,double *k=NULL,double kp=0);

double DMBosonCrossSection(const FourVector &e_in,  const FourVector &e_out,
                           const FourVector &q_out, const double mA);

class eventGenerator {
protected:
  reaction *Reaction;
  double Ex;              // Excitation Energy
  simDetectorBase *sime, *sim1, *sim2, *simD1, *simD2;
  double targetmass, dcte, dphie, dpcms;
  double Gamma;           // virtual photon flux
  double weight;          // weight for each event
  double scalefactor;     // weight for each event

  double q2_o, s_o, eps_o;

  void   generateLabAngles(Particle *P, double p, double theta0, double phi0, 
			 double dcostheta, double dphi);
  double eout_p_Q2_Phi(double setq2, double dpcms, double dq2);
  double eout_W_Q2_Phi(double setq2, double dW,    double dq2);
  double eout_Q2_Theta_Phi(double setq2, double dq2);
  double p_cms(double s, double m1, double m2) {
    return sqrt((s - (m1+m2)*(m1+m2)) * (s - (m1-m2)*(m1-m2)) / 4 / s); 
  }
  double calc_dphie() {
    Momentum dir;  
    dir.initPolar(1.0, M_PI/2, sime->getDphi()); // upper corner for limits
    dir.rot_theta(fabs(sime->getAngle()) - sime->getDtheta() - M_PI/2);       
    return dir.phi();
  }

public: 
  char *Label;
  char *Unit;
  virtual ~eventGenerator() { ; }
  virtual double generateEvent(double helicity) = 0;
  virtual double integrationVolume() = 0;
  virtual double crosssectionOfThisEvent() { 
    cerr <<" Event generator broken "<<endl;
    exit(-1);
    return 0;
  }
  double getGamma() { return Gamma;};
  double getScale() { return scalefactor;};
  double generateElectron(double setq2, double dW, double dq2, double dpcms) {
    if (dW    != 0) return eout_W_Q2_Phi(setq2,dW,dq2);
    if (dpcms != 0) return eout_p_Q2_Phi(setq2,dpcms,dq2);
    return eout_Q2_Theta_Phi(setq2,dq2);
  };
  double eoutIntegrationVolume(double dW, double dq2, double dpcms);
};

class generateElastic : public eventGenerator {
public:
  virtual ~generateElastic() { ; }
  generateElastic(reaction *r, simDetectorBase *Sime, simDetectorBase *Sim1,
		  SIMUL *rundb)
  { 
    Label      = "d[W]'";
    Unit       = "sr";
    Reaction   = r; sime=Sime; sim1=Sim1;
    targetmass = Reaction->getTarget()->getMass();
    dcte       = sin(sime->getDtheta());
    dphie      = calc_dphie();
    Ex         = rundb->excitation;
    Gamma = 1;
  };
  
  double generateEvent(double helicity);
  double integrationVolume() { return 4 * sime->getDphi() * dcte; }
};

static const double mup    = 2.79285;

class generateDMBoson : public eventGenerator {
  double m_A, Z, A, LOG, CH, ctmin, ctmax, phmin, phmax, Emax, Emin;
  // (A18) of J. Bjorken et al., PRD 80, 075018 (2009)
  inline double G2el(double t) {
    double a = 111.0 * pow(Z,-1.0/3)/m_electron;
    double d = 0.164 * pow(A,-2.0/3);
    return pow(a*a*t/(1+a*a*t)/(1+t/d)*Z, 2);
  }
  // (A19) of J. Bjorken et al., PRD 80, 075018 (2009)
  inline double G2in(double t) {
    double a = 773.0 * pow(Z, -2.0/3)/m_electron;
    return pow(a*a*t/(1+a*a*t)
	       *((1+t/4/pow(m_proton,2)*(mup*mup-1))/pow(1+t/0.71,4)),2)*Z;
  }
  // (A17) of J. Bjorken et al., PRD 80, 075018 (2009)
  inline double Chi(double m, double E0) {
    double res=0, tmin=pow(m * m / 2 / E0, 2), tmax = m * m;
    int nx=10000;
    for (int i=0;i<=nx;i++) {
      double t = tmin + i * (tmax - tmin) / nx;
      res += (t - tmin) / t / t * (G2el(t) + G2in(t)) / (i==0||i==nx ? 2 : 1);
    }
    return res * (tmax - tmin) / nx;
  }
  // dsigma/(dx), (8) of J. Bjorken et al., PRD 80, 075018 (2009)
  double crosssection(double E0, double E, double m, double eps) {
    double x = E / E0;
    if (LOG==0) LOG = Chi(m, E0) / Z / Z;
    return 8 * pow(Z * eps / m, 2) * pow(alpha,3)*x*LOG*(1+x*x/3/(1-x))*mubarn;
  }
  // dsigma/(dx dcostheta), (5) of J. Bjorken et al., PRD 80, 075018 (2009)
  double crosssection(double E0, double E, double theta, double m, double eps)
  {
    double x = E / E0;
    if (LOG==0) LOG = Chi(m, E0) / Z / Z;
    double U = E0*E0*x*theta*theta + m*m*(1-x)/x + pow(m_electron,2)*x;
    return 8*pow(Z*eps*E0/U,2)*pow(alpha,3)*x*LOG
      * ((1-x+x*x/2)-x*(1-x)*m*m*(E0*E0*x*theta*theta)/U/U)*mubarn;
  }
  // Bethe-Heitler cs, (C5) von J. Bjorken et al., PRD 80, 075018 (2009)
  double BHBackground(double E0, double E, double theta, 
		      double m, double dm, double c) 
  { 
    double x=E/E0;
    if (CH == 0) CH = 2 * pow(alpha,4) * Chi(m, E0) / M_PI;
    return CH * dm / (m*m*m) * (1+(1-x)*(1-x))/(theta*theta*x)
      * ((1+c)/(1-c)+ (1-c)/(1+c))*mubarn;
  }
  // Form Factor of a homogenous sphere with radius R
  double F_sphere(double q, double R) {
    if (q<1e-4) return 1;
    return 3/pow(q*R,2)*(sin(q*R)/(q*R)-cos(q*R));
  }
  // Coherent cross section for A' production from heavy target
  double DMHeavyCS(const FourVector &e_in,  const FourVector &e_out,
		   const FourVector &q_out, const double mA)
  {
    FourVector p_in    = FourVector(targetmass, 0, 0, 0);  
    FourVector CM      = p_in + e_in - e_out;
    FourVector e_inCM  = e_in.Lorentz(-CM).rotate(CM);
    FourVector e_outCM = e_out.Lorentz(-CM).rotate(CM);
    FourVector q_outCM = q_out.Lorentz(-CM).rotate(CM);
    FourVector p_inCM  = p_in.Lorentz(-CM).rotate(CM);
    FourVector p_outCM = e_inCM - e_outCM + p_inCM - q_outCM;
    FourVector q_p     = e_inCM - e_outCM;
    FourVector q       = p_outCM - p_inCM;
    FourVector ne1     = e_inCM  - q_outCM;           // internal electron line
    FourVector ne2     = e_outCM + q_outCM;           // internal electron line

    double q2          = q.square();                    // virtual photon q^2 
    Complex i_e3_q2    =  i_e3 / q2;                    // electron charge -
    double M_square = 0;
    
    Spinor ei[2] = {Spinor(e_inCM,  0.5), Spinor(e_inCM,  -0.5)};
    Spinor eo[2] = {Spinor(e_outCM, 0.5), Spinor(e_outCM, -0.5)};
    Spinor po[2] = {Spinor(p_outCM, 0.5), Spinor(p_outCM, -0.5)};
    Spinor pi[2] = {Spinor(p_inCM,  0.5), Spinor(p_inCM,  -0.5)};
    
    Tensor t1 = (dag(ne1) + eID) / (ne1*ne1 - m_electron*m_electron);
    Tensor t2 = (dag(ne2) + eID) / (ne2*ne2 - m_electron*m_electron);
    Tensor ee2gmunu[16], gNNnu[4];
    
    for(int mu=0; mu<4; mu++) {
      gNNnu[mu]=Z*F_sphere(sqrt(-q.square()),1.21/0.19733*pow(A,1.0/3))*gam[mu];
      for(int nu=0; nu<4; nu++)
	ee2gmunu[nu+(mu<<2)]=gam[nu]*t1*gam[mu] + gam[mu]*t2*gam[nu];
    }
    
    FourVector v1(0,1,0,0), v2(0,0,1,0), 
    v3(q_outCM.momentum()/mA,0,0,q_outCM.energy()/mA);

    for(int spin_e_in  = 0; spin_e_in  <= 1; spin_e_in++ )
      for(int spin_p_in  = 0; spin_p_in  <= 1; spin_p_in++ )
	for(int spin_e_out = 0; spin_e_out <= 1; spin_e_out++) {
	  for(int spin_p_out = 0; spin_p_out <= 1; spin_p_out++) {
	    Complex M_if[4] = {0,0,0,0};
	    for(int mu=0; mu<4; mu++) {
	      for (int nu=0; nu<4; nu++) {
		Tensor ee2g = ee2gmunu[nu+(mu<<2)];
		Complex el_part = !eo[spin_e_out] * ee2g      * ei[spin_e_in];
		Complex nu_part = !po[spin_p_out] * gNNnu[nu] * pi[spin_p_in];
		M_if[mu] += i_e3_q2 * gmn[nu] * el_part * nu_part;
	      }  
	    }
	    for(int mu=0; mu<4; mu++) 
	      for (int nu=0; nu<4; nu++) 
		M_square += real(conj(M_if[mu])*M_if[nu]
				 *(-g_mu_nu[mu][nu] +gmn[mu]*gmn[nu]
				   *q_outCM[mu]*q_outCM[nu]/mA/mA));
	  }    
	}
    const double c = mubarn/64/pow(2*M_PI, 5)/targetmass; 
    double s = CM.square();
    return c*(s-targetmass*targetmass)
      /s*e_out.momentum()/e_in.momentum()*M_square/4;
  }
  // the same in the varibles of Bjorken
  double heavycrosssection(double E0, double E, 
			   double theta, double m, double eps)
  {
    if (E<m) return 0;
    FourVector e_in(energy(m_electron,E0),0,0,E0);
    FourVector q_out(E, 0, momentum(E,m)*sin(theta), momentum(E,m)*cos(theta));
    FourVector e_out,p_in(targetmass,0,0,0);
    
    FourVector epsystem = e_in + p_in - q_out;// = e_out + p_out;
    double s_ep = epsystem * epsystem;
    if (s_ep < pow(m_electron+targetmass, 2)) return 0;
    double p = sqrt((s_ep - pow(m_electron + targetmass,2)) *
		    (s_ep - pow(m_electron - targetmass,2))/4/s_ep); 
    static PseudoRandom rndm;
    e_out = Polar(energy(m_electron,p), p, 
		  acos(rndm()*2-1), rndm()*2*M_PI).Lorentz(epsystem);
    return DMHeavyCS(e_in,  e_out, q_out, m)*4*M_PI;
  }

public:
  virtual ~generateDMBoson() { ; }
  generateDMBoson(reaction *r, SIMUL *rdb)
  { 
    Label      = "dE_A'!";
    Unit       = "GeV";
    Z = r->getTarget()->getCharge();
    A = (int) (r->getTarget()->getMass()/(0.93827231-0.007)+0.5);
    LOG = CH = 0;

    Reaction   = r;
    m_A = rdb->massA/1000;
    Emin = rdb->BHmin/1000;
    ctmin  = cos(rdb->CMSTheta[1] * rad);
    ctmax  = cos(rdb->CMSTheta[0] * rad);
    phmin  = rdb->CMSPhi[0] * rad;
    phmax  = rdb->CMSPhi[1] * rad;
    targetmass = Reaction->getTarget()->getMass();
    Gamma = 1;
  }; 

  double generateEvent(double);
  double integrationVolume() { 
    return (rundb.Ebeam/1000 - Emin)                 // dE_A'
      * (m_A<0 ? 0.100 : 1)
      * (ctmax - ctmin) * (phmax - phmin) / 4 / M_PI; 
  }
};

class generateDMQEDBackground : public eventGenerator {
  double m_A, Z, ctmin, ctmax, phmin, phmax, Emax, Emin;
  FourVector A;
  double mA, mtheta, mphi, jacobian;
public:
  virtual ~generateDMQEDBackground() { ; }
  generateDMQEDBackground(reaction *r, SIMUL *rdb)
  { 
    Label      = "dE'!d[W]'d[W]_ee!dm_ee!d[W]_decay!";
    Unit       = "GeV^2!sr^3!";
    Z = r->getTarget()->getCharge();
    //    A = (int) (r->getTarget()->getMass()/(0.93827231-0.007)+0.5);

    Reaction   = r;
    m_A = rdb->massA/1000;
    Emin = rdb->BHmin/1000;
    ctmin  = cos(rdb->CMSTheta[1] * rad);
    ctmax  = cos(rdb->CMSTheta[0] * rad);
    phmin  = rdb->CMSPhi[0] * rad;
    phmax  = rdb->CMSPhi[1] * rad;
    targetmass = Reaction->getTarget()->getMass();
    Gamma = 1;
  }; 

  double generateEvent(double);
  double crosssectionOfThisEvent() {
     return QEDBackground(Reaction->electronIn, 
			  Reaction->electronOut, 
			  A, mA, mtheta, mphi) * jacobian;
  }
  double integrationVolume() { 
    return (rundb.Ebeam/1000 - Emin)          // dE'
      * (4 * M_PI)                            // dΩ_E'
      * 0.35*fabs(m_A)                        // dm_ee
      * (4 * M_PI)                            // dΩ_ee
      * (ctmax - ctmin) * (phmax - phmin);    // dΩ_Decay
  }
};

class generateElasticCollider : public eventGenerator {
public:
  virtual ~generateElasticCollider() { ; }
  generateElasticCollider(reaction *r, simDetectorBase *Sime,
			  simDetectorBase *Sim1, SIMUL *rundb)
  { 
    Label      = "d[W]'";
    Unit       = "sr";
    Reaction   = r; sime=Sime; sim1=Sim1;
    targetmass = Reaction->getTarget()->getMass();
    dcte       = sin(sime->getDtheta());
    dphie      = calc_dphie();
    Ex         = rundb->excitation;
    Gamma = 1;
  };
  
  double generateEvent(double helicity);
  double integrationVolume() { return 4 * sime->getDphi() * dcte; }
};

class generateBremsstrahlung : public eventGenerator {
  // for a test
  PseudoRandom prndm;
protected:
  double BHmin, BHmax, RadCutOff;
public:
  virtual ~generateBremsstrahlung() { ; }
  generateBremsstrahlung(reaction *r, simDetectorBase *Sime, 
			 simDetectorBase *Sim1, SIMUL *rundb) { 
    Label      = "d[W]'";
    Unit       = "[m]b";
    Reaction   = r; sime=Sime; sim1=Sim1;
    targetmass = Reaction->getTarget()->getMass();
    dcte       = sin(sime->getDtheta());
    dphie      = calc_dphie();
    Ex         = rundb->excitation;
    Gamma      = 1;
    BHmin      = rundb->BHmin;
    BHmax      = rundb->BHmax;
    RadCutOff  = rundb->RadiationCutOff;
  };
  
  double generateEvent(double helicity);
  double integrationVolume() { return 4 * sime->getDphi() * dcte; }
};

class generateElasticProton : public eventGenerator {
public:
  virtual ~generateElasticProton() { ; }
  generateElasticProton(reaction *r, simDetectorBase *Sime, 
			simDetectorBase *Sim1, SIMUL *rundb) { 
    Label      = "d[W]'";
    Unit       = "[m]b";
    Reaction   = r; sime=Sime; sim1=Sim1;
    targetmass = Reaction->getTarget()->getMass();
    dcte       = sin(sime->getDtheta());
    dphie      = calc_dphie();
    Gamma      = 1;
  };
  
  double generateEvent(double helicity);
  double integrationVolume() { return 4 * sime->getDphi() * dcte; }
};

class generateTwoBodyPWIA : public eventGenerator {
  class PWIA *pwia;
public:
  virtual ~generateTwoBodyPWIA() { ; }
  generateTwoBodyPWIA(reaction *r, simDetectorBase *Sime, 
	       simDetectorBase *Sim1, SIMUL *rundb, modeltype ModelType) { 
    Label      = "dE'd[W]'d[W]_1!";
    Unit       = "GeV sr^2!";
    Reaction   = r; sime=Sime; sim1=Sim1;
    targetmass = Reaction->getTarget()->getMass();
    dcte       = sin(sime->getDtheta());
    dphie      = calc_dphie();
    Gamma      = 1;
    pwia       = NULL;
    if (ModelType == TwoBodyHe3) {
      pwia = new PWIA_He3_2Body();
      Label      = "d[s]";
      Unit       = "[m]barn";
    }
  };
  
  double generateEvent(double helicity);
  double integrationVolume() { return 1; //volume hidden in generate;
  }
};



class generateIsotropic6D : public eventGenerator {
public:
  virtual ~generateIsotropic6D() { ; }
  generateIsotropic6D(reaction *r, simDetectorBase *Sime, 
		      simDetectorBase *Sim1, SIMUL *rundb,
		      modeltype ModelType) { 
    Reaction   = r;
    sime       = Sime;
    sim1       = Sim1;
    targetmass = Reaction->getTarget()->getMass();
    dcte       = sin(sime->getDtheta());
    dphie      = calc_dphie();
    Gamma      = 1;
    Label      = "d[s]";
    Unit       = "[m]barn";
  };
  
  double generateEvent(double helicity);
  double integrationVolume() {
    return(1);                    // volume hidden in generate;
  }
};

class generateDeuteronBreakup : public eventGenerator {
public:
  virtual ~generateDeuteronBreakup() { ; }
  generateDeuteronBreakup(reaction *r, simDetectorBase *Sime, 
			  simDetectorBase *Sim1, SIMUL *rundb,
			  modeltype ModelType) { 
    Reaction   = r;
    sime       = Sime;
    sim1       = Sim1;
    targetmass = Reaction->getTarget()->getMass();
    dcte       = sin(sime->getDtheta());
    dphie      = calc_dphie();
    Gamma      = 1;
    Label      = "d[s]";
    Unit       = "[m]barn";
  };
  
  double generateEvent(double helicity);
  double integrationVolume() {
    return(1);                    // volume hidden in generate;
  }
};


class generateThreeBodyPWIA : public eventGenerator {
  class PWIA *pwia;
public:
  virtual ~generateThreeBodyPWIA() { ; }
  generateThreeBodyPWIA(reaction *r, simDetectorBase *Sime, 
			simDetectorBase *Sim1, SIMUL *rundb, 
			modeltype ModelType) { 
    Label      = "dE'd[W]'dE_1!d[W]_1!";
    Unit       = "GeV^2! sr^2!";
    Reaction   = r; sime=Sime; sim1=Sim1;
    targetmass = Reaction->getTarget()->getMass();
    dcte       = sin(sime->getDtheta());
    dphie      = calc_dphie();
    Gamma      = 1;
    pwia       = NULL;
    if (ModelType == ThreeBodyHe3) {
      pwia = new PWIA_He3_3Body();
      Label      = "d[s]";
      Unit       = "[m]barn";
    }
  };
  
  double generateEvent(double helicity);
  double integrationVolume() { return 1; //volume hidden in generate;
  }
};

class generateTwoBody : public eventGenerator {
private:
  double coscmsmin;
  double coscmsmax;
  double phicmsmin;
  double phicmsmax;
  double solidangleCMS;

  int scale;
  class model *m;
  modeltype ModelType;
public:
  virtual ~generateTwoBody() { ; }
  generateTwoBody(reaction *r, simDetectorBase *Sime, 
		  simDetectorBase *Sim1, simDetectorBase *Sim2,
		  SIMUL *rundb, modeltype ModelType, int scaling=0) 
  { 
    generateTwoBody::ModelType = ModelType;

    if (scaling != 0) scale=scaling; else scale = 0;

    const double parameter[7] = {1.8, -1.5, 15.0, -15.0, 19.0, 4, -.5};
    m = NULL;
    if (ModelType == Pi0Threshold) m = new SPwaves(parameter);   
    if (ModelType == ChPTh)        m = new ChPT(); 
    if (ModelType == EtaMaid)      m = new etaMaid("maid.dat"); 
    if (ModelType == Maid) {
      if (*r->getM1()==P_proton && *r->getM2()==P_pi0)     m = new maid2000(1);
      if (*r->getM1()==P_piplus && *r->getM2()==P_neutron) m = new maid2000(3);
      if (!m) std::cerr << "Wrong reaction "<<r->getName()<<" for MAID model\n";
    }
    if (ModelType == kMaid) {
      if (*r->getM1()==P_kplus && *r->getM2()==P_Lambda)  m = new kmaid(1);
      if (*r->getM1()==P_kplus && *r->getM2()==P_Sigma0)  m = new kmaid(3);
      if (!m) std::cerr << "Wrong reaction "<<r->getName()<<" for K-MAID model\n";
    }
    
    if (m) {
      Label = "d[s]";
      Unit  = "[m]b";
    } else {
      Label = "d[W]_2!^*!";
      Unit  = "sr";
    }

    Reaction   = r; sime=Sime; sim1=Sim1; sim2=Sim2;
    targetmass = Reaction->getTarget()->getMass();
    dcte       = sin(sime->getDtheta());
    dphie      = calc_dphie();
    coscmsmin  = cos(rundb->CMSTheta[1] * rad);
    coscmsmax  = cos(rundb->CMSTheta[0] * rad);
    phicmsmin  = rundb->CMSPhi[0] * rad;
    phicmsmax  = rundb->CMSPhi[1] * rad;
    solidangleCMS = (coscmsmax-coscmsmin) * (phicmsmax-phicmsmin);
  };
  double generateEvent(double helicity);
  double integrationVolume() { return solidangleCMS; }
};

class generateDMProton : public generateTwoBody {
private:
  double mA;
public:
  generateDMProton(reaction *r, simDetectorBase *Sime, 
		   simDetectorBase *Sim1, simDetectorBase *Sim2,
		   SIMUL *rundb) 
    : generateTwoBody(r, Sime, Sim1, Sim2,rundb, Isotropic) 
  { Label = "d[s]";
    Unit  = "[m]b";
    mA = rundb->massA/1000;
    r->m2 = new Particle("A'",0,mA);
    r->Out2 = *(r->m2);//Particle("A'",0,mA);
    r->threshold=m_proton+mA;
  };
  double generateEvent(double helicity) {
    double res = generateTwoBody::generateEvent(helicity);
   
    /* Momentum P_Spin, P_Spin_CM;    double Msqr,k=-1,kp=0;
       FourVector qout(Reaction->Out2.momentum(),Reaction->Out2);
       double cs1=ElasticCrossSection(Reaction->electronIn,
       Reaction->electronOut, qout,1,0,0, &P_Spin, &P_Spin_CM, &Msqr, &k, kp)/k;
       cout <<cs1<< " " <<DMBosonCrossSection(Reaction->electronIn, 
       Reaction->electronOut, qout, 0.0)/cs1<<endl; */
    // cout <<"DMProton: "<< endl;
    // cout <<Reaction->electronIn<<endl;
    // cout <<Reaction->electronOut<<endl;
    // cout <<Reaction->Out2<<endl;
    return res / Gamma * 
      DMBosonCrossSection(Reaction->electronIn, Reaction->electronOut,
			  Reaction->Out2, mA);
  };
};

class generateTripleCMS : public eventGenerator {
private:
  double coscmsmin;
  double coscmsmax;
  double phicmsmin;
  double phicmsmax;
  double decay_m_min;
  double decay_m_max;
  double cosDmin;
  double cosDmax;
  double phiDmin;
  double phiDmax;
public:
  virtual ~generateTripleCMS() { ; }
  generateTripleCMS(reaction *r, simDetectorBase *Sime, 
		    simDetectorBase *Sim1, 
		    simDetectorBase *SimD1, simDetectorBase *SimD2,
		    SIMUL *rundb) 
  { 
    Label = "dm_12!d[W]_12!d[W]_d!";
    Unit  = "GeV/c^2! sr^2!";
    Reaction   = r;
    sime = Sime; sim1 = Sim1; simD1 = SimD1; simD2 = SimD2;
    targetmass  = Reaction->getTarget()->getMass();
    dcte        = sin(sime->getDtheta());
    dphie       = calc_dphie();
    coscmsmin   = cos(rundb->CMSTheta[1] * rad);
    coscmsmax   = cos(rundb->CMSTheta[0] * rad);
    phicmsmin   = rundb->CMSPhi[0] * rad;
    phicmsmax   = rundb->CMSPhi[1] * rad;
    decay_m_min = rundb->DecayMass[0];
    decay_m_max = rundb->DecayMass[1];
    cosDmin     = cos(rundb->DecayTheta[1] * rad);
    cosDmax     = cos(rundb->DecayTheta[0] * rad);
    phiDmin     = rundb->DecayPhi[0] * rad;
    phiDmax     = rundb->DecayPhi[1] * rad; 
  };
  double generateEvent(double helicity);
  double integrationVolume() {
    return 1; // variable volume hidden in generateEvent (dirty!)
  }
};

class generateResonance : public eventGenerator {
private:
  double coscmsmin;
  double coscmsmax;
  double phicmsmin;
  double phicmsmax;
  double decay_m_min;
  double decay_m_max;
  double cosDmin;
  double cosDmax;
  double phiDmin;
  double phiDmax;
  double solidangleCMS;
  double solidangleD;
  modeltype ModelType;
public:
  virtual ~generateResonance() { ; }
  generateResonance(reaction *r, simDetectorBase *Sime, 
		    simDetectorBase *Sim1,
		    simDetectorBase *SimD1, simDetectorBase *SimD2,
		    SIMUL *rundb, modeltype Model) 
  { 
    Label = "d[W]"; //???
    Unit  = "sr";   //???
    sime = Sime; sim1 = Sim1; simD1 = SimD1; simD2 = SimD2;
    Reaction    = r;
    ModelType   = Model;
    targetmass  = Reaction->getTarget()->getMass();
    dcte        = sin(sime->getDtheta());
    dphie       = calc_dphie();
    coscmsmin   = cos(rundb->CMSTheta[1] * rad);
    coscmsmax   = cos(rundb->CMSTheta[0] * rad);
    phicmsmin   = rundb->CMSPhi[0] * rad;
    phicmsmax   = rundb->CMSPhi[1] * rad;
    decay_m_min = rundb->DecayMass[0];
    decay_m_max = rundb->DecayMass[1];
    cosDmin     = cos(rundb->DecayTheta[1] * rad);
    cosDmax     = cos(rundb->DecayTheta[0] * rad);
    phiDmin     = rundb->DecayPhi[0] * rad;
    phiDmax     = rundb->DecayPhi[1] * rad;
    solidangleCMS = (coscmsmax-coscmsmin) * (phicmsmax-phicmsmin);
    solidangleD = (cosDmax-cosDmin) * (phiDmax-phiDmin);
  };
  double generateEvent(double helicity);
  double integrationVolume() {return solidangleCMS * solidangleD;}
  //normalization probably wrong! One dimension of volume is missing...
};

class He3eepnBase {
private:
  double  AA, BB, CC, Eb, thee, al;
  int ps;
  int t1mode;
public:
  He3eepnBase();
  void Photon_rotation(FourVector* Photon, FourVector* Proton,FourVector* Neutron);
  double T1max_nr(double q, double omega_, double thg,
		  double th1, double ph1, double th2, double ph2);
  void setup(double Ebeam, double omega, double q, double the1, double phi1, double the2, double phi2);
  double calculate(double p1, double targetmass, double omega, double q);
  double bisec(double start1, double start2, int trace, double targetmass, double omega, double q);
  double getWeight(FourVector EIn, FourVector EOut,
		   FourVector Proton, FourVector Neutron,
		   double targetmass, gridinterpolation* CS,int match,double T1rel);
};

class generateTriple : public eventGenerator, He3eepnBase {
private:
  double p1_0, p2_0, dcte, dphie, dct1, dphi1, dp1, dct2, dphi2, dp2;
  modeltype ModelType;
  gridinterpolation *CrossSection;
  int phasespace;
  double Eprime(double E0, double m0, double mm, double the, double phe,
		double m1, double p1, double th1, double ph1, 
		double m2, double p2, double th2, double ph2);
  double Jacobian(double E0, double m0, double mm, double the, double phe, 
		  double m1, double p1, double th1, double ph1, 
		  double m2, double p2, double th2, double ph2);
  void   Transform(double in[8], double out[8]);
  double JacobiDet(double the, double phie,
		   double p1, double th1, double phi1,
		   double p2, double th2, double phi2);
public:
  virtual ~generateTriple() { ; }
  generateTriple(reaction *r, simDetectorBase *Sime, 
		 simDetectorBase *SimD1, simDetectorBase *SimD2,
		 SIMUL *rundb, modeltype Model);
  double generateEvent(double helicity);

  double integrationVolume() { return 8 * dct1 * dphi1 * dp1 // d\Omega_1 dp_1
                                    * 8 * dct2 * dphi2 * dp2 // d\Omega_2 dp_2
                                    * 4 * dcte * dphie;      // d\Omega_e;
  };
};




class generateHe3eepn : public eventGenerator,  He3eepnBase{
private:
  double domega, dq, omega, q, omega_0, q_0,
    p1_0, dp1, dct1, dphi1, dct2, dphi2, dp2,
    AA, BB, CC, Eb, thee, al;
  int phasespace;
  int counter;
  double m;
  double Ebind;
  modeltype ModelType;
  gridinterpolation *CrossSection;
  FourVector Target;
  bool makefile;
  bool readfile;
  double th1min, th1max, th1stepsize, th1steps, 
         ph1min, ph1max, ph1stepsize, ph1steps, 
         th2min, th2max, th2stepsize, th2steps, 
         ph2min, ph2max, ph2stepsize, ph2steps, 
         gridT1min, T1stepsize, T1steps,
         gridomega, gridq;
//   double T1max_nr(double q, double omega_, double thg,
//                 double th1, double ph1, double th2, double ph2);
  double ftest(double w, double p1, double p2, double p3);
  double fp1guess(double C, double D, double E, double p3);
  double fp2root(double C, double D, double E, double p1, double p3);
  int tstsort(double *p1, double *tst, int n, int usefabs);
  double fp2(double A, double B, double SQRT, double p1);
  double fp2test(double w, double A, double B, double SQRT, double p1, double p3);
  double cap(double x);
  double p3toT1(double q, double w, double thg,
		double th1, double ph1,
		double th2, double ph2,
		double p3, double T1guess, double p2guess);
  //  void setup(double the1, double phi1, double the2, double phi2);
  //  double calculate(double p1);
  //  double bisec(double start1, double start2, int trace);




public:
  generateHe3eepn(reaction *r, simDetectorBase *Sime, 
		  simDetectorBase *SimD1, simDetectorBase *SimD2,
		  SIMUL *rundb, modeltype Model);
  double generateEvent(double helicity);
  double integrationVolume() {
    return 8 * domega * dq * dphie
      * 4 * dct1 * dphi1 * 2 *dp1
    * 4 * dct2 * dphi2;
  };
  double debugOutput();
  virtual ~generateHe3eepn() {if (!phasespace) delete CrossSection;}
};

class generateHe3nr : public eventGenerator{
private:
  double omega_0,domega;
  double domdq;
  double q_0,dq;
  double dct1, dphi1, dct2, dphi2;
  double Ebind;
  double M;
  double dt;
  bool dodE;
  double consta,cpx,cpy,cpz,calpha,salpha,qx,qy,qz;
  void setupbisec(double t,double w,double qx_,double qy_,double qz_,
		    double p1x,double p1y,double p1z,
		    double p2x,double p2y,double p2z);
  double bisecfunc(double l);
public:
  generateHe3nr(reaction *r, simDetectorBase *Sime, 
		  simDetectorBase *SimD1, simDetectorBase *SimD2,
		  SIMUL *rundb, modeltype Model);
  double generateEvent(double helicity);
  double integrationVolume() {
    return  domdq 
      * 2 * dphie
      * 4 * dct1 * dphi1 *2 * M_PI //* dp1 *2
    * 4 * dct2 * dphi2;

  };
  double debugOutput();
  double fixnrenergy(double mass,double p);
  double bisec(double a,double b);
  virtual ~generateHe3nr() {}
};


class generateHe3fast : public eventGenerator,  He3eepnBase{
private:
  double domega, dq, omega, q, omega_0, q_0,
    p1_0, dp1, dct1, dphi1, dct2, dphi2, dp2;
  int phasespace;
  double m;
  double Ebind;
  modeltype ModelType;
  gridinterpolation *CrossSection;
  FourVector Target;
  double         gridomega, gridq;
  int energycode;
  int jacob;

public:
  generateHe3fast(reaction *r, simDetectorBase *Sime, 
		  simDetectorBase *SimD1, simDetectorBase *SimD2,
		  SIMUL *rundb, modeltype Model);
  double generateEvent(double helicity);
  double integrationVolume() {
    return 8 * domega * dq * dphie
      * 4 * dct1 * dphi1  //* dp *2
    * 4 * dct2 * dphi2;
  };
  double debugOutput();
  double fixnrenergy(double mass,double p);
  virtual ~generateHe3fast() {if (!phasespace) delete CrossSection;}
};


class generateBetheHeitlerPeak {
  double a, Norm, a1, a2, norm1, norm2, ct, E, Egen, fweight;
  FourVector in, out;
  PseudoRandom prndm;

  double f(double ct) {  return (1-ct*ct)/(a-ct)/(a-ct)/Norm; };

  double F(double ct) {
    return ((1-a*a)/(a-ct)-ct-2*a*log((a-ct)/(a+1))-2+a)/Norm;
  };
  double transform(double y) {
    double fx, x, x1=-1, x2=1;
    int i=0;
    while((x2-x1)>2e-8) {
      fx = F(x = (x1+x2)/2);
      if (fx<y) x1=x; else x2=x;
      i++;
    }
    // for(int i=0;i<100;i++) x -= (F(x)-y)/f(x); // Newton formula to improve
    return x;
  };

public:
  generateBetheHeitlerPeak(FourVector ein, FourVector eout);
  double BetheHeitlerCrossSection(FourVector in, FourVector out, FourVector k);
  class FourVector generate (double &weight,bool *lowerlimit, double *EBH,double *helicity = 0,
		      class Momentum *P_Spin = NULL);
  class FourVector generate();
};

class generateInclusive : public eventGenerator {
private:
  double dpe;
public:
  generateInclusive(reaction *r, simDetectorBase *Sime, SIMUL *rundb);
  virtual ~generateInclusive() { ; }
  double generateEvent(double helicity);
  double integrationVolume();
};

#endif /* __Generate_h__ */
