//                                                                    -*-C++-*-
 
#include <cmath>
#include <iostream>

using namespace std;

static const long double alpha         =  1/137.035989561;
static const long double hc            = 0.1973269631;       // GeV fm
static const long double mubarn        = hc * hc * 10000;    // GeV² μbarn
static const long double e             = sqrt(alpha * 4 * M_PI);
static const long double m_electron    =  0.00051099907;     //GeV/c²
static const long double m_muon        =  0.1134289267;      //GeV/c²
//static const long double m_heavytarget = 168.514725563;      //GeV/c² for  181Ta
static const long double m_heavytarget = 0.93827231;         //GeV/c² for proton

#include "FourVector.h"
#include "Spinor.h"
#include <stdlib.h>

//Tantalum 181:
//long double Z=73;
//long double A=181;

//proton:
long double Z=1;
long double A=1;

/// form factor of solid sphere with radius R
inline long double F_sphere(long double q, long double R) 
{
  return q<1e-4 ? 1 : 3/pow(q*R,2)*(sin(q*R)/(q*R)-cos(q*R));
}

/*
Calculating the A' production cross section on heavy target

(1)             (2)     	   
      	 .                     .
      	. A'                  . A'
       .                     .	  
  e   .	      e'    e	    .	e'
  ---o---o----->    ---o---o----->	   
        (   	      (  	  	   
         )	       )          
  ------o------>    --o---------->
  Z           Z'    Z           Z'
*/

// Parameter:
//   e_in:      incoming electron e
//   e_out:     outgoing electron e'
//   q_out:     outgoing A'
//   mA:        mass of A'
// Output
//   d^5 sigma / ((d Omega_e' * d E')_LAB * (d Omega A')_CMS)
//   in  mubarn / (sr^2 GeV)

long double DMHeavyCS(const FourVector &e_in,  const FourVector &e_out,
		 const FourVector &q_out, const long double mA)
{
  FourVector p_in    = FourVector(m_heavytarget, 0, 0, 0);  
  FourVector p_out   = e_in - e_out + p_in - q_out;
  FourVector q       = p_out - p_in;
  FourVector ne1     = e_in  - q_out;             // internal electron line
  FourVector ne2     = e_out + q_out;             // internal electron line

  long double q2          = q.square();                // virtual photon q^2 
  Complex i_e3_q2    =  i * e*e*e / q2;           // electron charge -
  long double M_square = 0;

  Spinor ei[2] = {Spinor(e_in,  0.5), Spinor(e_in,  -0.5)};
  Spinor eo[2] = {Spinor(e_out, 0.5), Spinor(e_out, -0.5)};
  Spinor po[2] = {Spinor(p_out, 0.5), Spinor(p_out, -0.5)};
  Spinor pi[2] = {Spinor(p_in,  0.5), Spinor(p_in,  -0.5)};
 
  Tensor t1 = (dag(ne1) + ID * m_electron) / (ne1*ne1 - m_electron*m_electron);
  Tensor t2 = (dag(ne2) + ID * m_electron) / (ne2*ne2 - m_electron*m_electron);
  Tensor ee2gmunu[16];
  Tensor gNNnu[4];

  for(int mu=0; mu<4; mu++) {
    gNNnu[mu]= Z*F_sphere(sqrt(-q.square()),1.21/0.19733*pow(A,1.0/3))* gam[mu];
    for(int nu=0; nu<4; nu++)
      ee2gmunu[nu+(mu<<2)]=gam[nu]*t1*gam[mu] + gam[mu]*t2*gam[nu];
  }

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
	       *(-g_mu_nu[mu][nu]
		 +gmn[mu]*gmn[nu]*q_out[mu]*q_out[nu]/mA/mA));
	}    
      }
  
  const long double c = mubarn/64/pow(2*M_PI, 5)/m_heavytarget; 
  long double s = (p_in + e_in - e_out).square();
  return c 
    * sqrt((s - pow(m_heavytarget + mA,2)) * (s - pow(m_heavytarget - mA,2)))/s
    * e_out.momentum()/e_in.momentum() * M_square/4;
}

/*
Calculating the four QED Background graphs:

(1)   |e+          (2)       |e+      (3)	        (4)		   
      |	  e-	             |  e-      -----o----->      ------o----->
      o--->	             o--->      e     )   e'	  e     )   e'
     (	                    (	             (    e+	        o   /e+
  e   )	      e'     e	     )	e'            o-----            |\ /
  ---o---o----->  +  ---o---o----->  +        |       +	        | X
        (   	       (  	              o----> 	        |/ \	   
         )	        )                    (    e-	        o   \e-
  ------o------>     --o---------->     Z     )	 Z'       Z    (     Z'
  Z           Z'     Z           Z'     -----o----->      ------o----->

*/

// Parameter:
//   e_in:      incoming electron e
//   e_out:     outgoing electron e'
//   q_out:     outgoing A', i.e. the sum of the leptons: A' = e+ + e'
//   mA:        mass of A',  i.e. sqrt((e+ + e')^2)
//   theta_e12:
//   phi_e12:   A' decay angles in A' rest frame
// Output
//   d⁸σ / dΩ_e' dE' dΩ_A' dΩ_decay dm_A'
//   in  mubarn / (sr^3 GeV^2)
//   dΩ_e', dE' in LAB System
//   dΩ_A'      in photon-target CM system
//   dΩ_decay   in A' rest system
//   use "export TRIDENT=0" to disable trident graphs.
//   use "export TRIDENT=1" or unset TRIDENT for all four graphs.
//   use "export TRIDENT=2" for trident graphs only.
//   use "export TRIDENT=3" for interference term only.

long double QEDBackground(const FourVector &e_in,  const FourVector &e_out,
		     const FourVector &q_out, const long double mA, 
		     const long double theta_e12, const long double phi_e12)
{


  static int muon   = getenv("MUON")    ? atoi(getenv("MUON"))    : 0;
  static int graphs = getenv("TRIDENT") ? atoi(getenv("TRIDENT")) : 1;
  static int asymm  = getenv("ANTISYM") ? atoi(getenv("ANTISYM")) : 1;
  
  long double lepton = muon ? m_muon : m_electron;

  FourVector p_in  = FourVector(m_heavytarget, 0, 0, 0);  
  FourVector p_out = e_in - e_out + p_in - q_out;
  FourVector q_p   = e_in - e_out;
  FourVector q     = p_out - p_in;

  const long double peq = sqrt((mA*mA - pow(2*lepton,2))/4);
  FourVector e1outq = Polar(energy(lepton,peq), peq, 
			    theta_e12, phi_e12).rotateTo(q_out);
  FourVector e1out = ( e1outq).Lorentz(q_out);
  FourVector e2out = (-e1outq).Lorentz(q_out);

  FourVector ne1   = e_in  - q_out;  // virtual electron lines
  FourVector ne2   = e_out + q_out;
  FourVector ne3   = q_p - e1out;    // virtual lepton lines
  FourVector ne4   = e2out - q_p;

  Tensor t1 = (dag(ne1) + ID * m_electron) / (ne1*ne1 - m_electron*m_electron);
  Tensor t2 = (dag(ne2) + ID * m_electron) / (ne2*ne2 - m_electron*m_electron);
  Tensor t3 = (dag(ne3) + ID * lepton) / (ne3*ne3 - lepton*lepton);
  Tensor t4 = (dag(ne4) + ID * lepton) / (ne4*ne4 - lepton*lepton);

  long double q2    = q.square();         // virtual photon lines
  long double qq2   = q_out.square();
  long double qe2   = q_p.square();

  Spinor ei[2] = {Spinor(e_in,  0.5),       Spinor(e_in,  -0.5)};
  Spinor eo[2] = {Spinor(e_out, 0.5),       Spinor(e_out, -0.5)};
  Spinor e1[2] = {Antiparticle(e1out, 0.5), Antiparticle(e1out,-0.5)}; 
  Spinor e2[2] = {Spinor(e2out, 0.5),       Spinor(e2out, -0.5)};

  long double F = Z * F_sphere(sqrt(-q.square()), 1.21/hc * pow(A, 1.0/3));

  long double M_square = 0;

  Tensor A12[4][4], A34[4][4];
  for(int mu=0; mu<4; mu++) 
    for(int nu=0; nu<4; nu++) {
      A12[mu][nu] = gam[nu] * t1 * gam[mu] + gam[mu] * t2 * gam[nu];
      A34[mu][nu] = gam[nu] * t3 * gam[mu] + gam[mu] * t4 * gam[nu];
    }

  // for antisymmetrization we exchanged e2out <-> e_out
  FourVector ASq_out = e1out + e_out;  
  FourVector ASq_p   = e_in  - e2out;
  FourVector ASne1   = e_in  - ASq_out;  // virtual electron lines
  FourVector ASne2   = e2out + ASq_out;
  FourVector ASne3   = ASq_p - e1out;
  FourVector ASne4   = e_out - ASq_p;

  Tensor ASt1 = (dag(ASne1) + ID * m_electron) / (ASne1*ASne1 - m_electron*m_electron);
  Tensor ASt2 = (dag(ASne2) + ID * m_electron) / (ASne2*ASne2 - m_electron*m_electron);
  Tensor ASt3 = (dag(ASne3) + ID * m_electron) / (ASne3*ASne3 - m_electron*m_electron);
  Tensor ASt4 = (dag(ASne4) + ID * m_electron) / (ASne4*ASne4 - m_electron*m_electron);
  long double ASqq2 = ASq_out.square();
  long double ASqe2 = ASq_p.square();

  Tensor ASA12[4][4], ASA34[4][4];
  for(int mu=0; mu<4; mu++) 
    for(int nu=0; nu<4; nu++) {
      ASA12[mu][nu] = gam[nu] * ASt1 * gam[mu] + gam[mu] * ASt2 * gam[nu];
      ASA34[mu][nu] = gam[nu] * ASt3 * gam[mu] + gam[mu] * ASt4 * gam[nu];
    }
 
  int spin_e_in=0;
  // for(int spin_e_in=0; spin_e_in<=1; spin_e_in++) // skipped by parity cons
  for(int spin_e_out       = 0; spin_e_out <= 1; spin_e_out++) 
    for(int spin_e1        = 0; spin_e1    <= 1; spin_e1++   )
      for(int spin_e2      = 0; spin_e2    <= 1; spin_e2++   ) {
	Complex M_if = 0, M_e1, M_e2, M_if1 = 0, M_if2 = 0;
	for(int mu=0; mu<4; mu++) 
	  for(int nu=0; nu<4; nu++) {
	    Complex el_part = 0;
	    if (graphs != 2)
	      el_part = M_e1 = (!eo[spin_e_out] * A12[mu][nu] * ei[spin_e_in])
		* gmn[mu]/qq2 *  (!e2[spin_e2] * gam[mu] * e1[spin_e1]); 
	    
	    if (graphs != 0)
	      el_part += M_e2 = (!e2[spin_e2] * A34[mu][nu] * e1[spin_e1])
		* gmn[mu]/qe2 * (!eo[spin_e_out] * gam[mu]* ei[spin_e_in]);
	    if (asymm != 0 && !muon) {
	      el_part += -(!e2[spin_e2] * ASA12[mu][nu] * ei[spin_e_in])
		* gmn[mu]/ASqq2 *  (!eo[spin_e_out] * gam[mu] * e1[spin_e1]); 
	      el_part += -(!eo[spin_e_out] * ASA34[mu][nu] * e1[spin_e1])
	      	* gmn[mu]/ASqe2 * (!e2[spin_e2] * gam[mu]* ei[spin_e_in]);
	    }
	    Complex nu_part =  F * (p_out+p_in)[nu];
	    
	    if (graphs!=3) 
	      M_if += i * pow(e,4)/q2 * el_part * gmn[nu] * nu_part;
	    else {
	      M_if1 +=  i * pow(e,4)/q2 * M_e1 * gmn[nu] * nu_part;
	      M_if2 +=  i * pow(e,4)/q2 * M_e2 * gmn[nu] * nu_part;
	    }
	  }
	if (graphs!=3) 
	  M_square += norm(M_if); 
	else
	  M_square += 2*(real(M_if1)*real(M_if2) +imag(M_if1)*imag(M_if2));
      }

  const long double c = mubarn/64/pow(2*M_PI, 8)/m_heavytarget; 
  long double s = (p_in + e_in - e_out).square();
  long double res = c
    * sqrt((s - pow(m_heavytarget + mA,2)) * (s - pow(m_heavytarget - mA,2)))/s
    * e_out.momentum()/e_in.momentum()*M_square/2 * peq;

  return res;
}

//0.000038030287518954825260776575
