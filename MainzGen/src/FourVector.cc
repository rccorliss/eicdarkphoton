//
// Created 1993-2002 by	Michael O. Distler,
//                      Harald Merkel
//			Institut fuer Kernphysik,
//			Johannes Gutenberg-Universitaet Mainz
//
// $Header: /tmp/cvsroot/Cola/FourVector/FourVector.cc,v 3.9 2002-11-26 18:11:50 distler Exp $
//
// Implementation der Vierervektorenklasse
//

#include "FourVector.h"

double
Momentum::specPhi(double angle) const
{
  Momentum part = *this;
  part.rot_theta(-angle);

  return part.specPhi()*1000;
}

double
Momentum::specTheta(double angle) const
{
  Momentum part = *this;
  part.rot_theta(-angle);

  return part.specTheta()*1000;
}

int
Momentum::specCheck(double pMin, double pMax,
		      double th0, double dth, double dph) const
{
  double p = square();
  if (p < pMin*pMin) return 0;
  if (pMax*pMax < p) return 0;

  Momentum part = *this;
  part.rot_theta(-th0);

  if (dth < fabs(part.specPhi())) return 0;
  if (dph < fabs(part.specTheta())) return 0;

  return 1;
}

int
Momentum::specCheck(double pMin, double pMax, double th0,
		      double thMin, double thMax,
		      double phMin, double phMax) const
{
  double p = square();
  if (p < pMin*pMin) return 0;
  if (pMax*pMax < p) return 0;

  Momentum part = *this;
  part.rot_theta(-th0);

  double thetaTmp = part.specPhi();
  if ((thetaTmp < thMin) || (thMax < thetaTmp)) return 0;

  double phiTmp = part.specTheta();
  if ((phiTmp < phMin) || (phMax < phiTmp)) return 0;

  return 1;
}

FourVector 
FourVector::rotate(const Momentum& direction) const
{ FourVector n = *this;  
  n.rot_phi(   -direction.phi());
  n.rot_theta( -direction.theta());
  return n;
}

FourVector 
FourVector::rotateTo(const Momentum& direction) const
{ FourVector n = *this;  
  n.rot_theta( direction.theta());
  n.rot_phi(   direction.phi());
  return n;
}


FourVector Polar(double E, double p, double theta, double phi)
{ 
  double pst = p * sin(theta);
  return FourVector(E, pst * cos(phi), pst * sin(phi), p * cos(theta));
}

FourVector rotate(FourVector a, Momentum direction) 
{ 
  return a.rotate(direction);
}

FourVector rotate(FourVector a, FourVector direction) 
{ 
  return a.rotate(direction);
}

FourVector rotateTo(FourVector a, Momentum direction) 
{ 
  return a.rotateTo(direction);
}

FourVector rotateTo(FourVector a, FourVector direction) 
{ 
  return a.rotateTo(direction);
}

void
FourVector::boost(double gamma)
{
  double p0 = E;
  double beta = sqrt(1-1/(gamma * gamma));

  E    = (p0 + p[2] * beta) * gamma;
  p[2] = (p0 * beta + p[2]) * gamma;
}

void
FourVector::boost(double gamma, double theta, double phi)
{
  rot_phi(-phi);
  rot_theta(-theta);
  boost(gamma);
  rot_theta(theta);
  rot_phi(phi);
}

FourVector
FourVector::Lorentz(const FourVector& reference) const
{ 
  double   g = reference.gamma();
  Momentum p = *this;
  Momentum beta = reference.beta();

  double beta_p = beta * p;
  double factor = g * (beta_p * g / (1 + g) + E);
  return FourVector(g * (E + beta_p), p + beta * factor);
}

FourVector 
Lorentz(FourVector a, FourVector reference) 
{ 
  return a.Lorentz(reference);
}

double
angle(const Momentum &a, const Momentum &b)
{
  double zwerg = a.abs() * b.abs();
  return zwerg == 0 ? 0 : acos(a*b/zwerg);
}

FourVector
FourVector::operator+(const double dE) const
{
  FourVector sum = *this;

  double massSqr = square();
  double p = momentum();

  if ((Ekin() < -dE) || (massSqr < 0.0) || (p == 0.0)) return sum;

  sum.E += dE;
  double factor =  ::momentum(sum.E, sqrt(massSqr>0.0 ? massSqr : 0.0))/p;
  sum.p[0] *= factor;
  sum.p[1] *= factor;
  sum.p[2] *= factor;

  return sum;
}

FourVector
FourVector::operator-(const double dE) const
{
  FourVector sum = *this;

  double massSqr = square();
  double p = momentum();

  if ((Ekin() < dE) || (massSqr < 0.0) || (p == 0.0)) return sum;

  sum.E -= dE;
  double factor = ::momentum(sum.E, sqrt(massSqr>0 ? massSqr : 0))/p;
  sum.p[0] *= factor;
  sum.p[1] *= factor;
  sum.p[2] *= factor;

  return sum;
}

FourVector&
FourVector::operator+=(const double dE)
{
  double massSqr = square();
  double mom = momentum();

  if ((Ekin() < -dE) || (massSqr < 0) || (mom == 0)) return *this;

  E += dE;
  double factor = sqrt(E*E > massSqr ? E * E - massSqr : 0)/mom;
  p[0] *= factor;
  p[1] *= factor;
  p[2] *= factor;

  return *this;
}

FourVector&
FourVector::operator-=(const double dE)
{
  double massSqr = square();
  double mom = momentum();

  if (Ekin() < dE || massSqr < 0 ) return *this;

  E -= dE;
  double factor = sqrt(E*E > massSqr ? E * E - massSqr : 0)/mom;
  p[0] *= factor;
  p[1] *= factor;
  p[2] *= factor;

  return *this;
}

std::istream& operator>>(std::istream& in, Momentum& m)
{
  return in >> m.p[0] >> m.p[1] >> m.p[2];
}

std::istream& operator>>(std::istream& in, FourVector& f)
{
  return in >> f.E >> f.p[0] >> f.p[1] >> f.p[2];
}

std::ostream& operator<<(std::ostream& out, Momentum m) 
{
  out.setf(std::ios::fixed, std::ios::floatfield);
  int p = out.precision(6);
  out << std::setw(10) << m[0]
      << std::setw(10) << m[1] << std::setw(10) << m[2];
  out.setf(std::ios::right, std::ios::adjustfield);
  out.precision(p);

  return out;
}

std::ostream& operator<<(std::ostream& out, FourVector f) 
{
  out.setf(std::ios::fixed, std::ios::floatfield);
  int p = out.precision(6);
  out << std::setw(10) << f.energy() << std::setw(10) << f[1]
      << std::setw(10) << f[2] << std::setw(10) << f[3]
      << std::setw(10) << f.square() << std::setw(10) << f.Ekin()
      << std::setw(10) << f.momentum();
  out.precision(p);
  out.setf(std::ios::right, std::ios::adjustfield);

  return out;
}

void
FourVector::print(char *comment) const
{
  std::cout << comment << '\t' << *this << std::endl;
}

double
FourVector::phaseSpace() const
{
  double Pxy = p[0] * p[0] + p[1] * p[1];
  return sqrt(Pxy * (Pxy + p[2] * p[2])) / (2 * E);
}

FourVector
particle(double m, double p, double th0, double th, double ph)
{
  double tan_th = tan(th);
  double tan_ph = tan(ph);
  double p_z = p / sqrt(1 + square(tan_th) + square(tan_ph));

  FourVector part(energy(p,m), p_z * tan_th, p_z * tan_ph, p_z);
  part.rot_theta(th0);

  return part;
}

FourVector
particle(double m, double p,
	 double scat0, double oop0, double scat, double oop)
{
  double tan_sca = tan(scat);
  double tan_oop = tan(oop);
  double p_z = p / sqrt(1 + square(tan_sca) + square(tan_oop));

  FourVector part(energy(p,m), p_z * tan_sca, p_z * tan_oop, p_z);
  part.rot_theta(M_PI_2);
  part.rot_phi(oop0);
  part.rot_theta(scat0-M_PI_2);

  return part;
}

double
epsilon(const FourVector &in, const FourVector &out)
{
  FourVector photon = in - out;
  double q2 = photon.square();

  return 1/(1 - 0.5 * q2 * (square(photon.energy()) - q2)/
	      square(in.energy() * out.energy() * sin(out.theta())));
}
