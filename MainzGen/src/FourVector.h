//								      -*-c++-*-
// Institut fuer Kernphysik,
// Johannes Gutenberg-Universitaet Mainz
//
// $Id: FourVector.h 2216 2008-06-13 21:13:47Z distler $
//
// Header file of the four vector library
//

#ifndef __FourVector__
#define __FourVector__

#include <cmath>
#include <iostream>
#include <istream>
#include <ostream>
#include <iomanip>
#include <iosfwd>

/** @name FourVector 

    @memo A relativistic Four Momentum Library
*/
//@{
/**  
  This class represents a cartesian momentum vector with 3 components.
  The standard arithmetics are defined, including some special geometric
  member functions to handle a spectrometer coordinate system.

  We use a right handed coordinate system \TEX{$\vec{z}=\vec{x}\times\vec{y}$}.
  The polar coordinates are defined in respect to the z-axis, i.e.

  \TEX{ 
  \begin{eqnarray*}
    x & = & r \cdot \sin\theta \cos\phi\\
    y & = & r \cdot \sin\theta \sin\phi\\
    z & = & r \cdot \cos\theta\\[3mm]
    r      & = & \sqrt{x^2+y^2+z^2}\\
    \theta & = & acos \frac{z}{r}\\
    \phi   & = & atan \frac{y}{x}\\
  \end{eqnarray*}
  }
*/

class Momentum {
protected:
  double p[3];
public:
  /// Constructor
  Momentum(double x=0, double y=0, double z=0) { 
    p[0] = x; p[1] = y; p[2] = z; }
  /// Construct a Momentum from an array of floats 
  Momentum(float v[])    { p[0] = v[0]; p[1] = v[1]; p[2] = v[2];}
  /// Construct a Momentum from an array of doubles 
  Momentum(double v[])   { p[0] = v[0]; p[1] = v[1]; p[2] = v[2];}
  /// Dump Momentum into an array of floats
  inline void fill(float v[])  const { v[0] = p[0]; v[1] = p[1]; v[2] = p[2];}
  /// Dump Momentum into an array of doubles
  inline void fill(double v[]) const { v[0] = p[0]; v[1] = p[1]; v[2] = p[2];}
  /// Set explicit the components of Momentum
  inline void set(double x, double y=0, double z=0) { p[0]=x; p[1]=y; p[2]=z;}
  /// Read access to the components of Momentum
  inline double operator[] (int i) const { return p[i]; }
  /// Add two vectors  \TEX{$\vec{a}+\vec{b}$}
  inline Momentum operator+(const Momentum &m) const {
    return Momentum(p[0]+m.p[0], p[1]+m.p[1], p[2]+m.p[2]); }
  /// Subtract two vectors  \TEX{$\vec{a} - \vec{b}$}
  inline Momentum operator - (const Momentum &m) const {
    return Momentum(p[0]-m.p[0], p[1]-m.p[1], p[2]-m.p[2]); }
  /// Negative vector  \TEX{$-\vec{a}$}
  inline Momentum operator - () const {
    return Momentum(-p[0], -p[1], -p[2]); }
  /// Comparison between Momenta
  inline int operator == (const Momentum &m) {
    return p[0]==m.p[0] && p[1]==m.p[1] && p[2]==m.p[2]; }
  /// Comparison between Momenta
  inline int operator != (const Momentum &m) {
    return p[0]!=m.p[0] || p[1]!=m.p[1] || p[2]!=m.p[2]; }
  /// Increment operator
  inline Momentum& operator += (const Momentum &m) {
    p[0]+=m.p[0]; p[1]+=m.p[1]; p[2]+=m.p[2]; return *this;}
  /// Decrement operator
  inline Momentum& operator -= (const Momentum &m) {
    p[0]-=m.p[0]; p[1]-=m.p[1]; p[2]-=m.p[2]; return *this;}
  /// Scalar product \TEX{$\vec{a}\cdot\vec{b}$}
  inline double operator*(const Momentum &m) const {
    return p[0]*m.p[0] + p[1]*m.p[1] + p[2]*m.p[2]; }
  /// Multiplication with scalar constant \TEX{$\vec{a}\cdot c$}
  inline Momentum operator*(double d) const {
    return Momentum(p[0] * d, p[1] * d, p[2] * d); }
  /// Multiplication with scalar constant \TEX{$c \cdot \vec{a}$}
  friend inline Momentum operator*(double d, Momentum m) {
    return m * d; }
  /// \TEX{$\vec{a} = c \cdot \vec{a}$}
  inline Momentum& operator*=(double d)  {
    p[0] *= d; p[1] *= d; p[2] *= d; return *this; }
  /// \TEX{$\vec{a} = \vec{a} / c$}
  inline Momentum& operator/=(double d)  {
    p[0] /= d; p[1] /= d; p[2] /= d; return *this; }
  /// \TEX{$\vec{a} / c$}
  inline Momentum operator/(double d) const {
    return Momentum(p[0]/d, p[1]/d, p[2]/d); }
  /// Vector product of two vectors  \TEX{$\vec{a} \times \vec{b}$}
  inline Momentum mult(Momentum &a, Momentum &b) {
    return Momentum(a[1]*b[2]-a[2]*b[1],
		    a[2]*b[0]-a[0]*b[2],
		    a[0]*b[1]-a[1]*b[0]); }
  /// Square of the absolute value 
  inline double square() const {
    return p[0] * p[0] + p[1] * p[1] + p[2] * p[2]; }
  /// Absolute Value
  inline double abs() const {
    return sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]); }
  /// Angle between third component and Momentum 
  inline double theta() const {
    double mom = abs(); return mom == 0 ? 0 : acos(p[2] / mom); }
  /// Polar angle between x component and Momentum around z-axis 
  inline double phi() const {
    return p[0] == 0 && p[1] == 0 ? 0 : atan2(p[1], p[0]); }
  /// Azimuth angle defined in the range \TEX{$\left[-\pi/2,\pi/2\right]$}
  inline double varPhi() const {
    return phi() + (p[0]>0 ? 0 : p[1]>0 ? - M_PI : M_PI); }
  /// Polar angle corresponding to vartheta \TEX{$\left[-\pi,\pi\right]$}
  inline double varTheta() const {
    return p[0]<0 ? -theta() : theta(); }
  /// Initalizes a Momentum with polar coordinates
  void initPolar(double m, double th, double ph) {
    p[0]=p[1]=m*sin(th); p[0]*=cos(ph); p[1]*=sin(ph); p[2]=m*cos(th); }

  /// Rotation z-axis 
 inline void rot_phi(double phi) {
    set(p[0] * cos(phi) - p[1] * sin(phi), p[1] * cos(phi) + p[0] * sin(phi), p[2]); }
  /// Rotation y-axis
  inline void rot_theta(double theta) {
    set(p[0] * cos(theta) + p[2] * sin(theta), p[1], p[2] * cos(theta) - p[0] * sin(theta)); }

  /// Deviation from central momentum in \TEX{\%}
  inline double specDelta(double centMom) const {
    return abs()*100/centMom-100; }
  /// Cartesian angle in plane for spectrometer coordinates
  double specPhi() const {
    return (p[0]==0 && p[2]==0 ? 0 : atan2(p[0], p[2])); }
  /// Cartesian angle out of plane for spectrometer coordinates
  double specTheta() const {
    return (p[2]==0 && p[1]==0 ? 0 : atan2(p[1], p[2])); }
  /// Spectrometer ccordinates
  double specPhi(double angle) const;
  /// Spectrometer ccordinates
  double specTheta(double angle) const;
  /// Test of spetrometer acceptance
  int specCheck(double pMin, double pMax,
		double th0, double dth, double dph) const;
  /// Test of spetrometer acceptance (long target)
  int specCheck(double pMin, double pMax, double th0,
		double thMin, double thMax, double phMin, double phMax) const;
  /// IO-Stream Inputoperator
  friend std::istream& operator>>(std::istream&, Momentum&);
  /// IO-Stream Outputoperator
  friend std::ostream& operator<<(std::ostream&, Momentum);
  /// Angle between two vectors
  friend double angle(const Momentum &, const Momentum &);
};

/** This class represents a relativistic four vector.
 *  The three momentum arithmetics are inherited form Momentum,
 *  in addition the energy is taken into account.
 *
 */
class FourVector : public Momentum {
protected:
  double E;
public:
  /// Constructor 
  FourVector(double En, Momentum m) { 
    E = En; p[0] = m[0]; p[1] = m[1]; p[2] = m[2]; }
  /// Constructor with four components
  FourVector(double En=0, double px=0, double py=0, double pz=0) {
    E=En; p[0]=px; p[1]=py; p[2]=pz; }
  /// Constructor with an array of floats
  FourVector(float v[])      { E = v[0]; p[0]=v[1]; p[1]=v[2]; p[2]=v[3]; }
  /// Constructor with an array of doubles
  FourVector(double v[])     { E = v[0]; p[0]=v[1]; p[1]=v[2]; p[2]=v[3]; }
  /// Initalizes a FourVector with polar coordinates
  void initPolar(double energy, double momentum, double theta, double phi) {
    E=energy; Momentum::initPolar(momentum, theta, phi);}
  /// Returns a FourVector, initialized by polar coordinates 
  friend FourVector Polar(double E, double p, double theta, double phi);
  /// Dump FourVector into array of floats
  void fill(float *v) const  { v[0] = E; Momentum::fill(&v[1]); }
  /// Dump FourVector into array of floats
  void fill(double *v) const { v[0] = E; Momentum::fill(&v[1]); }
  /// Set the FourVector Energy
  void setE(const double En) { E = En ;}
  /// Set the FourVector's Momentum 
  void setP(const Momentum & m) {
    p[0]=m[0]; p[1]=m[1]; p[2]=m[2];}
  /// FourVector's Momentum
  inline Momentum getP(void) {return (Momentum) *this;}
  /// Mass square
  inline double square() const      { return E*E - Momentum::square(); }
  /// Invariant mass 
  inline double mass() const        { return sqrt(E*E - Momentum::square()); }
  /// Lorentz factor \TEX{$\gamma$}
  inline double gamma() const       { return E / mass(); }
  /// Just the energy
  inline double energy() const      { return E; }
  /// Just the absolute value of the momentum
  inline double momentum() const    { return Momentum::abs(); }
  /// same as square
  inline double momentumSqr() const { return Momentum::square(); }
  /// Relativistic velocity \TEX{$\vec\beta$}
  inline Momentum beta() const      { return (Momentum) *this/E; } 
  /// Kinetic energy
  double Ekin() const        { double e2=square(); return e2>0?E-sqrt(e2):0; }
  /// Access to components [0..3]
  double operator[](int i) const { return i ? p[i-1] : E; }
  /// Phase space factor
  double phaseSpace() const;
  /// Print fourvector with comment string to cout
  void print(char *comment) const;
  /// Rotate a FourVector into a frame where the given vector defines de Z axis
  FourVector rotate(const Momentum& direction) const;
  /// Rotate a FourVector in the direction given by another vector
  FourVector rotateTo(const Momentum& direction) const;
  /// Lorentz boost in z-Direction. Obsolete, use method Lorentz()
  void boost(double gamma);
  /// Lorentz boost. Obsolete, use method Lorentz
  void boost(double gamma, double theta, double phi);
  /// Lorentz transformation
  FourVector Lorentz(const FourVector&) const;
  /// Comparison
  inline int operator == (const FourVector &b) const {
    return E==b.E && (Momentum) *this == (Momentum) b; }
  ///Comparison
  inline int operator != (const FourVector &b) const {
    return E!=b.E || (Momentum) *this != (Momentum) b; }
  /// Increment operator
  inline FourVector& operator += (const FourVector &b) {
    E+= b.E; p[0]+=b.p[0]; p[1]+=b.p[1]; p[2]+=b.p[2]; return *this; }
  /// Decrement operator
  inline FourVector& operator -= (const FourVector &b) {
    E-= b.E; p[0]-=b.p[0]; p[1]-=b.p[1]; p[2]-=b.p[2]; return *this; }
  /// a+b
  inline FourVector operator + (const FourVector &b) const {
    return FourVector(E + b.E, p[0]+b.p[0], p[1]+b.p[1], p[2]+b.p[2]); }
  /// a-b
  inline FourVector operator - (const FourVector &b) const {
    return FourVector(E - b.E, p[0]-b.p[0], p[1]-b.p[1], p[2]-b.p[2]); }
  /// -a, only the direction is inverted!
  inline FourVector operator - () const {
    return FourVector(E, -p[0], -p[1], -p[2]); } 
  /// gauge invariant product a*b
  inline double     operator * (const FourVector &b) const {
    return E*b.E - p[0]*b.p[0] - p[1]*b.p[1] - p[2]*b.p[2]; }
  /// scalar multiplication
  friend inline FourVector operator * (double d, FourVector f)  {
    return f * d; }
  /// scalar multiplication
  inline FourVector operator * (double d) const {
    return FourVector(E*d, p[0]*d, p[1]*d, p[2]*d); }
  /// scalar division
  inline FourVector operator / (double d) const {
    return FourVector(E/d, p[0]/d, p[1]/d, p[2]/d); }
  /// Increment operator for energylosscorrection
  FourVector& operator+=(const double dE);
  /// Decrement operator for energylosscorrection
  FourVector& operator-=(const double dE);
  /// Addition operator for energylosscorrection 
  FourVector  operator+(const double) const;
  /// Subtraction operator for energylosscorrection
  FourVector  operator-(const double) const;  
  /// IOStream input operator
  friend std::istream& operator>>(std::istream&, FourVector&);
  /// IOStream output operator
  friend std::ostream& operator<<(std::ostream&, FourVector);
  /// Initalize FourVector from spectrometer coordinates
  friend FourVector particle(double m, double p,
			     double th0, double th, double ph);
  friend FourVector particle(double m, double p, double scat0,
			     double oop0, double scat, double oop);
  /// Lorentz transformation with FV reference 
  friend FourVector Lorentz(FourVector a, FourVector reference);
  /// Rotate a FourVector
  friend FourVector rotate(FourVector a, Momentum b);
  /// Rotate a FourVector
  friend FourVector rotate(FourVector a, FourVector b);
  /// Rotate a FourVector
  friend FourVector rotateTo(FourVector a, Momentum b);
  /// Rotate a FourVector
  friend FourVector rotateTo(FourVector a, FourVector b);
  /** Virtual photon polarization 
    @return \TEX{$\epsilon$}
    @param in incoming electron
    @param out scattered electron
  */
  friend double epsilon(const FourVector &in, const FourVector &out);
};

class RNG;
class Landau;
class Multiple;
class Normal;
class Uniform;

/**
  This class is not used by Cola, but by former "plain Cindy++" programs
 */
class SimFourVector : public FourVector {
  Landau *dEnergy;
  Normal *dAngle;
  static Uniform *dPhi;
public:
  ///
  SimFourVector() { dEnergy = NULL; dAngle = NULL; }
  ///
  SimFourVector(double dE, double dA, RNG *gen);
  ///
  ~SimFourVector();
  ///
  void init(FourVector v);
  ///
  void initPolar(double E, double P, double theta, double phi);
};

/**
  This class is not used by Cola, but by former "plain Cindy++" programs
 */
class SimulFourVector : public FourVector {
  static Landau *dEnergy;
  static Normal *dAngle;
  static Uniform *dPhi;
public:
  SimulFourVector(RNG *gen);
  void init(FourVector v, double dE, double dA);
  void initPolar(double E, double P, double theta,
		 double phi, double dE, double dA);
  void scatter(double dE, double dA);
};

/**
  This class is not used by Cola, but by former "plain Cindy++" programs
 */
class FourVectorEloss : public FourVector {
    static Landau   *dEnergy;
    static Multiple *dAngle;
    static Uniform  *dPhi;
public:
  ///
  FourVectorEloss(RNG *gen);
  ///
  void init(FourVector v, double dA, double dEmp, double dEmean);
  ///
  void initPolar(double E, double P, double theta, double phi,
		 double dA, double dEmp, double dEmean);
  ///
  void scatter(double dA, double dEmp, double dEmean);
};

inline double square(double x) {  return x*x; }
inline double cubic(double x)  {  return x*x*x;}
FourVector Polar(double E, double p, double theta, double phi);
FourVector particle(double m, double p,
		    double th0, double th, double ph);
FourVector particle(double m, double p, double scat0,
		    double oop0, double scat, double oop);
FourVector Lorentz(FourVector a, FourVector reference);
FourVector rotate(FourVector a, Momentum b);
FourVector rotate(FourVector a, FourVector b);
FourVector rotateTo(FourVector a, Momentum b);
FourVector rotateTo(FourVector a, FourVector b);
double epsilon(const FourVector &in, const FourVector &out);

/// relativistic momentum
inline double
momentum(double energy, double mass) 
{
  return sqrt(energy*energy-mass*mass);
}

/// relativistic energy
inline double
energy(double momentum, double mass) 
{
  return sqrt(momentum*momentum+mass*mass);
}

inline double
lambda(double x, double y, double z)
{
  return square(x-y-z)-4.0*y*z;
}

/// velocity \TEX{$|\beta| = |p|/E$}
inline double
beta(double momentum, double mass) 
{
  return momentum/sqrt(momentum*momentum+mass*mass);
}

/// \TEX{$|\beta|^2$}
inline double
betaSquare(double momentum, double mass) 
{
  return momentum*momentum/(momentum*momentum+mass*mass);
}

/// Lorentz gamma \TEX{$\gamma$}
inline double
gamma(double momentum, double mass) 
{
  return sqrt((momentum*momentum)/(mass*mass) + 1);
}

/// Lorentz gamma \TEX{$\gamma^2$}
inline double
gammaSquare(double momentum, double mass) 
{
  return (momentum*momentum)/(mass*mass) + 1;
}
//@}
#endif

