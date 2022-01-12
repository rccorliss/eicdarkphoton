//                                                                    -*-c++-*-

#ifndef __Spinor_h__
#define __Spinor_h__
#include <math.h> 
#include "FourVector.h"

class Complex {
  double Re,Im;
public:
  Complex(double R=0, double I=0) : Re(R), Im(I) {};

  inline Complex  operator*(const Complex& b) const {
    return Complex(Re*b.Re-Im*b.Im, Re*b.Im+Im*b.Re);
  };
  inline Complex  operator*(const double b) const {
    return Complex(Re*b, Im*b);
  };
  //???????????????????????
  inline friend Complex operator*(const double a, const Complex& b) {
    return Complex(a*b.Re, a*b.Im);
  };
  //???????????????????????
  inline Complex& operator*=(const Complex& b) {
    double f=Re*b.Im+Im*b.Re; Re=Re*b.Re-Im*b.Im; Im=f; return *this;
  };
  inline Complex operator/(const double b) const {
    return Complex(Re/b, Im/b);
  }; 
  inline Complex  operator/(const Complex& b) const { 
    double nenner = b.Re*b.Re + b.Im*b.Im;
    return Complex((Re*b.Re+Im*b.Im)/nenner, (-Re*b.Im+Im*b.Re)/nenner);
  };
  inline Complex& operator+=(const Complex& b) {
    Re+=b.Re; Im+=b.Im; return *this;
  };

  inline Complex operator-() const { return Complex(-Re, -Im); };
  inline Complex operator+(const Complex& b) const {
    return Complex(Re+b.Re, Im+b.Im);
  }; 
  inline friend Complex operator+(const double a, const Complex& b) {
    return Complex(a+b.Re, b.Im);
  };
  inline friend Complex operator+(const Complex& b, const double a) {
    return Complex(a+b.Re, b.Im);
  };
  inline        Complex operator-(const Complex& b) const {
    return Complex(Re-b.Re, Im-b.Im);
  };
  inline friend Complex operator-(const double a, const Complex& b) {
    return Complex(a-b.Re,-b.Im);
  };
  inline friend Complex operator-(const Complex& b, const double a) {
    return Complex(b.Re-a, b.Im);
  };
  inline friend Complex conj(const Complex& a) {
    return Complex(a.Re, -a.Im);
  }; 
  inline friend double real(const Complex& a) { return a.Re; }; 
  inline friend double imag(const Complex& a) { return a.Im; }; 
  inline friend double norm(const Complex& a) { return a.Re*a.Re+a.Im*a.Im; };
};

static const Complex i(0,1);

class Tensor;
class Spinor;

class Tensor { 
private: 
  Complex m[4][4];
public:
  Tensor(const Complex& n00=0, const Complex& n01=0,
	 const Complex& n02=0, const Complex& n03=0,
	 const Complex& n10=0, const Complex& n11=0,
	 const Complex& n12=0, const Complex& n13=0,
	 const Complex& n20=0, const Complex& n21=0,
	 const Complex& n22=0, const Complex& n23=0,
	 const Complex& n30=0, const Complex& n31=0,
	 const Complex& n32=0, const Complex& n33=0) { 
    m[0][0]=n00; m[0][1]=n01; m[0][2]=n02; m[0][3]=n03;
    m[1][0]=n10; m[1][1]=n11; m[1][2]=n12; m[1][3]=n13;
    m[2][0]=n20; m[2][1]=n21; m[2][2]=n22; m[2][3]=n23;
    m[3][0]=n30; m[3][1]=n31; m[3][2]=n32; m[3][3]=n33;
  };
  
  const Complex *operator[](const int i        ) const { return m[i]; };
  Tensor         operator+ (const Tensor &     ) const;
  Tensor         operator- (const Tensor &     ) const;
  Tensor         operator* (const Tensor &     ) const;
  Tensor         operator* (const Complex &    ) const;
  Tensor         operator* (const double       ) const;
  Tensor         operator/ (const Complex &    ) const;
  Spinor         operator* (const class Spinor&) const;
  Tensor         operator+=(const Tensor &     );

  friend Tensor dag(const FourVector &);    // Feynman dagger
  friend Tensor operator*(const Complex &, const Tensor &);
  friend std::ostream &operator<<(std::ostream &os, const Tensor &t);
};

extern Tensor sigma[4][4], sigma_gmn[4][4];
extern const Tensor ID, gam[5], g_mu_nu, Tensor0;
extern const Complex    gmn[4];
extern const double  re_gmn[4];

class Spinor { 
private:
  Complex s[4];
public:
  Spinor(const Complex& s0=0, const Complex& s1=0, 
	 const Complex& s2=0, const Complex& s3=0) { 
    s[0]=s0; s[1]=s1; s[2]=s2; s[3]=s3;
  };

  Spinor(const FourVector& a, const double pol);

  Complex operator[](const int i    ) const { return s[i];};
  Spinor  operator! (               ) const;    // that's the adjoint Spinor
  Complex operator* (const Spinor & ) const;
  Spinor  operator* (const Complex &k) const;
  Spinor  operator+ (const Spinor &a) const {
    return Spinor(s[0]+a.s[0], s[1]+a.s[1], s[2]+a.s[2], s[3]+a.s[3]);
  }
  Spinor  operator- (const Spinor &a) const {
    return Spinor(s[0]-a.s[0], s[1]-a.s[1], s[2]-a.s[2], s[3]-a.s[3]);
  }
  Spinor  operator* (const Tensor & ) const;
  Spinor  operator/ (const Complex &) const;


  friend Spinor Tensor::operator*(const Spinor &) const;
  friend std::ostream &operator<<(std::ostream &, const Spinor &);
  friend Spinor Antiparticle(const FourVector& a, const double pol);

};

inline Spinor 
Spinor::operator!() const
{
  return Spinor(conj(s[0]), conj(s[1]), -conj(s[2]), -conj(s[3]));
}

inline Tensor dag(const FourVector &pp)
{
  return Tensor(pp[0],                     0,         -pp[3], -pp[1]+i*pp[2],
                    0,                 pp[0], -pp[1]-i*pp[2],          pp[3],
                pp[3],         pp[1]-i*pp[2],         -pp[0],              0,
                pp[1]+i*pp[2],         -pp[3],             0,         -pp[0]);
}

inline Spinor 
Spinor::operator*(const Complex &k) const 
{
  return Spinor(s[0]*k, s[1]*k, s[2]*k, s[3]*k);
}

inline Complex 
Spinor::operator*(const Spinor &sp) const 
{
  return s[0] * sp.s[0] + s[1] * sp.s[1] + s[2] * sp.s[2] + s[3] * sp.s[3];
}

inline Spinor Spinor::operator/(const Complex &k) const
{
  return Spinor(s[0]/k, s[1]/k, s[2]/k, s[3]/k);
}

inline
Spinor::Spinor(const FourVector& a, const double pol) 
{                                         // at the moment polarization is
    double n = a.energy() + a.mass();     // only possible in direction of p[2]
    Momentum p(a);
    if (pol>0) *this=Spinor(1, 0, p[2]/n,   Complex(p[0]/n, p[1]/n))*sqrt(n);
    if (pol<0) *this=Spinor(0, 1, Complex(p[0]/n, -p[1]/n), -p[2]/n)*sqrt(n);
}

inline
Spinor Antiparticle(const FourVector& a, const double pol) 
{                                         // at the moment polarization is
    double n = a.energy() + a.mass();     // only possible in direction of p[2]
    Momentum p(a);
    if (pol>0) return Spinor(Complex(p[0]/n, -p[1]/n),-p[2]/n, 0, 1)*sqrt(n);
    else       return Spinor(p[2]/n, Complex(p[0]/n, p[1]/n),  1, 0)*sqrt(n);
}

#endif /** __Spinor_h__ **/













