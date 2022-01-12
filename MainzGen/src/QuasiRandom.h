#ifndef __QuasiRandom_h__
#define __QuasiRandom_h__

/** @name Random Generators.
 *
 *  These classes provides random generators for different distributions.
 * There are only quasi random numbers and no pseudo random numbers, since
 * we need this for simulation! See Numerical Recipies for the difference.
 *
 * For Sobol' Sequence see e.g. \\
 *    "Numerical Recipes in C"\\
 *    by W.H.Press,S.A.Teukolsky,W.T.Vetterling,B.P.Flannery\\
 *
 * usage:\\
 *    SobolSequence sobol(6);     *initializes a new sobol sequence in 6 dim\\
 *                                *maximum dimension = 6 at the moment \\
 *    sobol(double x[n])          *Next n-dim element of Sobol Sequence\\
 *       or\\
 *    sobol.nextValues();\\
 *    x1=sobol();\\
 *    x2=sobol();\\
 *      ...\\
 */
//@{

#include <math.h>

const int MAXDIM = 18;
const int MAXBIT = 53;

static const int    lengthUL  = sizeof(unsigned long) * 8;
static const double highvalue = pow(2, (int) (sizeof(unsigned long) * 8));

/// internal class to handle bit operations on 8 byte integers
class lint { 
private:
  unsigned long high, low;
public:
  lint() {};
  lint(unsigned long l, unsigned long h=0) { low=l; high=h; };
  lint(double d) { 
    low  = (unsigned long) fmod(d, highvalue); 
    high = (unsigned long) floor(d/highvalue); 
  };
  lint operator^(const lint &b) const {return lint(low^b.low, high^b.high); };
  lint operator&(const lint &b) const {return lint(low&b.low, high&b.high); };
  lint operator|(const lint &b) const {return lint(low|b.low, high|b.high); };
  lint operator ^= (const lint &b) { low^=b.low; high^=b.high; return *this;};
  lint operator &= (const lint &b) { low&=b.low; high&=b.high; return *this;};
  lint operator |= (const lint &b) { low|=b.low; high|=b.high; return *this;};
  lint operator <<=(const int s)   { return *this = *this << s; };
  lint operator >>= (const int s)  { return *this = *this >> s; };
  lint operator >> (const int s) const {
    if(s==0) return *this;
    return (s < lengthUL ?
	    lint((low  >> s) | (high << (lengthUL-s)), high >> s) : 
	    lint(high  >> (s-lengthUL) , 0));
  };
  lint operator << (const int s) const {
    if(s==0) return *this;
    return (s<lengthUL ?
	    lint(low << s, (high<<s) | (low >> (lengthUL-s))) :
	    lint(0,low <<(s-lengthUL)));
  };
  operator double() {
    const double highvalue = pow(2,lengthUL);
    return highvalue * high + low;
  }
  //  friend ostream& operator<<(ostream& out, lint l) {
  //  for (int i=lengthUL-1; i>=0; i--) out << ((l.high >> i) & 1);
  //  for (int i=lengthUL-1; i>=0; i--) out << ((l.low  >> i) & 1);
  //  return out;
  // }
};

/// virtual class of Random Generators
class RandomGenerator {
public:
  virtual ~RandomGenerator() { ; }
  /// returns the next number of the random sequence 
  virtual double operator()() = 0;
};

/// virtual class for uniformly distributed random generators
class Uniform : public RandomGenerator {
public:
  virtual ~Uniform() { ; }
  virtual double operator()() = 0;
};

/** Pseudo Random Generator.
 *  Random numbers between 0, 1. See Numerical recipes.
 */
class PseudoRandom : public Uniform {
public:
  PseudoRandom();
  virtual ~PseudoRandom() { ; }
  double operator()();
};
/** Sobol-Antonov-Saleev sequences.
 * This class provides Sobol-Antonov-Saleev sequences upto
 * 18 dimensions and with a period of \TEX{$2^{53}$}.
 */
class SobolSequence : public Uniform {
private:
  short dim, mdeg[MAXDIM+1];
  double in;
  lint ix[MAXDIM+1],*iu[MAXBIT+1], ip[19], iv[MAXDIM*MAXBIT+1];
  short  access;
  double *xSet;
public:
  virtual ~SobolSequence() { ; }
  /// Constructor
  SobolSequence(const int d = MAXDIM);
  /// fill the next elements in an array. Use this XOR nextValues() + Op() 
  void   operator()(double x[]);
  /// Calulate next elements and store them internal
  void   nextValues() { access = 0;  (*this)(xSet); }
  /// Return the next internal stored element
  double operator()();
  /// Skip the first (start) values of the sequence
  void   init(double start);
};
  
#endif /* __QuasiRandom_h__ */
















