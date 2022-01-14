#include <math.h>
#include <iostream>
#include "Spinor.h"

Tensor 
Tensor::operator*(const Tensor &b) const
{
  return Tensor(m[0][0]*b.m[0][0] + m[0][1]*b.m[1][0] + m[0][2]*b.m[2][0] +m[0][3]*b.m[3][0],
		m[0][0]*b.m[0][1] + m[0][1]*b.m[1][1] + m[0][2]*b.m[2][1] +m[0][3]*b.m[3][1],
		m[0][0]*b.m[0][2] + m[0][1]*b.m[1][2] + m[0][2]*b.m[2][2] +m[0][3]*b.m[3][2],
		m[0][0]*b.m[0][3] + m[0][1]*b.m[1][3] + m[0][2]*b.m[2][3] +m[0][3]*b.m[3][3],
		m[1][0]*b.m[0][0] + m[1][1]*b.m[1][0] + m[1][2]*b.m[2][0] +m[1][3]*b.m[3][0],
		m[1][0]*b.m[0][1] + m[1][1]*b.m[1][1] + m[1][2]*b.m[2][1] +m[1][3]*b.m[3][1],
		m[1][0]*b.m[0][2] + m[1][1]*b.m[1][2] + m[1][2]*b.m[2][2] +m[1][3]*b.m[3][2],
		m[1][0]*b.m[0][3] + m[1][1]*b.m[1][3] + m[1][2]*b.m[2][3] +m[1][3]*b.m[3][3],
		m[2][0]*b.m[0][0] + m[2][1]*b.m[1][0] + m[2][2]*b.m[2][0] +m[2][3]*b.m[3][0],
		m[2][0]*b.m[0][1] + m[2][1]*b.m[1][1] + m[2][2]*b.m[2][1] +m[2][3]*b.m[3][1],
		m[2][0]*b.m[0][2] + m[2][1]*b.m[1][2] + m[2][2]*b.m[2][2] +m[2][3]*b.m[3][2],
		m[2][0]*b.m[0][3] + m[2][1]*b.m[1][3] + m[2][2]*b.m[2][3] +m[2][3]*b.m[3][3],
		m[3][0]*b.m[0][0] + m[3][1]*b.m[1][0] + m[3][2]*b.m[2][0] +m[3][3]*b.m[3][0],
		m[3][0]*b.m[0][1] + m[3][1]*b.m[1][1] + m[3][2]*b.m[2][1] +m[3][3]*b.m[3][1],
		m[3][0]*b.m[0][2] + m[3][1]*b.m[1][2] + m[3][2]*b.m[2][2] +m[3][3]*b.m[3][2],
		m[3][0]*b.m[0][3] + m[3][1]*b.m[1][3] + m[3][2]*b.m[2][3] +m[3][3]*b.m[3][3]);
}

Tensor 
Tensor::operator*(const Complex &b) const
{
  return Tensor(m[0][0]*b, m[0][1]*b, m[0][2]*b, m[0][3]*b,
		m[1][0]*b, m[1][1]*b, m[1][2]*b, m[1][3]*b,
		m[2][0]*b, m[2][1]*b, m[2][2]*b, m[2][3]*b,
		m[3][0]*b, m[3][1]*b, m[3][2]*b, m[3][3]*b);
}

Tensor 
Tensor::operator*(const long double b) const
{
  return Tensor(m[0][0]*b, m[0][1]*b, m[0][2]*b, m[0][3]*b,
		m[1][0]*b, m[1][1]*b, m[1][2]*b, m[1][3]*b,
		m[2][0]*b, m[2][1]*b, m[2][2]*b, m[2][3]*b,
		m[3][0]*b, m[3][1]*b, m[3][2]*b, m[3][3]*b);
}

Tensor 
Tensor::operator/(const Complex &b) const
{
  Tensor result;
  for(int mu=0;mu<4;mu++)
    for(int nu=0;nu<4;nu++)
      result.m[mu][nu] = m[mu][nu]/b;
  return result;
}

Tensor 
Tensor::operator+(const Tensor &b) const
{
  Tensor result;
  for(int mu=0; mu<4; mu++)
    for(int nu=0; nu<4; nu++)
      result.m[mu][nu] = m[mu][nu]+b[mu][nu];
  return result;
}

Tensor 
Tensor::operator+=(const Tensor &b)
{
  for(int mu=0; mu<4; mu++)
    for(int nu=0; nu<4; nu++) 
      m[mu][nu] += b[mu][nu];
  return *this;
}

Tensor 
Tensor::operator-(const Tensor &b) const
{
  Tensor result;
  for(int mu=0; mu<4; mu++)
    for(int nu=0; nu<4; nu++)
      result.m[mu][nu] = m[mu][nu] - b[mu][nu];
  return result;
}

Spinor 
Tensor::operator*(const Spinor &s) const
{
  Spinor result(0,0,0,0);
  for(int mu=0; mu<4; mu++)
    for(int nu=0; nu<4; nu++)
      result.s[mu] += m[mu][nu] * s[nu];
  return result;
}

Tensor 
operator*(const Complex &b, const Tensor &t)
{
  Tensor result;
  for(int mu=0; mu<4; mu++)
    for(int nu=0; nu<4; nu++)
      result.m[mu][nu] = t.m[mu][nu] * b;
  return result;
}

/*
ostream &operator<<(ostream &os, const Tensor &t)
{
  os << "(" << t.m[0][0] << t.m[0][1] << t.m[0][2] << t.m[0][3]<< ")" <<endl
     << "(" << t.m[1][0] << t.m[1][1] << t.m[1][2] << t.m[1][3]<< ")" <<endl
     << "(" << t.m[2][0] << t.m[2][1] << t.m[2][2] << t.m[2][3]<< ")" <<endl
     << "(" << t.m[3][0] << t.m[3][1] << t.m[3][2] << t.m[3][3]<< ")" <<endl
     <<endl;
  return os;
}
*/

const Tensor 

  Tensor0 ,

  ID(              1,  0,  0,  0,  // unity
		   0,  1,  0,  0,  
		   0,  0,  1,  0,  
		   0,  0,  0,  1),

  g_mu_nu(         1,  0,  0,  0,  // metric tensor
		   0, -1,  0,  0,  
		   0,  0, -1,  0,  
		   0,  0,  0, -1),
  
  gam[5] = {Tensor(  1,  0,  0,  0,  //  gamma matrices 
		       0,  1,  0,  0,  
		       0,  0, -1,  0,  
		       0,  0,  0, -1),
	      
	      Tensor(  0,  0,  0,  1,  
		       0,  0,  1,  0,  
		       0, -1,  0,  0,  
		      -1,  0,  0,  0),
	      
	      Tensor(  0,  0,  0, -i,  
		       0,  0,  i,  0,  
		       0,  i,  0,  0,  
		      -i,  0,  0,  0),
	      
	      Tensor(  0,  0,  1,  0,  
		       0,  0,  0, -1,  
		      -1,  0,  0,  0,  
		       0,  1,  0,  0),  
	      
	      Tensor(  0,  0,  1,  0,  
		       0,  0,  0,  1,  
		       1,  0,  0,  0,  
		       0,  1,  0,  0)};


const Complex   gmn[4] = { g_mu_nu[0][0], g_mu_nu[1][1], 
			   g_mu_nu[2][2], g_mu_nu[3][3] }; 

Tensor sigma[4][4] = {
  {Tensor0,    i*gam[0]*gam[1], i*gam[0]*gam[2], i*gam[0]*gam[3] },
  {i*gam[1]*gam[0],Tensor0,     i*gam[1]*gam[2], i*gam[1]*gam[3] },
  {i*gam[2]*gam[0],i*gam[2]*gam[1], Tensor0,     i*gam[2]*gam[3] },
  {i*gam[3]*gam[0],i*gam[3]*gam[1], i*gam[3]*gam[2], Tensor0     }
};

Tensor sigma_gmn[4][4] = {
  {sigma[0][0]*gmn[0], sigma[0][1]*gmn[1], sigma[0][2]*gmn[2], sigma[0][3]*gmn[3]},
  {sigma[1][0]*gmn[0], sigma[1][1]*gmn[1], sigma[1][2]*gmn[2], sigma[1][3]*gmn[3]},
  {sigma[2][0]*gmn[0], sigma[2][1]*gmn[1], sigma[2][2]*gmn[2], sigma[2][3]*gmn[3]},
  {sigma[3][0]*gmn[0], sigma[3][1]*gmn[1], sigma[3][2]*gmn[2], sigma[3][3]*gmn[3]}};

const long double re_gmn[4] = { real(gmn[0]), real(gmn[1]), 
			   real(gmn[2]), real(gmn[3]) }; 

//////////////////////////////////////////////////////////////////////////////
Spinor 
Spinor::operator*(const Tensor &t) const
{
  Spinor result;
  for(int mu=0; mu<4; mu++)
    for(int nu=0; nu<4; nu++)
      result.s[nu] += t[mu][nu] * s[mu];
  return result;
} 
     
/*
ostream &operator<<(ostream &os, const Spinor &s)
{ 
  return os << "(" << s[0] << "," << s[1] << "," << s[2] << "," << s[3] << ")";
}

*/

















