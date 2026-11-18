#ifndef COMBINEDLIMIT_HZG_EXP_POLY_ERF_MATH_H
#define COMBINEDLIMIT_HZG_EXP_POLY_ERF_MATH_H
#include <cmath>
#include <limits>

namespace hzg {
// c0=0 in both native factories. Accumulate in extended precision to reduce
// cancellation; the public PDF still returns double and uses the same exp/CDF.
inline double expPolyExponent(int mode, int n, double value, double low,
                              double high, const double* c) {
 const long double u=(static_cast<long double>(value)-low)/(static_cast<long double>(high)-low);
 if(mode==0) {
  long double b[7]={0.};
  for(int i=1;i<=n;++i)b[i]=c[i-1];
  // de Casteljau is stable throughout the construction-time Bernstein domain.
  for(int k=n;k>0;--k)for(int i=0;i<k;++i)b[i]=(1-u)*b[i]+u*b[i+1];
  return static_cast<double>(b[0]);
 }
 const long double z=2*u-1;
 if(z < -1 || z > 1)return std::numeric_limits<double>::quiet_NaN();
 long double b1=0,b2=0;
 for(int i=n;i>=1;--i) {
  const long double b0=2*z*b1-b2+c[i-1];b2=b1;b1=b0;
 }
 return static_cast<double>(z*b1-b2);
}
}
#endif
