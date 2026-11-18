#include "HiggsAnalysis/CombinedLimit/interface/RooErfProductPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
using Real=long double;
const Real pi=acosl(-1.L), sqrt2=sqrtl(2.L);
Real phi(Real z) { return expl(-z*z/2)/sqrtl(2*pi); }
Real Phi(Real z) { return erfcl(-z/sqrt2)/2; }
Real cdfDifference(Real a,Real b) {
 return a>0 ? (erfcl(a/sqrt2)-erfcl(b/sqrt2))/2 : Phi(b)-Phi(a);
}
Real choose(int n,int j) { Real v=1;for(int i=1;i<=j;++i)v=v*(n-i+1)/i;return v; }
// All moments through degree n share endpoints and the Gaussian recurrence.
// The largest caller is the degree-eight small-slope expansion below.
std::array<Real,9> erfMoments(int n,Real a,Real b,Real t,Real s) {
 const Real za=(a-t)/s,zb=(b-t)/s,pa=phi(za),pb=phi(zb);
 const Real cdfa=Phi(za),cdfb=Phi(zb);
 std::array<Real,10> normal{},tp{},sp{};
 std::array<Real,9> moments{};
 normal[0]=cdfDifference(za,zb);normal[1]=pa-pb;
 for(int k=2;k<=n+1;++k)normal[k]=(k-1)*normal[k-2]+powl(za,k-1)*pa-powl(zb,k-1)*pb;
 // Retain powl and the original summation order to avoid changing rounding
 // through a new power recurrence or a rearranged polynomial expression.
 for(int k=0;k<=n+1;++k) {tp[k]=powl(t,k);sp[k]=powl(s,k);}
 for(int j=0;j<=n;++j) {
  Real integral=0;
  for(int k=0;k<=j+1;++k)integral+=choose(j+1,k)*tp[j+1-k]*sp[k]*normal[k];
  moments[j]=(powl(b,j+1)*cdfb-powl(a,j+1)*cdfa-integral)/(j+1);
 }
 return moments;
}
// Integrate exp(k*(x-anchor))*Phi((x-t)/s); use a convergent moment series near k=0.
Real expIntegral(Real lo,Real hi,Real t,Real s,Real k,Real anchor,Real span) {
 if(fabsl(k*span)<1e-3L) {
  Real sum=0,factor=1;
  const auto moments=erfMoments(8,(lo-anchor)/span,(hi-anchor)/span,(t-anchor)/span,s/span);
  for(int j=0;j<=8;++j) {
   sum+=factor*moments[j];
   factor*=k*span/(j+1);
  }
  return span*sum;
 }
 const Real za=(lo-t)/s,zb=(hi-t)/s;
 // Difference the CDFs before scaling. Separate primitives can be huge and
 // nearly equal for a steep falling exponential with a broad turn-on.
 return (expl(k*(hi-anchor))*Phi(zb)-expl(k*(lo-anchor))*Phi(za)
         -expl(k*(t-anchor)+k*k*s*s/2)*cdfDifference(za-k*s,zb-k*s))/k;
}
Real expNorm(Real k,Real span) { return k==0 ? span : (k<0 ? expm1l(k*span)/k : -expm1l(-k*span)/k); }
Real get(const RooListProxy& list,int i) { return static_cast<const RooAbsReal&>(list[i]).getVal(); }
std::vector<Real> powerCoefficients(const RooListProxy& p) {
 int n=p.size();std::vector<Real> a(n+1,0);
 for(int i=0;i<=n;++i)for(int j=i;j<=n;++j)
  a[j]+=(i==0?1:get(p,i-1))*choose(n,j)*choose(j,i)*((j-i)%2?-1:1);
 return a;
}
}
ClassImp(RooErfProductPdf)
RooErfProductPdf::RooErfProductPdf(const char* n,const char* title,RooAbsReal& v,RooAbsReal& t,RooAbsReal& s,const RooArgList& p,const RooArgList& f,int m):
 RooAbsPdf(n,title),x("x","mass",this,v),step("step","threshold",this,t),sigma("sigma","Gaussian width",this,s),pars("pars","tail parameters",this),fractions("fractions","normalized component fractions",this),mode(m) {
 pars.add(p);fractions.add(f);
 if((m==0&&(p.size()<1||p.size()>5||f.size()!=0))||(m==1&&(p.size()<1||f.size()!=p.size()-1))||(m!=0&&m!=1))throw std::invalid_argument("unsupported erf-product specification");
}
RooErfProductPdf::RooErfProductPdf(const RooErfProductPdf& o,const char* n):RooAbsPdf(o,n),x("x",this,o.x),step("step",this,o.step),sigma("sigma",this,o.sigma),pars("pars",this,o.pars),fractions("fractions",this,o.fractions),mode(o.mode) {}
void RooErfProductPdf::updateTailCache() const {
 const double L=x.min(),H=x.max();
 bool changed=cacheKey.size()!=pars.size()||L!=cacheLow||H!=cacheHigh;
 if(cacheKey.size()!=pars.size())cacheKey.resize(pars.size());
 for(size_t i=0;i<pars.size();++i) {
  const double value=static_cast<const RooAbsReal&>(pars[i]).getVal();
  if(value!=cacheKey[i])changed=true;
  cacheKey[i]=value;
 }
 if(!changed)return;
 cacheLow=L;cacheHigh=H;
 if(mode==0)cachePower=powerCoefficients(pars);
 else {
  cacheInvNorm.resize(pars.size());
  for(size_t i=0;i<pars.size();++i)cacheInvNorm[i]=static_cast<double>(1/expNorm(cacheKey[i],H-L));
 }
}
Double_t RooErfProductPdf::evaluate() const {
 const Real v=x,t=step,s=sigma,L=x.min(),H=x.max(),D=H-L;
 if(!(s>0&&D>0)||!std::isfinite(v)||!std::isfinite(t)||!std::isfinite(s))return 0.;
 updateTailCache();
 Real tail=0;
 if(mode==0) {
  const Real u=(v-L)/D;
  for(auto it=cachePower.rbegin();it!=cachePower.rend();++it)tail=tail*u+*it;
 } else {
  Real used=0;
  for(size_t i=0;i<pars.size();++i) {
   const Real f=i<fractions.size()?get(fractions,i):1-used;used+=f;
   const double k=cacheKey[i],anchor=k<=0?cacheLow:cacheHigh;
   tail+=f*std::exp(k*(static_cast<double>(v)-anchor))*cacheInvNorm[i];
  }
 }
 const double turn=std::max(0.5*std::erfc(-(static_cast<double>(v)-static_cast<double>(t))/(std::sqrt(2.)*static_cast<double>(s))),std::numeric_limits<double>::min());
 return static_cast<double>(tail*turn);
}
Int_t RooErfProductPdf::getAnalyticalIntegral(RooArgSet& all,RooArgSet& integrated,const char*) const { return matchArgs(all,integrated,x)?1:0; }
Double_t RooErfProductPdf::analyticalIntegral(Int_t code,const char* range) const {
 if(code!=1)return 0.;
 if(range&&std::string(range).find(',')!=std::string::npos) {
  std::string r(range);Real total=0;size_t begin=0;
  while(begin<r.size()) {auto end=r.find(',',begin);auto token=r.substr(begin,end-begin);if(!token.empty())total+=analyticalIntegral(1,token.c_str());if(end==std::string::npos)break;begin=end+1;}
  return static_cast<double>(total);
 }
 const Real L=x.min(),H=x.max(),D=H-L,lo=x.min(range),hi=x.max(range),t=step,s=sigma;
 if(!(D>0&&hi>lo&&s>0)||!std::isfinite(t)||!std::isfinite(s))return 0.;
 Real result=0;
 if(mode==0) {
  updateTailCache();
  const auto moments=erfMoments(pars.size(),(lo-L)/D,(hi-L)/D,(t-L)/D,s/D);
  for(size_t i=0;i<cachePower.size();++i)result+=cachePower[i]*moments[i];
  result*=D;
 } else {
  Real used=0;
  for(size_t i=0;i<pars.size();++i) {
   const Real f=i<fractions.size()?get(fractions,i):1-used;used+=f;
   const Real k=get(pars,i),anchor=k<=0?L:H;
   result+=f*expIntegral(lo,hi,t,s,k,anchor,D)/expNorm(k,D);
  }
 }
 return static_cast<double>(result);
}
