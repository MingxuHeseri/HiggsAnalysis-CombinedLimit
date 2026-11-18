#include "HiggsAnalysis/CombinedLimit/interface/RooPowerErfPdf.h"
#include "RooArgList.h"
#include "RooNumIntConfig.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
namespace {
double invNorm(double p,double low,double high) {
 const long double a=static_cast<long double>(p)+1;
 const long double span=logl(static_cast<long double>(high)/low);
 const long double norm=static_cast<long double>(low)*(a==0?span:expm1l(a*span)/a);
 return static_cast<double>(1/norm);
}
double turnOn(double value,double step,double sigma) {
 return std::max(.5*std::erfc(-(value-step)/(std::sqrt(2.)*sigma)),std::numeric_limits<double>::min());
}
double get(const RooListProxy& p,int i){return static_cast<const RooAbsReal&>(p[i]).getVal();}
}
ClassImp(RooPowerErfPdf)
RooPowerErfPdf::RooPowerErfPdf(const char* name,const char* title,RooAbsReal& mass,
 RooAbsReal& threshold,RooAbsReal& width,const RooArgList& p,const RooArgList& f):
 RooAbsPdf(name,title),x("x","mass",this,mass),step("step","threshold",this,threshold),
 sigma("sigma","Gaussian width",this,width),powers("powers","component powers",this),
 fractions("fractions","normalized component fractions",this) {
 if(p.size()<1||p.size()>6||f.size()!=p.size()-1)throw std::invalid_argument("unsupported power-erf mixture");
 powers.add(p);fractions.add(f);
 // Retain relative accuracy for small normalization integrals.
 auto* config = specialIntegratorConfig(true);
 config->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
 config->setEpsAbs(1e-50);
 config->setEpsRel(1e-10);
}
RooPowerErfPdf::RooPowerErfPdf(const RooPowerErfPdf& p,const char* name):
 RooAbsPdf(p,name),x("x",this,p.x),step("step",this,p.step),sigma("sigma",this,p.sigma),
 powers("powers",this,p.powers),fractions("fractions",this,p.fractions) {}
void RooPowerErfPdf::updateCache(double low, double high) const {
 const auto count = powers.size();
 const bool newRange = cachedPowers.size() != count || low != cachedLow || high != cachedHigh;
 if (cachedPowers.size() != count) {
   cachedPowers.resize(count);
   cachedInvNorm.resize(count);
 }
 for (std::size_t i = 0; i < count; ++i) {
   const double power = get(powers, i);
   if (newRange || power != cachedPowers[i]) {
     cachedPowers[i] = power;
     cachedInvNorm[i] = invNorm(power, low, high);
   }
 }
 cachedLow = low;
 cachedHigh = high;
}

Double_t RooPowerErfPdf::evaluate() const {
 const double value=x,t=step,s=sigma,low=x.min(),high=x.max();
 if(!(value>0&&low>0&&high>low&&s>0)||!std::isfinite(value)||!std::isfinite(t)||!std::isfinite(s))return 0.;
 updateCache(low, high);

 // Positive mass and lower bound make zero an unused initial key.
 constexpr std::size_t slots = 1024;
 if (pointFactors.empty()) pointFactors.resize(slots);
 std::uint64_t bits;
 static_assert(sizeof(bits) == sizeof(value));
 std::memcpy(&bits, &value, sizeof(bits));
 bits ^= bits >> 33;
 bits *= UINT64_C(0xff51afd7ed558ccd);
 auto& point = pointFactors[bits >> 54];
 const bool sameMass = point.mass == value;
 const bool sameBasis = sameMass && point.low == low;
 if (!sameBasis) {
   point.mass = value;
   point.low = low;
   point.logRatio = std::log(value / low);
 }
 if (!sameMass || point.threshold != t || point.width != s) {
   point.threshold = t;
   point.width = s;
   point.turn = turnOn(value, t, s);
 }
 double tail = 0., used = 0.;
 for (std::size_t i = 0; i < cachedPowers.size(); ++i) {
   if (!sameBasis || point.exponents[i] != cachedPowers[i]) {
     point.exponents[i] = cachedPowers[i];
     point.basis[i] = std::exp(cachedPowers[i] * point.logRatio);
   }
   const double f = i < fractions.size() ? get(fractions, i) : 1 - used;
   used += f;
   tail += f * point.basis[i] * cachedInvNorm[i];
 }
 return tail * point.turn;
}
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
void RooPowerErfPdf::doEval(RooFit::EvalContext& ctx) const {
 const auto xs=ctx.at(x),ts=ctx.at(step),ss=ctx.at(sigma);
 const size_t n=powers.size();const double low=x.min(),high=x.max();
 std::array<decltype(ctx.at(x)),6> ps,fs;std::array<double,6> norms{};
 for(size_t i=0;i<n;++i){ps[i]=ctx.at(&powers[i]);if(ps[i].size()==1)norms[i]=invNorm(ps[i][0],low,high);}
 for(size_t i=0;i<fractions.size();++i)fs[i]=ctx.at(&fractions[i]);
 auto at=[](const auto& a,size_t i){return a[a.size()==1?0:i];};auto output=ctx.output();
 for(size_t j=0;j<output.size();++j){
  double value=at(xs,j),t=at(ts,j),s=at(ss,j);
  if(!(value>0&&low>0&&high>low&&s>0)||!std::isfinite(value)||!std::isfinite(t)||!std::isfinite(s)){output[j]=0.;continue;}
  const double logx=std::log(value/low);double tail=0,used=0;
  for(size_t i=0;i<n;++i){double p=at(ps[i],j),f=i<fractions.size()?at(fs[i],j):1-used;used+=f;tail+=f*std::exp(p*logx)*(ps[i].size()==1?norms[i]:invNorm(p,low,high));}
  output[j]=tail*turnOn(value,t,s);
 }
}
#endif
