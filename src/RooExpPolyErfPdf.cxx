#include "HiggsAnalysis/CombinedLimit/interface/RooExpPolyErfPdf.h"
#include "RooArgList.h"
#include "HZGExpPolyErfMath.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

ClassImp(RooExpPolyErfPdf)
RooExpPolyErfPdf::RooExpPolyErfPdf(const char* name,const char* title,
 RooAbsReal& mass,RooAbsReal& threshold,RooAbsReal& width,const RooArgList& pars,
 double lower,double upper,int basis):
 RooAbsPdf(name,title),x("x","mass",this,mass),step("step","threshold",this,threshold),
 sigma("sigma","Gaussian width",this,width),coefficients("coefficients","exponent",this),
 low(lower),high(upper),mode(basis) {
 const int n=pars.size();
 if(!std::isfinite(low)||!std::isfinite(high)||!(high>low)||n<1||
    (mode==0&&n>6)||(mode==1&&n>8)||(mode!=0&&mode!=1))
  throw std::invalid_argument("unsupported exponential-polynomial erf specification");
 coefficients.add(pars);
}
RooExpPolyErfPdf::RooExpPolyErfPdf(const RooExpPolyErfPdf& other,const char* name):
 RooAbsPdf(other,name),x("x",this,other.x),step("step",this,other.step),
 sigma("sigma",this,other.sigma),coefficients("coefficients",this,other.coefficients),
 low(other.low),high(other.high),mode(other.mode) {}
Double_t RooExpPolyErfPdf::evaluate() const {
 const double value=x,t=step,s=sigma;
 if(!(s>0.)||!std::isfinite(value)||!std::isfinite(t)||!std::isfinite(s))return 0.;
 const int n=coefficients.size();
 double cs[8];
 for(int i=0;i<n;++i)cs[i]=static_cast<const RooAbsReal&>(coefficients[i]).getVal();
 const double exponent=hzg::expPolyExponent(mode,n,value,low,high,cs);
 const double argument=(value-t)/(std::sqrt(2.)*s);
 const double turn=std::isfinite(argument)
     ? std::min(std::max(.5*std::erfc(-argument),std::numeric_limits<double>::min()),1.)
     : (argument>0. ? 1. : 0.);
 return std::exp(exponent)*turn;
}

// ROOT selects this method via RooFit::EvalBackend("cpu") at NLL creation.
// Parameter spans may be scalar or event-wise; use the same arithmetic as evaluate().
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
void RooExpPolyErfPdf::doEval(RooFit::EvalContext& ctx) const {
 const auto values=ctx.at(x),steps=ctx.at(step),widths=ctx.at(sigma);
 const int n=coefficients.size();
 std::array<decltype(ctx.at(x)),8> cs;
 for(int i=0;i<n;++i)cs[i]=ctx.at(&coefficients[i]);
 auto output=ctx.output();
 auto at=[](const auto& a,size_t i) {return a[a.size()==1?0:i];};
 for(size_t j=0;j<output.size();++j) {
  const double value=at(values,j),t=at(steps,j),s=at(widths,j);
  if(!(s>0.)||!std::isfinite(value)||!std::isfinite(t)||!std::isfinite(s)) {output[j]=0.;continue;}
  double coeff[8];
  for(int i=0;i<n;++i)coeff[i]=at(cs[i],j);
  const double exponent=hzg::expPolyExponent(mode,n,value,low,high,coeff);
  const double argument=(value-t)/(std::sqrt(2.)*s);
  const double turn=std::isfinite(argument)
      ? std::min(std::max(.5*std::erfc(-argument),std::numeric_limits<double>::min()),1.)
      : (argument>0. ? 1. : 0.);
  output[j]=std::exp(exponent)*turn;
 }
}
#endif
