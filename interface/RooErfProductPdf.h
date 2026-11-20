#ifndef HIGGSANALYSIS_COMBINEDLIMIT_ROO_ERF_PRODUCT_PDF_H
#define HIGGSANALYSIS_COMBINEDLIMIT_ROO_ERF_PRODUCT_PDF_H
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include <vector>
#include <limits>
class RooErfProductPdf : public RooAbsPdf {
public:
  RooErfProductPdf() = default;
  // mode=0: Bernstein ratios b1..bn (b0=1).
  // mode=1: slopes and all normalized fractions, or n-1 plus an implicit remainder.
  RooErfProductPdf(const char*,const char*,RooAbsReal&,RooAbsReal&,RooAbsReal&,
                   const RooArgList&,const RooArgList&,int);
  RooErfProductPdf(const RooErfProductPdf&,const char* = nullptr);
  TObject* clone(const char* name) const override { return new RooErfProductPdf(*this,name); }
  Int_t getAnalyticalIntegral(RooArgSet&,RooArgSet&,const char* = nullptr) const override;
  Double_t analyticalIntegral(Int_t,const char* = nullptr) const override;
protected:
  RooRealProxy x,step,sigma;
  RooListProxy pars,fractions;
  Int_t mode=0;
  void updateTailCache() const;
  mutable std::vector<double> cacheKey; //! exact tail parameters, transient
  mutable std::vector<long double> cachePower; //! derived Bernstein coefficients
  mutable std::vector<double> cacheInvNorm; //! derived exponential normalizations
  mutable double cacheLow=std::numeric_limits<double>::quiet_NaN(); //!
  mutable double cacheHigh=std::numeric_limits<double>::quiet_NaN(); //!
  Double_t evaluate() const override;
  ClassDefOverride(RooErfProductPdf,1)
};
#endif
