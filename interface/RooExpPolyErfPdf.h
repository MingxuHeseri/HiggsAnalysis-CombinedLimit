#ifndef HIGGSANALYSIS_COMBINEDLIMIT_ROO_EXP_POLY_ERF_PDF_H
#define HIGGSANALYSIS_COMBINEDLIMIT_ROO_EXP_POLY_ERF_PDF_H
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

// Complete exp(polynomial) * Gaussian-CDF PDF. RooFit supplies numerical
// normalization. Coordinates retain the native factory's construction bounds.
class RooExpPolyErfPdf : public RooAbsPdf {
public:
  RooExpPolyErfPdf() = default;
  // mode 0: Bernstein exponent c1..cn with c0=0; mode 1: Chebyshev c1..cn.
  RooExpPolyErfPdf(const char*,const char*,RooAbsReal&,RooAbsReal&,RooAbsReal&,
                  const RooArgList&,double,double,int);
  RooExpPolyErfPdf(const RooExpPolyErfPdf&,const char* = nullptr);
  TObject* clone(const char* name) const override {return new RooExpPolyErfPdf(*this,name);}
protected:
  RooRealProxy x,step,sigma;
  RooListProxy coefficients;
  Double_t low=0.,high=1.;
  Int_t mode=0;
  Double_t evaluate() const override;
 #if ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
 void doEval(RooFit::EvalContext&) const override;
#endif
  ClassDefOverride(RooExpPolyErfPdf,1)
};
#endif
