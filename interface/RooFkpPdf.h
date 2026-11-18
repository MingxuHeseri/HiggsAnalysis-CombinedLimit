#ifndef HIGGSANALYSIS_COMBINEDLIMIT_ROOFKPPDF_H
#define HIGGSANALYSIS_COMBINEDLIMIT_ROOFKPPDF_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooAbsReal;

// Generalised ATLAS-style Fk background function:
//
//   (1 - m / scale)^b * (m / scale)^
//       (a0 + a1 * log(m / scale) + a2 * log(m / scale)^2)
//
// The supported orders retain the first order + 1 terms in the exponent:
// Fk0 uses b,a0; Fk1 uses b,a0,a1; and Fk2 uses b,a0,a1,a2.  The a1/a2
// proxies are retained for a single serialisable class; callers should pass
// constant zero RooAbsReals for terms above the selected order.
class RooFkpPdf : public RooAbsPdf {
public:
  RooFkpPdf() : order_(0) {}
  RooFkpPdf(const char* name, const char* title, Int_t order,
            RooAbsReal& _observable, RooAbsReal& _scale, RooAbsReal& _b,
            RooAbsReal& _a0, RooAbsReal& _a1, RooAbsReal& _a2);
  RooFkpPdf(const RooFkpPdf& other, const char* name = nullptr);
  TObject* clone(const char* newname) const override
  {
    return new RooFkpPdf(*this, newname);
  }
  ~RooFkpPdf() override = default;

  Int_t order() const { return order_; }

protected:
  RooRealProxy observable;
  RooRealProxy scale;
  RooRealProxy b;
  RooRealProxy a0;
  RooRealProxy a1;
  RooRealProxy a2;

  Double_t evaluate() const override;

private:
  Int_t order_;

  ClassDefOverride(RooFkpPdf, 1)
};

#endif
