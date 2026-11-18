#ifndef HIGGSANALYSIS_COMBINEDLIMIT_ROO_POWER_ERF_PDF_H
#define HIGGSANALYSIS_COMBINEDLIMIT_ROO_POWER_ERF_PDF_H
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include <array>
#include <vector>

// Complete mixture of separately normalized powers times Gaussian CDF.
// RooFit retains numerical integration for the complete product.
class RooPowerErfPdf : public RooAbsPdf {
public:
 RooPowerErfPdf()=default;
 RooPowerErfPdf(const char*,const char*,RooAbsReal&,RooAbsReal&,RooAbsReal&,
                const RooArgList&,const RooArgList&);
 RooPowerErfPdf(const RooPowerErfPdf&,const char* =nullptr);
 TObject* clone(const char* name) const override {return new RooPowerErfPdf(*this,name);}
protected:
 RooRealProxy x,step,sigma;
 RooListProxy powers,fractions;
 mutable std::vector<double> cachedPowers; //!
 mutable std::vector<double> cachedInvNorm; //!
 mutable double cachedLow=0.,cachedHigh=0.; //!
 // Exact scalar factors. The hash selects a slot; equality decides reuse.
 struct PointFactors {
   double mass = 0., low = 0., logRatio = 0.;
   std::array<double,6> exponents{}, basis{};
   double threshold = 0., width = 0., turn = 0.;
 };
 mutable std::vector<PointFactors> pointFactors; //!
 void updateCache(double low, double high) const;
 Double_t evaluate() const override;
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
 void doEval(RooFit::EvalContext&) const override;
#endif
 ClassDefOverride(RooPowerErfPdf,1)
};
#endif
