#ifndef HIGGSANALYSIS_COMBINEDLIMIT_HZGTHRESHOLDPDFS_H
#define HIGGSANALYSIS_COMBINEDLIMIT_HZGTHRESHOLDPDFS_H

#include "RooAbsPdf.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"

class RooAbsReal;
class RooArgList;

// Gaussian-smeared step written directly as its error-function CDF:
//   0.5 * (1 + erf((x - x0) / (sqrt(2) * sigma))).
// Here sigma has the same standard-deviation convention as the Gaussian in
// an explicit Step (x) Gaussian convolution.
class RooHZGErfTurnOnPdf : public RooAbsPdf {
public:
  RooHZGErfTurnOnPdf() = default;
  RooHZGErfTurnOnPdf(const char* name, const char* title, RooAbsReal& _x,
                     RooAbsReal& _x0, RooAbsReal& _sigma);
  RooHZGErfTurnOnPdf(const RooHZGErfTurnOnPdf& other,
                     const char* name = nullptr);
  TObject* clone(const char* newname) const override
  {
    return new RooHZGErfTurnOnPdf(*this, newname);
  }
  ~RooHZGErfTurnOnPdf() override = default;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName = nullptr) const override;
  Double_t analyticalIntegral(Int_t code,
                              const char* rangeName = nullptr) const override;

protected:
  RooRealProxy x;
  RooRealProxy x0;
  RooRealProxy sigma;

  Double_t evaluate() const override;

private:
  ClassDefOverride(RooHZGErfTurnOnPdf, 1)
};

// Smooth error-function turn-on multiplied by a falling exponential:
//   0.5 * (1 + erf((x - x0) / sigma)) * exp(-lambda * x)
class RooHZGErfExpPdf : public RooAbsPdf {
public:
  RooHZGErfExpPdf() = default;
  RooHZGErfExpPdf(const char* name, const char* title, RooAbsReal& _x,
                  RooAbsReal& _x0, RooAbsReal& _sigma,
                  RooAbsReal& _lambda);
  RooHZGErfExpPdf(const RooHZGErfExpPdf& other,
                  const char* name = nullptr);
  TObject* clone(const char* newname) const override
  {
    return new RooHZGErfExpPdf(*this, newname);
  }
  ~RooHZGErfExpPdf() override = default;

protected:
  RooRealProxy x;
  RooRealProxy x0;
  RooRealProxy sigma;
  RooRealProxy lambda;

  Double_t evaluate() const override;

private:
  ClassDefOverride(RooHZGErfExpPdf, 1)
};

// Thresholded power law multiplied by a falling exponential:
//   (x - x0)^alpha * exp(-lambda * x), x > x0
// and zero otherwise.
class RooHZGPowerExpPdf : public RooAbsPdf {
public:
  RooHZGPowerExpPdf() = default;
  RooHZGPowerExpPdf(const char* name, const char* title, RooAbsReal& _x,
                    RooAbsReal& _x0, RooAbsReal& _alpha,
                    RooAbsReal& _lambda);
  RooHZGPowerExpPdf(const RooHZGPowerExpPdf& other,
                    const char* name = nullptr);
  TObject* clone(const char* newname) const override
  {
    return new RooHZGPowerExpPdf(*this, newname);
  }
  ~RooHZGPowerExpPdf() override = default;

protected:
  RooRealProxy x;
  RooRealProxy x0;
  RooRealProxy alpha;
  RooRealProxy lambda;

  Double_t evaluate() const override;

private:
  ClassDefOverride(RooHZGPowerExpPdf, 1)
};

// Error-function turn-on times a Bernstein polynomial. The coefficient list
// must contain exactly three (order 2) or four (order 3) RooAbsReal objects.
// Non-negative coefficients make the Bernstein component non-negative over
// the observable range; invalid coefficient values evaluate to zero.
class RooErfBernsteinPdf : public RooAbsPdf {
public:
  RooErfBernsteinPdf() = default;
  RooErfBernsteinPdf(const char* name, const char* title, RooAbsReal& _x,
                     RooAbsReal& _x0, RooAbsReal& _sigma,
                     const RooArgList& _coefficients);
  RooErfBernsteinPdf(const RooErfBernsteinPdf& other,
                     const char* name = nullptr);
  TObject* clone(const char* newname) const override
  {
    return new RooErfBernsteinPdf(*this, newname);
  }
  ~RooErfBernsteinPdf() override = default;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName = nullptr) const override;
  Double_t analyticalIntegral(Int_t code,
                              const char* rangeName = nullptr) const override;

protected:
  RooRealProxy x;
  RooRealProxy x0;
  RooRealProxy sigma;
  RooListProxy coefficients;

  Double_t evaluate() const override;

private:
  ClassDefOverride(RooErfBernsteinPdf, 1)
};

#endif
