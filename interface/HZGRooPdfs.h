/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: HZGRooPdfs.h,v 1.1 2013/02/18 11:52:16 gpetrucc Exp $
 * Authors:                                                                  *
 *   Kyle Cranmer (L. Gray for step function piece)
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_STEPBERNSTEIN
#define ROO_STEPBERNSTEIN

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;
class RooArgList ;

class RooStepBernstein : public RooAbsPdf {
public:

  RooStepBernstein() ;
  RooStepBernstein(const char *name, const char *title,
		   RooAbsReal& _x, RooAbsReal& _stepThresh,
		   const RooArgList& _coefList) ;

  RooStepBernstein(const RooStepBernstein& other, const char* name = 0);
  TObject* clone(const char* newname) const override { return new RooStepBernstein(*this, newname); }
  inline ~RooStepBernstein() override { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override ;

private:

  RooRealProxy _x;
  RooRealProxy _stepThresh;
  RooListProxy _coefList ;

  Double_t evaluate() const override;

  ClassDefOverride(RooStepBernstein,1) // Bernstein polynomial PDF with step function
};

class RooGaussStepBernstein : public RooAbsPdf {
public:

  RooGaussStepBernstein() ;
  RooGaussStepBernstein(const char *name, const char *title,
			RooAbsReal& _x, RooAbsReal& _mean,
			RooAbsReal& _sigma, RooAbsReal& _stepVal,
			const RooArgList& _coefList) ;

  RooGaussStepBernstein(const RooGaussStepBernstein& other,
const char* name = 0);
  TObject* clone(const char* newname) const override
  { return new RooGaussStepBernstein(*this, newname); }
  inline ~RooGaussStepBernstein() override { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override ;

private:

  RooRealProxy _x;
  RooRealProxy _mean,_sigma,_stepVal;
  RooListProxy _coefList ;

  Double_t evaluate() const override;

  ClassDefOverride(RooGaussStepBernstein,1) // Bernstein polynomial PDF with step function convoluted with gaussian
};

// Factorized Bernstein step kernel used with an explicit RooFFTConvPdf.
// This is the optimized counterpart of RooGaussStepBernstein; the suffix is
// intentional so both constructions can coexist in validation workspaces.
class RooGaussStepBernsteinFast : public RooAbsPdf {
public:

  RooGaussStepBernsteinFast();
  RooGaussStepBernsteinFast(const char *name, const char *title,
                            RooAbsReal& _x, RooAbsReal& _stepThresh,
                            const RooArgList& _coefList,
                            Double_t _xmax);

  RooGaussStepBernsteinFast(const RooGaussStepBernsteinFast& other,
                            const char* name = 0);
  TObject* clone(const char* newname) const override
  { return new RooGaussStepBernsteinFast(*this, newname); }
  inline ~RooGaussStepBernsteinFast() override { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override;

private:

  RooRealProxy _x;
  RooRealProxy _stepThresh;
  RooListProxy _coefList;
  Double_t _xmax;

  Double_t evaluate() const override;

  ClassDefOverride(RooGaussStepBernsteinFast,1)
};

// Compiled HZG step kernels for the Gaussian-convolved background families.
// The Gaussian convolution remains represented explicitly by RooFFTConvPdf in
// PdfModelBuilder, while these classes replace the RooGenericPdf formula used
// for the step component.  The names follow RooGaussStepBernstein and the
// public HZG model names (*StepxGau).
class RooGaussStepPowerLaw : public RooAbsPdf {
public:

  RooGaussStepPowerLaw();
  RooGaussStepPowerLaw(const char *name, const char *title,
                       RooAbsReal& _x, RooAbsReal& _stepThresh,
                       const RooArgList& _powerList,
                       const RooArgList& _coefList,
                       Double_t _xmax);

  RooGaussStepPowerLaw(const RooGaussStepPowerLaw& other,
                       const char* name = 0);
  TObject* clone(const char* newname) const override
  { return new RooGaussStepPowerLaw(*this, newname); }
  inline ~RooGaussStepPowerLaw() override { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override;

private:

  RooRealProxy _x;
  RooRealProxy _stepThresh;
  RooListProxy _powerList;
  RooListProxy _coefList;
  Double_t _xmax;

  Double_t evaluate() const override;

  ClassDefOverride(RooGaussStepPowerLaw,1)
};

class RooGaussStepExponential : public RooAbsPdf {
public:

  RooGaussStepExponential();
  RooGaussStepExponential(const char *name, const char *title,
                          RooAbsReal& _x, RooAbsReal& _stepThresh,
                          const RooArgList& _slopeList,
                          const RooArgList& _coefList,
                          Double_t _xmax);

  RooGaussStepExponential(const RooGaussStepExponential& other,
                          const char* name = 0);
  TObject* clone(const char* newname) const override
  { return new RooGaussStepExponential(*this, newname); }
  inline ~RooGaussStepExponential() override { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override;

private:

  RooRealProxy _x;
  RooRealProxy _stepThresh;
  RooListProxy _slopeList;
  RooListProxy _coefList;
  Double_t _xmax;

  Double_t evaluate() const override;

  ClassDefOverride(RooGaussStepExponential,1)
};

class RooGaussStepLaurent : public RooAbsPdf {
public:

  RooGaussStepLaurent();
  RooGaussStepLaurent(const char *name, const char *title,
                      RooAbsReal& _x, RooAbsReal& _stepThresh,
                      const RooArgList& _powerList,
                      const RooArgList& _coefList,
                      Double_t _xmax);

  RooGaussStepLaurent(const RooGaussStepLaurent& other,
                      const char* name = 0);
  TObject* clone(const char* newname) const override
  { return new RooGaussStepLaurent(*this, newname); }
  inline ~RooGaussStepLaurent() override { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override;

private:

  RooRealProxy _x;
  RooRealProxy _stepThresh;
  RooListProxy _powerList;
  RooListProxy _coefList;
  Double_t _xmax;

  Double_t evaluate() const override;

  ClassDefOverride(RooGaussStepLaurent,1)
};

// Combine-owned implementations of the two non-step HZG background families.
// The legacy package-local ModGaus remains available for old workspace
// deserialization, while new workspaces use RooModGaus.
class RooModGaus : public RooAbsPdf {
public:

  RooModGaus();
  RooModGaus(const char *name, const char *title,
             RooAbsReal& _x, RooAbsReal& _m0,
             RooAbsReal& _vl, RooAbsReal& _vr,
             RooAbsReal& _s0, RooAbsReal& _sl,
             RooAbsReal& _sh, Double_t _xlow,
             Double_t _xhigh);

  RooModGaus(const RooModGaus& other, const char* name = 0);
  TObject* clone(const char* newname) const override
  { return new RooModGaus(*this, newname); }
  inline ~RooModGaus() override { }

  bool selfNormalized() const override { return true; }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code,
                              const char* rangeName=0) const override;

private:

  RooRealProxy _x;
  RooRealProxy _m0;
  RooRealProxy _vl;
  RooRealProxy _vr;
  RooRealProxy _s0;
  RooRealProxy _sl;
  RooRealProxy _sh;
  Double_t _xlow;
  Double_t _xhigh;

  mutable Double_t _norm;
  mutable Double_t _lastM0;
  mutable Double_t _lastVl;
  mutable Double_t _lastVr;
  mutable Double_t _lastS0;
  mutable Double_t _lastSl;
  mutable Double_t _lastSh;
  mutable bool _normValid;

  void updateNormCache() const;
  Double_t evaluate() const override;

  ClassDefOverride(RooModGaus,1)
};

// Exponentially modified Gaussian (exGaussian) PDF with an analytic
// finite-range integral, replacing the legacy RooGenericPdf formula.
class RooExpModGauss : public RooAbsPdf {
public:

  RooExpModGauss();
  RooExpModGauss(const char *name, const char *title,
                 RooAbsReal& _x, RooAbsReal& _mu,
                 RooAbsReal& _sigma, RooAbsReal& _lambda);

  RooExpModGauss(const RooExpModGauss& other, const char* name = 0);
  TObject* clone(const char* newname) const override
  { return new RooExpModGauss(*this, newname); }
  inline ~RooExpModGauss() override { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code,
                              const char* rangeName=0) const override;

private:

  RooRealProxy _x;
  RooRealProxy _mu;
  RooRealProxy _sigma;
  RooRealProxy _lambda;

  Double_t evaluate() const override;

  ClassDefOverride(RooExpModGauss,1)
};



#endif
