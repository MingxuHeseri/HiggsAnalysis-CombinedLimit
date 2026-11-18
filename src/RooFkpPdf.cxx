#include "../interface/RooFkpPdf.h"

#include "RooAbsReal.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

Int_t supportedOrder(const Int_t order)
{
  return std::max(0, std::min(2, order));
}

Double_t positiveExpFromLog(const Double_t logValue)
{
  const Double_t largest = std::numeric_limits<Double_t>::max();
  const Double_t smallest = std::numeric_limits<Double_t>::min();
  const Double_t logLargest = std::log(largest);
  const Double_t logSmallest = std::log(smallest);

  if (std::isnan(logValue) || logValue <= logSmallest) {
    return 0.;
  }
  if (logValue >= logLargest) {
    return largest;
  }
  return std::exp(logValue);
}

}  // namespace

ClassImp(RooFkpPdf)

RooFkpPdf::RooFkpPdf(const char* name, const char* title, const Int_t order,
                     RooAbsReal& _observable, RooAbsReal& _scale,
                     RooAbsReal& _b, RooAbsReal& _a0, RooAbsReal& _a1,
                     RooAbsReal& _a2)
    : RooAbsPdf(name, title),
      observable("observable", "Observable", this, _observable),
      scale("scale", "Dimensionful Fk scale", this, _scale),
      b("b", "Endpoint exponent", this, _b),
      a0("a0", "Leading logarithmic exponent", this, _a0),
      a1("a1", "Linear logarithmic exponent", this, _a1),
      a2("a2", "Quadratic logarithmic exponent", this, _a2),
      order_(supportedOrder(order))
{
}

RooFkpPdf::RooFkpPdf(const RooFkpPdf& other, const char* name)
    : RooAbsPdf(other, name),
      observable("observable", this, other.observable),
      scale("scale", this, other.scale),
      b("b", this, other.b),
      a0("a0", this, other.a0),
      a1("a1", this, other.a1),
      a2("a2", this, other.a2),
      order_(other.order_)
{
}

Double_t RooFkpPdf::evaluate() const
{
  const Double_t observableValue = observable;
  const Double_t scaleValue = scale;
  const Double_t endpointExponent = b;
  const Double_t exponent0 = a0;
  const Double_t exponent1 = a1;
  const Double_t exponent2 = a2;
  if (!std::isfinite(observableValue) || !std::isfinite(scaleValue) ||
      !std::isfinite(endpointExponent) || !std::isfinite(exponent0) ||
      !std::isfinite(exponent1) || !std::isfinite(exponent2) ||
      !(scaleValue > 0.)) {
    return 0.;
  }

  const Double_t x = observableValue / scaleValue;
  if (!std::isfinite(x) || !(x > 0.) || !(x < 1.)) {
    return 0.;
  }

  const Double_t logX = std::log(x);
  const Double_t logOneMinusX = std::log1p(-x);
  if (!std::isfinite(logX) || !std::isfinite(logOneMinusX)) {
    return 0.;
  }

  Double_t exponent = exponent0;
  if (order_ >= 1) {
    exponent += exponent1 * logX;
  }
  if (order_ >= 2) {
    exponent += exponent2 * logX * logX;
  }
  if (!std::isfinite(exponent)) {
    return 0.;
  }

  const Double_t logValue = endpointExponent * logOneMinusX + exponent * logX;
  return positiveExpFromLog(logValue);
}
