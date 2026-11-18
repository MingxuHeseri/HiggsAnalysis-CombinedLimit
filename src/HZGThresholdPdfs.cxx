#include "../interface/HZGThresholdPdfs.h"

#include "RooAbsReal.h"
#include "RooArgList.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <stdexcept>

namespace {

Double_t gaussianErfTurnOn(const Double_t value, const Double_t threshold,
                          const Double_t sigma)
{
  if (!std::isfinite(value) || !std::isfinite(threshold) ||
      !std::isfinite(sigma) || !(sigma > 0.)) {
    return 0.;
  }

  const Double_t argument =
      (value - threshold) / (std::sqrt(2.) * sigma);
  if (!std::isfinite(argument)) {
    return argument > 0. ? 1. : 0.;
  }

  return std::min(
      std::max(0.5 * std::erfc(-argument),
               std::numeric_limits<Double_t>::min()),
      1.);
}

Double_t gaussianErfTurnOnPrimitive(const Double_t value,
                                    const Double_t threshold,
                                    const Double_t sigma)
{
  const Double_t displacement = value - threshold;
  const Double_t z = displacement / sigma;
  const Double_t cdf = gaussianErfTurnOn(value, threshold, sigma);
  const Double_t density =
      std::exp(-0.5 * z * z) / std::sqrt(2. * std::acos(-1.));
  return displacement * cdf + sigma * density;
}

Double_t erfTurnOn(const Double_t value, const Double_t threshold,
                  const Double_t width)
{
  if (!std::isfinite(value) || !std::isfinite(threshold) ||
      !std::isfinite(width) || !(width > 0.)) {
    return 0.;
  }

  const Double_t argument = (value - threshold) / width;
  if (!std::isfinite(argument)) {
    return argument > 0. ? 1. : 0.;
  }

  // erfc avoids the catastrophic cancellation in 1 + erf(argument) for a
  // negative turn-on argument.  Keep an infinitesimal positive floor so a
  // mathematically non-zero tail cannot become an invalid zero-valued PDF.
  const Double_t turnOn = std::max(
      0.5 * std::erfc(-argument), std::numeric_limits<Double_t>::min());
  return turnOn < 1. ? turnOn : 1.;
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

Double_t coefficientValue(const RooListProxy& coefficients, const Int_t index)
{
  return static_cast<const RooAbsReal&>(coefficients[index]).getVal();
}

// Primitive of y^power * 0.5 * (1 + erf(y)), for power <= 3.
Double_t erfMomentPrimitive(const Int_t power, const Double_t y)
{
  const Double_t turnOn = 0.5 * std::erfc(-y);
  const Double_t gaussian = std::exp(-y * y);
  const Double_t invSqrtPi = 1. / std::sqrt(std::acos(-1.));
  switch (power) {
    case 0:
      return y * turnOn + 0.5 * gaussian * invSqrtPi;
    case 1:
      return 0.5 * y * y * turnOn - 0.125 * std::erf(y) +
             0.25 * y * gaussian * invSqrtPi;
    case 2:
      return y * y * y * turnOn / 3. +
             (y * y + 1.) * gaussian * invSqrtPi / 6.;
    case 3:
      return y * y * y * y * turnOn / 4. -
             3. * std::erf(y) / 32. +
             (y * y * y / 8. + 3. * y / 16.) * gaussian * invSqrtPi;
    default:
      return 0.;
  }
}

// Primitive of t^power * 0.5 * (1 + erf((t-threshold)/width)).
Double_t scaledErfMomentPrimitive(const Int_t power, const Double_t t,
                                  const Double_t threshold,
                                  const Double_t width)
{
  const Double_t y = (t - threshold) / width;
  const Double_t i0 = erfMomentPrimitive(0, y);
  if (power == 0) {
    return width * i0;
  }
  const Double_t i1 = erfMomentPrimitive(1, y);
  if (power == 1) {
    return width * (threshold * i0 + width * i1);
  }
  const Double_t i2 = erfMomentPrimitive(2, y);
  if (power == 2) {
    return width * (threshold * threshold * i0 +
                    2. * threshold * width * i1 + width * width * i2);
  }
  const Double_t i3 = erfMomentPrimitive(3, y);
  return width * (threshold * threshold * threshold * i0 +
                  3. * threshold * threshold * width * i1 +
                  3. * threshold * width * width * i2 +
                  width * width * width * i3);
}

}  // namespace

ClassImp(RooHZGErfTurnOnPdf)
ClassImp(RooHZGErfExpPdf)
ClassImp(RooHZGPowerExpPdf)
ClassImp(RooErfBernsteinPdf)

RooHZGErfTurnOnPdf::RooHZGErfTurnOnPdf(
    const char* name, const char* title, RooAbsReal& _x, RooAbsReal& _x0,
    RooAbsReal& _sigma)
    : RooAbsPdf(name, title),
      x("x", "Observable", this, _x),
      x0("x0", "Turn-on threshold", this, _x0),
      sigma("sigma", "Gaussian standard deviation", this, _sigma)
{
}

RooHZGErfTurnOnPdf::RooHZGErfTurnOnPdf(
    const RooHZGErfTurnOnPdf& other, const char* name)
    : RooAbsPdf(other, name),
      x("x", this, other.x),
      x0("x0", this, other.x0),
      sigma("sigma", this, other.sigma)
{
}

Double_t RooHZGErfTurnOnPdf::evaluate() const
{
  return gaussianErfTurnOn(x, x0, sigma);
}

Int_t RooHZGErfTurnOnPdf::getAnalyticalIntegral(
    RooArgSet& allVars, RooArgSet& analVars, const char*) const
{
  return matchArgs(allVars, analVars, x) ? 1 : 0;
}

Double_t RooHZGErfTurnOnPdf::analyticalIntegral(
    const Int_t code, const char* rangeName) const
{
  if (code != 1) {
    return 0.;
  }

  // RooFit can pass a disjoint sideband domain as a comma-separated range
  // name. Sum the same closed-form primitive over each component range.
  if (rangeName != nullptr) {
    const std::string ranges(rangeName);
    const std::size_t comma = ranges.find(',');
    if (comma != std::string::npos) {
      Double_t sum = 0.;
      std::size_t begin = 0;
      while (begin < ranges.size()) {
        const std::size_t end = ranges.find(',', begin);
        const std::string token = ranges.substr(begin, end - begin);
        if (!token.empty()) {
          sum += analyticalIntegral(code, token.c_str());
        }
        if (end == std::string::npos) {
          break;
        }
        begin = end + 1;
      }
      return std::isfinite(sum) && sum > 0. ? sum : 0.;
    }
  }

  const Double_t threshold = x0;
  const Double_t width = sigma;
  const Double_t lower = x.min(rangeName);
  const Double_t upper = x.max(rangeName);
  if (!std::isfinite(threshold) || !std::isfinite(width) || !(width > 0.) ||
      !std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
    return 0.;
  }

  const Double_t integral =
      gaussianErfTurnOnPrimitive(upper, threshold, width) -
      gaussianErfTurnOnPrimitive(lower, threshold, width);
  return std::isfinite(integral) && integral > 0. ? integral : 0.;
}

RooHZGErfExpPdf::RooHZGErfExpPdf(const char* name, const char* title,
                                  RooAbsReal& _x, RooAbsReal& _x0,
                                  RooAbsReal& _sigma, RooAbsReal& _lambda)
    : RooAbsPdf(name, title),
      x("x", "Observable", this, _x),
      x0("x0", "Turn-on threshold", this, _x0),
      sigma("sigma", "Turn-on width", this, _sigma),
      lambda("lambda", "Exponential slope", this, _lambda)
{
}

RooHZGErfExpPdf::RooHZGErfExpPdf(const RooHZGErfExpPdf& other,
                                  const char* name)
    : RooAbsPdf(other, name),
      x("x", this, other.x),
      x0("x0", this, other.x0),
      sigma("sigma", this, other.sigma),
      lambda("lambda", this, other.lambda)
{
}

Double_t RooHZGErfExpPdf::evaluate() const
{
  const Double_t value = x;
  const Double_t threshold = x0;
  const Double_t width = sigma;
  const Double_t slope = lambda;
  if (!std::isfinite(slope)) {
    return 0.;
  }

  const Double_t turnOn = erfTurnOn(value, threshold, width);
  if (!(turnOn > 0.)) {
    return 0.;
  }

  const Double_t result = turnOn * positiveExpFromLog(-slope * value);
  return std::isfinite(result) && result >= 0. ? result : 0.;
}

RooHZGPowerExpPdf::RooHZGPowerExpPdf(
    const char* name, const char* title, RooAbsReal& _x, RooAbsReal& _x0,
    RooAbsReal& _alpha, RooAbsReal& _lambda)
    : RooAbsPdf(name, title),
      x("x", "Observable", this, _x),
      x0("x0", "Threshold", this, _x0),
      alpha("alpha", "Power exponent", this, _alpha),
      lambda("lambda", "Exponential slope", this, _lambda)
{
}

RooHZGPowerExpPdf::RooHZGPowerExpPdf(const RooHZGPowerExpPdf& other,
                                      const char* name)
    : RooAbsPdf(other, name),
      x("x", this, other.x),
      x0("x0", this, other.x0),
      alpha("alpha", this, other.alpha),
      lambda("lambda", this, other.lambda)
{
}

Double_t RooHZGPowerExpPdf::evaluate() const
{
  const Double_t value = x;
  const Double_t threshold = x0;
  const Double_t power = alpha;
  const Double_t slope = lambda;
  if (!std::isfinite(value) || !std::isfinite(threshold) ||
      !std::isfinite(power) || !std::isfinite(slope) ||
      !(value > threshold)) {
    return 0.;
  }

  const Double_t distance = value - threshold;
  if (!(distance > 0.) || !std::isfinite(distance)) {
    return 0.;
  }

  const Double_t logValue = power * std::log(distance) - slope * value;
  return positiveExpFromLog(logValue);
}

RooErfBernsteinPdf::RooErfBernsteinPdf(
    const char* name, const char* title, RooAbsReal& _x, RooAbsReal& _x0,
    RooAbsReal& _sigma, const RooArgList& _coefficients)
    : RooAbsPdf(name, title),
      x("x", "Observable", this, _x),
      x0("x0", "Turn-on threshold", this, _x0),
      sigma("sigma", "Turn-on width", this, _sigma),
      coefficients("coefficients", "Bernstein coefficients", this)
{
  if (_coefficients.getSize() != 3 && _coefficients.getSize() != 4) {
    throw std::invalid_argument(
        "RooErfBernsteinPdf requires exactly 3 or 4 Bernstein coefficients");
  }

  for (RooAbsArg* coefficient : _coefficients) {
    RooAbsReal* realCoefficient = dynamic_cast<RooAbsReal*>(coefficient);
    if (realCoefficient == nullptr) {
      throw std::invalid_argument(
          "RooErfBernsteinPdf coefficients must be RooAbsReal objects");
    }
    coefficients.add(*realCoefficient);
  }
}

RooErfBernsteinPdf::RooErfBernsteinPdf(const RooErfBernsteinPdf& other,
                                       const char* name)
    : RooAbsPdf(other, name),
      x("x", this, other.x),
      x0("x0", this, other.x0),
      sigma("sigma", this, other.sigma),
      coefficients("coefficients", this, other.coefficients)
{
}

Double_t RooErfBernsteinPdf::evaluate() const
{
  const Int_t degree = coefficients.getSize() - 1;
  if (degree != 2 && degree != 3) {
    return 0.;
  }

  const Double_t value = x;
  const Double_t minimum = x.min();
  const Double_t maximum = x.max();
  if (!std::isfinite(value) || !std::isfinite(minimum) ||
      !std::isfinite(maximum) || !(maximum > minimum) ||
      value < minimum || value > maximum) {
    return 0.;
  }

  const Double_t turnOn = erfTurnOn(value, x0, sigma);
  if (!(turnOn > 0.)) {
    return 0.;
  }

  const Double_t t = (value - minimum) / (maximum - minimum);
  const Double_t u = 1. - t;
  const Double_t c0 = coefficientValue(coefficients, 0);
  const Double_t c1 = coefficientValue(coefficients, 1);
  const Double_t c2 = coefficientValue(coefficients, 2);
  if (!std::isfinite(c0) || !std::isfinite(c1) || !std::isfinite(c2) ||
      c0 < 0. || c1 < 0. || c2 < 0.) {
    return 0.;
  }

  Double_t bernstein = c0 * u * u + 2. * c1 * t * u + c2 * t * t;
  if (degree == 3) {
    const Double_t c3 = coefficientValue(coefficients, 3);
    if (!std::isfinite(c3) || c3 < 0.) {
      return 0.;
    }
    bernstein = c0 * u * u * u + 3. * c1 * t * u * u +
                3. * c2 * t * t * u + c3 * t * t * t;
  }

  const Double_t result = turnOn * bernstein;
  return std::isfinite(result) && result >= 0. ? result : 0.;
}

Int_t RooErfBernsteinPdf::getAnalyticalIntegral(
    RooArgSet& allVars, RooArgSet& analVars, const char*) const
{
  return matchArgs(allVars, analVars, x) ? 1 : 0;
}

Double_t RooErfBernsteinPdf::analyticalIntegral(
    const Int_t code, const char* rangeName) const
{
  if (code != 1) {
    return 0.;
  }

  // RooFit represents a disjoint normalization domain such as the blinded
  // low and high sidebands as a comma-separated list of named ranges.
  if (rangeName != nullptr) {
    const std::string ranges(rangeName);
    const std::size_t comma = ranges.find(',');
    if (comma != std::string::npos) {
      Double_t sum = 0.;
      std::size_t begin = 0;
      while (begin < ranges.size()) {
        const std::size_t end = ranges.find(',', begin);
        const std::string token = ranges.substr(begin, end - begin);
        if (!token.empty()) {
          sum += analyticalIntegral(code, token.c_str());
        }
        if (end == std::string::npos) {
          break;
        }
        begin = end + 1;
      }
      return std::isfinite(sum) && sum > 0. ? sum : 0.;
    }
  }

  const Int_t degree = coefficients.getSize() - 1;
  const Double_t minimum = x.min();
  const Double_t maximum = x.max();
  const Double_t lower = std::max(x.min(rangeName), minimum);
  const Double_t upper = std::min(x.max(rangeName), maximum);
  const Double_t span = maximum - minimum;
  const Double_t width = static_cast<Double_t>(sigma) / span;
  const Double_t threshold = (static_cast<Double_t>(x0) - minimum) / span;
  if ((degree != 2 && degree != 3) || !std::isfinite(span) || !(span > 0.) ||
      !std::isfinite(width) || !(width > 0.) || !(upper > lower)) {
    return 0.;
  }

  const Double_t tLower = (lower - minimum) / span;
  const Double_t tUpper = (upper - minimum) / span;
  Double_t moments[4] = {0., 0., 0., 0.};
  for (Int_t power = 0; power <= degree; ++power) {
    moments[power] =
        scaledErfMomentPrimitive(power, tUpper, threshold, width) -
        scaledErfMomentPrimitive(power, tLower, threshold, width);
  }

  const Double_t c0 = coefficientValue(coefficients, 0);
  const Double_t c1 = coefficientValue(coefficients, 1);
  const Double_t c2 = coefficientValue(coefficients, 2);
  if (!std::isfinite(c0) || !std::isfinite(c1) || !std::isfinite(c2) ||
      c0 < 0. || c1 < 0. || c2 < 0.) {
    return 0.;
  }

  Double_t integral = c0 * moments[0] +
                      (-2. * c0 + 2. * c1) * moments[1] +
                      (c0 - 2. * c1 + c2) * moments[2];
  if (degree == 3) {
    const Double_t c3 = coefficientValue(coefficients, 3);
    if (!std::isfinite(c3) || c3 < 0.) {
      return 0.;
    }
    integral = c0 * moments[0] +
               (-3. * c0 + 3. * c1) * moments[1] +
               (3. * c0 - 6. * c1 + 3. * c2) * moments[2] +
               (-c0 + 3. * c1 - 3. * c2 + c3) * moments[3];
  }
  integral *= span;
  return std::isfinite(integral) && integral > 0. ? integral : 0.;
}
