/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: HZGRooPdfs.cxx,v 1.2 2013/02/18 22:33:08 gpetrucc Exp $
 * Authors:                                                                  *
 *   Kyle Cranmer
 *                                                                           *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Bernstein basis polynomials are positive-definite in the range [0,1].
// In this implementation, we extend [0,1] to be the range of the parameter.
// There are n+1 Bernstein basis polynomials of degree n.
// Thus, by providing N coefficients that are positive-definite, there 
// is a natural way to have well bahaved polynomail PDFs.
// For any n, the n+1 basis polynomials 'form a partition of unity', eg.
//  they sum to one for all values of x. See
// http://www.idav.ucdavis.edu/education/CAGDNotes/Bernstein-Polynomials.pdf
// Step function piece means that for x < value f(x) = 0 and for
// x >= value: f(x) = RooBernstein
// END_HTML
//

#include "RooFit.h"

#include "Riostream.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <math.h>
#include "TMath.h"
#include "Math/SpecFunc.h"
#include "../interface/HZGRooPdfs.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"

using namespace std;

ClassImp(RooStepBernstein)


//_____________________________________________________________________________
RooStepBernstein::RooStepBernstein()
{
}


//_____________________________________________________________________________
RooStepBernstein::RooStepBernstein(const char* name, const char* title, 
				   RooAbsReal& x, RooAbsReal& step_thresh,
				   const RooArgList& coefList): 
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _stepThresh("stepThresh","value where step function changes",this,step_thresh),
  _coefList("coefficients","List of coefficients",this)
{
  // Constructor
  for (RooAbsArg* coef : coefList) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      cout << "RooStepBernstein::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() 
	   << " is not of type RooAbsReal" << endl ;
      assert(0) ;
    }
    _coefList.add(*coef) ;
  }
}



//_____________________________________________________________________________
RooStepBernstein::RooStepBernstein(const RooStepBernstein& other, const char* name) :
  RooAbsPdf(other, name), 
  _x("x", this, other._x), 
  _stepThresh("stepThresh",this,other._stepThresh),
  _coefList("coefList",this,other._coefList)
{
}


//_____________________________________________________________________________
Double_t RooStepBernstein::evaluate() const 
{
  // now xmin is stepThresh and we return when less than stepThresh
  // warning! this means we are only useful as a PDF when we are convoluted with
  // something that's positive definite!!!!!!
  Double_t xmin = _x.min(); // old  

  //Double_t xmin = _stepThresh;
  Double_t x = (_x - xmin) / (_x.max() - xmin); // rescale to [0,1] where 0 is on stepThresh
  if( x < _stepThresh) return 0; // lazy return

  x = (x - _stepThresh)/(1.0-_stepThresh); // now we get the full polynomial  

  Int_t degree = _coefList.getSize() - 1; // n+1 polys of degree n

  if(degree == 0) {

    return ((RooAbsReal &)_coefList[0]).getVal();

  } else if(degree == 1) {

    Double_t a0 = ((RooAbsReal &)_coefList[0]).getVal(); // c0
    Double_t a1 = ((RooAbsReal &)_coefList[1]).getVal() - a0; // c1 - c0
    return a1 * x + a0;

  } else if(degree == 2) {

    Double_t a0 = ((RooAbsReal &)_coefList[0]).getVal(); // c0
    Double_t a1 = 2 * (((RooAbsReal &)_coefList[1]).getVal() - a0); // 2 * (c1 - c0)
    Double_t a2 = ((RooAbsReal &)_coefList[2]).getVal() - a1 - a0; // c0 - 2 * c1 + c2
    return (a2 * x + a1) * x + a0;

  } else if(degree > 2) {

    Double_t t = x;
    Double_t s = 1 - x;

    Double_t result = ((RooAbsReal &)_coefList[0]).getVal() * s;    
    for(Int_t i = 1; i < degree; i++) {
      result = (result + t * TMath::Binomial(degree, i) * ((RooAbsReal &)_coefList[i]).getVal()) * s;
      t *= x;
    }
    result += t * ((RooAbsReal &)_coefList[degree]).getVal(); 

    return result;
  }

  // in case list of arguments passed is empty
  return TMath::SignalingNaN();
}


//_____________________________________________________________________________
Int_t RooStepBernstein::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  // No analytical calculation available (yet) of integrals over subranges
  if (rangeName && strlen(rangeName)) {
    return 0 ;
  }

  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}


//_____________________________________________________________________________
Double_t RooStepBernstein::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(code==1) ;
  // since we are basically squishing the bernstein polynomial to after some
  // threshold the normalization only need be modified by the range in which
  // the function is zero
  Double_t xmin = _x.min(rangeName);
  Double_t step = _stepThresh;
  Double_t xmax = _x.max(rangeName);
  Int_t degree= _coefList.getSize()-1; // n+1 polys of degree n
  Double_t norm(0) ;

  for (int i=0; i<=degree; ++i){
    // for each of the i Bernstein basis polynomials
    // represent it in the 'power basis' (the naive polynomial basis)
    // where the integral is straight forward.
    Double_t temp = 0;
    for (int j=i; j<=degree; ++j){ // power basis≈ß
      temp += pow(-1.,j-i) * TMath::Binomial(degree, j) * TMath::Binomial(j,i) / (j+1);
    }
    temp *= ((RooAbsReal&)_coefList[i]).getVal(); // include coeff
    norm += temp; // add this basis's contribution to total
  }

  norm *= xmax-(1.0+step)*xmin;
  return norm;
}

ClassImp(RooGaussStepBernstein)

namespace {
  double zeroth(const double x, const double m, const double s, 
		const double a, const double b) {
    const double root2 = sqrt(2.0);    
    const double erf_parm_a = (a+m-x)/(root2*s);
    const double erf_parm_b = (b+m-x)/(root2*s);

    return ( sqrt(M_PI/2.0)*s*
	     ( erf(erf_parm_b) - erf(erf_parm_a) ));
  }
  double zeroth_integral(const double x, const double m, const double s,
			 const double a, const double b) {
    const double root2 = sqrt(2.0);    
    const double root2overpi = sqrt(2.0/M_PI);    
    const double rootpiover2 = sqrt(M_PI/2.0);  
    const double gaus_a = exp(-0.5*(a+m-x)*(a+m-x)/(s*s));
    const double gaus_b = exp(-0.5*(b+m-x)*(b+m-x)/(s*s));
    const double erf_a = erf((a+m-x)/(root2*s));
    const double erf_b = erf((b+m-x)/(root2*s));

    return rootpiover2*s*( gaus_a*s*root2overpi - 
			   gaus_b*s*root2overpi +
			   (a+m-x)*erf_a - 
			   (b+m-x)*erf_b           );
  }
  double first(const double x, const double m, const double s,
	       const double a, const double b) {
    const double root2 = sqrt(2);
    const double rootpiover2 = sqrt(M_PI/2.0);
    const double erf_parm_a = (a+m-x)/(root2*s);
    const double erf_parm_b = (b+m-x)/(root2*s);
    const double gaus_a = exp(-0.5*(a+m-x)*(a+m-x)/(s*s));
    const double gaus_b = exp(-0.5*(b+m-x)*(b+m-x)/(s*s));
    return ( s*s*(gaus_b-gaus_a) + 
	     rootpiover2*s*(a+m-x)*
	     (erf(erf_parm_b)-erf(erf_parm_a)) )/(a-b);
  }
  double first_integral(const double x, const double m, const double s,
			const double a, const double b) {
    const double root2 = sqrt(2.0);    
    const double root2overpi = sqrt(2.0/M_PI);    
    const double rootpiover2 = sqrt(M_PI/2.0);  
    const double gaus_a = exp(-0.5*(a+m-x)*(a+m-x)/(s*s));
    const double gaus_b = exp(-0.5*(b+m-x)*(b+m-x)/(s*s));
    const double erf_a = erf((a+m-x)/(root2*s));
    const double erf_b = erf((b+m-x)/(root2*s));

    return rootpiover2*s*( root2overpi*s*( gaus_a*(a+m-x) +
					   gaus_b*(-2*a+b-m+x) ) +
			   (s*s + (a+m-x)*(a+m-x))*erf_a - 
			   (s*s + (2*a-b+m-x)*(b+m-x))*erf_b       )/(2*(a-b));
  }
  double second(const double x, const double m, const double s,
		const double a, const double b) {
    const double root2 = sqrt(2);
    const double amb = pow(a-b,2.0);
    const double rootpiover2 = sqrt(M_PI/2.0);
    const double erf_a = erf((a+m-x)/(root2*s));
    const double erf_b = erf((b+m-x)/(root2*s));
    const double gaus_a = exp(-0.5*(a+m-x)*(a+m-x)/(s*s));
    const double gaus_b = exp(-0.5*(b+m-x)*(b+m-x)/(s*s));
    return s*( gaus_b*s*(2*a-b+m-x) -
	       gaus_a*s*(a+m-x) +
	       rootpiover2*(s*s + (a+m-x)*(a+m-x))*(erf_b - erf_a) )/amb;
  }
  double second_integral(const double x, const double m, const double s,
			 const double a, const double b) {
    const double root2 = sqrt(2);
    const double amb = 3*pow(a-b,2.0);
    const double rootpiover2 = sqrt(M_PI/2.0);
    const double erf_a = erf((a+m-x)/(root2*s));
    const double erf_b = erf((b+m-x)/(root2*s));
    const double gaus_a = exp(-0.5*(a+m-x)*(a+m-x)/(s*s));
    const double gaus_b = exp(-0.5*(b+m-x)*(b+m-x)/(s*s));


    const double s2 = s*s;

    const double a2 = a*a;

    const double b2 = b*b;
    const double b3 = b*b2;

    const double mmx = m-x;
    const double mmx2 = mmx*mmx;

    const double apmmx = a+m-x;
    const double apmmx2 = apmmx*apmmx;

    
    const double poly1 = (2*s2 + apmmx2);
    const double poly2 = (3*a2 + b2 + 2*s2 + mmx2 - b*mmx - 3*a*(b-m+x));
    const double poly3 = (3*s2 + apmmx2)*apmmx;
    const double poly4 = (b3 + 3*a*(-b2+s2+mmx2) + (3*s2+mmx2)*mmx +
			  3*a2*(b+m-x));

    return s*( gaus_a*s*poly1 -
	       gaus_b*s*poly2 +
	       rootpiover2*(poly3*erf_a - poly4*erf_b ) )/amb;
  }
  double third(const double x, const double m, const double s,
	       const double a, const double b) {
    const double root2 = sqrt(2);
    const double amb   = pow(a-b,3.0);
    const double rootpiover2 = sqrt(M_PI/2.0);
    const double erf_a = erf((a+m-x)/(root2*s));
    const double erf_b = erf((b+m-x)/(root2*s));
    const double gaus_a = exp(-0.5*(a+m-x)*(a+m-x)/(s*s));
    const double gaus_b = exp(-0.5*(b+m-x)*(b+m-x)/(s*s));   
    const double s2 = s*s;
    const double apmmx  = (a+m-x);
    const double apmmx2 = apmmx*apmmx;
    const double poly1 = (3*a*a + b*b + 2*s2 + (m-x)*(m-x) +
			  b*(x-m) - 3*a*(b-m+x));
    const double poly2 = (2*s2 + apmmx2);
    const double poly3 = (3*s2 + apmmx2)*apmmx;
    return s*( gaus_b*s*poly1 -
	       gaus_a*s*poly2 +
	       rootpiover2*poly3*(erf_b - erf_a) )/amb;
  }
  double third_integral(const double x, const double m, const double s,
			const double a, const double b) {
    const double root2 = sqrt(2);
    const double amb   = 4*pow(a-b,3.0);
    const double rootpiover2 = sqrt(M_PI/2.0);
    const double erf_a = erf((a+m-x)/(root2*s));
    const double erf_b = erf((b+m-x)/(root2*s));
    const double gaus_a = exp(-0.5*(a+m-x)*(a+m-x)/(s*s));
    const double gaus_b = exp(-0.5*(b+m-x)*(b+m-x)/(s*s));
    
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x*x3;

    const double s2 = s*s;
    const double s3 = s*s2;
    const double s4 = s*s3;

    const double a2 = a*a;
    const double a3 = a*a2;
    const double a4 = a*a3;

    const double b2 = b*b;
    const double b3 = b*b2;
    const double b4 = b*b3;

    const double m2 = m*m;
    const double m3 = m*m2;
    const double m4 = m*m3;

    const double mmx = m-x;
    const double mmx2 = mmx*mmx;

    const double apm = a+m;
    const double apm2 = apm*apm;

    const double apmmx = a+m-x;
    const double apmmx2 = apmmx*apmmx;

    
    const double poly1 = ( 5*s2 + apmmx2 )*apmmx;
    const double poly2 = ( -4*a3 + b3 + b*(3*s2 + mmx2) - (5*s2 + mmx2)*mmx - b2*mmx + 
			   6*a2*(b-m+x) - 4*a*(b2 + 2*s2 + mmx2 - b*mmx) );
    const double poly3 = ( a4 + 4*a3*m + 6*a2*m2 + 4*a*m3 + m4 + 6*a2*s2 + 12*a*m*s2 +
			   6*m2*s2 + 3*s4 - 4*apm*(apm2 + 3*s2)*x + 6*(apm2 + s2)*x2 -
			   4*apm*x3 + x4 );
    const double poly4 = ( -b4 + m4 + 6*m2*s2 + 3*s4 + 4*a*(b3 +(3*s2 + mmx2)*mmx) +
			   6*a2*(-b2 + s2 + mmx2) + 4*a3*(b+m-x) - 4*m*(m2 + 3*s2)*x +
			   6*(m2+s2)*x2 - 4*m*x3 + x4);
    return s*( gaus_a*s*poly1 +
	       gaus_b*s*poly2 +
	       rootpiover2*(poly3*erf_a - poly4*erf_b ) )/amb;    
  }
  double fourth(const double x, const double m, const double s, 
		const double a, const double b) {
    const double root2 = sqrt(2);
    const double amb   = pow(a-b,4.0);
    const double rootpiover2 = sqrt(M_PI/2.0);
    const double erf_a = erf((a+m-x)/(root2*s));
    const double erf_b = erf((b+m-x)/(root2*s));
    const double gaus_a = exp(-0.5*(a+m-x)*(a+m-x)/(s*s));
    const double gaus_b = exp(-0.5*(b+m-x)*(b+m-x)/(s*s));
    
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x*x3;
    
    const double m2 = m*m;
    const double m3 = m*m2;
    const double m4 = m*m3;
    
    const double a2 = a*a;
    const double a3 = a*a2;
    const double a4 = a*a3;

    const double b2 = b*b;
    const double b3 = b*b2;
    
    const double s2 = s*s;
    const double s3 = s2*s;
    const double s4 = s3*s;

    const double mmx = m-x;
    const double mmx2 = mmx*mmx;
    const double apm = a+m;
    const double apm2 = apm*apm;
    const double apmmx  = (a+m-x);
    const double apmmx2 = apmmx*apmmx;

    const double poly1 = ( -4*a3 + b3 + b*(3*s2 + mmx2) - (5*s2 + mmx2)*mmx -
			   b2*mmx + 6*a2*(b-m+x) - 
			   4*a*(b2 + 2*s2 + mmx2 - b*mmx) );
    const double poly2 = (5*s2 +  apmmx2)*apmmx;
    const double poly3 = ( a4 + 4*a3*m + 6*a2*m2 + 4*a*m3 + m4 +
			   6*a2*s2 + 12*a*m*s2 + 6*m2*s2 + 3*s4 -
			   4*apm*(apm2 + 3*s2)*x + 6*(apm2 + s2)*x2 -
			   4*apm*x3 + x4 );

    return s*( -gaus_b*s*poly1 -
	        gaus_a*s*poly2 +
	       rootpiover2*poly3*(erf_b - erf_a) )/amb;
  }
  double fourth_integral(const double x, const double m, const double s,
			 const double a, const double b) {
    const double root2 = sqrt(2);
    const double amb   = 5*pow(a-b,4.0);
    const double rootpiover2 = sqrt(M_PI/2.0);
    const double erf_a = erf((a+m-x)/(root2*s));
    const double erf_b = erf((b+m-x)/(root2*s));
    const double gaus_a = exp(-0.5*(a+m-x)*(a+m-x)/(s*s));
    const double gaus_b = exp(-0.5*(b+m-x)*(b+m-x)/(s*s));
    
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x*x3;
    const double x5 = x*x4;
    
    const double s2 = s*s;
    const double s3 = s*s2;
    const double s4 = s*s3;

    const double a2 = a*a;
    const double a3 = a*a2;
    const double a4 = a*a3;

    const double b2 = b*b;
    const double b3 = b*b2;
    const double b4 = b*b3;
    const double b5 = b*b4;

    const double m2 = m*m;
    const double m3 = m*m2;
    const double m4 = m*m3;
    const double m5 = m*m4;

    const double mmx = m-x;
    const double mmx2 = mmx*mmx;

    const double apm = a+m;
    const double apm2 = apm*apm;
    const double apm3 = apm*apm2;
    const double apm4 = apm*apm3;

    const double apmmx = a+m-x;
    const double apmmx2 = apmmx*apmmx;

    
    const double poly1 = (s2 + apmmx2)*(8*s2 + apmmx2);
    const double poly2 = ( 5*a4 + b4 + b2*(4*s2 + mmx2) + (s2 + mmx2)*(8*s2 +mmx2) - 
			   b*(7*s2 +mmx2)*mmx - b3*(mmx) - 10*a3*(b-m+x) + 
			   10*a2*(b2 + 2*s2 + mmx2 -b*mmx) - 
			   5*a*(b3 + b*(3*s2 + mmx2) - (5*s2 +mmx2)*mmx - b2*mmx) );
    const double poly3 = ( apmmx*(a4 + 4*a3*m + 6*a2*m2 + 4*a*m3 + m4 + 10*a2*s2 + 20*a*m*s2 +
				  10*m2*s2 + 15*s4 - 4*apm*(apm2 + 5*s2)*x + 2*(3*apm2 + 5*s2)*x2 -
				  4*apm*x3 +x4) );
    const double poly4 = ( 5*a4*b - 10*a3*b2 + 10*a2*b3 - 5*a*b4 + b5 +5*a4*m + 10*a3*m2 +
			   10*a2*m3 + 5*a*m4 + m5 + 10*a3*s2 + 30*a2*m*s2 + 30*a*m2*s2 +
			   10*m3*s2 + 15*a*s4 + 15*m*s4 - 5*(apm4 + 6*apm2*s2 + 3*s4)*x +
			   10*apm*(apm2 + 3*s2)*x2 - 10*(apm2 + s2)*x3 + 5*apm*x4 -x5);
    return s*( gaus_a*s*poly1 -
	       gaus_b*s*poly2 +
	       rootpiover2*(poly3*erf_a - poly4*erf_b ) )/amb;    
  }
  double fifth(const double x, const double m, const double s,
	       const double a, const double b) {
    const double root2 = sqrt(2);
    const double amb   = pow(a-b,5.0);
    const double rootpiover2 = sqrt(M_PI/2.0);
    const double erf_a = erf((a+m-x)/(root2*s));
    const double erf_b = erf((b+m-x)/(root2*s));
    const double gaus_a = exp(-0.5*(a+m-x)*(a+m-x)/(s*s));
    const double gaus_b = exp(-0.5*(b+m-x)*(b+m-x)/(s*s));
    
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x*x3;
    
    const double m2 = m*m;
    const double m3 = m*m2;
    const double m4 = m*m3;
    
    const double a2 = a*a;
    const double a3 = a*a2;
    const double a4 = a*a3;

    const double b2 = b*b;
    const double b3 = b*b2;
    const double b4 = b*b3;

    const double s2 = s*s;
    const double s3 = s2*s;
    const double s4 = s3*s;

    const double mmx = m-x;
    const double mmx2 = mmx*mmx;
    const double apm = a+m;
    const double apm2 = apm*apm;
    const double apmmx  = (a+m-x);
    const double apmmx2 = apmmx*apmmx;

    const double poly1 = ( 5*a4 + b4 + b2*(4*s2 + mmx2) +
			   (s2 + mmx2)*(8*s2 + mmx2) - b*(7*s2 + mmx2)*mmx -
			   b3*mmx - 10*a3*(b-m+x) + 
			   10*a2*(b2+2*s2+mmx2-b*mmx) - 
			   5*a*(b3 + b*(3*s2+mmx2) - (5*s2+mmx2)*mmx - b2*mmx) 
			   );
    const double poly2 = (s2 + apmmx2)*(8*s2 + apmmx2);
    const double poly3 = apmmx*( a4 + 4*a3*m + 6*a2*m2 +4*a*m3 + m4 + 
				 10*a2*s2 + 20*a*m*s2 + 10*m2*s2 + 15*s4 -
				 4*apm*(apm2 + 5*s2)*x + 
				 2*(3*apm2 + 5*s2)*x2 - 4*apm*x3 + x4 );

    return s*( gaus_b*s*poly1 -
	        gaus_a*s*poly2 +
	       rootpiover2*poly3*(erf_b - erf_a) )/amb;
  }
  double fifth_integral(const double x, const double m, const double s,
			const double a, const double b) {
    const double root2 = sqrt(2);
    const double amb   = 6*pow(a-b,5.0);
    const double rootpiover2 = sqrt(M_PI/2.0);
    const double erf_a = erf((a+m-x)/(root2*s));
    const double erf_b = erf((b+m-x)/(root2*s));
    const double gaus_a = exp(-0.5*(a+m-x)*(a+m-x)/(s*s));
    const double gaus_b = exp(-0.5*(b+m-x)*(b+m-x)/(s*s));
    
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x*x3;
    const double x5 = x*x4;
    const double x6 = x*x5;
    
    const double s2 = s*s;
    const double s3 = s*s2;
    const double s4 = s*s3;
    const double s5 = s*s4;
    const double s6 = s*s5;

    const double a2 = a*a;
    const double a3 = a*a2;
    const double a4 = a*a3;
    const double a5 = a*a4;
    const double a6 = a*a5;

    const double b2 = b*b;
    const double b3 = b*b2;
    const double b4 = b*b3;
    const double b5 = b*b4;
    const double b6 = b*b5;

    const double m2 = m*m;
    const double m3 = m*m2;
    const double m4 = m*m3;
    const double m5 = m*m4;
    const double m6 = m*m5;

    const double mmx = m-x;
    const double mmx2 = mmx*mmx;

    const double apm = a+m;
    const double apm2 = apm*apm;
    const double apm3 = apm*apm2;
    const double apm4 = apm*apm3;

    const double apmmx = a+m-x;
    const double apmmx2 = apmmx*apmmx;

    
    const double poly1 = (3*s2 + apmmx2)*(11*s2 + apmmx2)*apmmx;
    const double poly2 = ( -6*a5 + b5 + b3*(5*s2 + mmx2) - b2*(9*s2 + mmx2)*mmx -
			   (3*s2 + mmx2)*(11*s2 + mmx2)*mmx - b4*mmx + 
			   15*a4*(b-m+x) + b*(m4 + 12*m2*s2 + 15*s4 - 4*m*(m2+6*s2)*x + 6*(m2 + 2*s2)*x2
					      - 4*m*x3 + x4) -
			   20*a3*(b2 + 2*s2 + mmx2 - b*mmx) + 
			   15*a2*(b3+b*(3*s2+mmx2) - (5*s2 + mmx2)*mmx - b2*mmx) -
			   6*a*(b4 + b2*(4*s2 + mmx2) + (s2 + mmx2)*(8*s2 + mmx2) - b*(7*s2 + mmx2)*mmx -
				b3*mmx) );
    const double poly3 = ( a6 + 6*a5*m + 15*a4*m2 + 20*a3*m3 + 15*a2*m4 + 6*a*m5 + m6 + 15*a4*s2 +
			   60*a3*m*s2 + 90*a2*m2*s2 + 60*a*m3*s2 + 15*m4*s2 + 45*a2*s4 + 90*a*m*s4 +
			   45*m2*s4 + 15*s6 - 6*apm*(apm4 + 10*apm2*s2 +15*s4)*x +
			   15*(apm4 + 6*apm2*s2 + 3*s4)*x2 - 20*apm*(apm2 + 3*s2)*x3 + 
			   15*(apm2 + s2)*x4 - 6*apm*x5 + x6);
    const double poly4 = ( 6*a5*b - 15*a4*b2 + 20*a3*b3 - 15*a2*b4 + 6*a*b5 - b6 + 6*a5*m +15*a4*m2 +
			   20*a3*m3 + 15*a2*m4 + 6*a*m5 + m6 + 15*a4*s2 + 60*a3*m*s2 + 90*a2*m2*s2 +
			   60*a*m3*s2 + 15*m4*s2 + 45*a2*s4 + 90*a*m*s4 + 45*m2*s4 *15*s6 - 
			   6*apm*(apm4 + 10*apm2*s2 +15*s4)*x + 15*(apm4 + 6*apm2*s2 + 3*s4)*x2 -
			   20*apm*(apm2 + 3*s2)*x3 + 15*(apm2 + s2)*x4 - 6*apm*x5 + x6 );
    return s*( gaus_a*s*poly1 +
	       gaus_b*s*poly2 +
	       rootpiover2*(poly3*erf_a - poly4*erf_b ) )/amb;    
  }
  /*
  double sixth(const double x, const double m, const double s) {
    const double root2 = sqrt(2);
    const double rootpiover2 = sqrt(M_PI/2.0);
    const double erf_parm = (m-x)/(root2*s);
    const double gaus = exp(-0.5*(m-x)*(m-x)/(s*s));
    
    const double x2 = pow(x,2);
    const double x3 = x2*x;
    const double x4 = x3*x;
    const double x5 = x4*x;
    const double x6 = x5*x;

    const double mmx2 = pow(m-x,2);
    //const double mmx3 = mmx2*(m-x);
    //const double mmx4 = mmx3*(m-x);
    //const double mmx5 = mmx4*(m-x);

    const double m2 = pow(m,2);
    const double m3 = m2*m;
    const double m4 = m3*m;
    const double m5 = m4*m;
    const double m6 = m5*m;

    const double s2 = pow(s,2);
    const double s3 = s2*s;
    const double s4 = s3*s;
    const double s5 = s4*s;
    const double s6 = s5*s;
    
    const double poly1 =
      ( (3.0*s2 + mmx2)*(11.0*s2 + mmx2)*(m-x) );
    const double poly2 =
      ( m6 +15.0*m4*s2 + 45.0*m2*s4 + 15.0*s6 - 6.0*m*(m4 + 10.0*m2*s2 + 15.0*s4)*x +
15.0*(m4 + 6.0*m2*s2 +3.0*s4)*x2 - 20.0*m*(m2 + 3.0*s2)*x3 +
15.0*(m2+s2)*x4 - 6.0*m*x5 + x6 );
    const double poly3 =
      ( m6 + 15.0*s6 - 6.0*m5*x + 45.0*s4*x2 + 15.0*s2*x4 + x6 +15.0*m4*(s2 + x2) -
20.0*m3*(3.0*s2*x + x3) +15*m2*(3.0*s4 + 6.0*s2*x2 + x4) -
6.0*m*(15.0*s4*x + 10.0*s2*x3 + x5) );
    return ( s*gaus*( -s*poly1 +
gaus*rootpiover2*(m-x)*poly2 -
gaus*rootpiover2*poly3*erf(erf_parm) ) );
  }
  */

  double poly_conv(const double x, const double mean,
		   const double sigma, const double a,
		   const double b, const int i) {
    switch( i ) {
    case 0:
      return zeroth(x,mean,sigma,a,b);
      break;
    case 1:
      return first(x,mean,sigma,a,b);
      break;
    case 2:
      return second(x,mean,sigma,a,b);
      break;
    case 3:
      return third(x,mean,sigma,a,b);
      break;
    case 4:
      return fourth(x,mean,sigma,a,b);
      break;
    case 5:
      return fifth(x,mean,sigma,a,b);
      break;
    case 6:
      //return sixth(x,mean,sigma);
      //break;
    default:
      assert(0 && "You requested a convolution that we haven't calculated yet!");
      return -1;
    }
    return -1;
  }

  double poly_conv_integral(const double x, const double mean,
		   const double sigma, const double a,
		   const double b, const int i) {
    switch( i ) {
    case 0:
      return zeroth_integral(x,mean,sigma,a,b);
      break;
    case 1:
      return first_integral(x,mean,sigma,a,b);
      break;
    case 2:
      return second_integral(x,mean,sigma,a,b);
      break;
    case 3:
      return third_integral(x,mean,sigma,a,b);
      break;
    case 4:
      return fourth_integral(x,mean,sigma,a,b);
      break;
    case 5:
      return fifth_integral(x,mean,sigma,a,b);
      break;
    case 6:
      //return sixth(x,mean,sigma);
      //break;
    default:
      assert(0 && "You requested a convolution that we haven't calculated yet!");
      return -1;
    }
    return -1;
  }
}


//_____________________________________________________________________________
RooGaussStepBernstein::RooGaussStepBernstein()
{
}


//_____________________________________________________________________________
RooGaussStepBernstein::RooGaussStepBernstein(const char* name,
					     const char* title,
					     RooAbsReal& x,
					     RooAbsReal& mean,
					     RooAbsReal& sigma,
					     RooAbsReal& stepVal,
					     const RooArgList& coefList):
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _mean("mean","mean of gaussian convolved with polyn",this,mean),
  _sigma("sigma","sigma of gaussian convolved with polyn",this,sigma),
  _stepVal("stepVal","position of the step in absolute scale",this,stepVal),
  _coefList("coefficients","List of coefficients",this)
{
  // Constructor
  for (RooAbsArg *coef : coefList) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      cout << "RooGaussStepBernstein::ctor(" << GetName()
	   << ") ERROR: coefficient " << coef->GetName()
	   << " is not of type RooAbsReal" << endl ;
      assert(0) ;
    }    
    _coefList.add(*coef) ;
  }
}



//_____________________________________________________________________________
RooGaussStepBernstein::RooGaussStepBernstein(const RooGaussStepBernstein& other, 
					     const char* name) :
  RooAbsPdf(other, name),
  _x("x", this, other._x),
  _mean("mean",this,other._mean),
  _sigma("sigma",this,other._sigma),
  _stepVal("stepVal",this,other._stepVal),
  _coefList("coefList",this,other._coefList)
{
}


//_____________________________________________________________________________
Double_t RooGaussStepBernstein::evaluate() const
{  
  const Int_t degree = _coefList.getSize() - 1; // n+1 polys of degree n
  
  double poly_vals[degree+1];
  for( int i = 0; i <= degree; ++i ) {
    // we stretch the upper bound by 10% to remove the turn-off of 
    // covolution
    poly_vals[i] = poly_conv(_x,_mean,_sigma,_stepVal,_x.max()*1.10, i);
  }

  double beta = 0.0;
  // iterate through each
  Double_t result = 0.0;
  RooAbsReal* coef = NULL;
  for(Int_t i = 0; i <= degree; ++i) {
    // coef is the bernstein polynomial coefficient
    coef = &(RooAbsReal &)_coefList[i]; 
    // calculate the coefficient in the 'power basis'
    // i.e. the naive polynomial basis
    beta = 0.0;
    for(Int_t k=i; k <= degree; ++k) {      
      beta += (pow(-1.,k-i)*
	       TMath::Binomial(degree,k)*
	       TMath::Binomial(k,i)*
	       poly_vals[k]);
    }
    beta *= coef->getVal();     
    result += beta;
  }

  return result;
}


//_____________________________________________________________________________
Int_t RooGaussStepBernstein::getAnalyticalIntegral(RooArgSet& allVars,
						   RooArgSet& analVars,
						   const char* rangeName) const
{  
  if (matchArgs(allVars, analVars, _x) && _coefList.getSize() <= 5) return 1;
  
  return 0;
}


//_____________________________________________________________________________
Double_t RooGaussStepBernstein::analyticalIntegral(Int_t code,
						   const char* rangeName) const
{
  assert(code==1) ;

  const Int_t degree = _coefList.getSize() - 1; // n+1 polys of degree n
  
  double integral_vals_hi[degree+1],integral_vals_lo[degree+1];
  for( int i = 0; i <= degree; ++i ) {
    // we stretch the upper bound by 10% to remove the turn-off of 
    // covolution
    integral_vals_lo[i] = poly_conv_integral(_x.min(rangeName),_mean,_sigma,_stepVal,
					     _x.max()*1.10, i);
    integral_vals_hi[i] = poly_conv_integral(_x.max(rangeName),_mean,_sigma,_stepVal,
					     _x.max()*1.10, i);
  }

  double beta = 0.0;
  // iterate through each
  Double_t result = 0.0;
  for(Int_t i = 0; i <= degree; ++i) {
    // coef is the bernstein polynomial coefficient
    RooAbsReal *coef = &(RooAbsReal &)_coefList[i]; 
    // calculate the coefficient in the 'power basis'
    // i.e. the naive polynomial basis
    beta = 0.0;
    for(Int_t k=i; k <= degree; ++k) {      
      beta += (pow(-1.,k-i)*
	       TMath::Binomial(degree,k)*
	       TMath::Binomial(k,i)*
	       (integral_vals_hi[k] - integral_vals_lo[k]));
    }
    beta *= coef->getVal();     
    result += beta;
  }

  //cout << result << endl;

  return result;
}

namespace {

  constexpr double kHZGStepSlope = 153.85;
  constexpr double kHZGIntegralEpsilon = 1.e-10;

  double positiveLowerLimit(double value) {
    return std::max(value, std::numeric_limits<double>::min());
  }

  double stepWeight(double x, double step) {
    return std::min(std::max((x - step) * kHZGStepSlope, 0.0), 1.0);
  }

  void copyRealList(RooListProxy& destination, const RooArgList& source,
                    const char* className, const char* listName,
                    const char* ownerName) {
    for (RooAbsArg* argument : source) {
      RooAbsReal* realArgument = dynamic_cast<RooAbsReal*>(argument);
      if (!realArgument) {
        std::cerr << className << "::ctor(" << ownerName << ") ERROR: "
                  << listName << " entry is not a RooAbsReal" << std::endl;
        assert(0);
      }
      destination.add(*realArgument);
    }
  }

  double rawPowerIntegral(double lo, double hi, double power) {
    if (hi <= lo || hi <= 0.) return 0.;
    lo = positiveLowerLimit(lo);
    if (std::fabs(power + 1.) < kHZGIntegralEpsilon) {
      return std::log(hi / lo);
    }
    return (std::pow(hi, power + 1.) - std::pow(lo, power + 1.)) /
           (power + 1.);
  }

  double powerDensity(double x, double power, double step, double xmax) {
    if (x <= 0. || step <= 0. || xmax <= step) return 0.;
    const double alpha = power + 1.;
    if (std::fabs(alpha) < kHZGIntegralEpsilon) {
      return 1. / (x * std::log(xmax / step));
    }
    const double denominator = std::pow(xmax, alpha) - std::pow(step, alpha);
    if (denominator == 0.) return 0.;
    return std::pow(x, power) * alpha / denominator;
  }

  double powerIntegral(double lo, double hi, double power,
                       double step, double xmax) {
    if (hi <= lo || hi <= 0. || step <= 0. || xmax <= step) return 0.;
    lo = std::max(lo, step);
    if (hi <= lo) return 0.;
    const double alpha = power + 1.;
    if (std::fabs(alpha) < kHZGIntegralEpsilon) {
      return std::log(hi / lo) / std::log(xmax / step);
    }
    const double denominator = std::pow(xmax, alpha) - std::pow(step, alpha);
    if (denominator == 0.) return 0.;
    return (std::pow(hi, alpha) - std::pow(lo, alpha)) / denominator;
  }

  double powerFirstIntegral(double lo, double hi, double power,
                            double step, double xmax) {
    if (hi <= lo || hi <= 0. || step <= 0. || xmax <= step) return 0.;
    lo = std::max(lo, step);
    if (hi <= lo) return 0.;
    const double alpha = power + 1.;
    if (std::fabs(alpha) < kHZGIntegralEpsilon) {
      return ((hi - lo) - step * std::log(hi / lo)) /
             std::log(xmax / step);
    }
    const double denominator = std::pow(xmax, alpha) - std::pow(step, alpha);
    if (denominator == 0.) return 0.;
    return alpha / denominator *
           (rawPowerIntegral(lo, hi, power + 1.) -
            step * rawPowerIntegral(lo, hi, power));
  }

  double steppedPowerIntegral(double lo, double hi, double power,
                              double step, double xmax) {
    if (hi <= lo) return 0.;
    const double rampLo = std::max(lo, step);
    const double rampHi = std::min(hi, step + 1. / kHZGStepSlope);
    double result = 0.;
    if (rampHi > rampLo) {
      result += kHZGStepSlope *
                powerFirstIntegral(rampLo, rampHi, power, step, xmax);
    }
    const double tailLo = std::max(lo, step + 1. / kHZGStepSlope);
    if (hi > tailLo) {
      result += powerIntegral(tailLo, hi, power, step, xmax);
    }
    return result;
  }

  double rawExponentialIntegral(double lo, double hi, double slope,
                                double step) {
    if (hi <= lo) return 0.;
    const double shiftedLo = lo - step;
    const double shiftedWidth = hi - lo;
    if (std::fabs(slope) < kHZGIntegralEpsilon) return shiftedWidth;
    return std::exp(slope * shiftedLo) *
           std::expm1(slope * shiftedWidth) / slope;
  }

  double rawExponentialFirstIntegral(double lo, double hi, double slope,
                                     double step) {
    if (hi <= lo) return 0.;
    const double shiftedLo = lo - step;
    const double shiftedHi = hi - step;
    const double scaledWidth = slope * (shiftedHi - shiftedLo);

    // The series avoids cancellation in the primitive for slopes close to 0.
    if (std::fabs(scaledWidth) < 1.e-3) {
      double result = 0.;
      double factorial = 1.;
      for (int order = 0; order <= 8; ++order) {
        if (order > 0) factorial *= order;
        result += std::pow(slope, order) / factorial *
                  (std::pow(shiftedHi, order + 2) -
                   std::pow(shiftedLo, order + 2)) / (order + 2.);
      }
      return result;
    }

    const auto primitive = [slope](double value) {
      return std::exp(slope * value) *
             (value / slope - 1. / (slope * slope));
    };
    return primitive(shiftedHi) - primitive(shiftedLo);
  }

  double exponentialDensity(double x, double slope, double step,
                            double xmax) {
    if (x < step || xmax <= step) return 0.;
    const double width = xmax - step;
    const double denominator = std::fabs(slope) < kHZGIntegralEpsilon
                                 ? width
                                 : std::expm1(slope * width) / slope;
    if (denominator == 0.) return 0.;
    return std::exp(slope * (x - step)) / denominator;
  }

  double exponentialIntegral(double lo, double hi, double slope,
                             double step, double xmax) {
    if (hi <= lo || xmax <= step) return 0.;
    lo = std::max(lo, step);
    if (hi <= lo) return 0.;
    const double denominator = std::fabs(slope) < kHZGIntegralEpsilon
                                 ? xmax - step
                                 : std::expm1(slope * (xmax - step)) / slope;
    if (denominator == 0.) return 0.;
    return rawExponentialIntegral(lo, hi, slope, step) / denominator;
  }

  double exponentialFirstIntegral(double lo, double hi, double slope,
                                  double step, double xmax) {
    if (hi <= lo || xmax <= step) return 0.;
    lo = std::max(lo, step);
    if (hi <= lo) return 0.;
    const double denominator = std::fabs(slope) < kHZGIntegralEpsilon
                                 ? xmax - step
                                 : std::expm1(slope * (xmax - step)) / slope;
    if (denominator == 0.) return 0.;
    return rawExponentialFirstIntegral(lo, hi, slope, step) / denominator;
  }

  double steppedExponentialIntegral(double lo, double hi, double slope,
                                    double step, double xmax) {
    if (hi <= lo) return 0.;
    const double rampLo = std::max(lo, step);
    const double rampHi = std::min(hi, step + 1. / kHZGStepSlope);
    double result = 0.;
    if (rampHi > rampLo) {
      result += kHZGStepSlope *
                exponentialFirstIntegral(rampLo, rampHi, slope, step, xmax);
    }
    const double tailLo = std::max(lo, step + 1. / kHZGStepSlope);
    if (hi > tailLo) {
      result += exponentialIntegral(tailLo, hi, slope, step, xmax);
    }
    return result;
  }

  double bernsteinIntegral(double lo, double hi, double step, double xmax,
                           int degree, int index) {
    if (degree < 0 || index < 0 || index > degree || hi <= lo || xmax <= step) {
      return 0.;
    }
    const double width = xmax - step;
    const double clippedLo = std::max(lo, step);
    const double clippedHi = std::min(hi, xmax);
    if (clippedHi <= clippedLo) return 0.;

    const double zLo = std::min(std::max((clippedLo - step) / width, 0.), 1.);
    const double zHi = std::min(std::max((clippedHi - step) / width, 0.), 1.);
    const double cdfLo = zLo <= 0. ? 0. :
                         zLo >= 1. ? 1. :
                         TMath::BetaIncomplete(zLo, index + 1., degree - index + 1.);
    const double cdfHi = zHi <= 0. ? 0. :
                         zHi >= 1. ? 1. :
                         TMath::BetaIncomplete(zHi, index + 1., degree - index + 1.);
    return width * (cdfHi - cdfLo) / (degree + 1.);
  }

}

ClassImp(RooGaussStepBernsteinFast)

RooGaussStepBernsteinFast::RooGaussStepBernsteinFast() : _xmax(0.)
{
}

RooGaussStepBernsteinFast::RooGaussStepBernsteinFast(
    const char* name, const char* title, RooAbsReal& x,
    RooAbsReal& stepThresh, const RooArgList& coefList, Double_t xmax):
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _stepThresh("stepThresh", "value where the step starts", this, stepThresh),
  _coefList("coefficients", "Bernstein coefficients", this),
  _xmax(xmax)
{
  copyRealList(_coefList, coefList, "RooGaussStepBernsteinFast",
               "coefficient", GetName());
}

RooGaussStepBernsteinFast::RooGaussStepBernsteinFast(
    const RooGaussStepBernsteinFast& other, const char* name):
  RooAbsPdf(other, name),
  _x("x", this, other._x),
  _stepThresh("stepThresh", this, other._stepThresh),
  _coefList("coefficients", this, other._coefList),
  _xmax(other._xmax)
{
}

Double_t RooGaussStepBernsteinFast::evaluate() const
{
  const double x = _x;
  const double step = _stepThresh;
  if (_coefList.getSize() == 0 || _xmax <= step || x < step || x > _xmax) {
    return 0.;
  }

  const double z = (x - step) / (_xmax - step);
  const int degree = _coefList.getSize() - 1;
  double result = 0.;
  if (z <= 0.) {
    return static_cast<RooAbsReal&>(_coefList[0]).getVal();
  }
  if (z >= 1.) {
    return static_cast<RooAbsReal&>(_coefList[degree]).getVal();
  }

  // Build all basis values from the numerically safer endpoint.  This is
  // O(n), avoids repeated Binomial()/pow() calls, and remains stable close
  // to either endpoint for the orders used by the HZG fits.
  if (z <= 0.5) {
    double basis = std::pow(1. - z, degree);  // B_0,n
    for (int i = 0; i <= degree; ++i) {
      result += static_cast<RooAbsReal&>(_coefList[i]).getVal() * basis;
      if (i < degree) {
        basis *= (static_cast<double>(degree - i) / (i + 1)) *
                 z / (1. - z);
      }
    }
  } else {
    double basis = std::pow(z, degree);  // B_n,n
    for (int i = degree; i >= 0; --i) {
      result += static_cast<RooAbsReal&>(_coefList[i]).getVal() * basis;
      if (i > 0) {
        basis *= (static_cast<double>(i) / (degree - i + 1)) *
                 (1. - z) / z;
      }
    }
  }
  return result;
}

Int_t RooGaussStepBernsteinFast::getAnalyticalIntegral(
    RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}

Double_t RooGaussStepBernsteinFast::analyticalIntegral(
    Int_t code, const char* rangeName) const
{
  if (code != 1 || _coefList.getSize() == 0) return 0.;
  const double lo = _x.min(rangeName);
  const double hi = _x.max(rangeName);
  const int degree = _coefList.getSize() - 1;
  double result = 0.;
  for (int i = 0; i <= degree; ++i) {
    const double coefficient =
        static_cast<RooAbsReal&>(_coefList[i]).getVal();
    result += coefficient *
              bernsteinIntegral(lo, hi, _stepThresh, _xmax, degree, i);
  }
  return result;
}

ClassImp(RooGaussStepPowerLaw)
ClassImp(RooGaussStepExponential)
ClassImp(RooGaussStepLaurent)

RooGaussStepPowerLaw::RooGaussStepPowerLaw() : _xmax(0.)
{
}

RooGaussStepPowerLaw::RooGaussStepPowerLaw(
    const char* name, const char* title, RooAbsReal& x,
    RooAbsReal& stepThresh, const RooArgList& powerList,
    const RooArgList& coefList, Double_t xmax):
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _stepThresh("stepThresh", "value where the step starts", this, stepThresh),
  _powerList("powers", "Power-law exponents", this),
  _coefList("coefficients", "Power-law coefficients", this),
  _xmax(xmax)
{
  copyRealList(_powerList, powerList, "RooGaussStepPowerLaw", "power", GetName());
  copyRealList(_coefList, coefList, "RooGaussStepPowerLaw", "coefficient", GetName());
}

RooGaussStepPowerLaw::RooGaussStepPowerLaw(
    const RooGaussStepPowerLaw& other, const char* name):
  RooAbsPdf(other, name),
  _x("x", this, other._x),
  _stepThresh("stepThresh", this, other._stepThresh),
  _powerList("powers", this, other._powerList),
  _coefList("coefficients", this, other._coefList),
  _xmax(other._xmax)
{
}

Double_t RooGaussStepPowerLaw::evaluate() const
{
  const double x = _x;
  const double step = _stepThresh;
  if (stepWeight(x, step) <= 0. || _xmax <= step) return 0.;

  const int nTerms = std::min(_powerList.getSize(), _coefList.getSize());
  double result = 0.;
  for (int i = 0; i < nTerms; ++i) {
    const double power = static_cast<RooAbsReal&>(_powerList[i]).getVal();
    const double coefficient = static_cast<RooAbsReal&>(_coefList[i]).getVal();
    result += coefficient * powerDensity(x, power, step, _xmax);
  }
  return stepWeight(x, step) * result;
}

Int_t RooGaussStepPowerLaw::getAnalyticalIntegral(
    RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}

Double_t RooGaussStepPowerLaw::analyticalIntegral(
    Int_t code, const char* rangeName) const
{
  if (code != 1) return 0.;
  const double lo = _x.min(rangeName);
  const double hi = _x.max(rangeName);
  const int nTerms = std::min(_powerList.getSize(), _coefList.getSize());
  double result = 0.;
  for (int i = 0; i < nTerms; ++i) {
    const double power = static_cast<RooAbsReal&>(_powerList[i]).getVal();
    const double coefficient = static_cast<RooAbsReal&>(_coefList[i]).getVal();
    result += coefficient * steppedPowerIntegral(lo, hi, power, _stepThresh, _xmax);
  }
  return result;
}

RooGaussStepExponential::RooGaussStepExponential() : _xmax(0.)
{
}

RooGaussStepExponential::RooGaussStepExponential(
    const char* name, const char* title, RooAbsReal& x,
    RooAbsReal& stepThresh, const RooArgList& slopeList,
    const RooArgList& coefList, Double_t xmax):
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _stepThresh("stepThresh", "value where the step starts", this, stepThresh),
  _slopeList("slopes", "Exponential slopes", this),
  _coefList("coefficients", "Exponential coefficients", this),
  _xmax(xmax)
{
  copyRealList(_slopeList, slopeList, "RooGaussStepExponential", "slope", GetName());
  copyRealList(_coefList, coefList, "RooGaussStepExponential", "coefficient", GetName());
}

RooGaussStepExponential::RooGaussStepExponential(
    const RooGaussStepExponential& other, const char* name):
  RooAbsPdf(other, name),
  _x("x", this, other._x),
  _stepThresh("stepThresh", this, other._stepThresh),
  _slopeList("slopes", this, other._slopeList),
  _coefList("coefficients", this, other._coefList),
  _xmax(other._xmax)
{
}

Double_t RooGaussStepExponential::evaluate() const
{
  const double x = _x;
  const double step = _stepThresh;
  if (stepWeight(x, step) <= 0. || _xmax <= step) return 0.;

  const int nTerms = std::min(_slopeList.getSize(), _coefList.getSize());
  double result = 0.;
  for (int i = 0; i < nTerms; ++i) {
    const double slope = static_cast<RooAbsReal&>(_slopeList[i]).getVal();
    const double coefficient = static_cast<RooAbsReal&>(_coefList[i]).getVal();
    result += coefficient * exponentialDensity(x, slope, step, _xmax);
  }
  return stepWeight(x, step) * result;
}

Int_t RooGaussStepExponential::getAnalyticalIntegral(
    RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}

Double_t RooGaussStepExponential::analyticalIntegral(
    Int_t code, const char* rangeName) const
{
  if (code != 1) return 0.;
  const double lo = _x.min(rangeName);
  const double hi = _x.max(rangeName);
  const int nTerms = std::min(_slopeList.getSize(), _coefList.getSize());
  double result = 0.;
  for (int i = 0; i < nTerms; ++i) {
    const double slope = static_cast<RooAbsReal&>(_slopeList[i]).getVal();
    const double coefficient = static_cast<RooAbsReal&>(_coefList[i]).getVal();
    result += coefficient * steppedExponentialIntegral(lo, hi, slope,
                                                       _stepThresh, _xmax);
  }
  return result;
}

RooGaussStepLaurent::RooGaussStepLaurent() : _xmax(0.)
{
}

RooGaussStepLaurent::RooGaussStepLaurent(
    const char* name, const char* title, RooAbsReal& x,
    RooAbsReal& stepThresh, const RooArgList& powerList,
    const RooArgList& coefList, Double_t xmax):
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _stepThresh("stepThresh", "value where the step starts", this, stepThresh),
  _powerList("powers", "Laurent exponents", this),
  _coefList("coefficients", "Laurent coefficients", this),
  _xmax(xmax)
{
  copyRealList(_powerList, powerList, "RooGaussStepLaurent", "power", GetName());
  copyRealList(_coefList, coefList, "RooGaussStepLaurent", "coefficient", GetName());
}

RooGaussStepLaurent::RooGaussStepLaurent(
    const RooGaussStepLaurent& other, const char* name):
  RooAbsPdf(other, name),
  _x("x", this, other._x),
  _stepThresh("stepThresh", this, other._stepThresh),
  _powerList("powers", this, other._powerList),
  _coefList("coefficients", this, other._coefList),
  _xmax(other._xmax)
{
}

Double_t RooGaussStepLaurent::evaluate() const
{
  const double x = _x;
  const double step = _stepThresh;
  if (stepWeight(x, step) <= 0. || _xmax <= step) return 0.;

  const int nTerms = std::min(_powerList.getSize(), _coefList.getSize());
  double result = 0.;
  for (int i = 0; i < nTerms; ++i) {
    const double power = static_cast<RooAbsReal&>(_powerList[i]).getVal();
    const double coefficient = static_cast<RooAbsReal&>(_coefList[i]).getVal();
    result += coefficient * powerDensity(x, power, step, _xmax);
  }
  return stepWeight(x, step) * result;
}

Int_t RooGaussStepLaurent::getAnalyticalIntegral(
    RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}

Double_t RooGaussStepLaurent::analyticalIntegral(
    Int_t code, const char* rangeName) const
{
  if (code != 1) return 0.;
  const double lo = _x.min(rangeName);
  const double hi = _x.max(rangeName);
  const int nTerms = std::min(_powerList.getSize(), _coefList.getSize());
  double result = 0.;
  for (int i = 0; i < nTerms; ++i) {
    const double power = static_cast<RooAbsReal&>(_powerList[i]).getVal();
    const double coefficient = static_cast<RooAbsReal&>(_coefList[i]).getVal();
    result += coefficient * steppedPowerIntegral(lo, hi, power, _stepThresh, _xmax);
  }
  return result;
}

namespace {

constexpr int kRooModGausNormBins = 256;

double modGausRawValue(double x, double m0, double vl, double vr,
                       double s0, double sl, double sh,
                       double xlow, double xhigh) {
  if (!(xhigh > xlow) || x < xlow || x > xhigh) return 0.;

  const double power = vl + vr * (x - xlow) / (xhigh - xlow);
  if (!std::isfinite(power)) return 0.;

  // The historical parameterization allows s0=0.  At the mode the two
  // One-sided limits are finite in that case, whereas the literal 0/0
  // expression is not.  Average the available one-sided limits at the mode.
  if (x == m0 && s0 <= 0.) {
    double limitSum = 0.;
    int nLimits = 0;
    if (m0 > xlow && sl > s0) {
      const double leftRatio = (m0 - xlow) / (sl - s0);
      limitSum += std::exp(-std::pow(leftRatio, power));
      ++nLimits;
    }
    if (xhigh > m0 && sh > s0) {
      const double rightRatio = (xhigh - m0) / (sh - s0);
      limitSum += std::exp(-std::pow(rightRatio, power));
      ++nLimits;
    }
    return nLimits > 0 ? limitSum / nLimits : 0.;
  }

  double width = 0.;
  if (x <= m0 && m0 > xlow) {
    width = sl + (s0 - sl) * (x - xlow) / (m0 - xlow);
  } else if (x >= m0 && xhigh > m0) {
    width = s0 + (sh - s0) * (x - m0) / (xhigh - m0);
  }
  if (!(width > 0.) || !std::isfinite(width)) return 0.;

  const double ratio = std::fabs(x - m0) / width;
  if (ratio == 0.) {
    if (power < 0.) return 0.;
    return power == 0. ? std::exp(-1.) : 1.;
  }

  const double poweredRatio = std::pow(ratio, power);
  if (!std::isfinite(poweredRatio)) return poweredRatio > 0. ? 0. : 1.;
  return std::exp(-poweredRatio);
}

double midpointIntegral(double lo, double hi, double m0, double vl,
                        double vr, double s0, double sl, double sh,
                        double xlow, double xhigh) {
  if (!(hi > lo)) return 0.;
  const double width = (hi - lo) / kRooModGausNormBins;
  double result = 0.;
  for (int i = 0; i < kRooModGausNormBins; ++i) {
    const double x = lo + (i + 0.5) * width;
    result += modGausRawValue(x, m0, vl, vr, s0, sl, sh, xlow, xhigh);
  }
  return result * width;
}

double modGausIntegral(double lo, double hi, double m0, double vl,
                       double vr, double s0, double sl, double sh,
                       double xlow, double xhigh) {
  const double clippedLo = std::max(lo, xlow);
  const double clippedHi = std::min(hi, xhigh);
  if (!(clippedHi > clippedLo)) return 0.;

  const double split = std::min(std::max(m0, clippedLo), clippedHi);
  return midpointIntegral(clippedLo, split, m0, vl, vr, s0, sl, sh,
                          xlow, xhigh) +
         midpointIntegral(split, clippedHi, m0, vl, vr, s0, sl, sh,
                          xlow, xhigh);
}

double normalCdf(double value) {
  return 0.5 * TMath::Erfc(-value / std::sqrt(2.));
}

double expModGaussCdf(double x, double mu, double sigma, double lambda) {
  if (!(sigma > 0.) || !(lambda > 0.)) return 0.;
  const double a = (x - mu) / sigma;
  const double b = a - lambda * sigma;
  const double exponent = -lambda * (x - mu) +
                          0.5 * lambda * lambda * sigma * sigma;
  const double tailCdf = normalCdf(b);
  if (!(tailCdf > 0.)) return normalCdf(a);
  const double logTail = exponent + std::log(tailCdf);
  if (!std::isfinite(logTail)) return logTail > 0. ? 0. : normalCdf(a);
  if (logTail >= std::log(std::numeric_limits<double>::max())) return 0.;
  const double tail = logTail <= std::log(std::numeric_limits<double>::min())
                          ? 0.
                          : std::exp(logTail);
  return std::min(std::max(normalCdf(a) - tail, 0.), 1.);
}

}  // namespace

ClassImp(RooModGaus)

RooModGaus::RooModGaus()
    : RooAbsPdf(),
      _xlow(0.), _xhigh(0.), _norm(1.),
      _lastM0(0.), _lastVl(0.), _lastVr(0.), _lastS0(0.),
      _lastSl(0.), _lastSh(0.), _normValid(false) {
}

RooModGaus::RooModGaus(const char* name, const char* title, RooAbsReal& x,
                       RooAbsReal& m0, RooAbsReal& vl, RooAbsReal& vr,
                       RooAbsReal& s0, RooAbsReal& sl, RooAbsReal& sh,
                       Double_t xlow, Double_t xhigh)
    : RooAbsPdf(name, title),
      _x("x", "Dependent", this, x),
      _m0("m0", "mode", this, m0),
      _vl("vl", "left power", this, vl),
      _vr("vr", "power slope", this, vr),
      _s0("s0", "mode width", this, s0),
      _sl("sl", "left width", this, sl),
      _sh("sh", "right width", this, sh),
      _xlow(xlow), _xhigh(xhigh), _norm(1.),
      _lastM0(0.), _lastVl(0.), _lastVr(0.), _lastS0(0.),
      _lastSl(0.), _lastSh(0.), _normValid(false) {
}

RooModGaus::RooModGaus(const RooModGaus& other, const char* name)
    : RooAbsPdf(other, name),
      _x("x", this, other._x),
      _m0("m0", this, other._m0),
      _vl("vl", this, other._vl),
      _vr("vr", this, other._vr),
      _s0("s0", this, other._s0),
      _sl("sl", this, other._sl),
      _sh("sh", this, other._sh),
      _xlow(other._xlow), _xhigh(other._xhigh), _norm(other._norm),
      _lastM0(other._lastM0), _lastVl(other._lastVl),
      _lastVr(other._lastVr), _lastS0(other._lastS0),
      _lastSl(other._lastSl), _lastSh(other._lastSh),
      _normValid(other._normValid) {
}

void RooModGaus::updateNormCache() const {
  const double m0 = _m0;
  const double vl = _vl;
  const double vr = _vr;
  const double s0 = _s0;
  const double sl = _sl;
  const double sh = _sh;
  if (!_normValid || m0 != _lastM0 || vl != _lastVl || vr != _lastVr ||
      s0 != _lastS0 || sl != _lastSl || sh != _lastSh) {
    _norm = modGausIntegral(_xlow, _xhigh, m0, vl, vr, s0, sl, sh,
                            _xlow, _xhigh);
    if (!std::isfinite(_norm) || !(_norm > 0.)) _norm = 1.;
    _lastM0 = m0;
    _lastVl = vl;
    _lastVr = vr;
    _lastS0 = s0;
    _lastSl = sl;
    _lastSh = sh;
    _normValid = true;
  }
}

Double_t RooModGaus::evaluate() const {
  updateNormCache();
  const double m0 = _m0;
  const double vl = _vl;
  const double vr = _vr;
  const double s0 = _s0;
  const double sl = _sl;
  const double sh = _sh;
  return modGausRawValue(_x, m0, vl, vr, s0, sl, sh, _xlow, _xhigh) /
         _norm;
}

Int_t RooModGaus::getAnalyticalIntegral(RooArgSet& allVars,
                                        RooArgSet& analVars,
                                        const char* /*rangeName*/) const {
  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}

Double_t RooModGaus::analyticalIntegral(Int_t code,
                                        const char* rangeName) const {
  if (code != 1) return 0.;
  updateNormCache();
  return modGausIntegral(_x.min(rangeName), _x.max(rangeName), _m0, _vl, _vr,
                         _s0, _sl, _sh, _xlow, _xhigh) / _norm;
}

ClassImp(RooExpModGauss)

RooExpModGauss::RooExpModGauss() : RooAbsPdf() {
}

RooExpModGauss::RooExpModGauss(const char* name, const char* title,
                               RooAbsReal& x, RooAbsReal& mu,
                               RooAbsReal& sigma, RooAbsReal& lambda)
    : RooAbsPdf(name, title),
      _x("x", "Dependent", this, x),
      _mu("mu", "Gaussian mean", this, mu),
      _sigma("sigma", "Gaussian width", this, sigma),
      _lambda("lambda", "Exponential rate", this, lambda) {
}

RooExpModGauss::RooExpModGauss(const RooExpModGauss& other,
                               const char* name)
    : RooAbsPdf(other, name),
      _x("x", this, other._x),
      _mu("mu", this, other._mu),
      _sigma("sigma", this, other._sigma),
      _lambda("lambda", this, other._lambda) {
}

Double_t RooExpModGauss::evaluate() const {
  const double sigma = _sigma;
  const double lambda = _lambda;
  if (!(sigma > 0.) || !(lambda > 0.)) return 0.;

  const double arg = (_mu + lambda * sigma * sigma - _x) /
                     (std::sqrt(2.) * sigma);
  const double erfcValue = TMath::Erfc(arg);
  if (!(erfcValue > 0.)) return 0.;

  const double logValue = std::log(lambda * 0.5) +
                          lambda * (_mu - _x) +
                          0.5 * lambda * lambda * sigma * sigma +
                          std::log(erfcValue);
  if (logValue < std::log(std::numeric_limits<double>::min())) return 0.;
  if (logValue > std::log(std::numeric_limits<double>::max())) return
      std::numeric_limits<double>::max();
  if (!std::isfinite(logValue)) return logValue > 0. ?
      std::numeric_limits<double>::max() : 0.;
  return std::exp(logValue);
}

Int_t RooExpModGauss::getAnalyticalIntegral(
    RooArgSet& allVars, RooArgSet& analVars,
    const char* /*rangeName*/) const {
  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}

Double_t RooExpModGauss::analyticalIntegral(Int_t code,
                                            const char* rangeName) const {
  if (code != 1) return 0.;
  const double lo = _x.min(rangeName);
  const double hi = _x.max(rangeName);
  return std::max(0., expModGaussCdf(hi, _mu, _sigma, _lambda) -
                         expModGaussCdf(lo, _mu, _sigma, _lambda));
}
