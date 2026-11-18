"""Compare the cached analytical integrals to quadrature of the raw density."""
import math

import ROOT


def check_ranges(pdf, mass):
    saved = mass.getVal()
    # Reuse a range name with different endpoints, interleaving two intervals.
    for lo, hi in [(98., 116.), (131., 156.), (98., 116.)]:
        mass.setRange("probe", lo, hi)
        analytical = pdf.analyticalIntegral(1, "probe")
        intervals = 512
        width = (hi - lo) / intervals
        weights = []
        for j in range(intervals + 1):
            mass.setVal(lo + j * width)
            coefficient = 1 if j in (0, intervals) else 4 if j % 2 else 2
            weights.append(coefficient * pdf.getVal())
        reference = math.fsum(weights) * width / 3
        assert math.isclose(analytical, reference, rel_tol=1e-9, abs_tol=1e-12), (
            pdf.GetName(), lo, hi, analytical, reference)
    mass.setVal(saved)


def main():
    ROOT.gROOT.SetBatch(True)
    assert ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit") >= 0
    mass = ROOT.RooRealVar("mass", "", 110., 97., 162.)
    step = ROOT.RooRealVar("step", "", 105., 95., 115.)
    sigma = ROOT.RooRealVar("sigma", "", 3., 1., 10.)
    a = ROOT.RooRealVar("a", "", .4, -1., 1.)
    b = ROOT.RooRealVar("b", "", .1, -1., 1.)
    f = ROOT.RooRealVar("f", "", .3, 0., 1.)
    slopes = ROOT.RooArgList(a, b)
    bern = ROOT.RooErfProductPdf("bern", "", mass, step, sigma, slopes, ROOT.RooArgList(), 0)
    exp = ROOT.RooErfProductPdf("exp", "", mass, step, sigma, slopes, ROOT.RooArgList(f), 1)
    for parameter, value in [(a,.4), (f,.8), (a,-.03), (b,1e-8), (step,100.), (sigma,7.)]:
        parameter.setVal(value)
        for pdf in [bern, exp]:
            check_ranges(pdf, mass)
    mass.setRange(96., 165.)
    for pdf in [bern, exp]:
        check_ranges(pdf, mass)


if __name__ == "__main__":
    main()
