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


def check_bernstein_endpoint():
    mass = ROOT.RooRealVar("endpoint_mass", "", 125., 100., 165.)
    step = ROOT.RooRealVar("endpoint_step", "", 100., 95., 115.)
    sigma = ROOT.RooRealVar("endpoint_sigma", "", 1., .5, 15.)
    for case, logs in enumerate(((5., -20., 5., -20., -20.), (-20., 5., -20., 5., -20.))):
        step.setVal(100.)
        sigma.setVal(1.)
        ratios = [ROOT.RooConstVar(f"c{case}_{i}", "", math.exp(v)) for i, v in enumerate(logs)]
        pdf = ROOT.RooErfProductPdf("endpoint", "", mass, step, sigma, ROOT.RooArgList(*ratios), ROOT.RooArgList(), 0)
        coefficients = [1., *(p.getVal() for p in ratios)]
        for m in (100., 132.5, 164.935, 164.9935, 164.99935, 165.):
            mass.setVal(m)
            u = (m - 100.) / 65.
            complement = (165. - m) / 65.
            bern = math.fsum(c * math.comb(5, i) * u**i * complement**(5-i) for i, c in enumerate(coefficients))
            turn = .5 * math.erfc((100. - m) / math.sqrt(2))
            assert math.isclose(pdf.getVal(), bern * turn, rel_tol=1e-12)
        for threshold, width in ((100., 1.), (115., 15.)):
            step.setVal(threshold)
            sigma.setVal(width)
            for lo in (164.5, 164.75, 164.99, 164.999):
                mass.setRange("endpoint_bin", lo, 165.)
                span = 165. - lo
                # Independent composite Simpson rule on the original basis.
                panels = 512
                terms = []
                for j in range(panels + 1):
                    z = j / panels
                    u = ((lo - 100.) + z * span) / 65.
                    v = (1. - z) * span / 65.
                    bern = math.fsum(c * math.comb(5, i) * u**i * v**(5-i)
                                     for i, c in enumerate(coefficients))
                    turn = .5 * math.erfc(((threshold - lo) - z * span) / (math.sqrt(2) * width))
                    terms.append((1 if j in (0, panels) else 4 if j % 2 else 2) * bern * turn)
                reference = span * math.fsum(terms) / (3 * panels)
                assert math.isclose(pdf.analyticalIntegral(1, "endpoint_bin"), reference, rel_tol=2e-12)


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
    check_bernstein_endpoint()


if __name__ == "__main__":
    main()
