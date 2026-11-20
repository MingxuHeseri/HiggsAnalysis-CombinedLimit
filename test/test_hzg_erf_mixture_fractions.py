"""Check explicit and legacy mixture fractions against independent power/exponential tails."""
import math
import ROOT


def tail_reference(family, lo, hi):
    """Integrals of two unit-area components on [1, 100]."""
    if family == "power":
        return ((lo ** -9 - hi ** -9) / (1 - 100 ** -9),
                math.log(hi / lo) / math.log(100))
    return tuple(math.exp(k * (lo - 1)) * math.expm1(k * (hi - lo)) / math.expm1(k * 99)
                 for k in (-0.5, -0.001))


def point_reference(family, x):
    if family == "power":
        return 9 * x ** -10 / (1 - 100 ** -9), 1 / (x * math.log(100))
    return tuple(k * math.exp(k * (x - 1)) / math.expm1(k * 99) for k in (-0.5, -0.001))


def check_mixture(family):
    x = ROOT.RooRealVar("mass", "", 50., 1., 100.)
    obs = ROOT.RooArgSet(x)
    # The CDF equals one in floating point throughout this test interval.
    step = ROOT.RooConstVar("step", "", -100.)
    width = ROOT.RooConstVar("width", "", 1.)
    values = (-10., -1.) if family == "power" else (-0.5, -0.001)
    p = ROOT.RooConstVar("p", "", values[0])
    q = ROOT.RooConstVar("q", "", values[1])
    slopes = ROOT.RooArgList(p, q)
    ratio = ROOT.RooRealVar("ratio", "", 0.25, 0., 1.)
    first = ROOT.RooFormulaVar("first", "1/(1+@0)", ROOT.RooArgList(ratio))
    last = ROOT.RooFormulaVar("last", "@0/(1+@0)", ROOT.RooArgList(ratio))
    both = ROOT.RooArgList(first, last)
    legacy = ROOT.RooArgList(first)
    if family == "power":
        pdf = ROOT.RooPowerErfPdf("explicit", "", x, step, width, slopes, both)
        old = ROOT.RooPowerErfPdf("complement", "", x, step, width, slopes, legacy)
    else:
        pdf = ROOT.RooErfProductPdf("explicit", "", x, step, width, slopes, both, 1)
        old = ROOT.RooErfProductPdf("complement", "", x, step, width, slopes, legacy, 1)
    locations = (2., 30., 95.5, 99.5)
    data = ROOT.RooDataSet("data", "", obs)
    for m in locations:
        x.setVal(m)
        data.add(obs)
    x.setRange("high_tail", 95., 100.)
    high = pdf.createIntegral(obs, ROOT.RooFit.Range("high_tail"))
    cpu = pdf.createNLL(data, ROOT.RooFit.EvalBackend("cpu"), ROOT.RooFit.Offset(False))
    for r in (0., 2.**-60, 1., .25, 2.**-60):
        ratio.setVal(r)
        w1 = 1 / (1 + r)
        w2 = r / (1 + r)
        logs = []
        for m in locations:
            x.setVal(m)
            a, b = point_reference(family, m)
            ref = w1 * a + w2 * b
            assert math.isclose(pdf.getVal(obs), ref, rel_tol=1e-11)
            assert math.isclose(old.getVal(obs), w1*a + (1-w1)*b, rel_tol=1e-11)
            logs.append(math.log(ref))
        a, b = tail_reference(family, 95., 100.)
        assert math.isclose(high.getVal(), w1*a + w2*b, rel_tol=1e-10)
        assert math.isclose(cpu.getVal(), -math.fsum(logs), rel_tol=1e-12)


def check_saved_single_power():
    mass = ROOT.RooRealVar("mass", "", 125., 100., 165.)
    p = ROOT.RooConstVar("p", "", -6.576494165546397)
    step = ROOT.RooConstVar("step", "", 101.28335112633306)
    sigma = ROOT.RooConstVar("sigma", "", 1.8332665486942525)
    pdf = ROOT.RooPowerErfPdf("powerf1", "", mass, step, sigma, ROOT.RooArgList(p), ROOT.RooArgList())
    obs = ROOT.RooArgSet(mass)
    mass.setRange("left_sb", 100., 120.)
    mass.setRange("right_sb", 130., 165.)
    probability = sum(pdf.createIntegral(obs, ROOT.RooFit.NormSet(obs), ROOT.RooFit.Range(name)).getVal()
                      for name in ("left_sb", "right_sb"))
    assert abs(probability - 0.8485698430031652) < 1e-10


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    assert ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit") >= 0
    check_mixture("power")
    check_mixture("exponential")
    check_saved_single_power()
