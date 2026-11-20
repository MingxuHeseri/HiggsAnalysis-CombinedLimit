"""Round-trip all six full HZG background families using only CombinedLimit."""
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

import ROOT


def values(workspace):
    x = workspace.var("mass")
    obs = ROOT.RooArgSet(x)
    result = []
    for name in ("bern", "exp", "power", "laurent", "expbern", "expcheb", "explicit_power", "explicit_exp"):
        pdf = workspace.pdf(name)
        if name in ("power", "laurent", "explicit_power"):
            config = pdf.getIntegratorConfig()
            assert config.method1D().getLabel() == "RooAdaptiveGaussKronrodIntegrator1D"
            assert (config.epsAbs(), config.epsRel()) == (1e-50, 1e-10)
        for mass in (97., 101., 125., 145., 162.):
            x.setVal(mass)
            value = pdf.getVal(obs)
            assert math.isfinite(value) and value > 0, (name, mass, value)
            result.append(value)
    for name in ("explicit_power", "explicit_exp"):
        result.append(workspace.function(name + "_last").getVal())
        assert result[-1] == 2.**-60
    return result


def write_workspace(path):
    x = ROOT.RooRealVar("mass", "", 110., 97., 162.)
    step = ROOT.RooRealVar("step", "", 101., 85., 120.)
    sigma = ROOT.RooRealVar("sigma", "", 3., .5, 10.)
    w = ROOT.RooWorkspace("w")
    empty = ROOT.RooArgList()
    def store(pdf):
        getattr(w, "import")(pdf, ROOT.RooFit.RecycleConflictNodes())

    for name, mode, par in (("bern", 0, .6), ("exp", 1, -.02)):
        value = ROOT.RooRealVar(name+"_c", "", par)
        store(ROOT.RooErfProductPdf(name, name, x, step, sigma, ROOT.RooArgList(value), empty, mode))
    for name in ("power", "laurent"):
        floating = name == "power"
        p = ROOT.RooRealVar(name+"_p", "", -4., -15., 0.)
        p.setConstant(not floating)
        q = ROOT.RooConstVar(name+"_q", "", -5.)
        fraction = ROOT.RooRealVar(name+"_f", "", .2, 0., 1.)
        store(ROOT.RooPowerErfPdf(name, name, x, step, sigma, ROOT.RooArgList(p, q), ROOT.RooArgList(fraction)))
    for name, mode in (("expbern", 0), ("expcheb", 1)):
        c = ROOT.RooRealVar(name+"_c", "", -.2)
        store(ROOT.RooExpPolyErfPdf(name, name, x, step, sigma, ROOT.RooArgList(c), 97., 162., mode))
    for name, slopes in (("explicit_power", (-4., -6.)), ("explicit_exp", (-.05, -.001))):
        a, b = [ROOT.RooConstVar(name+f"_p{i}", "", k) for i, k in enumerate(slopes)]
        coeff = ROOT.RooRealVar(name+"_c", "", 2.**-60, 0., 1.)
        total = ROOT.RooFormulaVar(name+"_sum", "1+@0", ROOT.RooArgList(coeff))
        first = ROOT.RooFormulaVar(name+"_first", "1/@0", ROOT.RooArgList(total))
        last = ROOT.RooFormulaVar(name+"_last", "@0/@1", ROOT.RooArgList(coeff,total))
        parts = (x, step, sigma, ROOT.RooArgList(a,b), ROOT.RooArgList(first,last))
        pdf = (ROOT.RooPowerErfPdf(name, name, *parts) if name == "explicit_power"
               else ROOT.RooErfProductPdf(name, name, *parts, 1))
        store(pdf)
    expected = values(w)
    w.writeToFile(str(path/"pdf.root"))
    (path/"values.json").write_text(json.dumps(expected))


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    assert ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit") >= 0
    if len(sys.argv) == 2:
        directory = Path(sys.argv[1])
        source = ROOT.TFile.Open(str(directory/"pdf.root"))
        assert values(source.Get("w")) == json.loads((directory/"values.json").read_text())
    else:
        with tempfile.TemporaryDirectory() as tmp:
            write_workspace(Path(tmp))
            subprocess.run([sys.executable, __file__, tmp], check=True)
