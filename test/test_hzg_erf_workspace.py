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
    for name in ("bern", "exp", "power", "laurent", "expbern", "expcheb"):
        pdf = workspace.pdf(name)
        for mass in (97., 101., 125., 145., 162.):
            x.setVal(mass)
            value = pdf.getVal(obs)
            assert math.isfinite(value) and value > 0, (name, mass, value)
            result.append(value)
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
