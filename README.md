HiggsAnalysis-CombinedLimit
===========================

### Official documentation

All documentation, including installation instructions, is hosted at
http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest

The source code of this documentation can be found in the `docs/` folder in this repository.

### Publication 

The `Combine` tool publication can be found [here](https://arxiv.org/abs/2404.06614). Please consider citing this reference if you use the `Combine` tool. 

### HZG background PDFs in this fork

Build this branch in CMSSW and include the PDF headers from
`HiggsAnalysis/CombinedLimit/interface/`. The class implementations and ROOT
dictionaries are provided by `libHiggsAnalysisCombinedLimit.so`.

| PDF class | HZG family |
| --- | --- |
| `RooErfProductPdf` | Bernstein × Gaussian CDF (mode 0), exponential mixture × Gaussian CDF (mode 1) |
| `RooPowerErfPdf` | PowerLaw and Laurent mixtures × Gaussian CDF |
| `RooExpPolyErfPdf` | Exponentiated Bernstein (mode 0) and Chebyshev (mode 1) polynomials × Gaussian CDF |

FinalFit constructs the components, fractions, parameter ranges and sideband
mask. Link its factory against this CombinedLimit library and remove the three
class objects and dictionary entries from the factory build. A workspace using
these PDF classes can then be opened directly by Combine. Other custom classes
in a workspace still require their own dictionaries.

Run `python3 test/test_hzg_erf_workspace.py` in the built environment to write
all six families and read their normalized densities in a fresh process.

On the `hzg/erf-integral-reuse` candidate branch, `RooErfProductPdf` reuses
Bernstein integration moments and the exponential components' integral
numerators and denominators. The bounded caches use the numerical inputs
of each calculation and are local to each thread. Run
`python3 test/test_hzg_erf_integrals.py` to compare the resulting integrals
to independent quadrature after changing the geometry and shape parameters.
