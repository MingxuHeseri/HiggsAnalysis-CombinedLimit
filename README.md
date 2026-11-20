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

This branch includes the integral-reuse and power-erf integration updates.
`RooErfProductPdf` reuses
Bernstein integration moments and the exponential components' integral
numerators and denominators. The bounded caches use the numerical inputs
of each calculation and are local to each thread. Newly constructed
`RooPowerErfPdf` objects select finite-range adaptive Gauss-Kronrod integration
with absolute and relative tolerances of `1e-50` and `1e-10`.

Exponential, PowerLaw and Laurent mixtures accept either all `n` normalized
fractions or the legacy `n-1` fractions with a final complement. To retain tiny
components, pass every weight explicitly as `c_i / sum(c_j)`; the PDF does not
normalize supplied fractions for the caller. Existing workspaces retain their
stored fraction convention. FinalFit factories must pass all fractions to use
this improvement.

Bernstein point evaluation uses the original Bernstein basis to avoid
cancellation near the upper endpoint. Narrow intervals above the turn-on use
fixed Gauss-Legendre quadrature; saturated intervals use a rule exact for the
supported polynomial degrees. Other intervals retain the analytic moments.

After building, run the numerical and fresh-process workspace checks:

```bash
python3 test/test_hzg_erf_integrals.py
python3 test/test_hzg_erf_mixture_fractions.py
python3 test/test_hzg_erf_workspace.py
```

These check independent integral references, explicit and legacy fractions,
CPU NLL evaluation, and serialization. They do not establish analysis-specific
bias or coverage acceptance. The numerical updates were checked with ROOT
6.36.07 on CERN Linux and ROOT 6.36.08 on macOS; use a compatible build environment.
