# From Forecast to Position

Miquel Noguer Alonso · 28 September 2026  
DOI: https://doi.org/10.5281/zenodo.23027604

## Manuscript

Edit `p3_forecast.tex` and its section inputs. `paper.tex` is the Overleaf entry point;
select pdfLaTeX. `paper.pdf` is the compiled manuscript. Figures, bibliography,
tables, mathematical verification scripts and computational results are included.

## Reproduce

```bash
python -m pip install -r verification/requirements.txt
python scripts/release.py --check
python verification/verify_theory.py
python verification/certified_improvement.py
python verification/realism_stress.py
python deployment/verify_production.py
python deployment/verify_additional.py
python scripts/release.py --build
```

The release script checks references, layout warnings, SHA-256 inventories and
archive contents. The guarantees depend on the stated models and assumptions.
The numerical experiments are synthetic and do not establish market profitability.
