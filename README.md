# From Forecast to Position

**Miquel Noguer Alonso**  
Artificial Intelligence Finance Institute (AIFI)

A theory of converting order-book forecasts into causal positions under transaction costs, market impact, constraints, latency, and partial observation.

- DOI: [10.5281/zenodo.22759534](https://doi.org/10.5281/zenodo.22759534)
- Overleaf: [editable project](https://www.overleaf.com/project/6aa88b2f836c121c0176c8c5)
- Manuscript: [`paper.pdf`](paper.pdf)
- LaTeX: [`paper.tex`](paper.tex)

## Repository contents

This private repository is part 3 of the *Market Microstructure Trilogy*. It contains the reviewed manuscript, its LaTeX source, and the corresponding source and verification archive. The manuscript uses author-year citations and includes a table of contents.

## Build

The manuscript was built with pdfLaTeX. For Papers II and III, run BibTeX between LaTeX passes.

```bash
latexmk -pdf paper.tex
```

## Verification status

The released PDF was reproduced from the included source on 15 September 2026. The Overleaf build completed with zero errors and zero warnings. Numerical and symbolic checks are contained in the accompanying verification archive.

## Scope

The guarantees in the paper are conditional on the declared models, information sets, and uncertainty bounds. Synthetic calculations verify the stated identities and certificates; they do not claim live-market profitability.
