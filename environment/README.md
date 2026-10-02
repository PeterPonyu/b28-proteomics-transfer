# Environment boundaries

No dependency installation was performed for this export. Version files are
evidence/specifications, not promises that an arbitrary environment is identical.

`fit-recorded.json` / `requirements-fit.txt` describe the **recorded corrected
model-fit environment**: Python 3.13.5, CPTAC 1.5.14, XGBoost 3.3.0, NumPy 2.5.0,
pandas 3.0.3, scikit-learn 1.9.0, and SciPy 1.18.0. The first six package/version
facts come from existing fit/environment provenance; no fit was run here.

`verification-observed.json` / `requirements-offline.txt` describe the **actual
offline export-verifier environment** (including NumPy 2.2.6, pandas 2.3.3,
scikit-learn 1.8.0, pytest 8.4.2 and PyMuPDF 1.27.1). Tests verify stored patient
folds by identity; they do not regenerate or relabel the 1.9.0 GroupKFold outputs.
CPTAC is not required or imported for default offline checks. XGBoost is imported
by the actual model-code unit tests, but those tests do not fit a model.

PDF/figure verification additionally needs the external Poppler `pdftoppm`
utility. An explicit `--build` needs `Rscript`, `latexmk`, `pdflatex`, BibTeX, a
working TeX Live installation including `achemso`, `newtxtext/newtxmath`, and the
packages used by the manuscript, plus separately licensed Arial regular and
bold. No class/font/package binaries are vendored. Missing exact Arial faces
fail closed. Supplied PDFs/PNGs can be inspected without installing those fonts.

Observed R renderer dependencies at export were ggplot2 4.0.3, dplyr 1.2.1,
tidyr 1.3.2, systemfonts 1.3.2, ragg 1.5.2, scales 1.4.0, cowplot 1.2.0 and
gtable 0.3.6. These renderer/runtime versions are not original-model dependencies.

Default known reproduction scripts write only disposable verification/raster
output. Optional explicit table/figure/TeX builds write inside this repository.
Neither default nor explicit build installs anything, downloads assay data, or
fits a model. The recorded scientific evidence and public lock cannot be silently
rebound by a successful build.
