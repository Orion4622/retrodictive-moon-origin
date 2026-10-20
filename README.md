# The origin of the Moon without a giant impact

**A retrodictive theory: a secondary circumterrestrial disk redistributed by the retrograde flyby of a rocky planet**

- Michel Debailleul
- geophysicist M.Sc.
- Université libre de Bruxelles
- ORCID 0009-0003-1222-1433]

- https://orcid.org/0009-0003-1222-1433

Version V_1.0 (October 2026) - [DOI 10.5281/zenodo.23391968](https://doi.org/10.5281/zenodo.23391968) - licence CC BY 4.0
Version française : [README_FR.md](README_FR.md)

## The question
The Moon is a datum, and the question is retrodictive: which physical pathways lead to it without a giant impact?

## The pathway examined
A young, hot proto-Earth, soft but never liquid, flattened and inflated by a 6 h rotation, is surrounded by a disk of material that has not yet formed a satellite. A second rocky planet, **planet B**, passes close to the Earth against its rotation, without touching it and without giving it any material: it acts by gravity alone.

- About **47 %** of the disk stays bound to the Earth (14 % for a flyby in the sense of the rotation).
- The bound material assembles into a body of about **one lunar mass** in every realization computed; it is kept in four realizations out of ten, when its orbit settles beyond the corotation radius.
- Since planet B supplies no material, the Moon is made of terrestrial material: the Earth–Moon isotopic identity is expected from the outset.
- The angular momentum budget closes at ±0.05 L_EM in about a quarter to a third of the trajectories explored (±0.08 L_EM in two fifths to two thirds), under two conditions stated in the text: a weak early terrestrial tide and a slow contraction of the proto-Earth. These shares describe the parameter grid, not a probability.

## What is stated as open
The results are necessary conditions or families of solutions, not facts. The retention of the assembled body depends on the merging threshold; the iron sorting of the disk and the rate of such an encounter in an accreting system are not demonstrated. The limits are stated section by section.

## Contents
| Folder | Content |
|---|---|
| [`EN/`](EN/) | English monograph (45 pages): PDF, LaTeX source, figures, `refs.bib` |
| [`FR/`](FR/) | French monograph (46 pages): PDF, LaTeX source, figures, `refs.bib` |
| [`supplementary/massB_scan/`](supplementary/massB_scan/) | Supplementary note and scripts: scan of the mass of planet B (not part of the Zenodo record) |

## Compiling
From the `EN/` or `FR/` folder: `pdflatex` three times with `bibtex` once between the first two passes (or `latexmk -pdf`).

## Note on earlier records
Earlier Zenodo records of the author (September 2026) describe the flyby as prograde with the opposite sense labels. In the present work the sense is defined relative to the rotation of the Earth, and these labels supersede the earlier ones.

## Citation
Debailleul, Michel (2026). *A retrodictive theory of the origin of the Moon: a secondary circumterrestrial disk redistributed by the retrograde flyby of a rocky planet*. Version V_1.0. Zenodo. https://doi.org/10.5281/zenodo.23391968

Objections on the calculations and the stated limits are welcome (GitHub issues).
