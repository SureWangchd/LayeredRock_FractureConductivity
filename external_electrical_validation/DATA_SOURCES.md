# Sources for the public electrical comparison

The calculation compares the 5 MPa observations in Han et al. (2022), Table 1 (printed p. 1075): **Pressure-dependent joint elastic–electrical properties in brine-saturated artificial sandstones with aligned penny-shaped cracks—Part I: Experimental results**, DOI https://doi.org/10.1093/gji/ggab386. The attributed numerical transcription is `data/Han2022_table1_electrical_porosity.csv`; all 60 rows and 120 conductivity/porosity values match `sources/han2022/table1_independently_extracted.csv`. The original verification result is retained. These are six cores repeatedly measured at ten pressures, not 60 independent cores. Only 5 MPa is modeled here.

I-0 and I-90 calibrate the background. I-45 is an intact holdout; F-0, F-45 and F-90 are fractured targets. Their measured conductivities were not used to tune this model, but the observations were already public: this is a retrospective comparison, not a prospective independent registration. The measurements are at 1 kHz, 22±1 C and 35 g/L NaCl. The steady conduction model compares the reported conductive response; the available data do not establish a DC limit. Conductivity ±2% and porosity ±0.0001 are propagated as deterministic tolerances, not standard deviations.

Geometry and background conventions are documented in Part I and **Yan et al. (2022), Part II: Theoretical modelling**, DOI https://doi.org/10.1093/gji/ggab387. Original manufacturing values are 3 mm crack diameter, 0.063 mm thickness and crack density 0.062. The nominal 40 mm core length, regular centered-cell arrangement and aperture/spacing alternatives are stated modeling assumptions. The full per-input source distinction is in `sources/han2022/PARAMETER_SPECIFICATION_CN.md` and `INPUT_PROTOCOL.json`. The protocol records the original three mesh levels; the actual fourth refinement is documented in `RUN_AMENDMENTS.md`.

The nominal brine conductivity 4.69 S/m is transferred from **Han (2018)** at the same nominal salinity, DOI https://doi.org/10.1002/2018GL077541; it is not a direct measurement for the 2022 brine batch. Han2018 Table S1 was not obtained and is not used as an experimental target. The separate 5.2 S/m scenario is supported in magnitude by official manufacturer conductivity tables and temperature coefficients:

- Hach, *sensION+ DL Portable Meter*, DOC022.98.80155, July 2017, Table 3: https://cdn.hach.com/7FYZVWYB/at/bhtbpxpqn4t39tbsqjm44n4j/DOC0229880155.pdf
- YSI, *Model 3100 Operations Manual*, Appendix B, printed p. 39: https://www.ysi.com/File%20Library/Documents/Manuals/031041-YSI-3100-Operations-Manual-RevB.pdf

The manuals and publisher text excerpts are not redistributed here. Follow the official links for original documents and their applicable terms. No choice between the nominal and sensitivity scenarios was made by minimizing the fractured-target error.

The unused backup dataset **Mohammed and Schmitt (2020)** is available under its original CC BY 4.0 terms at https://doi.org/10.17632/j979d4r7ww.1; the associated paper is https://doi.org/10.1016/j.jappgeo.2020.104197. Those data were retrieved locally but no corresponding model was run, and they are not included in this focused public derivative.

The reported model increment is the model's fractured value minus its own background. A measured fractured value minus a transferred model background is a different, background-dependent estimate. The experiment does not provide each F core's same-body pre-fracture measurement. It therefore does not directly test this project's strict matched-background A/B configuration effect, peripheral-layer rearrangement or modal prediction claims.

Project code is GPL-3.0-only. Project-generated numerical results and Figure S1 data are CC BY 4.0. Third-party measurements retain their original authorship and DOI attribution; no license to redistribute third-party full text is asserted.
