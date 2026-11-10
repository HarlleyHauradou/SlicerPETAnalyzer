# SlicerPETAnalyzer

PET Analyzer is an English-language 3D Slicer extension for quantitative image-quality analysis of the six hot spheres of the IEC/NEMA body phantom. It supports calibrated PET and SPECT emission volumes and one- or two-syringe filling. It is intended for research and QA method development.

## Run from source

1. Clone this repository or download its source ZIP.
2. In 3D Slicer, open **Edit → Application Settings → Modules**. Add the repository's `PETAnalyzer` directory under **Additional module paths**, then restart Slicer. Alternatively, launch Slicer with `--additional-module-paths /path/to/SlicerPETAnalyzer/PETAnalyzer`.
3. Open **Quantification → PET Analyzer**. Slicer supplies the runtime dependencies, including NumPy and SciPy.

## Analyze a phantom

1. Load the reconstructed emission series with Slicer's DICOM module. Select the emission volume in PET Analyzer. An already registered CT can be used as a visual background.
2. Choose the filling method. Enter net activities in MBq and the final liquid volumes in mL. For **one syringe**, the first activity and volume describe the homogeneous stock. The software subtracts 47.838 mL for the six spheres and any unrecovered stock before computing background activity. For **two syringes**, the second activity is the net activity transferred to the background. Correct both assays to the same reference before entering them. Set the image reference delay and physical half-life explicitly.
3. Confirm the voxel units and decay reference. Bq/mL is converted to kBq/mL. If starting from other units, provide a validated multiplicative factor to kBq/mL for the selected image and calibration.
4. Select the sphere plane and run detection. Confirm the size labels against the emission/CT views. If necessary, reposition the points in Markups or create six manual points in ascending diameter order. Automatic size assignment expects an IEC lid with sequential sphere sizes around the ring.
5. Place a single background point inside homogeneous liquid. Optionally place a lung insert point and a separate list of 12 planar background points. Review the radius, wall clearance, lung insert clearance and the five axial sampling planes.
6. Measure, review the overlays and result table, then export. The package contains printable `report.html`, `recovery.svg`, `contrast.svg`, `spheres.csv` and `analysis.json`. Use the browser's Print dialog to save the report as PDF.

### One-syringe worked example

With 33.3 MBq in 1,270 mL, no unrecovered stock, a final background volume of 9,800 mL and zero decay delay, the true sphere and background concentrations are **26.2205** and **3.2700 kBq/mL** respectively. Substitute the actual final tank volume and timing for your scan.

## Metrics and limits

- Nominal 3D mean, standard deviation, maximum and nominal-volume contrast recovery.
- Background-adapted A50 connected mean and volume. A truncated or unavailable contour is recorded as missing.
- Fractional-occupancy 1 mL spherical peak with candidate centers constrained to the nominal sphere.
- Background calibration bias and within-VOI coefficient of variation, plus optional lung residual.
- Optional planar contrast and background variability with 12 locations and five native slices closest to 0, ±10 and ±20 mm.
- Published FDG EARL1 and EARL2 comparison intervals. EARL2 peak values are historical and under revision in the source publication. Results are visual and numerical comparisons. The implementation does not replace official accreditation or the full NEMA NU 2 protocol.

See [specification](docs/specification.md), [validation](docs/validation.md) and [release preparation](docs/release.md). Source is Apache-2.0.
