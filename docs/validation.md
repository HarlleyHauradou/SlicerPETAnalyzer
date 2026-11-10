# Validation

## Analytical cases

Run `PYTHONPATH=PETAnalyzer python -m unittest discover -s tests -v` from the repository root after installing NumPy and SciPy. The tests cover one-syringe activity conservation and decay, two-syringe truth, invalid preparations, exact voxelized spheres, blurred detection, blank data, ROI overlap, five-plane sampling and report export.

## Measured PET DICOM

The original archive [DICOM_SBR_4_1.zip](https://zenodo.org/records/3373973) has MD5 `c3b6c489a5861e5c1c9b1b0f5cd9e31e`. The measured series is `05 Q.Clear 150 180s`, 53 slices, with `BQML` units and `START` decay correction. The [article methods](https://doi.org/10.1186/s40658-020-0270-y) report 12.5 and 3.1 kBq/mL at the first acquisition for the nominal 4:1 case. The chemical tracer is fluorine-18 fluoride. FDG comparison bands are disabled for this measured check. The benchmark enters an equivalent pair of activity/volume values to represent the published truth, since the original syringe volumes are not part of the benchmark.

```sh
python -m pip install numpy scipy pydicom
PYTHONPATH=PETAnalyzer python tests/validate_measured_pet.py DICOM_SBR_4_1.zip \
  --dicom-dir /tmp/pet-analyzer-dicom --report-dir /tmp/pet-analyzer-report --repeatability
```

This check verifies scaling, plausible recovery, automatic six-point localization and report export. The independently placed liquid ROI measures approximately 2.972 kBq/mL (−4.12% from 3.1 kBq/mL). Across the 30 reconstructions in this archive, automatic localization succeeded for all series. The largest per-label shift was 2.10 mm relative to the first Q.Clear series. The fixture lacks an independent CT-based segmentation reference, so the shift assesses consistency rather than positional accuracy. The measured test verifies pipeline operation and is not an accuracy study.

[Example report](examples/measured-pet-report.html)

## Slicer integration

From the root, set `PYTHONPATH` to the absolute `tests` directory and invoke:

```sh
Slicer --no-splash --no-main-window --additional-module-paths "$PWD/PETAnalyzer" \
  --python-script "$PWD/tests/slicer_smoke.py"
```

Run the same command with `PET_ANALYZER_DICOM=/tmp/pet-analyzer-dicom` for measured DICOM import. The smoke test checks Slicer loading, widget setup, point creation, overlay masks, chart rendering, results table and report generation.

## Scientific references

- EFOMP, [Quality controls in PET/CT and PET/MR](https://www.efomp.org/uploads/10fd24d1-354d-48d6-bf8f-e93e98cb7d81/EFOMP%E2%80%99S%20GUIDELINE%20QUALITY%20CONTROLS%20IN%20PETCT%20AND%20PETMR.pdf).
- EARL, [PET/CT accreditation manual v4.2](https://earl.eanm.org/wp-content/uploads/2023/05/EARL_Manual_4.2.pdf).
- Prenosil et al., [EARL compliance measurements on the Biograph Vision Quadra PET/CT](https://pmc.ncbi.nlm.nih.gov/articles/PMC8994003/), Table 1, historical bands.
