"""Reproduce the external measured PET check using the original Zenodo archive."""

import argparse
from collections import defaultdict
from io import BytesIO
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pydicom

from PETAnalyzerLib.analysis import concentrations, detect, measure
from PETAnalyzerLib.report import export

CHECKSUM = "c3b6c489a5861e5c1c9b1b0f5cd9e31e"
DESCRIPTION = "05 Q.Clear 150 180s"


def load_series(archive_path, dicom_dir=None):
    with archive_path.open("rb") as stream:
        if hashlib.file_digest(stream, "md5").hexdigest() != CHECKSUM:
            raise ValueError("The archive differs from Zenodo DICOM_SBR_4_1.zip.")
    slices = []
    with ZipFile(archive_path) as archive:
        for name in archive.namelist():
            if not name.startswith("05/A/Z"):
                continue
            raw = archive.read(name)
            ds = pydicom.dcmread(BytesIO(raw))
            if ds.SeriesDescription == DESCRIPTION:
                slices.append((ds, raw))
    slices.sort(key=lambda item: float(item[0].ImagePositionPatient[2]))
    ds0 = slices[0][0]
    assert len(slices) == 53 and len({ds.SeriesInstanceUID for ds, _ in slices}) == 1
    assert all(ds.Units == "BQML" and ds.DecayCorrection == "START" and ds.AcquisitionTime == "173203"
               for ds, _ in slices)
    z = [float(ds.ImagePositionPatient[2]) for ds, _ in slices]
    dz = np.diff(z)
    assert np.allclose(dz, np.median(dz), atol=.001)
    spacing = np.array([np.median(dz), *ds0.PixelSpacing], dtype=float)
    bqml = np.stack([ds.pixel_array.astype(float) * float(ds.RescaleSlope) + float(ds.RescaleIntercept)
                    for ds, _ in slices])
    if dicom_dir:
        dicom_dir.mkdir(exist_ok=True)
        for index, (_, raw) in enumerate(slices):
            (dicom_dir / f"slice-{index:03}.dcm").write_bytes(raw)
    return bqml, spacing


def repeatability(archive_path, reference_centers):
    """Check detector repeatability against the first Q.Clear reconstruction."""
    groups = defaultdict(list)
    with ZipFile(archive_path) as archive:
        for name in archive.namelist():
            if "/A/Z" in name:
                ds = pydicom.dcmread(archive.open(name))
                groups[ds.SeriesInstanceUID].append(ds)
    assert len(groups) == 30
    reference = np.array(list(reference_centers.values()))
    largest = 0.0
    for slices in groups.values():
        slices.sort(key=lambda ds: float(ds.ImagePositionPatient[2]))
        positions = [float(ds.ImagePositionPatient[2]) for ds in slices]
        spacing = np.array([np.median(np.diff(positions)), *slices[0].PixelSpacing], dtype=float)
        image = np.stack([ds.pixel_array.astype(float) * float(ds.RescaleSlope) + float(ds.RescaleIntercept)
                          for ds in slices]) / 1000
        centers = detect(image, spacing, 26)
        distances = np.linalg.norm((np.array(list(centers.values())) - reference) * spacing, axis=1)
        largest = max(largest, float(distances.max()))
    assert largest < 4, largest
    return dict(series_count=len(groups), largest_center_shift_mm=largest)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--dicom-dir", type=Path)
    parser.add_argument("--report-dir", type=Path)
    parser.add_argument("--repeatability", action="store_true")
    args = parser.parse_args()
    bqml, spacing = load_series(args.archive, args.dicom_dir)
    image = bqml / 1000
    centers = detect(image, spacing, 26)
    # Manually selected central lung and uniform liquid positions in this series.
    background, lung = (26, 127, 164), (26, 127, 128)
    preparation = dict(method="two", source_mbq=12.5, source_ml=1000,
                       tank_ml=9800, background_mbq=30.38)
    truth = concentrations(**preparation)
    result, _ = measure(image, spacing, centers, background, 10, truth, lung)
    assert abs(result["background_bias_pct"]) < 10
    assert len(result["spheres"]) == 6
    assert all(0 < row["rc_max"] < 3 and row["rc_a50"] is not None for row in result["spheres"])
    provenance = dict(dataset="10.5281/zenodo.3373973", archive_md5=CHECKSUM, series=DESCRIPTION,
                      tracer="18F-fluoride", preparation_equivalent=preparation,
                      true_sphere_kbq_ml=truth[0], true_background_kbq_ml=truth[1],
                      centers_kji={d: c.tolist() for d, c in centers.items()},
                      background_kji=background, lung_kji=lung, roi_radius_mm=10,
                      spacing_kji_mm=spacing.tolist())
    if args.repeatability:
        provenance["repeatability"] = repeatability(args.archive, centers)
    if args.report_dir:
        export(args.report_dir, result, provenance, None)
    print(json.dumps({"provenance": provenance, "results": result}, indent=2))


if __name__ == "__main__":
    main()
