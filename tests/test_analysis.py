import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

from PETAnalyzerLib.analysis import DIAMETERS, VOLUMES, concentrations, detect, measure, planar_background
from PETAnalyzerLib.report import annotate_bands, export


def fixture(spacing=(2.5, 2, 2)):
    spacing = np.asarray(spacing)
    shape = (43, 100, 100)
    image = np.full(shape, 3.1)
    grid = np.ogrid[tuple(slice(n) for n in shape)]
    centers = {}
    for d, angle in zip(DIAMETERS, np.arange(6) * np.pi / 3):
        center = (21, 50 + 29 * np.cos(angle), 50 + 29 * np.sin(angle))
        centers[d] = center
        mask = sum(((a - c) * s)**2 for a, c, s in zip(grid, center, spacing)) <= (d / 2)**2
        image[mask] = 12.5
    return image, spacing, centers


class PhantomTest(unittest.TestCase):
    def test_one_syringe_conserves_activity_and_decays(self):
        hot, bg = concentrations("one", 33.3, 1270, 9800, lost_stock_ml=10,
                                 elapsed_minutes=109.77)
        self.assertAlmostEqual(sum(VOLUMES.values()), 47.83755, places=4)
        self.assertAlmostEqual(hot, 1000 * 33.3 / 1270 / 2)
        self.assertAlmostEqual((bg * 9800 + hot * (sum(VOLUMES.values()) + 10)) / 1000, 33.3 / 2)

    def test_two_syringes(self):
        self.assertEqual(concentrations("two", 12.5, 1000, 9800, background_mbq=30.38),
                         (12.5, 3.1))

    def test_invalid_preparation(self):
        for args in [("one", 33.3, 40, 9800), ("one", 33.3, 1270, 900),
                     ("two", 33.3, 1270, 9800), ("invalid", 33.3, 1270, 9800)]:
            with self.subTest(args=args), self.assertRaises(ValueError):
                concentrations(*args)

    def test_exact_measurements(self):
        image, spacing, centers = fixture()
        result, overlays = measure(image, spacing, centers, (21, 50, 50), 10, (12.5, 3.1))
        self.assertEqual(len(overlays), 13)
        for row in result["spheres"]:
            self.assertAlmostEqual(row["rc_nominal_mean"], 1)
            self.assertAlmostEqual(row["rc_a50"], 1)
            self.assertAlmostEqual(row["rc_max"], 1)
            self.assertAlmostEqual(row["contrast_recovery_pct"], 100)
            if row["diameter_mm"] >= 17:
                self.assertAlmostEqual(row["rc_peak"], 1)
        self.assertLess(result["spheres"][0]["rc_peak"], 1)
        self.assertAlmostEqual(result["background_bias_pct"], 0)
        self.assertAlmostEqual(result["background_cv_pct"], 0)

    def test_detect_blurred_hot_spheres(self):
        image, spacing, centers = fixture()
        found = detect(gaussian_filter(image, 2.2 / spacing), spacing, 21)
        for d in DIAMETERS:
            self.assertLess(np.linalg.norm((found[d] - centers[d]) * spacing), 2.6)

    def test_blank_image_detection(self):
        image, spacing, _ = fixture()
        with self.assertRaises(ValueError):
            detect(np.zeros_like(image), spacing, 21)

    def test_reject_overlapping_rois(self):
        image, spacing, centers = fixture()
        with self.assertRaisesRegex(ValueError, "overlap"):
            measure(image, spacing, centers, centers[10], 10, (12.5, 3.1))

    def test_planar_background(self):
        image, spacing, centers = fixture()
        image = np.pad(image, ((0, 0), (30, 30), (30, 30)), constant_values=3.1)
        centers = {d: np.asarray(point) + (0, 30, 30) for d, point in centers.items()}
        points = [(21, 80 + 55 * np.cos(angle), 80 + 55 * np.sin(angle))
                  for angle in np.arange(12) * np.pi / 6]
        metrics, offsets = planar_background(image, spacing, points, centers, (12.5, 3.1))
        self.assertEqual(offsets, [-20, -10, 0, 10, 20])
        for row in metrics.values():
            self.assertAlmostEqual(row["planar_hot_contrast_pct"], 100)
            self.assertAlmostEqual(row["planar_background_variability_pct"], 0)
        with self.assertRaisesRegex(ValueError, "overlap"):
            planar_background(image, spacing, [points[0]] * 12, centers, (12.5, 3.1))

    def test_export(self):
        image, spacing, centers = fixture()
        result, _ = measure(image, spacing, centers, (21, 50, 50), 10, (12.5, 3.1))
        annotate_bands(result, "EARL2")
        self.assertEqual(result["spheres"][0]["rc_max_reference_comparison"], "Outside")
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp) / "export"
            report = export(folder, result, {"label": "<unsafe>"}, "EARL2")
            self.assertEqual({p.name for p in folder.iterdir()},
                             {"report.html", "recovery.svg", "contrast.svg", "spheres.csv", "analysis.json"})
            self.assertIn("&lt;unsafe&gt;", report.read_text())
            self.assertEqual(len(json.loads((folder / "analysis.json").read_text())["results"]["spheres"]), 6)
            with self.assertRaises(FileExistsError):
                export(folder, result, {}, "EARL2")


if __name__ == "__main__":
    unittest.main()
