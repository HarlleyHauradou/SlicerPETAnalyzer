"""Physical-space phantom measurements. Arrays and coordinates use KJI order."""

import itertools
import math

import numpy as np
from scipy import ndimage
from scipy.optimize import linear_sum_assignment

DIAMETERS = (10, 13, 17, 22, 28, 37)
VOLUMES = {d: math.pi * d**3 / 6000 for d in DIAMETERS}

# Historical published FDG reference intervals (PMC8994003, Table 1).
# The EARL2 peak intervals in this table were under revision.
BANDS = {
    "EARL1": dict(zip(DIAMETERS, (
        ((.27, .43), (.34, .57)), ((.44, .60), (.59, .85)),
        ((.57, .73), (.73, 1.01)), ((.63, .78), (.83, 1.09)),
        ((.72, .85), (.91, 1.13)), ((.76, .89), (.95, 1.16)),
    ))),
    "EARL2": dict(zip(DIAMETERS, (
        ((.39, .61), (.52, .88), (.27, .41)), ((.63, .86), (.85, 1.22), (.45, .70)),
        ((.76, .97), (1.00, 1.38), (.75, .99)), ((.80, .99), (1.01, 1.32), (.90, 1.10)),
        ((.82, .97), (1.01, 1.26), (.90, 1.10)), ((.85, 1.00), (1.05, 1.29), (.90, 1.10)),
    ))),
}


def concentrations(method, source_mbq, source_ml, tank_ml, *, background_mbq=0,
                   lost_stock_ml=0, elapsed_minutes=0, half_life_minutes=109.77):
    """Return true concentrations in kBq/mL at the image reference time."""
    values = (source_mbq, source_ml, tank_ml, background_mbq, lost_stock_ml,
              elapsed_minutes, half_life_minutes)
    if not all(map(math.isfinite, values)) or min(source_mbq, source_ml, tank_ml, half_life_minutes) <= 0:
        raise ValueError("Activity, volumes and half-life must be finite and positive.")
    if method == "one":
        transfer = source_ml - sum(VOLUMES.values()) - lost_stock_ml
        if lost_stock_ml < 0 or not 0 < transfer <= tank_ml:
            raise ValueError("Check the stock, unrecovered stock and final background volumes.")
        background_mbq = source_mbq * transfer / source_ml
    elif method == "two":
        if source_ml < sum(VOLUMES.values()) or background_mbq <= 0:
            raise ValueError("Enter a sufficient sphere dilution volume and positive background activity.")
    else:
        raise ValueError("Unknown filling method.")
    try:
        factor = 1000 * 2 ** (-elapsed_minutes / half_life_minutes)
    except OverflowError as error:
        raise ValueError("The decay interval is too large for this half-life.") from error
    if not math.isfinite(factor) or factor <= 0:
        raise ValueError("The decay interval is too large for this half-life.")
    hot, background = factor * source_mbq / source_ml, factor * background_mbq / tank_ml
    if not 0 < background < hot:
        raise ValueError("This analysis requires hot spheres and a positive, lower background concentration.")
    return hot, background


def _region(shape, spacing, center, radius):
    center = np.asarray(center)
    extent = radius / spacing
    if np.any(center - extent < -.5) or np.any(center + extent > np.asarray(shape) - .5):
        raise ValueError("A measurement ROI extends outside the emission image.")
    lo = np.maximum(0, np.floor(center - extent).astype(int))
    hi = np.minimum(shape, np.ceil(center + extent).astype(int) + 1)
    bounds = tuple(slice(a, b) for a, b in zip(lo, hi))
    distances = np.ogrid[bounds]
    squared = sum(((axis - c) * step)**2 for axis, c, step in zip(distances, center, spacing))
    return bounds, squared <= radius**2


def detect(image, spacing, axial_index, slab_mm=36):
    """Fit six approximately equally spaced maxima on a common circle."""
    spacing = np.asarray(spacing, dtype=float)
    image = np.asarray(image)
    if image.ndim != 3 or spacing.shape != (3,) or not np.all(np.isfinite(spacing)) or np.any(spacing <= 0):
        raise ValueError("Expected a three-dimensional image with positive spacing.")
    if not 0 <= axial_index < image.shape[0]:
        raise ValueError("The slab center lies outside the emission image.")
    pad = math.ceil(slab_mm / (2 * spacing[0]))
    start, stop = max(0, round(axial_index) - pad), min(image.shape[0], round(axial_index) + pad + 1)
    slab = ndimage.gaussian_filter(image[start:stop].astype(float), 2 / spacing)
    footprint = tuple(2 * max(1, round(9 / s)) + 1 for s in spacing)
    peaks = (slab == ndimage.maximum_filter(slab, size=footprint)) & (slab > np.percentile(slab, 90))
    candidates = np.argwhere(peaks)
    candidates = candidates[np.argsort(slab[tuple(candidates.T)])[::-1]][:14]
    if len(candidates) < 6:
        raise ValueError("Fewer than six candidates found. Use manual centers.")
    candidates[:, 0] += start
    points = candidates * spacing
    best = None
    for indices in itertools.combinations(range(len(points)), 3):
        triangle = points[list(indices), 1:]
        matrix = 2 * (triangle[1:] - triangle[0])
        if abs(np.linalg.det(matrix)) < 100:
            continue
        center = np.linalg.solve(matrix, np.sum(triangle[1:]**2 - triangle[0]**2, axis=1))
        radius = np.linalg.norm(triangle[0] - center)
        if not 45 <= radius <= 75:
            continue
        angle = math.atan2(*(triangle[0] - center)[::-1]) + np.arange(6) * math.pi / 3
        model = center + radius * np.column_stack((np.cos(angle), np.sin(angle)))
        distance = np.linalg.norm(model[:, None] - points[None, :, 1:], axis=2)
        _, matched = linear_sum_assignment(distance)
        errors = distance[np.arange(6), matched]
        spread = np.ptp(points[matched, 0])
        score = float(np.mean(errors**2) + spread**2 / 3)
        if np.max(errors) <= 10 and spread <= 14 and (best is None or score < best[0]):
            best = score, matched
    if best is None:
        raise ValueError("No credible six-sphere ring found. Adjust the slab or use manual centers.")
    ordered = candidates[best[1]].astype(float)
    # Score apparent widths in a smoothed local neighborhood. The
    # intensity at the center alone does not order the larger PSF spheres.
    refined, scores = [], []
    for point in ordered:
        bounds, search = _region(image.shape, spacing, point, 26)
        signal = ndimage.gaussian_filter(image[bounds].astype(float), 2 / spacing)
        origin = np.array([s.start for s in bounds])
        peak = signal[tuple((point - origin).astype(int))]
        baseline = np.percentile(signal[search], 20)
        labels, _ = ndimage.label(search & (signal > (peak + baseline) / 2))
        label = labels[tuple((point - origin).astype(int))]
        if not label:
            raise ValueError("A sphere candidate has insufficient local contrast.")
        component = np.argwhere(labels == label) + origin
        centroid = component.mean(axis=0)
        if np.linalg.norm((centroid - point) * spacing) > 14:
            raise ValueError("A sphere candidate is poorly isolated. Use manual centers.")
        refined.append(centroid)
        scores.append(len(component))
    # IEC lid designs place sphere sizes in order around the ring. Consider
    # both directions and all rotations, then rank the width evidence.
    rank = np.argsort(np.argsort(scores))
    sequences = [np.roll(order, offset) for order in (np.arange(6), np.arange(5, -1, -1))
                 for offset in range(6)]
    order = min(sequences, key=lambda seq: np.sum((rank[seq] - np.arange(6))**2))
    return {diameter: refined[index] for diameter, index in zip(DIAMETERS, order)}


def measure(image, spacing, centers, background_center, background_radius, truth,
            lung_center=None, lung_radius=10):
    """Measure positive hot-sphere contrast and return overlay masks and statistics."""
    image = np.asarray(image, dtype=float)
    spacing = np.asarray(spacing, dtype=float)
    if image.ndim != 3 or spacing.shape != (3,) or not np.all(np.isfinite(spacing)) or np.any(spacing <= 0):
        raise ValueError("Expected a three-dimensional image with positive spacing.")
    if not np.all(np.isfinite(image)) or set(centers) != set(DIAMETERS):
        raise ValueError("Expected finite intensities and six labeled sphere centers.")
    hot, background_true = truth
    if not 0 < background_true < hot:
        raise ValueError("Invalid reference concentrations.")
    def extract(point, radius):
        if not np.isfinite(radius) or radius <= 0 or np.shape(point) != (3,) or not np.all(np.isfinite(point)):
            raise ValueError("Each ROI requires a valid center and positive radius.")
        bounds, mask = _region(image.shape, spacing, point, radius)
        if np.count_nonzero(mask) < 2:
            raise ValueError("An ROI is too small for this voxel spacing.")
        return bounds, mask, image[bounds][mask]

    sphere_centers = [(np.asarray(centers[d]), d / 2) for d in DIAMETERS]
    regions = sphere_centers + [(np.asarray(background_center), background_radius)]
    if lung_center is not None:
        regions.append((np.asarray(lung_center), lung_radius))
    for (a, ra), (b, rb) in itertools.combinations(regions, 2):
        if np.linalg.norm((a - b) * spacing) < ra + rb:
            raise ValueError("Measurement ROIs overlap. Review the centers and radii.")

    bounds, mask, samples = extract(background_center, background_radius)
    measured_bg = float(samples.mean())
    if measured_bg <= 0:
        raise ValueError("The measured background must be positive.")
    overlays = {"Background": (bounds, mask)}
    output = {
        "background_mean_kbq_ml": measured_bg,
        "background_bias_pct": 100 * (measured_bg / background_true - 1),
        "background_cv_pct": float(100 * samples.std(ddof=1) / measured_bg),
        "spheres": [],
    }
    if lung_center is not None:
        bounds, mask, samples = extract(lung_center, lung_radius)
        overlays["Lung insert"] = (bounds, mask)
        output["lung_residual_pct"] = float(100 * samples.mean() / measured_bg)
    voxel_ml = np.prod(spacing) / 1000
    peak_radius = (3000 / (4 * math.pi))**(1 / 3)
    peak_extent = np.ceil(peak_radius / spacing).astype(int)
    peak_grid = np.ogrid[tuple(slice(-n, n + 1) for n in peak_extent)]
    # Supersample fractional voxel occupancy of the 1 mL spherical kernel.
    kernel = sum(sum(((axis + offset) * s)**2 for axis, offset, s in zip(peak_grid, offsets, spacing))
                 <= peak_radius**2 for offsets in itertools.product((-0.4, -0.2, 0, 0.2, 0.4), repeat=3))
    kernel = kernel / kernel.sum()
    for d in DIAMETERS:
        center = np.asarray(centers[d])
        bounds, mask, samples = extract(center, d / 2)
        maximum, mean = float(samples.max()), float(samples.mean())
        overlays[f"{d} mm nominal"] = (bounds, mask)
        seed = np.argwhere(mask)[np.argmax(samples)] + np.array([s.start for s in bounds])
        row = dict(diameter_mm=d, nominal_voxels=int(samples.size), nominal_mean_kbq_ml=mean,
                   nominal_sd_kbq_ml=float(samples.std(ddof=1)), max_kbq_ml=maximum,
                   rc_nominal_mean=mean / hot, rc_max=maximum / hot,
                   contrast_recovery_pct=100 * (mean / measured_bg - 1) / (hot / background_true - 1),
                   rc_a50=None, a50_ml=None, rc_peak=None)
        if maximum > measured_bg:
            try:
                crop, local_mask = _region(image.shape, spacing, center, d / 2 + 12)
            except ValueError:
                crop = None
            if crop is not None:
                signal = image[crop]
                labels, _ = ndimage.label((signal >= (maximum + measured_bg) / 2) & local_mask)
                label = labels[tuple(seed - [s.start for s in crop])]
                contour = labels == label
                if label and not np.any(contour & ~ndimage.binary_erosion(local_mask)):
                    row.update(rc_a50=float(signal[contour].mean() / hot), a50_ml=float(contour.sum() * voxel_ml))
                    overlays[f"{d} mm A50"] = (crop, contour)
        indices = np.argwhere(mask) + np.array([s.start for s in bounds])
        valid = np.all((indices >= peak_extent) & (indices + peak_extent < image.shape), axis=1)
        if valid.any():
            lo, hi = indices[valid].min(axis=0), indices[valid].max(axis=0)
            cut = tuple(slice(a - n, b + n + 1) for a, b, n in zip(lo, hi, peak_extent))
            averaged = ndimage.convolve(image[cut], kernel, mode="constant")
            grid = indices[valid] - (lo - peak_extent)
            row["rc_peak"] = float(averaged[tuple(grid.T)].max() / hot)
        output["spheres"].append(row)
    return output, overlays


def planar_background(image, spacing, template_centers, sphere_centers, truth):
    """Sample 12 template locations across the sphere plane and four offset planes.

    The closest native slices to 0, ±10 and ±20 mm are used. Placement of
    the 12 templates in uniform liquid is an operator responsibility.
    """
    spacing = np.asarray(spacing)
    template = np.asarray(template_centers)
    if template.shape != (12, 3) or not np.isfinite(template).all():
        raise ValueError("Provide twelve finite planar background points.")
    center_k = float(np.mean([point[0] for point in sphere_centers.values()]))
    if np.ptp(template[:, 0]) * spacing[0] > spacing[0] + 1e-6:
        raise ValueError("The twelve background points must be in one axial plane.")
    if abs(np.mean(template[:, 0]) - center_k) * spacing[0] > spacing[0] + 1e-6:
        raise ValueError("Place the twelve background points on the sphere plane.")
    planes = np.rint(center_k + np.arange(-20, 21, 10) / spacing[0]).astype(int)
    if len(set(planes)) < 5 or np.any(planes < 0) or np.any(planes >= image.shape[0]):
        raise ValueError("The selected axial range cannot sample five distinct background planes.")
    for a, b in itertools.combinations(template[:, 1:], 2):
        if np.linalg.norm((a - b) * spacing[1:]) < 37:
            raise ValueError("The 37 mm planar background ROIs must not overlap.")
    for point in template:
        for diameter, center in sphere_centers.items():
            if np.linalg.norm((point[1:] - np.asarray(center)[1:]) * spacing[1:]) < (37 + diameter) / 2:
                raise ValueError("A 37 mm planar background ROI overlaps a sphere.")

    hot, bg = truth
    metrics = {}
    for diameter in DIAMETERS:
        means = []
        for k in planes:
            for point in template:
                bounds, mask = _region(image.shape[1:], spacing[1:], point[1:], diameter / 2)
                pixels = image[k][bounds][mask]
                if pixels.size < 2:
                    raise ValueError("A planar ROI is too small for the image sampling.")
                means.append(float(pixels.mean()))
        sphere = np.asarray(sphere_centers[diameter])
        bounds, mask = _region(image.shape[1:], spacing[1:], sphere[1:], diameter / 2)
        pixels = image[round(sphere[0])][bounds][mask]
        if pixels.size < 2:
            raise ValueError("A planar sphere ROI is too small for the image sampling.")
        mean = float(np.mean(means))
        if mean <= 0:
            raise ValueError("The planar background average must be positive.")
        metrics[diameter] = dict(planar_background_kbq_ml=mean,
                                 planar_background_variability_pct=float(100 * np.std(means, ddof=1) / mean),
                                 planar_hot_contrast_pct=float(100 * (pixels.mean() / mean - 1) / (hot / bg - 1)))
    return metrics, ((planes - center_k) * spacing[0]).tolist()
