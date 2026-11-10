"""3D Slicer workflow for quantitative IEC hot-sphere phantom analysis."""

from datetime import datetime, timezone
from pathlib import Path

import ctk
import numpy as np
import qt
import slicer
import vtk
from slicer.ScriptedLoadableModule import ScriptedLoadableModule, ScriptedLoadableModuleTest, ScriptedLoadableModuleWidget

from PETAnalyzerLib.analysis import DIAMETERS, concentrations, detect, measure, planar_background
from PETAnalyzerLib.report import annotate_bands, export, plot


class PETAnalyzer(ScriptedLoadableModule):
    def __init__(self, parent):
        super().__init__(parent)
        parent.title = "PET Analyzer"
        parent.icon = qt.QIcon(str(Path(__file__).parent / "Resources" / "PETAnalyzer.png"))
        parent.categories = ["Quantification"]
        parent.contributors = ["Harlley Hauradou"]
        parent.helpText = ("Analyze quantitative IEC hot-sphere phantom images. "
                           "For research and scanner quality assurance development.")


class PETAnalyzerWidget(ScriptedLoadableModuleWidget):
    def setup(self):
        super().setup()
        self.review = None
        self.overlay = None
        self.forms = []
        self.section("1. Image and phantom preparation")
        self.method = qt.QComboBox()
        self.method.addItems(["One syringe · shared stock", "Two syringes · separate solutions"])
        self.form.addRow("Filling method", self.method)
        self.volume = self.selector("Emission volume", "vtkMRMLScalarVolumeNode")
        self.ct = self.selector("Registered CT for review (optional)", "vtkMRMLScalarVolumeNode")
        self.units = qt.QComboBox()
        self.units.addItems(["Bq/mL", "kBq/mL", "Other · apply calibration factor"])
        self.form.addRow("Loaded voxel units", self.units)
        self.factor = self.number("Calibration factor to kBq/mL", 1, .00000001, 1e8, 8)
        self.source = self.number("Net sphere stock activity (MBq)", 0)
        self.stock = self.number("Sphere stock volume (mL)", 0)
        self.tank = self.number("Final background volume (mL)", 0)
        self.loss = self.number("Unrecovered stock volume (mL)", 0)
        self.bg_activity = self.number("Net background activity (MBq)", 0)
        self.elapsed = self.number("Image reference minus assay (min)", 0, -1000000, 1000000)
        self.half_life = self.number("Physical half-life (min)", 109.77, .001, 1e9)
        self.bands = qt.QComboBox()
        self.bands.addItems(["None · general PET or quantitative SPECT", "EARL1", "EARL2"])
        self.form.addRow("FDG PET comparison profile", self.bands)
        self.note("Published comparison bands are for FDG PET only. EARL2 peak uses historical reference values.")
        self.reference = qt.QCheckBox("I verified image units, geometry and common assay reference")
        self.form.addRow(self.reference)
        self.concentration_label = self.note("")
        self.method.currentIndexChanged.connect(self.refresh)
        self.units.currentIndexChanged.connect(self.refresh)
        for field in (self.source, self.stock, self.tank, self.loss, self.bg_activity, self.elapsed, self.half_life):
            field.valueChanged.connect(self.refresh)

        self.section("2. Locate spheres and reference regions")
        self.note("Set the slab center near the sphere plane in native emission slices, then detect. "
                  "Review the diameter labels in the slice views. Use Edit Markups to reposition or reorder points.")
        self.axial = qt.QSpinBox()
        self.axial.maximum = 100000
        self.form.addRow("Sphere plane K index (zero-based)", self.axial)
        self.volume.currentNodeChanged.connect(self.on_volume)
        self.button("Use center of current axial view", self.use_axial_view)
        self.spheres = self.selector("Sphere centers", "vtkMRMLMarkupsFiducialNode")
        self.button("Detect six sphere centers", self.auto)
        self.button("Create manual center list", self.manual)
        self.note("For manual placement, add six points in order: 10, 13, 17, 22, 28, 37 mm. "
                  "Place one background point in uniform liquid and optionally one point within the lung insert.")
        self.background = self.selector("Background center", "vtkMRMLMarkupsFiducialNode")
        self.bg_radius = self.number("Background VOI radius (mm)", 15, 1, 100)
        self.lung = self.selector("Lung center (optional)", "vtkMRMLMarkupsFiducialNode")
        self.lung_radius = self.number("Lung VOI radius (mm)", 10, 1, 100)
        self.planar = self.selector("Planar background (optional, 12 points)", "vtkMRMLMarkupsFiducialNode")
        self.note("For planar contrast and background variability, place 12 well-spaced points "
                  "in uniform liquid on the sphere plane. Each must accommodate a 37 mm circle.")
        self.button("Edit Markups", lambda: slicer.util.selectModule("Markups"))

        self.section("3. Review and export")
        self.button("Measure and show overlays", self.run)
        self.summary = self.note("No measurements yet.")
        self.table = qt.QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(["mm", "Nominal", "A50", "Max", "Peak", "Contrast (%)"])
        self.table.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(qt.QHeaderView.Stretch)
        self.table.setMinimumHeight(220)
        self.form.addRow(self.table)
        self.chart = qt.QLabel()
        self.form.addRow(self.chart)
        self.confirm = qt.QCheckBox("I reviewed the overlay, diameter labels and ROI placement")
        self.form.addRow(self.confirm)
        self.export_button = self.button("Export report and data", self.save)
        self.export_button.enabled = False
        self.confirm.toggled.connect(self.update_export)
        self.layout.addStretch(1)
        self.refresh()
        self.on_volume()

    def section(self, title):
        panel = ctk.ctkCollapsibleButton()
        panel.text = title
        self.layout.addWidget(panel)
        self.form = qt.QFormLayout(panel)

    def selector(self, title, node_type):
        selector = slicer.qMRMLNodeComboBox()
        selector.nodeTypes = [node_type]
        selector.noneEnabled = True
        selector.addEnabled = node_type == "vtkMRMLMarkupsFiducialNode"
        selector.removeEnabled = False
        selector.setMRMLScene(slicer.mrmlScene)
        self.form.addRow(title, selector)
        if selector.addEnabled:
            placement = slicer.qSlicerMarkupsPlaceWidget()
            placement.setMRMLScene(slicer.mrmlScene)
            selector.currentNodeChanged.connect(placement.setCurrentNode)
            self.form.addRow("", placement)
        return selector

    def note(self, text):
        label = qt.QLabel(text)
        label.wordWrap = True
        self.form.addRow(label)
        return label

    def number(self, title, default, minimum=0, maximum=1e7, decimals=4):
        field = qt.QDoubleSpinBox()
        field.setDecimals(decimals)
        field.setRange(minimum, maximum)
        field.value = default
        self.form.addRow(title, field)
        return field

    def button(self, title, callback):
        button = qt.QPushButton(title)
        self.form.addRow(button)
        def call():
            with slicer.util.tryWithErrorDisplay("PET Analyzer could not complete the operation.", waitCursor=True):
                callback()
        button.clicked.connect(call)
        return button

    def preparation(self):
        return dict(method="one" if self.method.currentIndex == 0 else "two",
                    source_mbq=self.source.value, source_ml=self.stock.value, tank_ml=self.tank.value,
                    background_mbq=self.bg_activity.value if self.method.currentIndex else 0,
                    lost_stock_ml=self.loss.value if not self.method.currentIndex else 0,
                    elapsed_minutes=self.elapsed.value, half_life_minutes=self.half_life.value)

    def refresh(self, *_):
        one = self.method.currentIndex == 0
        self.loss.enabled = one
        self.bg_activity.enabled = not one
        self.factor.enabled = self.units.currentIndex == 2
        try:
            hot, bg = concentrations(**self.preparation())
            self.concentration_label.text = f"True sphere: {hot:.3f} kBq/mL · Background: {bg:.3f} kBq/mL · Ratio: {hot/bg:.2f}:1"
        except ValueError as error:
            self.concentration_label.text = str(error)

    def on_volume(self, *_):
        node = self.volume.currentNode()
        if node and node.GetImageData():
            self.axial.maximum = node.GetImageData().GetDimensions()[2] - 1
            self.axial.value = self.axial.maximum // 2

    def use_axial_view(self):
        _, _, _, matrix = self.image()
        red = slicer.app.layoutManager().sliceWidget("Red").mrmlSliceNode()
        width, height, _ = red.GetDimensions()
        point = red.GetXYToRAS().MultiplyPoint((width / 2, height / 2, 0, 1))
        inverse = vtk.vtkMatrix4x4()
        vtk.vtkMatrix4x4.Invert(matrix, inverse)
        k = round(inverse.MultiplyPoint(point)[2])
        if not 0 <= k <= self.axial.maximum:
            raise ValueError("The current axial view lies outside the emission image.")
        self.axial.value = k

    def image(self):
        node = self.volume.currentNode()
        if not node or not node.GetImageData() or node.GetParentTransformNode():
            raise ValueError("Select an emission volume with image data and no parent transform. Harden transforms first.")
        ijk_to_ras = vtk.vtkMatrix4x4()
        node.GetIJKToRASMatrix(ijk_to_ras)
        direction = np.array([[ijk_to_ras.GetElement(i, j) for j in range(3)] for i in range(3)])
        direction /= np.linalg.norm(direction, axis=0)
        if not np.allclose(direction.T @ direction, np.eye(3), atol=1e-4):
            raise ValueError("Resample the volume to an orthogonal grid before measuring.")
        return node, slicer.util.arrayFromVolume(node), np.array(node.GetSpacing()[::-1]), ijk_to_ras

    def manual(self):
        node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode", "Manual sphere centers (ascending diameter)")
        self.spheres.setCurrentNode(node)
        self.review = None
        self.confirm.checked = False
        self.update_export()

    def auto(self):
        volume, image, spacing, matrix = self.image()
        centers = detect(image, spacing, self.axial.value)
        node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode", "Detected sphere centers")
        for d, kji in centers.items():
            ras = matrix.MultiplyPoint((*map(float, kji[::-1]), 1))[:3]
            node.AddControlPoint(vtk.vtkVector3d(*ras), f"{d} mm")
        self.spheres.setCurrentNode(node)
        self.review = None
        self.confirm.checked = False
        self.update_export()
        self.show_volume(volume)

    def show_volume(self, volume):
        ct = self.ct.currentNode()
        if ct and ct != volume:
            slicer.util.setSliceViewerLayers(background=ct, foreground=volume, foregroundOpacity=.4)
        else:
            slicer.util.setSliceViewerLayers(background=volume, foreground=None)

    def points(self, selector, count, matrix, optional=False):
        node = selector.currentNode()
        if node is None and optional:
            return None
        if node is None or node.GetNumberOfControlPoints() != count or node.GetNumberOfDefinedControlPoints() != count:
            raise ValueError(f"Select a markup with exactly {count} defined point(s).")
        inverse = vtk.vtkMatrix4x4()
        vtk.vtkMatrix4x4.Invert(matrix, inverse)
        points = []
        for index in range(count):
            ras = [0.0, 0.0, 0.0]
            node.GetNthControlPointPositionWorld(index, ras)
            points.append(list(inverse.MultiplyPoint((*ras, 1))[2::-1]))
        return points

    def inputs(self):
        node, image, spacing, matrix = self.image()
        if not self.reference.checked:
            raise ValueError("Verify image units, geometry and the common assay reference.")
        preparation = self.preparation()
        hot, bg = concentrations(**preparation)
        centers = self.points(self.spheres, 6, matrix)
        labels = [self.spheres.currentNode().GetNthControlPointLabel(i) for i in range(6)]
        known = {f"{d} mm": d for d in DIAMETERS}
        if any(label in known for label in labels):
            if set(labels) != set(known):
                raise ValueError("Use six unique diameter labels or six manual points in ascending diameter order.")
            centers = {known[label]: point for label, point in zip(labels, centers)}
        else:
            centers = dict(zip(DIAMETERS, centers))
        background = self.points(self.background, 1, matrix)[0]
        lung = self.points(self.lung, 1, matrix, optional=True)
        planar = self.points(self.planar, 12, matrix, optional=True)
        metadata = dict(preparation=preparation, true_sphere_kbq_ml=hot, true_background_kbq_ml=bg,
                        image_units=self.units.currentText, calibration_to_kbq_ml=(.001, 1, self.factor.value)[self.units.currentIndex],
                        spacing_kji_mm=spacing.tolist(), image_shape_kji=list(image.shape),
                        ijk_to_ras=[[matrix.GetElement(i, j) for j in range(4)] for i in range(4)],
                        sphere_centers_kji=centers, background_center_kji=background,
                        background_radius_mm=self.bg_radius.value, lung_center_kji=None if lung is None else lung[0],
                        lung_radius_mm=self.lung_radius.value, planar_background_kji=planar, comparison=self.bands.currentText if self.bands.currentIndex else None)
        return node, image, spacing, metadata

    def run(self):
        self.review = None
        self.confirm.checked = False
        self.update_export()
        node, image, spacing, p = self.inputs()
        result, overlays = measure(image * p["calibration_to_kbq_ml"], spacing, p["sphere_centers_kji"],
                                   p["background_center_kji"], p["background_radius_mm"],
                                   (p["true_sphere_kbq_ml"], p["true_background_kbq_ml"]),
                                   p["lung_center_kji"], p["lung_radius_mm"])
        if p["planar_background_kji"] is not None:
            metrics, offsets = planar_background(image * p["calibration_to_kbq_ml"], spacing,
                                                 p["planar_background_kji"], p["sphere_centers_kji"],
                                                 (p["true_sphere_kbq_ml"], p["true_background_kbq_ml"]))
            result["planar_offsets_mm"] = offsets
            for row in result["spheres"]:
                row.update(metrics[row["diameter_mm"]])
        if p["comparison"]:
            annotate_bands(result, p["comparison"])
        self.draw_overlays(node, image.shape, overlays)
        pixmap = qt.QPixmap()
        pixmap.loadFromData(qt.QByteArray(plot(result, p["comparison"]).encode("utf-8")), "SVG")
        self.chart.setPixmap(pixmap.scaledToWidth(430, qt.Qt.SmoothTransformation))
        self.table.setRowCount(len(result["spheres"]))
        for i, row in enumerate(result["spheres"]):
            for j, key in enumerate(("diameter_mm", "rc_nominal_mean", "rc_a50", "rc_max",
                                     "rc_peak", "contrast_recovery_pct")):
                value = row[key]
                label = "—" if value is None else f"{value:.3f}" if j else str(value)
                item = qt.QTableWidgetItem(label)
                status = row.get(key + "_reference_comparison")
                if status:
                    item.setToolTip(f"{p['comparison']} published band: {status}")
                    if status == "Outside":
                        item.setBackground(qt.QBrush(qt.QColor("#ffcbb0")))
                self.table.setItem(i, j, item)
        self.review = (result, p, node.GetID(), node.GetImageData().GetMTime())
        self.summary.text = (f"Measured background: {result['background_mean_kbq_ml']:.3f} kBq/mL. "
                             "Review the segmentation in the slice views before export.")

    def draw_overlays(self, volume, shape, overlays):
        scene = slicer.mrmlScene
        if self.overlay and self.overlay.GetScene():
            scene.RemoveNode(self.overlay)
        self.overlay = scene.AddNewNodeByClass("vtkMRMLSegmentationNode", "PET Analyzer measurement VOIs")
        self.overlay.CreateDefaultDisplayNodes()
        self.overlay.SetReferenceImageGeometryParameterFromVolumeNode(volume)
        for name, (bounds, local_mask) in overlays.items():
            label = np.zeros(shape, np.uint8)
            label[bounds] = local_mask
            segment = self.overlay.GetSegmentation().AddEmptySegment("", name)
            slicer.util.updateSegmentBinaryLabelmapFromArray(label, self.overlay, segment, volume)
            if name.endswith("A50"):
                self.overlay.GetDisplayNode().SetSegmentVisibility(segment, False)
        self.show_volume(volume)

    def update_export(self, *_):
        self.export_button.enabled = bool(self.review and self.confirm.checked)

    def save(self):
        result, saved, node_id, modified = self.review
        node, _, _, current = self.inputs()
        if current != saved or node_id != node.GetID() or modified != node.GetImageData().GetMTime():
            self.review = None
            self.confirm.checked = False
            raise ValueError("The analysis inputs changed. Measure again and review the overlay before exporting.")
        parent = qt.QFileDialog.getExistingDirectory(None, "Select the output parent directory")
        if parent:
            time = datetime.now(timezone.utc)
            metadata = dict(saved, exported_utc=time.isoformat(), slicer_version=slicer.app.applicationVersion)
            path = export(Path(parent) / time.strftime("pet-analyzer-%Y%m%dT%H%M%S%fZ"),
                          result, metadata, metadata["comparison"])
            self.summary.text = f"Report saved to {path}"
            qt.QDesktopServices.openUrl(qt.QUrl.fromLocalFile(str(path)))


class PETAnalyzerTest(ScriptedLoadableModuleTest):
    def runTest(self):
        hot, bg = concentrations("one", 33.3, 1270, 9800)
        self.assertAlmostEqual(hot, 26.22047, places=4)
        self.assertAlmostEqual(bg, 3.26997, places=4)
