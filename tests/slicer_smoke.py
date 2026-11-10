"""Run inside Slicer with --python-script. Optional PET_ANALYZER_DICOM selects the measured fixture."""

import os
from pathlib import Path
import tempfile
import traceback

import numpy as np
import slicer
import vtk
from PETAnalyzerLib.report import export


def run():
    widget = slicer.modules.petanalyzer.widgetRepresentation().self()
    dicom = os.environ.get("PET_ANALYZER_DICOM")
    if dicom:
        from DICOMLib import DICOMUtils
        with DICOMUtils.TemporaryDICOMDatabase() as database:
            DICOMUtils.importDicom(dicom, database)
            ids = DICOMUtils.loadPatientByUID(database.patients()[0])
        volume = next(slicer.mrmlScene.GetNodeByID(node_id) for node_id in ids
                      if slicer.mrmlScene.GetNodeByID(node_id).IsA("vtkMRMLScalarVolumeNode"))
        background, lung = (26, 127, 164), (26, 127, 128)
        expected_background = 2.97232
    else:
        from test_analysis import fixture
        from scipy.ndimage import gaussian_filter
        image, spacing, _ = fixture()
        volume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "Synthetic phantom")
        volume.SetSpacing(*spacing[::-1])
        slicer.util.updateVolumeFromArray(volume, gaussian_filter(image, 2.2 / spacing) * 1000)
        background, lung = (21, 50, 50), None
        expected_background = 3.1
    widget.volume.setCurrentNode(volume)
    widget.method.currentIndex = 1
    widget.source.value = 12.5
    widget.stock.value = 1000
    widget.tank.value = 9800
    widget.bg_activity.value = 30.38
    widget.bg_radius.value = 10
    widget.reference.checked = True
    widget.bands.currentIndex = 0 if dicom else 2
    widget.auto()
    _, _, _, matrix = widget.image()
    for selector, point in ((widget.background, background), (widget.lung, lung)):
        if point is not None:
            node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
            ras = matrix.MultiplyPoint((*point[::-1], 1))[:3]
            node.AddControlPoint(vtk.vtkVector3d(*ras))
            selector.setCurrentNode(node)
    widget.run()
    widget.confirm.checked = True
    assert widget.export_button.enabled
    assert not widget.chart.pixmap.isNull()
    assert widget.table.rowCount == 6
    result, metadata, *_ = widget.review
    assert abs(result["background_mean_kbq_ml"] - expected_background) < .04, result
    assert widget.overlay.GetSegmentation().GetNumberOfSegments() >= 13
    with tempfile.TemporaryDirectory() as folder:
        path = export(Path(folder) / "report", result, metadata, metadata["comparison"])
        assert path.exists()
        assert "<svg" in path.read_text()
    print("PET ANALYZER SLICER SMOKE PASSED", "DICOM" if dicom else "synthetic", flush=True)


try:
    run()
except Exception:
    traceback.print_exc()
    slicer.util.exit(1)
else:
    slicer.util.exit(0)
