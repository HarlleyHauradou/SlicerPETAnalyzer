set(MODULE_NAME PETAnalyzer)
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_NAME}.py PETAnalyzerLib/__init__.py PETAnalyzerLib/analysis.py PETAnalyzerLib/report.py
  RESOURCES Resources/PETAnalyzer.png
  WITH_GENERIC_TESTS
)
