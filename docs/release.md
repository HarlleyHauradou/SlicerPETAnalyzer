# Distribution and publication

## Slicer Extensions Index

The source tree contains the Slicer CMake entry point, scripted module and icon. The tested host is 3D Slicer 5.10.0 on Linux. Before submitting a pull request to [Slicer Extensions Index](https://github.com/Slicer/ExtensionsIndex):

- Build and install the extension package against a supported Slicer build tree. Confirm Python package installation and generic CTest discovery.
- Run the workflow on Windows, macOS and Linux. Capture an English workflow screenshot, check appearance with light and dark themes, and verify the icon URL.
- Tag an immutable version. In the Extensions Index submission use the repository URL and the immutable revision. Follow the [official extension submission guide](https://slicer.readthedocs.io/en/latest/developer_guide/extensions.html).

## Archive DOI and paper

The repository contains `CITATION.cff`, the scientific scope, tests and a reproducible measured DICOM benchmark. A versioned archival DOI can be minted by activating the repository in Zenodo and archiving a tagged GitHub release. Verify the author list and metadata before publishing.

A software paper should report detection success and size-label accuracy against independent CT or manual centers, inter-operator agreement, robustness to reconstruction and low counts, and agreement with EFOMP/EARL tooling for known FDG measurements. Include confidence intervals and failure cases, and keep the full NEMA acquisition and accreditation claims separate from supported software measurements.
