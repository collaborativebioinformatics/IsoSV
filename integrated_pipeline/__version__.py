"""Package version information for the IsoSV integrated pipeline.

This module centralises the version string so it can be imported
consistently from both the library code and the CLI entry point.

The version is deliberately kept as a simple ``major.minor.patch`` style
string rather than being coupled to any specific packaging tool.  When
cutting a new release, update ``__version__`` here and, if applicable,
in any downstream packaging metadata.
"""

__all__ = ["__version__"]

__version__: str = "0.1.0"
