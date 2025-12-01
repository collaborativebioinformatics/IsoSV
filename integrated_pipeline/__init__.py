from .__version__ import __version__

# Re-export key functionality for a clean, discoverable public API.
from .bam_parser import (
		walk_cigar_indels_and_clips,
		parse_sa_tag,
		cigarstr_to_tuples,
)
from .clustering import (
		cluster_indels,
		cluster_clips,
		cluster_splits,
)
from .annotator import annotate_candidate
from . import variant_objects
from . import utils


__all__ = [
		"__version__",
		# Core parsing utilities
		"walk_cigar_indels_and_clips",
		"parse_sa_tag",
		"cigarstr_to_tuples",
		# Clustering helpers
		"cluster_indels",
		"cluster_clips",
		"cluster_splits",
		# Annotation utilities
		"annotate_candidate",
		# Data containers and shared utilities
		"variant_objects",
		"utils",
]
