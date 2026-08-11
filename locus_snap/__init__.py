"""LocusSnap: IGV-like genomic snapshots from an indexed BAM, built on pysam."""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("locus-snap")
except PackageNotFoundError:
    # Running from a source checkout that was never `pip install`-ed.
    __version__ = "0.0.0+unknown"

from locus_snap.cli import apply_config_preferences, build_parser, main
from locus_snap.molecule import (
    MOLECULE_TAGS,
    MoleculeBuildResult,
    build_consensus_read,
    build_molecule_consensus_reads,
    resolve_molecule_tag,
)
from locus_snap.rna import (
    CANONICAL_SPLICE_MOTIFS,
    RNA_STRANDNESS_MODES,
    FusionEvidence,
    JunctionEvidence,
    annotated_junctions,
    collect_fusion_evidence,
    collect_splice_junctions,
    rna_evidence_rows,
    splice_motif,
    transcript_strand,
    write_rna_evidence_tsv,
)
from locus_snap.track_plugin import (
    TRACK_PLUGIN_API_VERSION,
    TRACK_PLUGIN_ENTRY_POINT_GROUP,
    LoadedPluginTrack,
    PluginTrackSource,
    TrackCanvas,
    TrackPlugin,
    TrackPluginError,
    TrackRegion,
    available_track_plugins,
    load_track_plugin,
)

__all__ = [
    "__version__", "apply_config_preferences", "build_parser", "main",
    "TRACK_PLUGIN_API_VERSION", "TRACK_PLUGIN_ENTRY_POINT_GROUP",
    "TrackPlugin", "TrackCanvas", "TrackRegion", "TrackPluginError",
    "PluginTrackSource", "LoadedPluginTrack", "available_track_plugins",
    "load_track_plugin",
    "MOLECULE_TAGS", "MoleculeBuildResult", "resolve_molecule_tag",
    "build_consensus_read", "build_molecule_consensus_reads",
    "RNA_STRANDNESS_MODES", "CANONICAL_SPLICE_MOTIFS",
    "JunctionEvidence", "FusionEvidence", "collect_splice_junctions",
    "collect_fusion_evidence", "annotated_junctions", "splice_motif",
    "transcript_strand", "rna_evidence_rows", "write_rna_evidence_tsv",
]
