"""LocusSnap: IGV-like genomic snapshots from an indexed BAM, built on pysam."""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("locus-snap")
except PackageNotFoundError:
    # Running from a source checkout that was never `pip install`-ed.
    __version__ = "0.0.0+unknown"

from locus_snap.cli import apply_config_preferences, build_parser, main
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
]
