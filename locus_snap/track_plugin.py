"""Public API for third-party genomic track renderers.

Plugins are ordinary Python packages registered through the
``locus_snap.track_plugins.v1`` entry-point group.  They fetch any payload they
need for a genomic window, then draw it through :class:`TrackCanvas`.  The
canvas is intentionally smaller and more stable than Matplotlib's Axes API;
LocusSnap retains ownership of layout, genomic limits, grids, highlights, and
figure export.
"""
from __future__ import annotations

import importlib
import inspect
from abc import ABC, abstractmethod
from dataclasses import dataclass
from importlib import metadata
from types import MappingProxyType
from typing import Any, Dict, List, Mapping, Optional, Sequence

from matplotlib.colors import is_color_like
from matplotlib.patches import Rectangle


TRACK_PLUGIN_API_VERSION = 1
TRACK_PLUGIN_ENTRY_POINT_GROUP = "locus_snap.track_plugins.v1"

__all__ = [
    "TRACK_PLUGIN_API_VERSION", "TRACK_PLUGIN_ENTRY_POINT_GROUP",
    "TrackPlugin", "TrackCanvas", "TrackRegion", "TrackPluginError",
    "PluginTrackSource", "LoadedPluginTrack", "available_track_plugins",
    "load_track_plugin",
]


class TrackPluginError(ValueError):
    """A configuration, loading, fetching, or rendering plugin failure."""


@dataclass(frozen=True)
class TrackRegion:
    """A zero-based, half-open genomic window supplied to a track plugin."""

    chrom: str
    start: int
    end: int

    @property
    def span(self) -> int:
        return self.end - self.start


class TrackPlugin(ABC):
    """Base class implemented by external LocusSnap track plugins.

    ``options`` contains strings supplied after the plugin name on
    ``--plugin_track``. Plugins should validate their own option names and
    values in :meth:`fetch` and raise ``ValueError`` with an actionable message.
    """

    api_version = TRACK_PLUGIN_API_VERSION
    name = "track-plugin"
    default_label = "Plugin track"
    default_height_in = 0.85
    default_color = "#2c7fb8"

    @abstractmethod
    def fetch(self, region: TrackRegion, options: Mapping[str, str]) -> Any:
        """Return any payload needed to render ``region``."""

    @abstractmethod
    def render(
        self,
        canvas: "TrackCanvas",
        payload: Any,
        region: TrackRegion,
        options: Mapping[str, str],
    ) -> None:
        """Draw ``payload`` using only the stable canvas interface."""


@dataclass(frozen=True)
class LoadedPluginTrack:
    """A fetched plugin track consumed by the core renderer."""

    label: str
    height_in: float
    color: str
    plugin: TrackPlugin
    payload: Any
    region: TrackRegion
    options: Mapping[str, str]


class PluginTrackSource:
    """Configured plugin instance with the same ``fetch`` shape as core tracks."""

    def __init__(
        self,
        plugin: TrackPlugin,
        options: Optional[Mapping[str, str]] = None,
        label: Optional[str] = None,
        height_in: Optional[float] = None,
        color: Optional[str] = None,
    ):
        validate_track_plugin(plugin)
        selected_height = (
            plugin.default_height_in if height_in is None else height_in
        )
        if (
            isinstance(selected_height, bool)
            or not isinstance(selected_height, (int, float))
            or selected_height <= 0
        ):
            raise TrackPluginError("Plugin track height must be greater than zero.")
        selected_color = color or plugin.default_color
        if not isinstance(selected_color, str) or not is_color_like(selected_color):
            raise TrackPluginError(
                f"Invalid plugin track colour: {selected_color!r}."
            )
        selected_options = dict(options or {})
        for key, value in selected_options.items():
            if not isinstance(key, str) or not key:
                raise TrackPluginError("Plugin option names must be non-empty strings.")
            if not isinstance(value, str):
                raise TrackPluginError(
                    f"Plugin option '{key}' must have a string value."
                )
        selected_label = label or plugin.default_label or plugin.name
        if not isinstance(selected_label, str) or not selected_label.strip():
            raise TrackPluginError("Plugin track label cannot be empty.")
        self.plugin = plugin
        # Keep the source picklable for --batch_regions worker processes;
        # plugins still receive a read-only view at fetch/render time.
        self.options = selected_options
        self.label = selected_label
        self.height_in = float(selected_height)
        self.color = selected_color

    def fetch(self, chrom: str, start: int, end: int) -> LoadedPluginTrack:
        if end <= start:
            raise TrackPluginError("Plugin track region end must exceed its start.")
        region = TrackRegion(chrom, start, end)
        read_only_options = MappingProxyType(dict(self.options))
        try:
            payload = self.plugin.fetch(region, read_only_options)
        except Exception as exc:
            raise TrackPluginError(
                f"Track plugin '{self.plugin.name}' failed while fetching "
                f"{chrom}:{start + 1}-{end}: {exc}"
            ) from exc
        return LoadedPluginTrack(
            label=self.label,
            height_in=self.height_in,
            color=self.color,
            plugin=self.plugin,
            payload=payload,
            region=region,
            options=read_only_options,
        )


class TrackCanvas:
    """Versioned drawing facade passed to external track plugins.

    Coordinates supplied to drawing methods are data coordinates. The genomic
    x range is fixed by LocusSnap and cannot be changed through this API.
    Methods return their Matplotlib artist for optional plugin-side annotation.
    """

    api_version = TRACK_PLUGIN_API_VERSION

    def __init__(
        self,
        axis,
        region: TrackRegion,
        default_color: str,
        colors: Mapping[str, str],
    ):
        self._axis = axis
        self.region = region
        self.default_color = default_color
        self.colors = MappingProxyType(dict(colors))

    @property
    def pixels_per_base(self) -> float:
        return self._axis.get_window_extent().width / max(self.region.span, 1)

    def set_y_limits(self, minimum: float, maximum: float) -> None:
        if maximum <= minimum:
            raise TrackPluginError("Plugin y-axis maximum must exceed its minimum.")
        self._axis.set_ylim(minimum, maximum)

    def set_y_ticks(
        self, positions: Sequence[float], labels: Optional[Sequence[str]] = None,
    ) -> None:
        self._axis.set_yticks(list(positions))
        if labels is not None:
            if len(labels) != len(positions):
                raise TrackPluginError(
                    "Plugin y-tick labels must match the number of positions."
                )
            self._axis.set_yticklabels(list(labels), fontsize=6)
        self._axis.tick_params(
            left=True, labelleft=True, labelsize=6,
            colors=self.colors.get("axis", "#898781"), length=2,
        )
        self._axis.spines["left"].set_visible(True)
        self._axis.spines["left"].set_color(
            self.colors.get("axis", "#898781")
        )

    def line(
        self, x, y, *, color: Optional[str] = None, linewidth: float = 0.9,
        alpha: float = 1.0, linestyle: str = "-", zorder: float = 3,
    ):
        return self._axis.plot(
            x, y, color=color or self.default_color, linewidth=linewidth,
            alpha=alpha, linestyle=linestyle, zorder=zorder,
        )[0]

    def step(
        self, x, y, *, color: Optional[str] = None, linewidth: float = 0.9,
        alpha: float = 1.0, where: str = "mid", zorder: float = 3,
    ):
        return self._axis.step(
            x, y, color=color or self.default_color, linewidth=linewidth,
            alpha=alpha, where=where, zorder=zorder,
        )[0]

    def fill_between(
        self, x, y1, y2=0, *, color: Optional[str] = None,
        alpha: float = 0.30, zorder: float = 2,
    ):
        return self._axis.fill_between(
            x, y1, y2, color=color or self.default_color,
            alpha=alpha, linewidth=0, zorder=zorder,
        )

    def bars(
        self, x, height, *, width=1.0, bottom=0, color: Optional[str] = None,
        alpha: float = 0.85, zorder: float = 3,
    ):
        return self._axis.bar(
            x, height, width=width, bottom=bottom,
            color=color or self.default_color, alpha=alpha,
            linewidth=0, zorder=zorder,
        )

    def scatter(
        self, x, y, *, color: Optional[str] = None, size: float = 12.0,
        marker: str = "o", alpha: float = 0.9, zorder: float = 4,
    ):
        return self._axis.scatter(
            x, y, c=color or self.default_color, s=size, marker=marker,
            alpha=alpha, linewidths=0, zorder=zorder,
        )

    def rectangle(
        self, x: float, y: float, width: float, height: float, *,
        color: Optional[str] = None, alpha: float = 0.85,
        edgecolor: str = "none", linewidth: float = 0.0,
        zorder: float = 3,
    ):
        patch = Rectangle(
            (x, y), width, height, facecolor=color or self.default_color,
            edgecolor=edgecolor, linewidth=linewidth, alpha=alpha,
            zorder=zorder,
        )
        self._axis.add_patch(patch)
        return patch

    def span(
        self, start: float, end: float, *, color: Optional[str] = None,
        alpha: float = 0.20, zorder: float = 1,
    ):
        return self._axis.axvspan(
            start, end, facecolor=color or self.default_color,
            edgecolor="none", alpha=alpha, zorder=zorder,
        )

    def horizontal_line(
        self, y: float, *, color: Optional[str] = None,
        linewidth: float = 0.6, alpha: float = 1.0,
        linestyle: str = "-", zorder: float = 2,
    ):
        return self._axis.axhline(
            y, color=color or self.colors.get("axis", self.default_color),
            linewidth=linewidth, alpha=alpha,
            linestyle=linestyle, zorder=zorder,
        )

    def text(
        self, x: float, y: float, value: str, *,
        color: Optional[str] = None, fontsize: float = 6.5,
        horizontal_alignment: str = "left",
        vertical_alignment: str = "center", rotation: float = 0,
        zorder: float = 5,
    ):
        return self._axis.text(
            x, y, str(value), color=color or self.default_color,
            fontsize=fontsize, ha=horizontal_alignment,
            va=vertical_alignment, rotation=rotation,
            clip_on=True, zorder=zorder,
        )


def validate_track_plugin(plugin: Any) -> TrackPlugin:
    if not isinstance(plugin, TrackPlugin):
        raise TrackPluginError(
            "Track plugin must be an instance of locus_snap.TrackPlugin."
        )
    if plugin.api_version != TRACK_PLUGIN_API_VERSION:
        raise TrackPluginError(
            f"Track plugin '{plugin.name}' uses API version {plugin.api_version}; "
            f"LocusSnap requires version {TRACK_PLUGIN_API_VERSION}."
        )
    if not isinstance(plugin.name, str) or not plugin.name.strip():
        raise TrackPluginError("Track plugin name cannot be empty.")
    return plugin


def _entry_points_for_group():
    discovered = metadata.entry_points()
    if hasattr(discovered, "select"):
        return list(discovered.select(group=TRACK_PLUGIN_ENTRY_POINT_GROUP))
    return list(discovered.get(TRACK_PLUGIN_ENTRY_POINT_GROUP, []))


def available_track_plugins() -> Dict[str, str]:
    """Return installed entry-point names and their import targets."""
    return {
        entry_point.name: entry_point.value
        for entry_point in _entry_points_for_group()
    }


def load_track_plugin(identifier: str) -> TrackPlugin:
    """Load an installed plugin name or an explicit ``module:object`` target."""
    identifier = identifier.strip()
    if not identifier:
        raise TrackPluginError("Plugin identifier cannot be empty.")
    loaded = None
    if ":" in identifier:
        module_name, object_name = identifier.rsplit(":", 1)
        try:
            module = importlib.import_module(module_name)
            loaded = getattr(module, object_name)
        except (ImportError, AttributeError) as exc:
            raise TrackPluginError(
                f"Cannot import track plugin '{identifier}': {exc}"
            ) from exc
    else:
        matches = [
            entry_point for entry_point in _entry_points_for_group()
            if entry_point.name == identifier
        ]
        if not matches:
            installed = ", ".join(sorted(available_track_plugins())) or "none"
            raise TrackPluginError(
                f"Unknown track plugin '{identifier}'. Installed plugins: {installed}. "
                "Use an entry-point name or module:object target."
            )
        if len(matches) > 1:
            raise TrackPluginError(
                f"Multiple installed track plugins use the name '{identifier}'."
            )
        try:
            loaded = matches[0].load()
        except Exception as exc:
            raise TrackPluginError(
                f"Cannot load track plugin '{identifier}': {exc}"
            ) from exc

    if inspect.isclass(loaded):
        try:
            loaded = loaded()
        except Exception as exc:
            raise TrackPluginError(
                f"Cannot instantiate track plugin '{identifier}': {exc}"
            ) from exc
    elif callable(loaded) and not isinstance(loaded, TrackPlugin):
        try:
            loaded = loaded()
        except Exception as exc:
            raise TrackPluginError(
                f"Track plugin factory '{identifier}' failed: {exc}"
            ) from exc
    return validate_track_plugin(loaded)


def build_plugin_track_sources(
    specifications: Optional[Sequence[Sequence[str]]],
) -> List[PluginTrackSource]:
    """Parse repeatable ``PLUGIN key=value ...`` CLI definitions."""
    sources = []
    for specification in specifications or []:
        fields = list(specification)
        if not fields:
            raise TrackPluginError("--plugin_track requires a plugin name.")
        plugin = load_track_plugin(fields[0])
        parsed: Dict[str, str] = {}
        for field in fields[1:]:
            if "=" not in field:
                raise TrackPluginError(
                    "Plugin track options must use KEY=VALUE syntax: "
                    f"{field!r}."
                )
            key, value = field.split("=", 1)
            key = key.strip()
            if not key or not value:
                raise TrackPluginError(
                    f"Invalid plugin track option {field!r}; KEY and VALUE are required."
                )
            if key in parsed:
                raise TrackPluginError(f"Duplicate plugin track option: {key}.")
            parsed[key] = value

        label = parsed.pop("track_label", None)
        height_text = parsed.pop("track_height", None)
        color = parsed.pop("track_color", None)
        height = None
        if height_text is not None:
            try:
                height = float(height_text)
            except ValueError as exc:
                raise TrackPluginError(
                    f"Plugin track_height must be numeric: {height_text!r}."
                ) from exc
        sources.append(PluginTrackSource(
            plugin, parsed, label=label, height_in=height, color=color,
        ))
    return sources
