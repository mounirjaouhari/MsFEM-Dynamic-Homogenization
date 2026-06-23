"""Shared publication style: high-quality COLOUR figures (CMAME/Elsevier),
designed so that CURVES stay distinguishable even when the paper is printed in
black and white -- i.e. each curve carries a distinct COLOUR *and* a distinct
LINESTYLE *and* a distinct MARKER. Field maps keep attractive colour colormaps
(diverging RdBu for signed fields) with a thin zero-contour. Colours are the
colourblind-safe Okabe-Ito palette.

Fonts: Computer-Modern-like serif via mathtext; sizes tuned for a two-column page.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler

try:
    import scienceplots  # noqa: F401  (registers the 'science' style)
    _BASE = ["science", "no-latex"]
except Exception:
    _BASE = []

# attractive, colourblind-safe field colormap for SIGNED fields (correctors,
# modes, wave fields). A black zero-contour (added in the scripts) keeps the
# sign change visible if printed in grayscale.
FIELD_CMAP = "RdBu_r"

# Okabe-Ito colourblind-safe palette
_C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#000000"]

# per-curve styles: distinct COLOUR + LINESTYLE + MARKER  -> survives B&W print
BW = [dict(color=_C[0], ls="-",  marker="o"),
      dict(color=_C[1], ls="--", marker="s"),
      dict(color=_C[2], ls=":",  marker="^"),
      dict(color=_C[3], ls="-.", marker="D"),
      dict(color=_C[4], ls=(0, (3, 1, 1, 1)), marker="v")]

BW_CYCLE = (cycler("color", _C[:5])
            + cycler("linestyle", ["-", "--", ":", "-.", (0, (3, 1, 1, 1))])
            + cycler("marker", ["o", "s", "^", "D", "v"]))


def apply():
    if _BASE:
        try:
            plt.style.use(_BASE)
        except Exception:
            pass
    mpl.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "font.size": 9,
        "axes.titlesize": 9,
        "axes.labelsize": 9,
        "legend.fontsize": 7.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.7,
        "lines.linewidth": 1.3,
        "lines.markersize": 3.8,
        "patch.linewidth": 0.7,
        "savefig.dpi": 600,
        "figure.dpi": 120,
        "savefig.bbox": "tight",
        "image.cmap": FIELD_CMAP,
        "axes.prop_cycle": BW_CYCLE,
        "legend.frameon": True,
        "legend.framealpha": 0.9,
        "legend.edgecolor": "0.7",
    })
