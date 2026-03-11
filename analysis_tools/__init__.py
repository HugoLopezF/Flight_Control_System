import matplotlib as mpl

# Global plot style for analysis tools
mpl.rcParams.update(
    {
        # LaTeX-like look without requiring a TeX install
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif", "Times New Roman"],
        "mathtext.fontset": "stix",
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "axes.titleweight": "normal",
        "figure.titlesize": 14,
        "figure.titleweight": "normal",
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    }
)
