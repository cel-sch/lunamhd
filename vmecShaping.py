#!/Users/cellywelly/Dev/solis/bin/python3
"""
Analytic shaping coefficient extraction from VMEC wout Fourier coefficients.

Shaping coefficients are derived from the n=0 Fourier modes of R, Z, and B at
each flux surface.  For 3-D equilibria (stellarators) the n=0 modes represent
the toroidally-averaged shape; for axisymmetric cases they are exact.

Based on the parameterisation of Graves, PPCF 55 (2013) 074009, Appendix B1.
The VMEC Fourier harmonics are identified with analytical coefficients via:

    R(θ) = R₀(s) + R₁c cos θ + R₂c cos 2θ + …
    Z(θ) =         Z₁s sin θ + Z₂s sin 2θ + …
    B(θ) = B₀c    + B₁c cos θ + B₂c cos 2θ + …

where Rmc = rmnc[s, m, n=0], Zms = zmns[s, m, n=0], Bmc = bmnc[s, m, n=0].

Quantities computed  (Graves 2013 notation)
-------------------------------------------
R0(s)   major radius profile    rmnc[m=0, n=0]
a(s)    minor radius            rmnc[m=1, n=0]
ε(s)    inverse aspect ratio    a / R0
r(s)    mean semi-axis          (R₁c + Z₁s) / 2
S₂(s)   elongation amplitude    (R₁c − Z₁s) / 2
S₃(s)   triangularity amplitude rmnc[m=2, n=0]
κ(s)    elongation              (r − S₂) / (r + S₂) = Z₁s / R₁c
δ(s)    triangularity           4 S₃ / r
Δ(s)    Shafranov shift         R0[0] − R0(s)
F₂(s)   toroidal flux variation  bsubvmnc[m=0, n=0] / rbtor0 − 1
DI      Mercier shape factor     (3/4)(κ−1)(1 − 2δ/ε)

Sign conventions
----------------
δ > 0   inward triangularity (D-shape)
Δ > 0   outward shift from axis
F₂ < 0  typical (diamagnetic: finite-β reduces toroidal flux function)

Axis note
---------
At s = 0 the flux surface degenerates to a line so R₁c = Z₁s ≈ 0 and κ, δ, ζ
are set to NaN.  Boundary values (s = 1, index −1) are always well-defined.

Usage as a module
-----------------
    from vmecShaping import ShapingCoeffs
    from vmecReader import WoutReader

    w  = WoutReader('Output/KH/my_run/VMEC/wout/wout_my_run_0.nc')
    sc = ShapingCoeffs(w)
    print(sc.summary())

    kappa_profile = sc.kappa          # 1-D array, length ns
    delta_bdry    = sc.delta[-1]      # scalar at boundary

Usage as a script
-----------------
    python3 vmecShaping.py my_run                  # all wout files, print + save csv
    python3 vmecShaping.py my_run --index 0        # single file
    python3 vmecShaping.py path/to/wout.nc         # direct path
    python3 vmecShaping.py my_run --fmt npz --outdir ./shaping_out
    python3 vmecShaping.py my_run --no-save        # print summary only
"""

import sys
import numpy as np
from pathlib import Path

from vmecReader import WoutReader, find_wout, list_wouts, read_run_wouts


# ---------------------------------------------------------------------------
# Core class
# ---------------------------------------------------------------------------

class ShapingCoeffs:
    """
    Analytic shaping coefficients derived from the n=0 Fourier modes of a
    VMEC wout file.

    Parameters
    ----------
    wout : WoutReader
        A loaded WoutReader instance.

    Attributes (all 1-D arrays of length ns)
    -----------------------------------------
    s       normalised flux label
    R0      major radius profile
    a       minor radius (= R₁c)
    eps     inverse aspect ratio a/R0
    r       mean semi-axis (R₁c + Z₁s)/2
    S2      elongation amplitude (R₁c − Z₁s)/2
    S3      triangularity amplitude R₂c
    kappa   elongation κ = Z₁s/R₁c
    delta   triangularity δ = 4S₃/r
    shift   Shafranov shift Δ (positive = outward from axis)
    F2      toroidal flux variation bsubvmnc[m=0,n=0]/rbtor0 − 1
    di      Mercier shape factor DI
    """

    def __init__(self, wout):
        self.wout = wout
        self.s = wout.s.copy()
        self._extract_modes()
        self._compute()

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _get_n0(self, arr, m):
        """Return the radial profile for mode (m, n=0) from a 2-D coeff array (xm/xn grid)."""
        idx = np.where((self.wout.xm == m) & (self.wout.xn == 0))[0]
        if len(idx) == 0:
            return np.zeros(self.wout.ns)
        return arr[:, idx[0]].copy()

    def _get_n0_nyq(self, arr, m):
        """Return the radial profile for mode (m, n=0) from a 2-D coeff array (xm_nyq/xn_nyq grid)."""
        idx = np.where((self.wout.xm_nyq == m) & (self.wout.xn_nyq == 0))[0]
        if len(idx) == 0:
            return np.zeros(self.wout.ns)
        return arr[:, idx[0]].copy()

    def _extract_modes(self):
        """Cache the n=0 harmonic profiles needed for shaping."""
        self._R0c = self._get_n0(self.wout.rmnc, 0)
        self._R1c = self._get_n0(self.wout.rmnc, 1)
        self._R2c = self._get_n0(self.wout.rmnc, 2)
        self._R3c = self._get_n0(self.wout.rmnc, 3)
        self._Z1s = self._get_n0(self.wout.zmns, 1)
        self._Z2s = self._get_n0(self.wout.zmns, 2)
        self._Z3s = self._get_n0(self.wout.zmns, 3)
        # Toroidal flux function F = R·Bφ (Nyquist grid covariant toroidal B)
        # For axisymmetric case: B_ζ = ∂r/∂ζ · B = R·Bφ = F(ψ) exactly
        self._Fv = self._get_n0_nyq(self.wout.bsubvmnc, 0)

    def _safe_div(self, num, den):
        """Element-wise division; return NaN where den == 0 (axis singularity)."""
        with np.errstate(invalid='ignore', divide='ignore'):
            return np.where(np.abs(den) > 0, num / den, np.nan)

    # ------------------------------------------------------------------
    # Coefficient computation
    # ------------------------------------------------------------------

    # def _compute(self):
    #     # what claude produced
    #     R1c = self._R1c
    #     Z1s = self._Z1s

    #     self.R0    = self._R0c.copy()
    #     self.a     = R1c.copy()
    #     self.kappa = self._safe_div(Z1s, R1c)

    #     # Miller: δ = sin(arcsin(δ)) where arcsin(δ) = 2·R2c / R1c
    #     angle      = self._safe_div(2.0 * self._R2c, R1c)
    #     self.delta = np.where(np.isfinite(angle), np.sin(angle), np.nan)

    #     # Turnbull–Miller squareness: ζ = −Z3s / Z1s
    #     self.zeta  = self._safe_div(-self._Z3s, Z1s)

    #     # Shafranov shift relative to boundary (positive = outward)
    #     self.shift = self._R0c - self._R0c[-1]

    def _compute(self):
        # Based on PPCF 55 2013 074009, Appendix B1
        R0c = self._R0c.copy()
        R1c = self._R1c.copy()
        R2c = self._R2c.copy()
        Z1s = self._Z1s.copy()

        self.R0  = R0c
        self.a   = R1c
        self.eps = self._safe_div(self.a, self.R0)
        self.r   = (R1c + Z1s) / 2
        self.S2  = (R1c - Z1s) / 2   # elongation amplitude
        self.S3  = R2c                # triangularity amplitude

        # κ: elongation = Z₁s/R₁c = (r−S₂)/(r+S₂)
        self.kappa = self._safe_div(self.r - self.S2, self.r + self.S2)

        # δ: triangularity = 4S₃/r
        self.delta = self._safe_div(4 * self.S3, self.r)

        # Δ: Shafranov shift relative to axis (positive = outward)
        self.shift = self.R0[0] - R0c

        # F₂: fractional variation of toroidal flux function F = R·Bφ
        # F(r) = R₀B₀(1 + F₂), so F₂ = bsubvmnc[m=0,n=0] / rbtor0 − 1
        self.F2 = self._safe_div(self._Fv, self.wout.rbtor0) - 1.0

        # DI: Mercier shape factor
        self.di = (3/4) * (self.kappa - 1) * (1 - 2 * self._safe_div(self.delta, self.eps))

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def get_mode(self, coeff, m, n=0):
        """
        Return the radial profile of any Fourier mode (m, n) from the wout data.

        Parameters
        ----------
        coeff : 'rmnc' or 'zmns'
        m, n  : mode numbers

        Returns
        -------
        1-D array of length ns
        """
        arr = getattr(self.wout, coeff)
        idx = np.where((self.wout.xm == m) & (self.wout.xn == n))[0]
        if len(idx) == 0:
            raise ValueError(f"Mode ({coeff}, m={m}, n={n}) not found in wout data")
        return arr[:, idx[0]].copy()

    def as_dict(self):
        """Return all shaping profiles as a plain dict."""
        return {
            's':     self.s,
            'R0':    self.R0,
            'a':     self.a,
            'eps':   self.eps,
            'r':     self.r,
            'S2':    self.S2,
            'S3':    self.S3,
            'kappa': self.kappa,
            'delta': self.delta,
            'shift': self.shift,
            'F2':    self.F2,
            'DI':    self.di,
        }

    def summary(self, label=None):
        """Return a formatted string with axis and boundary values of each coefficient."""
        name = label or self.wout.path.name

        # Axis: use index 1 to avoid the s=0 singularity
        ax = 1

        def fmt(arr):
            ax_val = arr[ax]
            bd_val = arr[-1]
            ax_s = f'{ax_val:.4f}' if np.isfinite(ax_val) else '   NaN'
            bd_s = f'{bd_val:.4f}' if np.isfinite(bd_val) else '   NaN'
            return ax_s, bd_s

        rows = [
            ('R0  [m]',             *fmt(self.R0)),
            ('κ  (elongation)',     *fmt(self.kappa)),
            ('δ  (triangularity)',  *fmt(self.delta)),
            ('Δ  (Shafranov) [m]',  *fmt(self.shift)),
            ('F₂ (field 2nd harm)', *fmt(self.F2)),
            ('DI (Mercier mod)',    *fmt(self.di)),
        ]

        lines = [
            f"Shaping coefficients: {name}",
            f"  {'Quantity':<24} {'near-axis':>10} {'boundary':>10}",
            f"  {'-'*46}",
        ]
        for label_, ax_s, bd_s in rows:
            lines.append(f"  {label_:<24} {ax_s:>10} {bd_s:>10}")
        return '\n'.join(lines)

    def save(self, outpath=None, fmt='csv'):
        """
        Save shaping profiles to a file.

        Parameters
        ----------
        outpath : str or Path, optional
            Defaults to <wout_stem>_shaping.<fmt> in the wout's directory.
        fmt : 'csv' or 'npz'

        Returns
        -------
        Path of the saved file.
        """
        if outpath is None:
            outpath = self.wout.path.parent / f'{self.wout.path.stem}_shaping.{fmt}'
        outpath = Path(outpath)

        d = self.as_dict()

        if fmt == 'npz':
            np.savez(str(outpath), **d)
        elif fmt == 'csv':
            import csv
            keys = list(d.keys())
            with open(outpath, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(keys)
                for row in zip(*[d[k] for k in keys]):
                    writer.writerow(row)
        else:
            raise ValueError(f"Unknown format '{fmt}', choose 'csv' or 'npz'")

        return outpath


# ---------------------------------------------------------------------------
# Convenience functions (thin wrappers for interactive use)
# ---------------------------------------------------------------------------

def elongation(wout):
    """Return the κ(s) profile for a WoutReader."""
    return ShapingCoeffs(wout).kappa

def triangularity(wout):
    """Return the δ(s) profile for a WoutReader."""
    return ShapingCoeffs(wout).delta

def shafranov_shift(wout):
    """Return the Δ(s) profile for a WoutReader."""
    return ShapingCoeffs(wout).shift


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------

def compute_shaping(run_name_or_path, index=None, fmt='csv', outdir=None, save=True):
    """
    Compute and optionally save shaping coefficients for one or all wout files.

    Parameters
    ----------
    run_name_or_path : str or Path
        Path to a wout .nc file, or a run name (uses default output root).
    index : int or None
        Specific wout index. If None and a run name is given, all wout files
        in the run are processed.
    fmt : 'csv' or 'npz'
    outdir : str or Path, optional
        Output directory. Defaults to each wout file's own directory.
    save : bool
        Write output files. If False, only print summaries.

    Returns
    -------
    list of ShapingCoeffs objects.
    """
    path = Path(run_name_or_path)

    if path.suffix == '.nc':
        readers = [WoutReader(path)]
    elif index is not None:
        readers = [WoutReader(find_wout(run_name_or_path, index=index))]
    else:
        readers = read_run_wouts(run_name_or_path)

    results = []
    for w in readers:
        sc = ShapingCoeffs(w)
        print(sc.summary())
        if save:
            if outdir is not None:
                out = Path(outdir) / f'{w.path.stem}_shaping.{fmt}'
            else:
                out = None
            saved = sc.save(outpath=out, fmt=fmt)
            print(f"  → saved {saved}\n")
        results.append(sc)

    return results


# ---------------------------------------------------------------------------
# Script entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Compute analytic shaping coefficients from VMEC wout files.')
    parser.add_argument('target',
                        help='Path to wout .nc file, or run name (uses default output root)')
    parser.add_argument('--index', type=int, default=None,
                        help='Wout file index (default: all indices for a run name)')
    parser.add_argument('--fmt', choices=['csv', 'npz'], default='csv',
                        help='Output file format (default: csv)')
    parser.add_argument('--outdir', default=None,
                        help='Output directory for saved files')
    parser.add_argument('--no-save', action='store_true',
                        help='Print summary only, do not write files')
    args = parser.parse_args()

    compute_shaping(
        args.target,
        index=args.index,
        fmt=args.fmt,
        outdir=args.outdir,
        save=not args.no_save,
    )
