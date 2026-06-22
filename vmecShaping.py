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
    python3 vmecShaping.py my_run                  # all wout files, print + save npz
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
        self.r   = (R1c + Z1s) / 2
        # self.a = R1c?
        self.eps = self._safe_div(self.r, self.R0)
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

        # M(s): local Mach number profile, M_0 * omega(s)/omega_axis
        omega = self.wout.omega.copy()
        M0 = np.sqrt(max(self.wout.machsq, 0.0))
        self.mach = M0 * self._safe_div(omega, omega[0]) if omega[0] != 0.0 else np.full(len(self.s), M0)

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def reconstruct(self, theta, simple=True):
        """
        Reconstruct R(θ) and Z(θ) analytically from the Graves shaping coefficients.

        simple=True (default) — circular cross-section only, no elongation (S₂) or
            triangularity (S₃).  Matches the analytic parameterisation
                R = R₀ + r·cos θ
                Z = r·sin θ
            where R₀(s) already carries the Shafranov shift (it is the m=0
            Fourier component of R from VMEC, so R₀(s) = R₀_axis − Δ(s)).

        simple=False — full m=0,1,2 Graves expansion (Appendix B1):
            R(s,θ) = R₀ + (r + S₂)·cos θ + S₃·cos 2θ
            Z(s,θ) =      (r − S₂)·sin θ

        Parameters
        ----------
        theta  : 1-D array, shape (ntheta,)
        simple : bool

        Returns
        -------
        R, Z : 2-D arrays, shape (ns, ntheta)
        """
        R0 = self.R0[:, np.newaxis]
        r  = self.r[:, np.newaxis]

        if simple:
            R = R0 + r * np.cos(theta)
            Z = r  * np.sin(theta)
        else:
            S2 = self.S2[:, np.newaxis]
            S3 = self.S3[:, np.newaxis]
            R = R0 + (r + S2) * np.cos(theta) + S3 * np.cos(2 * theta)
            Z =      (r - S2) * np.sin(theta)
        return R, Z

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
            'eps':   self.eps,
            'r':     self.r,
            'S2':    self.S2,
            'S3':    self.S3,
            'kappa': self.kappa,
            'delta': self.delta,
            'shift': self.shift,
            'F2':    self.F2,
            'DI':    self.di,
            'mach':  self.mach,
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

    def save(self, outpath=None, fmt='npz'):
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

def compute_shaping(run_name_or_path, index=None, fmt='npz', outdir=None, save=True):
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
# Mach-scan fitting
# ---------------------------------------------------------------------------

#: Shaping quantities available for fitting and their display labels.
_FIT_QUANTITIES = {
    'kappa': r'$\kappa$',
    'delta': r'$\delta$',
    'shift': r'$\Delta$ [m]',
    'eps':   r'$\varepsilon$',
    'F2':    r'$F_2$',
    'DI':    r'$D_I$',
}


def fit_shaping_scan(npz_paths, s_index=1, degree=4, plot=True, outpath=None):
    """
    Fit shaping coefficients as polynomial functions of Mach number M.

    Loads a collection of shaping NPZ files (one per equilibrium), extracts
    the value of each coefficient at a chosen flux surface, and fits

        quantity(M) = a0 + a1·M + a2·M² + … + an·M^degree

    The x-axis Mach number is always taken from mach[0] (on-axis M₀ = √machsq),
    independent of s_index.

    Parameters
    ----------
    npz_paths : list of str or Path
        Shaping NPZ files to load (one per equilibrium/Mach number).
    s_index : int
        Flux-surface index at which to evaluate shaping quantities. Default 1
        (near-axis, avoiding the s=0 singularity).
    degree : int
        Polynomial degree in M (default 4).
    plot : bool
        Show a figure with data and fits for each quantity.
    outpath : str or Path, optional
        If given, save the fit coefficients to an NPZ file at this path.
        The file contains arrays '{qty}_coeffs' and '{qty}_perr' for each
        fitted quantity, readable by RealStability._load_shaping().

    Returns
    -------
    dict mapping quantity name → dict with keys 'coeffs', 'perr', 'mach', 'values'.
    """
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt

    # --- load all files ---
    npz_paths = [Path(p).expanduser() for p in npz_paths]
    records = []
    for p in sorted(npz_paths):
        d = np.load(str(p))
        mach_val = float(d['mach'][0]) if 'mach' in d else 0.0
        records.append((mach_val, d))

    if not records:
        print("fit_shaping_scan: no NPZ files loaded — check the path and that compute_shaping has been run.")
        return {}

    print(f"fit_shaping_scan: loaded {len(records)} files, M range [{min(r[0] for r in records):.3f}, {max(r[0] for r in records):.3f}]")

    records.sort(key=lambda x: x[0])
    machs = np.array([r[0] for r in records])

    results = {}

    # --- fit each quantity ---
    def poly_model(M, *coeffs):
        return sum(c * M**k for k, c in enumerate(coeffs))

    for qty, label in _FIT_QUANTITIES.items():
        vals = []
        for _, d in records:
            if qty in d:
                vals.append(float(d[qty][s_index]))
            else:
                vals.append(np.nan)
        vals = np.array(vals)

        mask = np.isfinite(vals) & np.isfinite(machs)
        if mask.sum() < degree + 1:
            print(f"  {qty}: skipped — only {mask.sum()} finite points, need {degree + 1} for degree-{degree} fit")
            continue

        p0 = np.zeros(degree + 1)
        p0[0] = np.nanmean(vals)
        try:
            popt, pcov = curve_fit(poly_model, machs[mask], vals[mask], p0=p0)
            perr = np.sqrt(np.diag(pcov))
        except RuntimeError:
            popt = p0
            perr = np.full_like(p0, np.nan)

        results[qty] = {'coeffs': popt, 'perr': perr, 'mach': machs, 'values': vals}

        terms = [f'{popt[0]:.4f}']
        for k in range(1, degree + 1):
            terms.append(f'({popt[k]:+.4f})·M^{k}')
        print(f"{qty:6s}(M) = {' '.join(terms)}")
        for k, (c, e) in enumerate(zip(popt, perr)):
            print(f"       a{k} = {c:.4f} ± {e:.4f}")

    # --- plot ---
    if plot and results:
        n = len(results)
        ncols = min(3, n)
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows), squeeze=False)

        M_dense = np.linspace(machs.min(), machs.max(), 200)

        for ax, (qty, res) in zip(axes.flat, results.items()):
            label = _FIT_QUANTITIES[qty]
            ax.scatter(res['mach'], res['values'], zorder=3, label='data')
            ax.plot(M_dense, poly_model(M_dense, *res['coeffs']), label='fit')
            ax.set_xlabel('M')
            ax.set_ylabel(label)
            ax.set_title(label)
            ax.legend(fontsize=8)

        for ax in axes.flat[len(results):]:
            ax.set_visible(False)

        fig.suptitle(f'Shaping vs Mach  (s-index {s_index})', y=1.01)
        fig.tight_layout()
        plt.show()

    if outpath is not None:
        save_dict = {}
        for qty, res in results.items():
            save_dict[f'{qty}_coeffs'] = res['coeffs']
            save_dict[f'{qty}_perr']   = res['perr']
        np.savez(str(outpath), **save_dict)
        print(f"  → fit coefficients saved to {outpath}")

    return results


# ---------------------------------------------------------------------------
# Shafranov shift vs Mach comparison plot (VMEC + VENUS-MHD)
# ---------------------------------------------------------------------------

def plot_shafranov_vs_mach(vmec_npz_paths=None, venus_h5_paths=None,
                           pluto_reader=None, pluto_scanparam=None, pluto_paramspecs=None,
                           s_value=None, ax=None):
    """
    Plot the Shafranov shift vs on-axis Mach number from VMEC, VENUS-MHD, and/or PlutoMHD.

    VMEC data is read from shaping NPZ files produced by compute_shaping().
    VENUS data is read directly from the stability h5 files written by Stability.Saveh5().
    PlutoMHD data is read via plutomhd.Reader.plutoread and the analytic shaf() formula.

    Both VMEC and VENUS plot Δ/R₀ (dimensionless Shafranov shift) and ε² (squared inverse
    aspect ratio) against on-axis Mach number.  PlutoMHD plots Δ/R₀ only (ε² is a fixed
    geometric quantity in the circular step model).

    The VENUS shift Δ/R₀(s) = R_mid(s) − R_mid(1), where R_mid = (R_max + R_min)/2
    at the chosen surface, normalised by R₀.  This gives the outward displacement of
    flux surface s relative to the plasma boundary (Δ = 0 at s = 1 by convention).
    The VENUS ε = (R_max − R_min)/(2 ⟨R⟩_axis) at the chosen surface.
    The VMEC shift is corrected analogously: Δ/R₀(s) = [shift(1) − shift(s)] / R₀,
    where shift(s) = R₀(0) − R₀(s) is stored in the shaping NPZ (zero at axis,
    maximum at boundary), and ε(s) = a(s)/R₀(s) from the Fourier m=0,1 coefficients.

    For PlutoMHD, the analytic Shafranov shift from plutorlstab.shaf() is evaluated at
    r = sqrt(s_value) (using s ≈ (r/a)² for circular surfaces) and normalised as
    Δ/R₀ = shaf(r) * eps_a, where shaf is in units of the minor radius a.

    VMEC and VENUS use different radial grids (different resolution and bunching),
    so a physical s value is used to locate the right point in each grid independently.

    Parameters
    ----------
    vmec_npz_paths : list of str or Path, optional
        Shaping NPZ files (one per VMEC equilibrium / Mach value).
    venus_h5_paths : list of str or Path, optional
        VENUS stability h5 files (one per run point).
    pluto_reader : plutoread object, str, or Path, optional
        A plutomhd.Reader.plutoread instance, or the full path to a PlutoMHD .npz
        data file (e.g. ``/path/to/Output/myrun/myrun.npz``).
    pluto_scanparam : str, optional
        Name of the scan parameter that sweeps Mach number (e.g. ``'omega0'`` or
        ``'omega_avg'``).  Auto-detected from the first entry of info['scanorder']
        if not given.
    pluto_paramspecs : dict, optional
        Fixed parameter values used to select a slice of a multi-dimensional
        PlutoMHD scan (e.g. ``{'omega_step': 0.0}``).
    s_value : float or None
        Normalised flux label s ∈ [0, 1] at which to evaluate the shift. Each
        grid is searched independently for its nearest point. Default None uses
        the outermost point in each grid (s → 1, i.e. the boundary).
    ax : matplotlib.axes.Axes, optional
        Axes to plot into. A new figure is created if None.

    Returns
    -------
    fig, ax
    """
    import h5py
    import matplotlib.pyplot as plt

    def _nearest_idx(s_arr, s_target):
        """Index of the grid point nearest to s_target."""
        return int(np.argmin(np.abs(s_arr - s_target)))

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4))
    else:
        fig = ax.get_figure()

    # Resolve the effective s-value: explicit > r0 from PlutoMHD > outermost
    s_eff = s_value
    r0_eff = None
    if s_eff is None and pluto_reader is not None:
        try:
            from plutomhd.Reader import plutoread
        except ImportError:
            pass
        else:
            _pr = pluto_reader if isinstance(pluto_reader, plutoread) else plutoread(
                Path(pluto_reader).expanduser().stem,
                filePath=Path(pluto_reader).expanduser().parent)
            _sp0 = list(_pr.info['scanparams'][
                _pr.info['scanorder'][0] if pluto_scanparam is None else pluto_scanparam])[0]
            _specs0 = {**(pluto_paramspecs or {}),
                       (_pr.info['scanorder'][0] if pluto_scanparam is None else pluto_scanparam): _sp0}
            _r0 = _pr('r0', _specs0)
            if _r0 is not None:
                r0_eff = float(_r0)
                s_eff = r0_eff ** 2

    if s_eff is not None:
        s_label = f'r₀ = {r0_eff:.3f}  (s = {s_eff:.3f})' if r0_eff is not None else f's = {s_eff:.3f}'
    else:
        s_label = 'outermost surface'

    # --- VMEC ---
    if vmec_npz_paths:
        vmec_npz_paths = [Path(p).expanduser() for p in vmec_npz_paths]
        machs_v, shifts_v, epssq_v = [], [], []
        for p in vmec_npz_paths:
            d = np.load(str(p))
            if not {'mach', 'shift', 's', 'R0', 'eps'}.issubset(d.files):
                print(f"  VMEC: skipping {p.name} — missing required key(s)")
                continue
            idx = _nearest_idx(d['s'], s_eff) if s_eff is not None else -1
            R0_axis = float(d['R0'][0])
            machs_v.append(float(d['mach'][0]))
            # shift(s) = R_axis − R_mid(s) = Δ_axis − Δ(s); subtract to get actual Δ(s)/R0
            shifts_v.append((float(d['shift'][-1]) - float(d['shift'][idx])) / R0_axis)
            epssq_v.append(float(d['eps'][idx])**2)
        if machs_v:
            order = np.argsort(machs_v)
            machs_v  = np.array(machs_v)[order]
            shifts_v = np.array(shifts_v)[order]
            epssq_v  = np.array(epssq_v)[order]
            ax.plot(machs_v, shifts_v, 'o-',  label=r'VMEC $\Delta/R_0$')
            ax.plot(machs_v, epssq_v,  'o--', label=r'VMEC $\varepsilon^2$')
            print(f"VMEC: {len(machs_v)} points, "
                  f"M ∈ [{machs_v.min():.3f}, {machs_v.max():.3f}], "
                  f"Δ/R₀ ∈ [{shifts_v.min():.4f}, {shifts_v.max():.4f}], "
                  f"ε² ∈ [{epssq_v.min():.4f}, {epssq_v.max():.4f}]")

    # --- VENUS-MHD ---
    if venus_h5_paths:
        venus_h5_paths = [Path(p).expanduser() for p in venus_h5_paths]
        machs_h, shifts_h, epssq_h = [], [], []
        for p in venus_h5_paths:
            with h5py.File(str(p), 'r') as f:
                M02  = float(f['normalisation']['M02'][()])
                s_v  = f['Grid']['S'][()]          # VENUS radial grid, shape (Nsurf,)
                R    = f['geometry']['R'][()]       # shape (Ntheta, Nsurf), normalised by R0
            idx = _nearest_idx(s_v, s_eff) if s_eff is not None else -1
            R_axis_n  = np.mean(R[:, 0])                              # normalised axis R
            R_col_n   = R[:, idx]                                     # normalised R at chosen surface
            R_center_n = (R_col_n.max() + R_col_n.min()) / 2
            R_edge_n   = (R[:, -1].max() + R[:, -1].min()) / 2       # boundary midpoint R
            eps_venus  = (R_col_n.max() - R_col_n.min()) / 2 / R_axis_n
            machs_h.append(np.sqrt(max(M02, 0.0)))
            # R_center_n − R_edge_n = Δ(s)/R0 (shift at s relative to boundary)
            shifts_h.append(R_center_n - R_edge_n)
            epssq_h.append(eps_venus**2)
        if machs_h:
            order = np.argsort(machs_h)
            machs_h  = np.array(machs_h)[order]
            shifts_h = np.array(shifts_h)[order]
            epssq_h  = np.array(epssq_h)[order]
            ax.plot(machs_h, shifts_h, 's-',  label=r'VENUS $\Delta/R_0$')
            ax.plot(machs_h, epssq_h,  's--', label=r'VENUS $\varepsilon^2$')
            print(f"VENUS: {len(machs_h)} points, "
                  f"M ∈ [{machs_h.min():.3f}, {machs_h.max():.3f}], "
                  f"Δ/R₀ ∈ [{shifts_h.min():.4f}, {shifts_h.max():.4f}], "
                  f"ε² ∈ [{epssq_h.min():.4f}, {epssq_h.max():.4f}]")

    # --- PlutoMHD analytic ---
    if pluto_reader is not None:
        try:
            from plutomhd.Reader import plutoread
        except ImportError:
            sys.exit("plutomhd package not found; ensure it is on the Python path")

        if not isinstance(pluto_reader, plutoread):
            p = Path(pluto_reader).expanduser()
            reader = plutoread(p.stem, filePath=p.parent)
        else:
            reader = pluto_reader

        if pluto_scanparam is None:
            pluto_scanparam = reader.info['scanorder'][0]

        spar_list = reader.info['scanparams'][pluto_scanparam]
        base_specs = dict(pluto_paramspecs or {})
        if r0_eff is not None:
            r_eval = np.array([r0_eff])
        elif s_value is not None:
            r_eval = np.array([np.sqrt(s_value)])
        else:
            r_eval = np.array([1.0])

        machs_p, shifts_p = [], []
        for spar in spar_list:
            paramSpecs = {**base_specs, pluto_scanparam: spar}
            shaf_arr, _, _ = reader.get_shaf(paramSpecs, r=r_eval)
            if shaf_arr is None:
                continue
            mach0 = reader('mach0', paramSpecs)
            eps_a = reader('eps_a', paramSpecs)
            if mach0 is None or eps_a is None:
                continue
            machs_p.append(float(mach0))
            shifts_p.append(-float(shaf_arr[0]) * float(eps_a)**3)
        if machs_p:
            order = np.argsort(machs_p)
            machs_p  = np.array(machs_p)[order]
            shifts_p = np.array(shifts_p)[order]
            ax.plot(machs_p, shifts_p, '^-', label=r'PlutoMHD analytic $\Delta/R_0$')
            print(f"PlutoMHD: {len(machs_p)} points, "
                  f"M ∈ [{machs_p.min():.3f}, {machs_p.max():.3f}], "
                  f"Δ/R₀ ∈ [{shifts_p.min():.4f}, {shifts_p.max():.4f}]")

    ax.set_xlabel(r'$\mathcal{M}$')
    ax.set_ylabel(r'$\Delta/R_0,\ \varepsilon^2$')
    ax.set_title(f'Shafranov shift vs Mach  ({s_label})')
    ax.legend()
    fig.tight_layout()
    plt.show()
    return fig, ax


# ---------------------------------------------------------------------------
# Script entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import argparse
    import glob

    parser = argparse.ArgumentParser(
        description='Compute analytic shaping coefficients from VMEC wout files.')
    parser.add_argument('target', nargs='?', default=None,
                        help='Path to wout .nc file, run name, or glob of shaping NPZ files (with --fit)')
    parser.add_argument('--index', type=int, default=None,
                        help='Wout file index (default: all indices for a run name)')
    parser.add_argument('--fmt', choices=['csv', 'npz'], default='npz',
                        help='Output file format (default: npz)')
    parser.add_argument('--outdir', default=None,
                        help='Output directory for saved files')
    parser.add_argument('--no-save', action='store_true',
                        help='Print summary only, do not write files')
    parser.add_argument('--fit', action='store_true',
                        help='Fit shaping coefficients vs Mach number across a scan')
    parser.add_argument('--s-value', type=float, default=None,
                        help='(--shafranov) Normalised flux label s ∈ [0,1] at which to '
                             'evaluate Δ. Each grid is searched independently. '
                             'Default: outermost point in each grid.')
    parser.add_argument('--s-index', type=int, default=1,
                        help='(--fit) Flux-surface array index for fit_shaping_scan '
                             '(default: 1, near-axis). Not used by --shafranov.')
    parser.add_argument('--degree', type=int, default=4,
                        help='Polynomial degree in M for fitting (default: 4)')
    parser.add_argument('--no-plot', action='store_true',
                        help='Skip the fit plot')
    parser.add_argument('--shafranov', action='store_true',
                        help='Plot Shafranov shift vs Mach number comparing VMEC, VENUS-MHD, and/or PlutoMHD')
    parser.add_argument('--vmec-npz', default=None,
                        help='Glob of VMEC shaping NPZ files for --shafranov (e.g. "run/wout/*_shaping.npz")')
    parser.add_argument('--venus-h5', default=None,
                        help='Glob of VENUS stability h5 files for --shafranov (e.g. "run/*.h5")')
    parser.add_argument('--pluto-npz', default=None,
                        help='Full path to a PlutoMHD .npz data file for --shafranov '
                             '(e.g. "Output/myrun/myrun.npz")')
    parser.add_argument('--pluto-scanparam', default=None,
                        help='PlutoMHD scan parameter name for the Mach axis '
                             '(default: first entry of info["scanorder"])')
    parser.add_argument('--pluto-fixed', nargs='*', default=None, metavar='KEY=VALUE',
                        help='Fixed PlutoMHD parameters for multi-dim scans, '
                             'e.g. --pluto-fixed omega_step=0.0 rho0=2.0')
    args = parser.parse_args()

    if args.shafranov:
        vmec_paths  = sorted(glob.glob(args.vmec_npz))  if args.vmec_npz  else []
        venus_paths = sorted(glob.glob(args.venus_h5))  if args.venus_h5  else []
        pluto_fixed = {}
        if args.pluto_fixed:
            for item in args.pluto_fixed:
                k, v = item.split('=', 1)
                try:
                    pluto_fixed[k] = float(v)
                except ValueError:
                    pluto_fixed[k] = v
        if not vmec_paths and not venus_paths and not args.pluto_npz:
            sys.exit("Provide at least one of --vmec-npz, --venus-h5, or --pluto-npz")
        plot_shafranov_vs_mach(
            vmec_npz_paths=vmec_paths or None,
            venus_h5_paths=venus_paths or None,
            pluto_reader=args.pluto_npz,
            pluto_scanparam=args.pluto_scanparam,
            pluto_paramspecs=pluto_fixed or None,
            s_value=args.s_value,
        )
    elif args.fit:
        if not args.target:
            sys.exit("Provide a glob target for --fit")
        npz_paths = sorted(glob.glob(args.target))
        if not npz_paths:
            sys.exit(f"No files matched: {args.target}")
        fit_shaping_scan(npz_paths, s_index=args.s_index, degree=args.degree, plot=not args.no_plot)  # noqa: E501
    else:
        if not args.target:
            sys.exit("Provide a target (wout path or run name)")
        compute_shaping(
            args.target,
            index=args.index,
            fmt=args.fmt,
            outdir=args.outdir,
            save=not args.no_save,
        )
