#!/Users/cellywelly/Dev/solis/bin/python3
"""
Plots VMEC equilibrium outputs from wout NetCDF files.

Since all current runs are axisymmetric (ntor=0), the Fourier reconstruction
collapses to a 1D sum over poloidal modes: R(s,θ) = Σ_m rmnc[s,m]*cos(m*θ),
and the result is the same at every toroidal angle.

Available plots
---------------
  flux_surfaces         R-Z cross-section with flux surface contours coloured by |B|
  b_profile             |B|(θ) on a chosen set of flux surfaces
  q_profile             safety factor and iota vs normalised flux label s
  flux_surfaces_slider  interactive slider over all Mach numbers for a run
  summary_slider        same but all three panels update together
  profiles_slider       Ω, P, q, ι profiles with Mach-number slider

Usage as a module
-----------------
  from vmecPlotter import VmecPlotter, flux_surfaces_slider, summary_slider, profiles_slider
  p = VmecPlotter('Output/KH/test_oldvmec/VMEC/wout/wout_test_oldvmec_0.nc')
  p.flux_surfaces()
  p.b_profile(s_indices=[0, 70, 140, 210, 280])
  p.q_profile()
  p.summary_plot()

  flux_surfaces_slider('QI_bstep0.7_omstep0.5')
  summary_slider('QI_bstep0.7_omstep0.5')

Usage as a script
-----------------
  python3 vmecPlotter.py <wout.nc> [flux_surfaces|b_profile|q_profile]
  python3 vmecPlotter.py test_oldvmec               # run name, index 0
  python3 vmecPlotter.py test_oldvmec 1 b_profile   # run name, index 1
  python3 vmecPlotter.py test_oldvmec slider        # Mach-number slider
  python3 vmecPlotter.py test_oldvmec summary_slider
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.widgets import Slider
from pathlib import Path

from vmecReader import WoutReader, find_wout, list_wouts


# ---------------------------------------------------------------------------
# Fourier reconstruction helpers
# These implement the same logic as PySTEL's cfunct/sfunct without needing
# the libstell.so compiled library.
# ---------------------------------------------------------------------------

def _reconstruct_cos(theta, fmnc, xm, xn):
    """
    Evaluate f(s,θ,ζ=0) = Σ_{m,n} fmnc(s,mn) * cos(m*θ - n*0)
    Returns array of shape (ns, ntheta).
    For axisymmetric cases (all n=0) this is exact at any ζ.
    Note: uses np.dot rather than @ to avoid a spurious BLAS warning on macOS ARM.
    """
    mt = np.outer(xm, theta)   # (nmn, ntheta)
    nz = np.outer(xn, np.zeros_like(theta))  # n*ζ = 0
    cos_arg = np.cos(mt - nz)  # (nmn, ntheta)
    return np.dot(fmnc, cos_arg)  # (ns, ntheta)


def _reconstruct_sin(theta, fmns, xm, xn):
    """
    Evaluate f(s,θ,ζ=0) = Σ_{m,n} fmns(s,mn) * sin(m*θ - n*0).
    Returns array of shape (ns, ntheta).
    """
    mt = np.outer(xm, theta)
    nz = np.outer(xn, np.zeros_like(theta))
    sin_arg = np.sin(mt - nz)
    return np.dot(fmns, sin_arg)  # (ns, ntheta)


# ---------------------------------------------------------------------------
# Main plotter class
# ---------------------------------------------------------------------------

class VmecPlotter:
    def __init__(self, wout_path, ntheta=256):
        self.w = WoutReader(wout_path)
        self.ntheta = ntheta
        self.theta = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
        self._build_geometry()

    def _build_geometry(self):
        w = self.w
        theta = self.theta

        # R(s, θ) and Z(s, θ) from rmnc/zmns (standard-grid modes)
        self.R = _reconstruct_cos(theta, w.rmnc, w.xm, w.xn)   # (ns, ntheta)
        self.Z = _reconstruct_sin(theta, w.zmns, w.xm, w.xn)   # (ns, ntheta)

        # |B|(s, θ) from bmnc (Nyquist-grid modes)
        self.B = _reconstruct_cos(theta, w.bmnc, w.xm_nyq, w.xn_nyq)  # (ns, ntheta)

    # -----------------------------------------------------------------------
    # Plot 1: flux surfaces in R-Z, coloured by |B| on each surface
    # -----------------------------------------------------------------------
    def flux_surfaces(self, n_surfaces=12, ax=None, show=True, title=None):
        """
        R-Z cross-section with flux surface contours.
        Surfaces are evenly spaced in s ∈ (0,1], skipping the axis.
        Contour colour = mean |B| on each surface.
        """
        w = self.w
        ns = w.ns

        # Evenly spaced flux surfaces (skip s=0 axis, include edge)
        indices = np.round(np.linspace(1, ns - 1, n_surfaces)).astype(int)

        B_mean = self.B.mean(axis=1)   # mean |B| per surface
        B_min, B_max = B_mean[indices].min(), B_mean[indices].max()
        cmap = cm.plasma

        if ax is None:
            fig, ax = plt.subplots(figsize=(5, 7))
        else:
            fig = ax.get_figure()

        for idx in indices:
            colour = cmap((B_mean[idx] - B_min) / (B_max - B_min + 1e-30))
            ax.plot(
                np.append(self.R[idx], self.R[idx, 0]),
                np.append(self.Z[idx], self.Z[idx, 0]),
                color=colour, lw=0.8
            )

        # Axis point (s=0 → single R,Z value from m=0 mode)
        R_ax = w.rmnc[0, 0]   # rmnc(s=0, m=0, n=0)
        Z_ax = w.zmns[0, 0]   # should be 0 for symmetric case
        ax.plot(R_ax, Z_ax, '+', color='black', ms=6)

        sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(B_min, B_max))
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label='⟨|B|⟩ [T]', fraction=0.046, pad=0.04)

        ax.set_xlabel('R [m]')
        ax.set_ylabel('Z [m]')
        ax.set_aspect('equal')
        ax.set_title(title or Path(self.w.path).name)
        ax.grid(True, alpha=0.3)

        if show:
            plt.tight_layout()
            plt.show()
        return ax

    # -----------------------------------------------------------------------
    # Plot 2: |B| as a function of poloidal angle θ on selected surfaces
    # -----------------------------------------------------------------------
    def b_profile(self, s_indices=None, ax=None, show=True, title=None):
        """
        |B|(θ) on a set of flux surfaces.
        s_indices: list of surface indices (0 = axis, ns-1 = edge).
                   Defaults to 5 evenly spaced surfaces.
        """
        w = self.w
        ns = w.ns

        if s_indices is None:
            s_indices = np.round(np.linspace(1, ns - 1, 5)).astype(int).tolist()

        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 4))

        theta_deg = np.degrees(self.theta)
        cmap = cm.viridis
        colours = cmap(np.linspace(0, 1, len(s_indices)))

        for colour, idx in zip(colours, s_indices):
            s_val = w.s[idx]
            ax.plot(theta_deg, self.B[idx], color=colour, label=f's={s_val:.2f}')

        ax.set_xlabel('θ [deg]')
        ax.set_ylabel('|B| [T]')
        ax.set_title(title or f'|B|(θ) — {Path(self.w.path).name}')
        ax.set_xlim(0, 360)
        ax.set_xticks(np.arange(0, 361, 60))
        ax.legend(fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3)

        if show:
            plt.tight_layout()
            plt.show()
        return ax

    # -----------------------------------------------------------------------
    # Plot 3: safety factor q(s) and rotational transform iota(s)
    # -----------------------------------------------------------------------
    def q_profile(self, ax=None, show=True, title=None):
        """
        Safety factor q(s) and rotational transform iota(s) vs normalised flux.
        """
        w = self.w

        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))
        ax2 = ax.twinx()

        ax.plot(w.s, w.q,    color='tab:blue',   label='q')
        ax2.plot(w.s, w.iotaf, color='tab:orange', label='ι', linestyle='--')

        ax.set_xlabel('s = √(Φ/Φ_edge)')
        ax.set_ylabel('q', color='tab:blue')
        ax2.set_ylabel('ι (rotational transform)', color='tab:orange')
        ax.set_title(title or f'q and ι — {Path(self.w.path).name}')

        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc='best')
        ax.grid(True, alpha=0.3)

        if show:
            plt.tight_layout()
            plt.show()
        return ax, ax2

    # -----------------------------------------------------------------------
    # Convenience: all three plots in a single figure
    # -----------------------------------------------------------------------
    def summary_plot(self, show=True):
        fig = plt.figure(figsize=(15, 5))
        ax1 = fig.add_subplot(1, 3, 1)
        ax2 = fig.add_subplot(1, 3, 2)
        ax3 = fig.add_subplot(1, 3, 3)

        self.flux_surfaces(ax=ax1, show=False,
                           title=f'Flux surfaces\n{Path(self.w.path).name}')
        self.b_profile(ax=ax2, show=False, title='|B|(θ)')
        self.q_profile(ax=ax3, show=False, title='q and ι')

        plt.tight_layout()
        if show:
            plt.show()
        return fig


# ---------------------------------------------------------------------------
# Interactive Mach-number slider
# ---------------------------------------------------------------------------

def flux_surfaces_slider(run_name, run_dir=None, n_surfaces=12):
    """
    Interactive R-Z flux surface plot with a slider to scrub through all
    Mach numbers in a run.  Wout files are sorted by M (ascending).

    Parameters
    ----------
    run_name : str
        Run directory name (e.g. 'QI_test_bstep0.7').
    run_dir : path-like, optional
        Override the default output root / KH / <run_name> location.
    n_surfaces : int
        Number of flux surfaces to draw per frame.

    Example
    -------
    from vmecPlotter import flux_surfaces_slider
    flux_surfaces_slider('QI_test_bstep0.7')
    """
    # --- load & sort ---
    paths = list_wouts(run_name, run_dir=run_dir)
    if not paths:
        raise FileNotFoundError(f"No wout files found for run '{run_name}'")

    print(f"Loading {len(paths)} wout files for '{run_name}'…")
    plotters = [VmecPlotter(str(p)) for p in paths]
    plotters.sort(key=lambda p: p.w.machsq)
    N = len(plotters)
    mach_vals = [np.sqrt(p.w.machsq) for p in plotters]

    # --- precompute surface indices (identical for all files) ---
    ns = plotters[0].w.ns
    s_indices = np.round(np.linspace(1, ns - 1, n_surfaces)).astype(int)

    # --- fixed R/Z axis limits (use outermost surface across all files) ---
    all_R = np.concatenate([p.R[s_indices[-1]] for p in plotters])
    all_Z = np.concatenate([p.Z[s_indices[-1]] for p in plotters])
    R_pad = (all_R.max() - all_R.min()) * 0.05
    Z_pad = (all_Z.max() - all_Z.min()) * 0.05
    R_lim = (all_R.min() - R_pad, all_R.max() + R_pad)
    Z_lim = (all_Z.min() - Z_pad, all_Z.max() + Z_pad)

    # --- figure layout ---
    fig = plt.figure(figsize=(6, 8))
    ax  = fig.add_axes([0.10, 0.15, 0.74, 0.78])
    cax = fig.add_axes([0.86, 0.15, 0.03, 0.78])
    ax_sl = fig.add_axes([0.15, 0.04, 0.70, 0.03])

    slider = Slider(ax_sl, 'M', 0, N - 1, valinit=0, valstep=1,
                    color='steelblue')

    # Tick the slider at integer positions, label with M values
    ax_sl.set_xticks(range(N))
    ax_sl.set_xticklabels([f'{m:.2f}' for m in mach_vals],
                           fontsize=7, rotation=45, ha='right')
    ax_sl.set_xlabel('Mach number  M', fontsize=8, labelpad=14)

    # --- initial draw: create all artists once, never clear the axes ---
    cmap   = cm.plasma
    p0     = plotters[0]
    B_mean = p0.B.mean(axis=1)
    B_min  = B_mean[s_indices].min()
    B_max  = B_mean[s_indices].max()
    norm0  = plt.Normalize(B_min, B_max)

    # one Line2D per flux surface + one for the magnetic axis marker
    lines = []
    for si in s_indices:
        R_closed = np.append(p0.R[si], p0.R[si, 0])
        Z_closed = np.append(p0.Z[si], p0.Z[si, 0])
        line, = ax.plot(R_closed, Z_closed, color=cmap(norm0(B_mean[si])), lw=0.9)
        lines.append(line)
    axis_marker, = ax.plot(p0.w.rmnc[0, 0], p0.w.zmns[0, 0],
                           '+', color='black', ms=7, mew=1.5)

    ax.set_xlim(*R_lim)
    ax.set_ylim(*Z_lim)
    ax.set_aspect('equal')
    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')
    ax.grid(True, alpha=0.3)

    sm   = cm.ScalarMappable(cmap=cmap, norm=norm0)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, label='⟨|B|⟩ [T]')

    title = ax.set_title(
        f'{run_name}\nM = {mach_vals[0]:.4f}   '
        f'M² = {p0.w.machsq:.4f}   '
        f'β = {p0.w.betatotal:.3e}',
        fontsize=9
    )

    def on_slide(val):
        idx    = int(round(val))
        p      = plotters[idx]
        B_mean = p.B.mean(axis=1)
        B_min  = B_mean[s_indices].min()
        B_max  = B_mean[s_indices].max()
        norm   = plt.Normalize(B_min, B_max)

        for line, si in zip(lines, s_indices):
            line.set_xdata(np.append(p.R[si], p.R[si, 0]))
            line.set_ydata(np.append(p.Z[si], p.Z[si, 0]))
            line.set_color(cmap(norm(B_mean[si])))

        axis_marker.set_xdata([p.w.rmnc[0, 0]])
        axis_marker.set_ydata([p.w.zmns[0, 0]])

        sm.set_clim(B_min, B_max)
        cbar.update_normal(sm)

        title.set_text(
            f'{run_name}\nM = {mach_vals[idx]:.4f}   '
            f'M² = {p.w.machsq:.4f}   '
            f'β = {p.w.betatotal:.3e}'
        )
        slider.valtext.set_text(f'{mach_vals[idx]:.4f}')
        fig.canvas.draw()

    slider.on_changed(on_slide)
    plt.show()
    return fig


def summary_slider(run_name, run_dir=None, n_surfaces=12, n_b_surfaces=5):
    """
    Three-panel interactive summary (flux surfaces, |B|(θ), q & ι) with a
    slider to scrub through all Mach numbers in a run.

    Parameters
    ----------
    run_name     : str   — run directory name
    run_dir      : path  — override default output root
    n_surfaces   : int   — flux surface contours in panel 1
    n_b_surfaces : int   — surfaces shown in |B|(θ) panel

    Example
    -------
    from vmecPlotter import summary_slider
    summary_slider('QI_bstep0.7_omstep0.5')
    """
    paths = list_wouts(run_name, run_dir=run_dir)
    if not paths:
        raise FileNotFoundError(f"No wout files found for run '{run_name}'")

    print(f"Loading {len(paths)} wout files for '{run_name}'…")
    plotters = [VmecPlotter(str(p)) for p in paths]
    plotters.sort(key=lambda p: p.w.machsq)
    N = len(plotters)
    mach_vals = [np.sqrt(p.w.machsq) for p in plotters]

    ns = plotters[0].w.ns
    surf_idx = np.round(np.linspace(1, ns - 1, n_surfaces)).astype(int)
    b_idx    = np.round(np.linspace(1, ns - 1, n_b_surfaces)).astype(int)

    # --- fixed axis limits computed across all runs ---
    all_R = np.concatenate([p.R[surf_idx[-1]] for p in plotters])
    all_Z = np.concatenate([p.Z[surf_idx[-1]] for p in plotters])
    R_pad = (all_R.max() - all_R.min()) * 0.05
    Z_pad = (all_Z.max() - all_Z.min()) * 0.05
    R_lim = (all_R.min() - R_pad, all_R.max() + R_pad)
    Z_lim = (all_Z.min() - Z_pad, all_Z.max() + Z_pad)

    B_all  = np.stack([p.B for p in plotters])          # (N, ns, ntheta)
    q_all  = np.stack([p.w.q for p in plotters])        # (N, ns)
    io_all = np.stack([p.w.iotaf for p in plotters])    # (N, ns)

    B_ylim  = (B_all[:, b_idx, :].min() * 0.995,
               B_all[:, b_idx, :].max() * 1.005)
    q_ylim  = (q_all.min() * 1.05,  q_all.max() * 0.95)   # q is negative
    io_ylim = (io_all.min() * 1.05, io_all.max() * 0.95)

    # --- figure layout ---
    fig = plt.figure(figsize=(15, 7))
    # main panels
    ax1 = fig.add_axes([0.05, 0.16, 0.22, 0.74])   # flux surfaces
    cax = fig.add_axes([0.28, 0.16, 0.01, 0.74])   # colorbar
    ax2 = fig.add_axes([0.38, 0.16, 0.25, 0.74])   # |B|(θ)
    ax3 = fig.add_axes([0.72, 0.16, 0.25, 0.74])   # q
    ax3b = ax3.twinx()                              # iota
    ax_sl = fig.add_axes([0.15, 0.04, 0.70, 0.025])

    slider = Slider(ax_sl, 'M', 0, N - 1, valinit=0, valstep=1,
                    color='steelblue')
    ax_sl.set_xticks(range(N))
    ax_sl.set_xticklabels([f'{m:.2f}' for m in mach_vals],
                           fontsize=7, rotation=45, ha='right')
    ax_sl.set_xlabel('Mach number  M', fontsize=8, labelpad=14)

    theta_deg = np.degrees(plotters[0].theta)
    cmap_fs  = cm.plasma
    cmap_b   = cm.viridis

    # --- panel 1: flux surfaces (initial) ---
    p0     = plotters[0]
    B_mean = p0.B.mean(axis=1)
    norm0  = plt.Normalize(B_mean[surf_idx].min(), B_mean[surf_idx].max())
    fs_lines = []
    for si in surf_idx:
        line, = ax1.plot(
            np.append(p0.R[si], p0.R[si, 0]),
            np.append(p0.Z[si], p0.Z[si, 0]),
            color=cmap_fs(norm0(B_mean[si])), lw=0.9
        )
        fs_lines.append(line)
    fs_axis, = ax1.plot(p0.w.rmnc[0, 0], p0.w.zmns[0, 0],
                        '+', color='black', ms=7, mew=1.5)
    ax1.set_xlim(*R_lim); ax1.set_ylim(*Z_lim)
    ax1.set_aspect('equal'); ax1.set_xlabel('R [m]'); ax1.set_ylabel('Z [m]')
    ax1.grid(True, alpha=0.3)

    sm   = cm.ScalarMappable(cmap=cmap_fs, norm=norm0)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, label='⟨|B|⟩ [T]')

    # --- panel 2: |B|(θ) (initial) ---
    b_colours = cmap_b(np.linspace(0, 1, n_b_surfaces))
    b_lines = []
    for si, colour in zip(b_idx, b_colours):
        line, = ax2.plot(theta_deg, p0.B[si], color=colour,
                         label=f's={p0.w.s[si]:.2f}')
        b_lines.append(line)
    ax2.set_xlim(0, 360); ax2.set_xticks(np.arange(0, 361, 60))
    ax2.set_ylim(*B_ylim)
    ax2.set_xlabel('θ [deg]'); ax2.set_ylabel('|B| [T]')
    ax2.set_title('|B|(θ)', fontsize=9); ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=7, ncol=2)

    # --- panel 3: q and ι (initial) ---
    q_line,  = ax3.plot(p0.w.s,  p0.w.q,     color='tab:blue',   label='q')
    io_line, = ax3b.plot(p0.w.s, p0.w.iotaf, color='tab:orange',
                          linestyle='--', label='ι')
    ax3.set_xlim(0, 1); ax3.set_ylim(*q_ylim)
    ax3b.set_ylim(*io_ylim)
    ax3.set_xlabel('s'); ax3.set_ylabel('q', color='tab:blue')
    ax3b.set_ylabel('ι', color='tab:orange')
    ax3.set_title('q and ι', fontsize=9); ax3.grid(True, alpha=0.3)
    lines_leg  = [q_line, io_line]
    labels_leg = ['q', 'ι']
    ax3.legend(lines_leg, labels_leg, fontsize=8)

    suptitle = fig.suptitle(
        f'{run_name}   M = {mach_vals[0]:.4f}   '
        f'M² = {p0.w.machsq:.4f}   β = {p0.w.betatotal:.3e}',
        fontsize=10, y=0.97
    )

    def on_slide(val):
        idx    = int(round(val))
        p      = plotters[idx]
        B_mean = p.B.mean(axis=1)
        B_min  = B_mean[surf_idx].min()
        B_max  = B_mean[surf_idx].max()
        norm   = plt.Normalize(B_min, B_max)

        # panel 1
        for line, si in zip(fs_lines, surf_idx):
            line.set_xdata(np.append(p.R[si], p.R[si, 0]))
            line.set_ydata(np.append(p.Z[si], p.Z[si, 0]))
            line.set_color(cmap_fs(norm(B_mean[si])))
        fs_axis.set_xdata([p.w.rmnc[0, 0]])
        fs_axis.set_ydata([p.w.zmns[0, 0]])
        sm.set_clim(B_min, B_max)
        cbar.update_normal(sm)

        # panel 2
        for line, si in zip(b_lines, b_idx):
            line.set_ydata(p.B[si])

        # panel 3
        q_line.set_ydata(p.w.q)
        io_line.set_ydata(p.w.iotaf)

        suptitle.set_text(
            f'{run_name}   M = {mach_vals[idx]:.4f}   '
            f'M² = {p.w.machsq:.4f}   β = {p.w.betatotal:.3e}'
        )
        slider.valtext.set_text(f'{mach_vals[idx]:.4f}')
        fig.canvas.draw()

    slider.on_changed(on_slide)
    plt.show()
    return fig


def profiles_slider(run_name, run_dir=None):
    """
    Four-panel profile plot (Ω, P, q, ι) with a slider to scrub through all
    Mach numbers in a run.  Y-limits are fixed across all files so changes
    are visible without rescaling.

    Panels
    ------
    top-left   : Ω(s)  — rotation profile
    top-right  : P(s)  — pressure profile
    bottom-left: q(s)  — safety factor
    bottom-right: ι(s) — rotational transform

    Example
    -------
    from vmecPlotter import profiles_slider
    profiles_slider('QI_bstep0.7_omstep0.5')
    """
    paths = list_wouts(run_name, run_dir=run_dir)
    if not paths:
        raise FileNotFoundError(f"No wout files found for run '{run_name}'")

    print(f"Loading {len(paths)} wout files for '{run_name}'…")
    plotters = [VmecPlotter(str(p)) for p in paths]
    plotters.sort(key=lambda p: p.w.machsq)
    N = len(plotters)
    mach_vals = [np.sqrt(p.w.machsq) for p in plotters]

    # --- fixed y-limits across all runs ---
    def _lim(arr, pad=0.08, ensure_zero_bottom=False):
        lo, hi = arr.min(), arr.max()
        span = hi - lo if hi != lo else abs(hi) * 0.1 or 0.1
        lo2, hi2 = lo - pad * span, hi + pad * span
        if ensure_zero_bottom:
            lo2 = min(lo2, -pad * span)
        return lo2, hi2

    all_omega = np.stack([p.w.omega for p in plotters])
    all_pres  = np.stack([p.w.presf for p in plotters])
    all_q     = np.stack([p.w.q     for p in plotters])
    all_iota  = np.stack([p.w.iotaf for p in plotters])

    lim_omega = _lim(all_omega, ensure_zero_bottom=True)
    lim_pres  = _lim(all_pres,  ensure_zero_bottom=True)
    lim_q     = _lim(all_q)
    lim_iota  = _lim(all_iota)

    s = plotters[0].w.s

    # --- figure layout: 2×2 panels + slider strip ---
    fig = plt.figure(figsize=(11, 8))
    ax_om  = fig.add_axes([0.08, 0.47, 0.38, 0.44])   # top-left:  Ω
    ax_pr  = fig.add_axes([0.57, 0.47, 0.38, 0.44])   # top-right: P
    ax_q   = fig.add_axes([0.08, 0.13, 0.38, 0.27])   # bot-left:  q
    ax_io  = fig.add_axes([0.57, 0.13, 0.38, 0.27])   # bot-right: ι
    ax_sl  = fig.add_axes([0.15, 0.03, 0.70, 0.025])

    slider = Slider(ax_sl, 'M', 0, N - 1, valinit=0, valstep=1,
                    color='steelblue')
    ax_sl.set_xticks(range(N))
    ax_sl.set_xticklabels([f'{m:.2f}' for m in mach_vals],
                           fontsize=7, rotation=45, ha='right')
    ax_sl.set_xlabel('Mach number  M', fontsize=8, labelpad=14)

    # --- initial draw ---
    p0 = plotters[0]

    om_line, = ax_om.plot(s, p0.w.omega, color='tab:blue',   lw=1.5)
    pr_line, = ax_pr.plot(s, p0.w.presf, color='tab:red',    lw=1.5)
    q_line,  = ax_q.plot( s, p0.w.q,    color='tab:blue',   lw=1.5)
    io_line, = ax_io.plot(s, p0.w.iotaf, color='tab:orange', lw=1.5)

    for ax, lim, ylabel, title in [
        (ax_om, lim_omega, 'Ω (normalised)', 'Rotation  Ω(s)'),
        (ax_pr, lim_pres,  'P [Pa]',         'Pressure  P(s)'),
        (ax_q,  lim_q,     'q',              'Safety factor  q(s)'),
        (ax_io, lim_iota,  'ι',              'Rotational transform  ι(s)'),
    ]:
        ax.set_xlim(0, 1)
        ax.set_ylim(*lim)
        ax.set_xlabel('s', fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(title, fontsize=9)
        ax.grid(True, alpha=0.3)

    suptitle = fig.suptitle(
        f'{run_name}   M = {mach_vals[0]:.4f}   '
        f'M² = {p0.w.machsq:.4f}   β = {p0.w.betatotal:.3e}',
        fontsize=10, y=0.99
    )

    def on_slide(val):
        idx = int(round(val))
        p   = plotters[idx]

        om_line.set_ydata(p.w.omega)
        pr_line.set_ydata(p.w.presf)
        q_line.set_ydata(p.w.q)
        io_line.set_ydata(p.w.iotaf)

        suptitle.set_text(
            f'{run_name}   M = {mach_vals[idx]:.4f}   '
            f'M² = {p.w.machsq:.4f}   β = {p.w.betatotal:.3e}'
        )
        slider.valtext.set_text(f'{mach_vals[idx]:.4f}')
        fig.canvas.draw()

    slider.on_changed(on_slide)
    plt.show()
    return fig


# ---------------------------------------------------------------------------
# Script entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)

    arg = sys.argv[1]
    path = Path(arg)

    # slider / summary_slider mode
    if len(sys.argv) > 2 and sys.argv[2] == 'slider':
        flux_surfaces_slider(arg)
        sys.exit(0)
    if len(sys.argv) > 2 and sys.argv[2] == 'summary_slider':
        summary_slider(arg)
        sys.exit(0)
    if len(sys.argv) > 2 and sys.argv[2] == 'profiles_slider':
        profiles_slider(arg)
        sys.exit(0)
    if len(sys.argv) == 2 and not path.suffix:
        # bare run name with no plot type → open summary_slider if multiple files exist
        candidate_paths = list_wouts(arg)
        if len(candidate_paths) > 1:
            summary_slider(arg)
            sys.exit(0)

    # Determine wout path for single-file plots
    if path.suffix == '.nc' and path.exists():
        wout_path = path
        plot_type = sys.argv[2] if len(sys.argv) > 2 else 'summary'
    else:
        idx = 0
        plot_type = 'summary'
        if len(sys.argv) > 2:
            try:
                idx = int(sys.argv[2])
                plot_type = sys.argv[3] if len(sys.argv) > 3 else 'summary'
            except ValueError:
                plot_type = sys.argv[2]
        wout_path = find_wout(arg, index=idx)
        if not wout_path.exists():
            print(f"ERROR: {wout_path} not found")
            sys.exit(1)

    p = VmecPlotter(str(wout_path))

    if plot_type == 'flux_surfaces':
        p.flux_surfaces()
    elif plot_type == 'b_profile':
        p.b_profile()
    elif plot_type == 'q_profile':
        p.q_profile()
    else:
        p.summary_plot()
