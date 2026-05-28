#!/Users/cellywelly/Dev/solis/bin/python3
"""
Reader for VMEC wout NetCDF output files.

Wout files live at:
    Output/KH/<run_name>/VMEC/wout/wout_<run_name>_<index>.nc

Usage as a module:
    from vmecReader import WoutReader
    w = WoutReader('Output/KH/test_oldvmec/VMEC/wout/wout_test_oldvmec_0.nc')
    print(w.q)          # safety factor profile
    print(w.summary())  # human-readable scalar summary

Usage as a script:
    python3 vmecReader.py Output/KH/test_oldvmec/VMEC/wout/wout_test_oldvmec_0.nc
"""

import sys
import socket
import numpy as np
from pathlib import Path
from scipy.io import netcdf_file


OUTPUT_ROOT_VIKING = Path('/users/cs2427/scratch/lunamhd-data')
OUTPUT_ROOT_LOCAL  = Path('/Users/cellywelly/Dev/lunamhd/Output')


def _output_root():
    return OUTPUT_ROOT_VIKING if 'viking' in socket.gethostname() else OUTPUT_ROOT_LOCAL


def find_wout(run_name, index=0, run_dir=None):
    """Return the Path to a specific wout file for a named run."""
    if run_dir is None:
        run_dir = _output_root() / 'KH' / run_name
    run_dir = Path(run_dir)
    return run_dir / 'VMEC' / 'wout' / f'wout_{run_name}_{index}.nc'


def list_wouts(run_name, run_dir=None):
    """Return a sorted list of all wout Paths for a named run."""
    if run_dir is None:
        run_dir = _output_root() / 'KH' / run_name
    wout_dir = Path(run_dir) / 'VMEC' / 'wout'
    return sorted(wout_dir.glob(f'wout_{run_name}_*.nc'))


class WoutReader:
    """
    Loads key physics quantities from a VMEC wout NetCDF file.

    Scalars
    -------
    ns, nfp, mpol, ntor, mnmax, mnmax_nyq
    aspect, Aminor, Rmajor, volume
    b0, rbtor0, betatotal, betapol, betator, betaxis
    machsq          (Mach number squared on axis; 0 if non-flow version)
    ier_flag        (0 = converged)

    1-D radial profiles  (length ns, on full-mesh flux surfaces)
    ------------------
    s               normalised flux label sqrt(phi/phi_edge)  [0, 1]
    iotaf           rotational transform
    q               safety factor
    presf           pressure profile
    phi             toroidal flux
    chi             poloidal flux
    omega           rotation profile (0 if non-flow VMEC)
    jcuru, jcurv    poloidal and toroidal current densities

    Mode arrays
    -----------
    xm, xn         poloidal/toroidal mode numbers  (mn_mode)
    xm_nyq, xn_nyq Nyquist mode numbers            (mn_mode_nyq)

    2-D Fourier coefficients  shape (ns, mn_mode)
    -------------------------
    rmnc, zmns, lmns

    2-D Nyquist Fourier coefficients  shape (ns, mn_mode_nyq)
    --------------------------------
    bmnc, gmnc
    bsubumnc, bsubvmnc, bsubsmns
    bsupumnc, bsupvmnc
    """

    def __init__(self, wout_path):
        self.path = Path(wout_path)
        self._load()

    def _load(self):
        f = netcdf_file(str(self.path), 'r', mmap=False)
        v = f.variables

        def get(name):
            # .copy() then astype to convert big-endian NetCDF bytes to native float64
            return v[name][()].copy().astype(float)

        # --- scalars ---
        self.ns       = int(get('ns'))
        self.nfp      = int(get('nfp'))
        self.mpol     = int(get('mpol'))
        self.ntor     = int(get('ntor'))
        self.mnmax    = int(get('mnmax'))
        self.mnmax_nyq = int(get('mnmax_nyq'))
        self.aspect   = float(get('aspect'))
        self.Aminor   = float(get('Aminor_p'))
        self.Rmajor   = float(get('Rmajor_p'))
        self.volume   = float(get('volume_p'))
        self.b0       = float(get('b0'))
        self.rbtor0   = float(get('rbtor0'))
        self.betatotal = float(get('betatotal'))
        self.betapol  = float(get('betapol'))
        self.betator  = float(get('betator'))
        self.betaxis  = float(get('betaxis'))
        self.ier_flag = int(get('ier_flag'))

        try:
            self.machsq = float(get('machsq'))
        except KeyError:
            self.machsq = 0.0

        # --- radial grid ---
        phi_edge = get('phi')[-1]
        self.phi  = get('phi')
        self.chi  = get('chi')
        # normalised flux label s = sqrt(phi / phi_edge)
        self.s = np.sqrt(np.abs(self.phi / phi_edge))

        # --- 1-D profiles ---
        self.iotaf  = get('iotaf')
        self.q      = get('q_factor')
        self.presf  = get('presf')
        self.jcuru  = get('jcuru')
        self.jcurv  = get('jcurv')

        try:
            self.omega = get('omega')
        except KeyError:
            self.omega = np.zeros(self.ns)

        # --- mode arrays ---
        self.xm     = get('xm')
        self.xn     = get('xn')
        self.xm_nyq = get('xm_nyq')
        self.xn_nyq = get('xn_nyq')

        # --- 2-D Fourier coefficients ---
        self.rmnc = get('rmnc')
        self.zmns = get('zmns')
        self.lmns = get('lmns')

        # --- 2-D Nyquist Fourier coefficients ---
        self.bmnc     = get('bmnc')
        self.gmnc     = get('gmnc')
        self.bsubumnc = get('bsubumnc')
        self.bsubvmnc = get('bsubvmnc')
        self.bsubsmns = get('bsubsmns')
        self.bsupumnc = get('bsupumnc')
        self.bsupvmnc = get('bsupvmnc')

        f.close()

    def summary(self):
        converged = 'YES' if self.ier_flag == 0 else f'NO (ier_flag={self.ier_flag})'
        lines = [
            f"VMEC wout: {self.path.name}",
            f"  Converged : {converged}",
            f"  ns={self.ns}  nfp={self.nfp}  mpol={self.mpol}  ntor={self.ntor}",
            f"  Aspect    : {self.aspect:.4f}",
            f"  Aminor    : {self.Aminor:.4f} m",
            f"  Rmajor    : {self.Rmajor:.4f} m",
            f"  Volume    : {self.volume:.4f} m³",
            f"  B0        : {self.b0:.4f} T",
            f"  beta_total: {self.betatotal:.4e}",
            f"  Mach² (M02): {self.machsq:.4e}",
            f"  q range   : [{self.q[1]:.3f}, {self.q[-1]:.3f}]",
            f"  iota range: [{self.iotaf[0]:.3f}, {self.iotaf[-1]:.3f}]",
        ]
        return '\n'.join(lines)

    def get_rz_coeffs(self):
        """Return R/Z Fourier coefficients as a structured dict.

        Returns
        -------
        dict with keys:
            s       : 1-D array (ns,)          normalised flux label sqrt(phi/phi_edge)
            m       : 1-D int array (mn_mode,)  poloidal mode numbers
            n       : 1-D int array (mn_mode,)  toroidal mode numbers
            rmnc    : 2-D array (ns, mn_mode)   R cosine coefficients [m]
            zmns    : 2-D array (ns, mn_mode)   Z sine coefficients   [m]
        """
        return {
            's':    self.s.copy(),
            'm':    self.xm.astype(int).copy(),
            'n':    self.xn.astype(int).copy(),
            'rmnc': self.rmnc.copy(),
            'zmns': self.zmns.copy(),
        }

    def save_rz_coeffs(self, outpath=None, fmt='csv'):
        """Save R/Z Fourier coefficients to a file.

        Parameters
        ----------
        outpath : str or Path, optional
            Output path. Defaults to <wout_stem>_rz_coeffs.<fmt> in the same directory.
        fmt : 'csv' or 'npz'
            'csv' writes a long-format table (s, m, n, rmnc, zmns).
            'npz' saves the arrays directly (load with np.load).

        Returns
        -------
        Path of the saved file.
        """
        if outpath is None:
            outpath = self.path.parent / f'{self.path.stem}_rz_coeffs.{fmt}'
        outpath = Path(outpath)

        coeffs = self.get_rz_coeffs()

        if fmt == 'npz':
            np.savez(str(outpath), **coeffs)
        elif fmt == 'csv':
            import csv
            with open(outpath, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['s', 'm', 'n', 'rmnc', 'zmns'])
                for i, s_val in enumerate(coeffs['s']):
                    for j, (m, n) in enumerate(zip(coeffs['m'], coeffs['n'])):
                        writer.writerow([s_val, int(m), int(n),
                                         coeffs['rmnc'][i, j], coeffs['zmns'][i, j]])
        else:
            raise ValueError(f"Unknown format '{fmt}', choose 'csv' or 'npz'")

        return outpath


# ---------------------------------------------------------------------------
# Convenience function: read all wout files for a run into a list
# ---------------------------------------------------------------------------

def read_run_wouts(run_name, run_dir=None):
    """Return a list of WoutReader objects for every wout file in a run."""
    paths = list_wouts(run_name, run_dir=run_dir)
    if not paths:
        raise FileNotFoundError(f"No wout files found for run '{run_name}'")
    return [WoutReader(p) for p in paths]


def extract_rz_coeffs(run_name_or_path, index=None, fmt='csv', outdir=None):
    """Extract and save R/Z Fourier coefficients for one or all wout files in a run.

    Parameters
    ----------
    run_name_or_path : str or Path
        Path to a wout .nc file, or a run name (uses the default output root).
    index : int or None
        Specific wout index to process. If None and a run name is given, all
        wout files in the run are processed.
    fmt : 'csv' or 'npz'
    outdir : str or Path, optional
        Directory for output files. Defaults to each wout file's own directory.

    Returns
    -------
    list of Path objects pointing to the saved files.
    """
    path = Path(run_name_or_path)

    if path.suffix == '.nc':
        readers = [WoutReader(path)]
    elif index is not None:
        readers = [WoutReader(find_wout(run_name_or_path, index=index))]
    else:
        readers = read_run_wouts(run_name_or_path)

    saved = []
    for w in readers:
        if outdir is not None:
            out = Path(outdir) / f'{w.path.stem}_rz_coeffs.{fmt}'
        else:
            out = None
        p = w.save_rz_coeffs(outpath=out, fmt=fmt)
        print(f"  Saved {p.name}  ({w.ns} surfaces, {w.mnmax} modes)")
        saved.append(p)

    return saved


# ---------------------------------------------------------------------------
# Script entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Read VMEC wout NetCDF files — print summary or extract R/Z coefficients.')
    parser.add_argument('target',
                        help='Path to wout .nc file, or run name (uses default output root)')
    parser.add_argument('--index', type=int, default=None,
                        help='Wout file index (default: 0 for summary; all indices for --coeffs)')
    parser.add_argument('--coeffs', action='store_true',
                        help='Extract R/Z Fourier coefficients to file')
    parser.add_argument('--fmt', choices=['csv', 'npz'], default='csv',
                        help='Output format when using --coeffs (default: csv)')
    parser.add_argument('--outdir', default=None,
                        help='Output directory for extracted coefficient files')
    args = parser.parse_args()

    if args.coeffs:
        print(f"Extracting R/Z coefficients ({args.fmt}) ...")
        extract_rz_coeffs(args.target, index=args.index, fmt=args.fmt, outdir=args.outdir)
    else:
        path = Path(args.target)
        if path.suffix == '.nc' and path.exists():
            w = WoutReader(path)
            print(w.summary())
        else:
            idx = args.index if args.index is not None else 0
            wout_path = find_wout(args.target, index=idx)
            if not wout_path.exists():
                print(f"ERROR: {wout_path} not found")
                sys.exit(1)
            w = WoutReader(wout_path)
            print(w.summary())
