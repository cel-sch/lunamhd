#README in progress

Running a 1D scan:
=================
### INPUTS ###
1) Specify which parameter is being scanned over in input file under 'scanparams'.
2) Specify scan array in the format min, max, nsteps.
3) Specify any other parameters.

### INITIALIZING A RUN ###
To initialize a run by creating all the necessary input files and subdirectories for each point in the run:
  
  run = lunaScan(runid = 'runname', inputfile = 'default.in', inputpath = 'path_to_input_file') # can specify inputpath, otherwise defaults to default input path
  run.init_run()

### RUNNING SCANS ###
If running the scan immediately in the same terminal:
  
  run.run(scan_saveloc = inputpath)

Else:

  scan = lunaScan(runid = scanid, inputfile = 'default.in', inputpath = 'path_to_input_file') # reads in the correct input file
  scan.run(scan_saveloc = 'path_to_input_file')

### SAVING A RUN ###
If saving the scan immediately in the same terminal:
  run.save_run()

Else:

  run = lunaScan(runid = runid, inputfile = 'path_to_input_file')
  run.scans = run._make_scan_list()
  run.save_run()

### NOTES ###
- When running an ND scan, the lowest level (i.e. first) scan parameter needs to go first in 'scanparams'. The order of the rest doesn't matter.

---

VMEC output tools
=================

Two modules for reading and plotting VMEC equilibrium outputs produced during a lunamhd run. VMEC wout files are NetCDF files located at:

```
Output/KH/<run_name>/VMEC/wout/wout_<run_name>_<index>.nc
```

Each index corresponds to a different Mach number in the scan.

**Python interpreter:** use `/Users/cellywelly/Dev/solis/bin/python3` (has scipy + numpy; system Python 3.9 does not).

---

### vmecReader.py

Reads a wout NetCDF file into a `WoutReader` object.

```python
from vmecReader import WoutReader, read_run_wouts, find_wout, list_wouts

# Single file
w = WoutReader('Output/KH/test_oldvmec/VMEC/wout/wout_test_oldvmec_0.nc')
print(w.summary())   # scalar overview (aspect, B0, beta, Mach², q range…)

# All wout files for a run, sorted by index
wouts = read_run_wouts('QI_bstep0.7_omstep0.5')

# Find / list paths
path  = find_wout('test_oldvmec', index=1)
paths = list_wouts('test_oldvmec')
```

**Key attributes on `WoutReader`:**

| Attribute | Shape | Description |
|---|---|---|
| `s` | `(ns,)` | Normalised flux label √(Φ/Φ_edge) |
| `q`, `iotaf` | `(ns,)` | Safety factor and rotational transform |
| `presf` | `(ns,)` | Pressure profile |
| `omega` | `(ns,)` | Rotation profile (0 for non-flow runs) |
| `jcuru`, `jcurv` | `(ns,)` | Poloidal / toroidal current density |
| `rmnc`, `zmns` | `(ns, mn_mode)` | Fourier coefficients of R and Z |
| `bmnc` | `(ns, mn_mode_nyq)` | Fourier coefficients of \|B\| |
| `machsq` | scalar | Mach number squared M² on axis |
| `b0`, `betatotal` | scalar | On-axis B field and total beta |

**As a script** (prints a scalar summary):
```bash
python3 vmecReader.py test_oldvmec           # run name, index 0
python3 vmecReader.py test_oldvmec 2         # run name, index 2
python3 vmecReader.py path/to/wout_foo_0.nc  # direct file path
```

---

### vmecPlotter.py

Plots VMEC equilibrium outputs. All current runs are axisymmetric (`ntor=0`), so the Fourier reconstruction is a 1D sum over poloidal modes and the result is the same at every toroidal angle. Field lines lie on the flux surfaces, so the R-Z cross-section is the relevant geometry.

#### Single-file plots (`VmecPlotter`)

```python
from vmecPlotter import VmecPlotter

p = VmecPlotter('Output/KH/test_oldvmec/VMEC/wout/wout_test_oldvmec_0.nc')

p.flux_surfaces()          # R-Z cross-section, surfaces coloured by ⟨|B|⟩
p.b_profile()              # |B|(θ) on 5 flux surfaces
p.q_profile()              # q(s) and ι(s) on twin axes
p.summary_plot()           # all three panels in one figure
```

#### Interactive Mach-number sliders

All slider functions sort wout files by ascending M, fix axis limits across all files, and update artist data in place (no axes clear) so they respond correctly in IPython.

```python
from vmecPlotter import flux_surfaces_slider, summary_slider, profiles_slider

# R-Z flux surfaces only
flux_surfaces_slider('QI_bstep0.7_omstep0.5')

# Flux surfaces + |B|(θ) + q/ι, all updating together
summary_slider('QI_bstep0.7_omstep0.5')

# Ω(s), P(s), q(s), ι(s) profiles — 2×2 layout
profiles_slider('QI_bstep0.7_omstep0.5')
```

> **Note:** The Shafranov shift is only clearly visible for runs reaching M > 1 (e.g. `QI_bstep0.7_omstep0.5`, M up to 2.6). For `QI_test_bstep0.7` (M up to 1.1) the axis shift is < 0.5 mm and not visible at the plot scale. The pressure profile does change visibly in the high-M case (centrifugal compression).

**As a script** — passing a run name with multiple wout files defaults to `summary_slider`:
```bash
python3 vmecPlotter.py QI_bstep0.7_omstep0.5                  # → summary_slider
python3 vmecPlotter.py QI_bstep0.7_omstep0.5 profiles_slider
python3 vmecPlotter.py QI_bstep0.7_omstep0.5 slider           # flux surfaces only
python3 vmecPlotter.py test_oldvmec 0 flux_surfaces            # single file
```
