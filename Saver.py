import argparse
from lunamhd import *

parser = argparse.ArgumentParser(description='Save run file by compiling scan files.')
parser.add_argument('inputfile', type=str, help='input file containing run info e.g. default.in')
parser.add_argument('runid', type=str, help='run name, used for final output file')

scan_args = parser.parse_args()

runid = scan_args.runid
run = lunaScan(runid = runid, inputfile = scan_args.inputfile)
run.scans = run._make_scan_list()
run.save_run()
print(f"Run {runid} saved.")

# Test plot
read = lunaRead(runid)
bp = read.basic_plot()
bp.open_plot()
bp.save_plot(f'/users/cs2427/scratch/lunamhd-data/KH/{runid}/{runid}_EVs.png')