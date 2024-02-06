import argparse
from AEmhd import *

parser = argparse.ArgumentParser(description='Save run file by compiling scan files.')
parser.add_argument('inputfile', type=str, help='input file containing run info e.g. default.in')
parser.add_argument('runid', type=str, help='run name, used for final output file')

scan_args = parser.parse_args()

runid = scan_args.runid
run = AEscan(runid = runid, inputfile = scan_args.inputfile)
run.scans = run._make_scan_list()
run.save_run()
print(f"Run {runid} saved.")
