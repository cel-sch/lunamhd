# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 15:17:03 2024

@author: celin
"""
import argparse
from lunamhd import lunaScan

parser = argparse.ArgumentParser(description='Read command line inputs to run AEscanner.')
parser.add_argument('inputfile', type=str, help='input file containing run info e.g. default.in')
parser.add_argument('runid', type=str, help='run name, used for final output file')

scan_args = parser.parse_args()

runid = scan_args.runid
run = lunaScan(runid = scan_args.runid, inputfile = scan_args.inputfile)
run.init_run()
print(f"Run {runid} initialised.")
