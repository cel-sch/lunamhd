import argparse
from lunamhd import lunaScan

parser = argparse.ArgumentParser(description='Read job ID to tell AEmhd scanner where to find input file')
parser.add_argument('input_path', type=str, help='path to input file')
parser.add_argument('inputfile', type=str, help='input file containing run info e.g. default.in')
parser.add_argument('scanid', type=str, help='scan name, used for final output file')

scan_args = parser.parse_args()

scanid = scan_args.scanid
scan = lunaScan(runid = scanid, inputfile = scan_args.inputfile, inputpath = scan_args.input_path) # reads in the correct input file
scan.run(scan_saveloc = f'{scan_args.input_path}')