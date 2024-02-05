import argparse
from lunamhd import lunaScan

parser = argparse.ArgumentParser(description='Read job ID to tell AEmhd scanner where to find input file')
parser.add_argument('input_path', type=str, help='path to input file')
parser.add_argument('inputfile', type=str, help='input file containing run info e.g. default.in')
parser.add_argument('scanid', type=str, help='scan name, used for final output file')

scan_args = parser.parse_args()

scan = lunaScan(inputfile = scan_args.inputfile, inputpath = scan_args.input_path) # reads in the correct input file
# don't need scanlist correctly added to object but it does need to be added in the final save
scan.init_scan(scan_args.scanid) # unclear what happens with the scanid, maybe just take the job id as the scanid?
# is there a way to get scanid to be the hex code?
# need to check naming so that the name inputs are correct for save_run function
scan.run(saveloc = f'{scan_args.input_path}')