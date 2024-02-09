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
