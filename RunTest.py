from lunamhd import lunaScan, lunaRead

runid = 'test_oldvmec'
run = lunaScan(runid = runid, inputfile = 'default_KH.in')
run.init_run()
run.run(scan_saveloc = f'/Users/cellywelly/Dev/lunamhd/Output/KH/{runid}')
run.save_run()

read = lunaRead(runid, filePath = f'/Users/cellywelly/Dev/lunamhd/Output/KH/{runid}')
# bp = read.basic_plot()
# bp.open_plot()
# bp.save_plot('/users/cs2427/scratch/lunamhd/KH/test/')