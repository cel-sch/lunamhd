from lunamhd import lunaScan, lunaRead

runid = 'test_oldvmec'
run = lunaScan(runid = runid, inputfile = 'default.in')
run.init_run()
run.run(scan_saveloc = f'/home/csch/VENUS-linux/lunamhd/Output/KH/{runid}')
run.save_run()

read = lunaRead(runid, filePath = f'/home/csch/VENUS-linux/lunamhd/Output/KH/{runid}')
# bp = read.basic_plot()
# bp.open_plot()
# bp.save_plot('/users/cs2427/scratch/lunamhd/KH/test/')