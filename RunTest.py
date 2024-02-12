from lunamhd import lunaScan, lunaRead

run = lunaScan(runid = 'test', inputfile = 'default.in')
run.init_run()
run.run(scan_saveloc = '/users/cs2427/scratch/lunamhd/KH/test')
run.save_run()

read = lunaRead('test', filePath = '/users/cs2427/scratch/lunamhd/KH/test')
bp = read.basic_plot()
bp.open_plot()
bp.save_plot('/users/cs2427/scratch/lunamhd/KH/test/')