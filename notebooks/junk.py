import numpy as np

piece=11
junk=np.load('/home/zajohns/dm-helium-flash/data/1.3_solar_mass/flash_scan-{0}.npy'.format(piece),allow_pickle=True, encoding='latin1')

print(junk)
