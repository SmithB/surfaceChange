import pointCollection as pc
import matplotlib.pyplot as plt
import numpy as np
import sys
import json
import os


pad=np.array([-1.e3, 1.e3])
tols=[0, 25, 50, 100]
tols.sort()

tile_file=sys.argv[1]
DEM_file=sys.argv[2]
dest_dir=sys.argv[3]

if not os.path.isdir(dest_dir):
	os.mkdir(dest_dir)
for tol in tols:
	tol_dir=os.path.join(dest_dir, str(tol))
	if not os.path.isdir(tol_dir):
		os.mkdir(tol_dir)

D=pc.data().from_h5(tile_file, group='data')
D=D[D.three_sigma_edit==1]
bounds=[np.c_[np.min(D.x), np.max(D.x)]+pad, np.c_[np.min(D.y), np.max(D.y)]+pad]
D.assign({'DEM':pc.grid.data().from_geotif(DEM_file, bounds=bounds).interp(D.x, D.y)})
r=np.abs(D.z-D.DEM)
report_file=os.path.basename(tile_file).replace('.h5','.DEM_deltas.json')
for tol in tols:
	if np.any(r>tol):
		this_tol=tol
		break
report_file=os.path.join(dest_dir, str(this_tol), report_file)
report={'tile_file':tile_file, 
				'DEM_file':DEM_file, 
				'tols':tols,
				'N':[np.abs(r>tol) for tol in tols],
				'sigma':np.nanstd(D.z-D.DEM),
				'std_r_data':np.nanstd(D.z-D.z_est)}
with open(report_file, 'w', encoding='utf-8') as f:
    json.dump(report, f, ensure_ascii=False, indent=4)


