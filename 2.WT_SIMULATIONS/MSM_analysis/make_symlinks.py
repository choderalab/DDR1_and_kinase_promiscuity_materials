"""
Make trajectory symlinks in local trajectories directory.
Looks for FAH data at the location given by environment variable
FAH_DATA_PATH
"""

import glob
import mdtraj as md
import os

DDR1_MUT_PROJECT = 'WT-11403'

RUNS = [1,2,3]

DATA_PATH = '../OSF'

PATH = "%s/%s/" % (DATA_PATH, DDR1_MUT_PROJECT)

filenames = [filename for filename in glob.glob('%s/*.h5'%PATH)]

try:
    os.mkdir("./WT-sims/")
except:
    pass

for filename in filenames:
    for RUN in RUNS:
        if RUN is None or "run%d-" % RUN in filename:
            base_filename = os.path.split(filename)[1]
            out_filename = "./WT-sims/%s" % (base_filename)
            if not os.path.exists(out_filename):
                os.symlink(filename, out_filename)
                
DDR1_MUT_PROJECT = 'WTnMUTS-11409'

RUNS = []
for i in range(40,80):
    RUNS.append(i)

PATH = "%s/%s/" % (DATA_PATH, DDR1_MUT_PROJECT)

filenames = [filename for filename in glob.glob('%s/*.h5'%PATH)]

for filename in filenames:
    for RUN in RUNS:
        if "run%d" % RUN in filename:
            base_filename = os.path.split(filename)[1]
            out_filename = "./WT-sims/%s" % (base_filename)
            if not os.path.exists(out_filename):
                os.symlink(filename, out_filename)


