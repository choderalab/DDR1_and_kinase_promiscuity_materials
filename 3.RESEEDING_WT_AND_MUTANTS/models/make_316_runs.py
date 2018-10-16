# BE VERY CAREFUL WHEN USING THIS SCRIPT IT RENAMES THINGS THAT ARE POTENTIALLY IMPORTANT

# this script makes target.txt files that correspond to the directory each RUN is in
# this script also renumbers runs from being 8 times 0-~39 to be one consecutive list of 0 to 316

from glob import glob
import re
import os

run_paths = glob('./packaged_models/fah-projects/*/*')

run_number = 0

for run_path in run_paths:

    # Define our target by the directory name
    target = run_path.split('/')[3]

    if run_path.split('.')[-1] == 'xml':
        #there are some xml files we don't want to include that come in the glob selection
        pass

    else:

       print run_path
       print target

       # write the target.txt file in the RUN* directory
       with open("%s/target.txt"%run_path, "w") as text_file:
           text_file.write("{0}".format(target))

       # Define our new run number

       print run_number

       # Rename our directory to be labeled by that new run number

       base_path = run_path.rsplit('/',1)[0]
 
       new_path = '%s/RUN%s'%(base_path,run_number)

       os.rename(run_path,new_path)

       # Iterate the run_number before it moves onto the next loop
       run_number += 1
