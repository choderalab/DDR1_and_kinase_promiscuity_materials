# Homology modeling by the automodel class
# With multiple templates.
#
from modeller import *
from modeller.automodel import * # Load the automodel class
 
log.verbose()
 
# Model:
env = environ()
# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']
 
a = automodel(env,
            alnfile = 'DDR1_dasa_4BKJ.ali' ,       # alignment filename
            knowns = ('DDR1_Dasatinib','4BKJ-loop'),             # codes of the templates
            sequence = 'DDR1')           # code of the target
           
a.starting_model= 1                      # index of the first model
a.ending_model = 10                      # index of the last model
                                        # (determines how many models to calculate)
a.make()                                 # do homology modeling
 

