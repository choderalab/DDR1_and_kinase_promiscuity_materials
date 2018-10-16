# This script generates the templates-resolved-seq.fa for all the 40 template structures.
# Each of these structures has the same WT DDR1 kinase domain sequence, but a different 
# name, which corresponds to it's PDB.

import shutil

# First let's copy these PDB's into the appropriate directory from 'picking-new-starting', 'multi-modeels', or 'final-models'.

destination_dir = 'templates/structures-resolved/'

source = {}

   # multi-model for DDR1_Dasatinib

source['DDR1_Dasatinib_model'] = '../multi-models/DDR1_Dasatinib/DDR1.B99990010.pdb'

   # previous models for DDR1_VX680, 4BKJ, and 3ZOS

source['DDR1_VX680_model'] = '../final-models/models/DDR1_HUMAN_D0/DDR1_VX680/model.pdb'  
source['DDR1_3ZOS_model']  = '../final-models/models/DDR1_HUMAN_D0/DDR1_HUMAN_3ZOS_A/model.pdb'
source['DDR1_4BKJ_model']  = '../final-models/models/DDR1_HUMAN_D0/DDR1_HUMAN_4BKJ_A/model.pdb'

   # new starting structures

from glob import glob
new_structures = glob('../picking_new_starting/*.pdb')

for structure in new_structures:

    short_filename = structure.split('/')[-1]
    short_name = short_filename.split('.')[0]
    source[short_name] = structure

   # copying everything in our source dictionary
for key in sorted(source):
    print 'copying %s' %key
    destination_pdb = '%s/%s.pdb'%(destination_dir,key)
    shutil.copy2(source[key],destination_pdb)

# Next let's make our templates-resolved-seq.fa using these

import mdtraj as md

for key in sorted(source):

    # using mdtraj allows us to get the sequence of the PDB to make the .fa file
    t = md.load(source[key])
    topol = t.topology

    # 'a' is append, make sure this file does not exist before you run this script.
    with open('templates/templates-resolved-seq.fa','a') as text_file:
        text_file.write('>%s\n'%key)
        text_file.write('%s\n'%topol.to_fasta()[0])

