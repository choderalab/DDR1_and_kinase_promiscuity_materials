mutant-models
=============

* Making 8 systems (targets) from 40 starting structures (templates). The 8 systems are WT, D671N, Y755A, and Y759A, each with the ASP of the DFG in a non-protonated or protonated state. 36 of the starting structures are from the script in `picking_new_starting`. The other four are from DDR1 crystal structures. The 4BKJ, 3ZOS, and DDR1_VX680 structures are as in `final-models`, but the DDR1_Dasatinib structure was regenerated using the `multi-models` script, since significant problems existed in the loops of the original set of simulations.
* Since all of the templates are manual the setup is a bit fancier than normal ensembler.
<pre>
ensembler init
python generate_templates-resolved.py
</pre>
* The `targets.fa` was copied from `final-models` and than adjusted manually for the other six targets.

* For the protonated DDR1 runs, need to explicitly define this in the `manual-overrides.yaml` file as, e.g.:

<pre>
refinement:
     ph: 8.0
     custom_residue_variants:
         DDR1_HUMAN_D0_PROTONATED:
             # keyed by 0-based residue index
             180: ASH
</pre>
* The rest just follows the [ensembler command line documentation](http://ensembler.readthedocs.org/en/latest/cli_docs.html) (run on cluster either using interactive session or [runscript](https://github.com/choderalab/dansu-dansu/tree/master/run_scripts)). Note: had to be using openmm version 6.3.0 for Folding@home compatibility: 

<pre>
ensembler align
ensembler build_models
ensembler cluster
ensembler refine_implicit
ensembler solvate
ensembler refine_explicit
ensembler package_models --package_for FAH
</pre>

* After packaging for FAH, the runs can be renumbered to go from 0 to 319 (or however many successful explicit-refined models were created) instead of from 0 to 39 for each target using:
<pre>
python make_316_runs.py
</pre>
This will also make a `target.txt` in each `RUN*` folder.

Additionally, in order to make sure each `RUN*` folder has a `system.xml` and an `integrator.xml`, these can be copied once inside the target folder, e.g. `packaged-models/fah-projects/DDR1_HUMAN_D0/` by:

<pre>
echo RUN* | xargs -n 1 cp integrator.xml
echo RUN* | xargs -n 1 cp system.xml
</pre>

The resulting simulations are available via the [Open Science Framework](https://osf.io/4r8x2/) in folder `WTnMUTS-11409`. Runs are defined as follows: 
<pre>
 * WT               =    RUN0-39
 * WT protonated    =    RUN40-79
 * D671N            =    RUN80-118
 * D671N protonated =    RUN119-158
 * Y755A            =    RUN159-198
 * Y755A protonated =    RUN199-238
 * Y759A            =    RUN239-276
 * Y759A protonated =    RUN277-316
</pre>
