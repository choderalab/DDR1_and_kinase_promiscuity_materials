final-models
=============

* A total of 36 models.
* Two sets of 18, one with D of DFG as ASP, the other as ASH (protonated ASP).
* Because we are extending beyond the kinase domain, a new ensembler run was done.
* The 16 templates in the pdb were gathered with:

<pre>
ensembler init
ensembler gather_templates --gather_from PDB --query 1IR3,1K3A,1Y57,2F4J,3DQW,4BKJ,3ZOS,2OIQ,1OPJ,2SRC,2G2I,2H8H,3U4W,4GT5,1LUF,1IRK
</pre>

* Note that manually deleted any B, C, or D chains. There is indeed a `--chainids` flag I probably should have used ([ref](https://github.com/choderalab/ensembler/blob/master/ensembler/cli_commands/gather_templates.py)).
* The other two templates were added manually to the `templates-resolved-seq.fa`, `templates-full-seq.fa`, and the `structures-resolved` directory.
* The are two identical target sequences for DDR1, since the only difference is the protonation state. They are named 
`DDR1_HUMAN_D0` and `DDR1_HUMAN_D0_PROTONATED`, and the sequence is

<pre>
DFPRSRLRFKEKLGEGQFGEVHLCEVDSPQDLVSLDFPLNVRKGHPLLVAVKILRPDATK
NARNDFLKEVKIMSRLKDPNIIRLLGVCVQDDPLCMITDYMENGDLNQFLSAHQLEDKAA
EGAPGDGQAAQGPTISYPMLLHVAAQIASGMRYLATLNFVHRDLATRNCLVGENFTIKIA
DFGMSRNLYAGDYYRVQGRAVLPIRWMAWECILMGKFTTASDVWAFGVTLWEVLMLCRAQ
PFGQLTDEQVIENAGEFFRDQGRQVYLSRPPACPQGLYELMLRCWSRESEQRPPFSQLHR
FLAEDALNTV
</pre>

* For the protonated DDR1, need to explicitly define this in the `manual-overrides.yaml` file as:

<pre>
refinement:
     ph: 8.0
     custom_residue_variants:
         DDR1_HUMAN_D0_PROTONATED:
             # keyed by 0-based residue index
             180: ASH
</pre>

* The rest is easy, just follow the [ensembler command line documentation](http://ensembler.readthedocs.org/en/latest/cli_docs.html) (run on cluster either using interactive session or [runscript](https://github.com/choderalab/dansu-dansu/tree/master/run_scripts)): 

<pre>
ensembler loopmodel
ensembler align
ensembler build_models
ensembler cluster
ensembler refine_implicit
ensembler solvate
ensembler refine_explicit
ensembler package_models --package_for FAH
</pre>

The resulting simulations are available via the [Open Science Framework](https://osf.io/4r8x2/) in folder `WT-11403`. All these simulations are for the protonated version of DDR1 and runs correspond to the following starting conformations: 
<pre>
 * run1 = `DDR1_VX680`
 * run2 = `DDR1_HUMAN_3ZOS_A`
 * run3 = `DDR1_HUMAN_4BKJ_A`
</pre>

