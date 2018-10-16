# Allosteric network results and figures

Results and figures used to generate Figures 4A, S5, and S6 relating to the Two Sample Logo Analysis. A description of the method can be found in [Vacic et al, 2006](https://doi.org/10.1093/bioinformatics/btl151), and the resource described can be found [here](http://www.twosamplelogo.org/cgi-bin/tsl/tsl.cgi).


## Manifest
* As inputs to this analysis we extracted the promiscuous kinases from the human kinase alignment presented in [Manning et al, 2002](https://doi.org/10.1126/science.1075762) and available at [kinase.com](http://kinase.com/human/kinome/groups/ePK.aln).
 * `eight-promiscuous-kinases.fasta`
 * `kinase-background.fasta`
* This resulted in two output files:
 * `results.txt`
 * `tsl-cgi-web-result.pdf`
* These results are summarized in Figures 4A, S5, and S6.
 * `TSL_network_w_Abl_res_muts.pse` - PyMol file used to generate `TSL-network.png` (Fig. 4A) and `TSL-ABL-resmut-overlap.png` (Fig. S6).
 * `all_prom_kin_allosteric_final.jvp` - JalView file used to generate sequence alignment presented in Figure S5.

