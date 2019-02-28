# unitig-counter
Uses a coloured de Bruijn graph (implemented in [GATB](https://github.com/GATB/gatb-core)) to count unitigs in bacterial populations.

## Details
This is a slightly modified version of the unitig and graph steps in [DBGWAS](https://gitlab.com/leoisl/dbgwas/) software, repurposed for input into [pyseer](https://pyseer.readthedocs.io/en/master/).

## Citation
If you use this, please cite the DBGWAS paper:

Jaillard M., Lima L. et al. A fast and agnostic method for bacterial genome-wide association studies: Bridging the gap between k-mers and genetic events. *PLOS Genetics*. **14**, e1007758 (2018). doi:[10.1371/journal.pgen.1007758](https://doi.org/10.1371/journal.pgen.1007758).

### List of changes

* Changes the format of the output from `step1` from bugwas matrix to pyseer input (Rtab or kmers).
* Removes all code for `step2` and `step3` in DBGWAS.
* Remove unused depencencies.
* Change installation procedure ready for bioconda.