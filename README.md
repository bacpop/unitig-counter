# unitig-counter
Uses a compressed de Bruijn graph (implemented in [GATB](https://github.com/GATB/gatb-core)) to count unitigs in bacterial populations.

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

## Usage
Run:
```
./unitig-counter -strains strain_list.txt -output output -nb-cores 4
```

Where `strain_list.txt` is a list of input files (assemblies) with a header, for example:
```
ID      Path
6925_1_49       assemblies/6925_1#49.contigs_velvet.fa
6925_1_50       assemblies/6925_1#50.contigs_velvet.fa
```

Output is in `output/unitigs.txt` and can be used with `--kmers` in pyseer. You can also test just the
unique patterns in `output/unitigs.unique_rows.txt` with the `--Rtab` option.

## Extending unitigs
Short unitigs can be extended by following paths in the graph to neightbouring nodes. This can help map
sequences which on their own are difficult to align in a specific manner.

Create a file `unitigs.txt` with the unitigs to extend (probably your significantly associated hits)
and run:
```
python unitig-graph/extend_hits.py --prefix output/graph --unitigs unitigs.txt > extended.txt
```

The output `extended.txt` will contain possible extensions, comma separated, with lines corresponding to unitigs
in the input. See the help for more options.