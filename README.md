# unitig-counter
[![Anaconda-Server Badge](https://anaconda.org/bioconda/unitig-counter/badges/version.svg)](https://anaconda.org/bioconda/unitig-counter)

Uses a compressed de Bruijn graph (implemented in [GATB](https://github.com/GATB/gatb-core)) to count unitigs in bacterial populations.

## Details
This is a slightly modified version of the unitig and graph steps in [DBGWAS](https://gitlab.com/leoisl/dbgwas/) software, repurposed for input into [pyseer](https://pyseer.readthedocs.io/en/master/).

**NB** We cannot offer support for unitig-counter, it is provided 'as-is'. Please consider using [unitig-caller](https://github.com/bacpop/unitig-caller) instead, which offers the same functionality.

## Citation
If you use this, please cite the DBGWAS paper:

Jaillard M., Lima L. et al. A fast and agnostic method for bacterial genome-wide association studies: Bridging the gap between k-mers and genetic events. *PLOS Genetics*. **14**, e1007758 (2018). doi:[10.1371/journal.pgen.1007758](https://doi.org/10.1371/journal.pgen.1007758).

### List of changes

* Changes the format of the output from `step1` from bugwas matrix to pyseer input (Rtab or kmers).
* Removes all code for `step2` and `step3` in DBGWAS.
* Remove unused depencencies.
* Change installation procedure ready for bioconda.

## Install

Recommended installation is through [conda](https://docs.conda.io/en/latest/miniconda.html):
```
conda install unitig-counter
```
If the package cannot be found, [ensure your channels are set up correctly for bioconda](http://bioconda.github.io/#set-up-channels).

For compilation from source, see `INSTALL.md`.

## Usage
Run:
```
unitig-counter -strains strain_list.txt -output output -nb-cores 4
```

Where `strain_list.txt` is a list of input files (assemblies) with a header, for example:
```
ID      Path
6925_1_49       assemblies/6925_1#49.contigs_velvet.fa
6925_1_50       assemblies/6925_1#50.contigs_velvet.fa
```

Output is in `output/unitigs.txt` and can be used with `--kmers` in pyseer. You can also test just the
unique patterns in `output/unitigs.unique_rows.txt` with the `--Rtab` option.

## Cleaning up output
Some unitigs in the output may span multiple input contigs. If you wish to restrict your unitig calls to those appearing in assembled contigs, you can either:

1. Run [unitig-caller](https://github.com/johnlees/unitig-caller) on the input genomes, using the unitig calls from your run.
2. Run the [script](https://github.com/GATB/bcalm/blob/master/scripts/split_unitigs.py) in the `gatb`/`bcalm` package, which will cut unitigs that span multiple contigs.

Thanks to @rchikhi and @apredeus for discovering and fixing this.

## Extracting distances
Two get the shortest sequence distance between two unitigs:
```
cdbg-ops dist --graph test_data/graph --source GTAATAAACAAA --target AAAAAAAAAAGTTAAAAAT
```

## Extending unitigs
Short unitigs can be extended by following paths in the graph to neightbouring nodes. This can help map
sequences which on their own are difficult to align in a specific manner.

Create a file `unitigs.txt` with the unitigs to extend (probably your significantly associated hits)
and run:
```
cdbg-ops extend --graph output/graph --unitigs unitigs.txt > extended.txt
```

The output `extended.txt` will contain possible extensions, comma separated, with lines corresponding to unitigs
in the input. See the help for more options.

### Python
A similar python script can be found in `unitig-graph`:
```
python unitig-graph/extend_hits.py --prefix output/graph --unitigs unitigs.txt > extended.txt
```
