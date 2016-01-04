## Parsing T- and B-cell receptor segments from IMGT database to a flexible plain-text format

### Getting raw sequences

Instructions for downloading raw IMGT files:

* Go to the [genedb page](http://www.imgt.org/genedb/) and click the first submit button
* Next, scroll down to the end of resulting page (loading can take a while) and mark ""*Select all genes*"
* Select ""*F+ORF+in-frame P nucleotide sequences with IMGT gaps*" format and click submit
* Copy-paste resulting FASTA records and use them as an input file

### Running the software

Get the compiled [binaries](https://github.com/antigenomics/imgtparser/releases) and run the software as
``java -jar imgtparser.jar [options] imgt_raw_file output_prefix``.

The following options can be selected:

* ``-n`` include non-functional segments into output (pseudogenes, etc)
* ``-m`` include minor alleles (segments with ``*02``, ``*03``, etc suffix)
* ``-s`` toggle species detalisation (e.g. BALB/c and C57Bl6 for MusMusculus)
* ``-b`` report IMGT records that cannot be parsed properly (missing conserved residues, etc)

ImgtParser generates a tab-delimited table with species name, gene and segment id, nucleotide sequence and
the reference point position: 0-based coordinate of first nucleotide after conserved Cys for Variable segments and
before first nucleotide before conserved Phe/Trp for Joining segments. The metadata table provided with results
lists all species and genes and tells if there are any V/D/J segments associated with them (1 in corresponding row).