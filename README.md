## Parsing T- and B-cell receptor segments from IMGT database to a flexible plain-text format

### Getting raw sequences

Instructions for downloading raw IMGT files:

* Go to the [genedb page](http://www.imgt.org/genedb/) and click the first submit button
* Next, scroll down to the end of resulting page (loading can take a while) and mark ""*Select all genes*"
* Select ""*F+ORF+in-frame P nucleotide sequences with IMGT gaps*" format and click submit
* Copy-paste resulting FASTA records and use them as an input file

### Running the software

Get the [binaries]() and run as ``java -jar imgtparser.jar [options] imgt_raw_file output_prefix``.

The following options can be selected:

* ``-n`` include non-functional segments into output (pseudogenes, etc)
* ``-m`` include minor alleles (segments with ``*02``, ``*03``, etc suffix)
* ``-b`` report IMGT records that cannot be parsed properly (missing conserved residues, etc)