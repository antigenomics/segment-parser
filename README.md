## Parsing T- and B-cell receptor segments from IMGT database to a flexible plain-text format

### Getting raw sequences

Instructions for downloading raw IMGT files:

* Go to the [genedb page](http://www.imgt.org/genedb/) and click the first submit button
* Next, scroll down to the end of resulting page (loading can take a while) and mark ``Select all genes``
* Select ``F+ORF+in-frame P nucleotide sequences with IMGT gaps`` format and click submit
* Copy-paste resulting FASTA records and use them as an input file

### Running the software

Get the compiled [binaries](https://github.com/antigenomics/imgtparser/releases) and run the software as
``java -jar segmentparser.jar [options] imgt_raw_file output_prefix``.

The following options can be selected:

* ``-n`` include non-functional segments into output (pseudogenes, etc)
* ``-m`` include minor alleles (segments with ``*02``, ``*03``, etc suffix)
* ``-s`` toggle species detalisation (e.g. BALBc and C57Bl6 for MusMusculus)
* ``-b`` report IMGT records that cannot be parsed properly (missing conserved residues, etc)

Output files include:

* A ``$output_prefix$.metadata.txt`` file with summary statistics.
* Files with erroneous/bad records: ``$output_prefix$.nojrefpoint.txt``, ``$output_prefix$.novrefpoint.txt``, ``$output_prefix$.othersegm.txt``.
* Output file containing sequences, CDR3 reference points and CDR1,2,2.5 coordinates: ``$output_prefix$.txt``.

> SegmentParser generates a tab-delimited table with species name, gene and segment id, nucleotide sequence and
the reference point position: 0-based coordinate of first nucleotide after conserved Cys for Variable segments and
before first nucleotide before conserved Phe/Trp for Joining segments. The metadata table provided with results
lists all species and genes and tells if there are any V/D/J segments associated with them (``0`` or ``1`` in corresponding row).

* A file with CDR1,2,2.5 nucleotide and amino acid sequences: ``$output_prefix$.txt`` (only includes V segments).

> Note that **CDR2.5** is a putative MHC-binding region of TCR V segment, defined in a recent work of Paul Thomas lab (Dash et al. Nature 2017).