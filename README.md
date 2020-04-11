# Ariadne: Barcoded Read Deconvolution

Linked-reads, however, come with their own challenges. Because multiple long fragments can be associated with the same 3' barcode, linked-read data is algorithmically difficult to deconvolve. The lack of a unique correspondence between a long fragment and a barcode, in conjunction with low sequencing depth, confounds the assignment of linkage between short-reads. 

Here, we introduce Ariadne, a novel assembly graph-based algorithm, that can be used to deconvolve a large metagenomic linked-read dataset. As demonstrated substantial increases in the largest alignments, contigs, and genomic contiguity, *de novo* assemblies from deconvolved reads represent advancements in terms of completeness and accuracy. Ariadne is intuitive, computationally efficient, and scalable to other large-scale linked-read problems, such as human genome phasing. 

The Ariadne manuscript will be publicly available soon. 

## Installation

Currently, Ariadne is implemented as a module of the SPAdes *de novo* assembly program (version 3.13.1). The following libraries to be pre-installed:

* g++ (version 5.3.1 or higher)
* cmake (version 2.8.12 or higher)
* zlib
* libbz2

From source: 
```
git clone <url>   
cd Ariadne
./spades_compile.sh
```
The installation directory can be set by `PREFIX=<destination_dir>` in the build step. 

In the future, Ariadne will be repackaged as a standalone program.

## Deconvolving Reads

Use the following command to run barcode deconvolution. `<fastq_name>` should be separate or interleaved fastq file where reads have a `BX` tag designating the barcode (this is the default output of [longranger basic](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/other-pipelines)). `<max_search_dist>` is the user-specified parameter for the maximum search distance, which should reflect the average length of a genomic fragment. Future versions of Ariadne will allow for the user to specify the barcode identifier. 
```
spades.py --only-assembler -1 <fastq_name>.read1.fq.gz -2 <fastq_name>.read2.fq.gz --barcode-distance <max_search_dist> -o /path/to/output_dir
```

For more SPAdes options, refer to the [Spades manual](http://cab.spbu.ru/files/release3.13.1/manual.html) or the command-line options.
```
spades.py
```

### Output

At the end of the deconvolution procedure, Ariadne outputs an interleaved fastq file with enhanced barcode assignments in the directory `/path/to/output/k55/<max_search_dist>enhanced.fastq`. The original barcode has been augmented with the enhanced grouping number. In this example, the first set of paired reads `@D00547:847:HYHNTBCXX:1:1101:10000:10626` have been assigned to the 13th group out of all reads with the barcode `BX:Z:CCTTCCCTCCTTCAAT`. This fastq can be directly provided as input to a *de novo* assembler or a read aligner.

```
$ head /path/to/output/k55/<max_search_dist>enhanced.fastq
@D00547:847:HYHNTBCXX:1:1101:10000:10626 BX:Z:CCTTCCCTCCTTCAAT-13
ATGCTGGGGTTTCCGCTGCAATTCTTTGTCCGGTTCTTTAAGAACCACGGCTTGCTGTCGATCAGCAACCGCCCACAGTGGTGCGTGATCGAAGGCGGCTCCAGCAGCTACATCGAGCCGCTGACCC
+
IIGGIIIGIIIIIGIIIGGGGGGIIIGIIIIIIGGGGIGGGIIIIIGIGGGGGGGIIGGGGGGGIIIIIIIIIIGGGIGIIIIGGIIGIIIGIIIAGGG.GGIIGGGGGGGGGGIIIGGIIIIGGGI
@D00547:847:HYHNTBCXX:1:1101:10000:10626 BX:Z:CCTTCCCTCCTTCAAT-13
TCATGTCGTAGGTAACGGCGGCCTGTGTCTGCGCATCGCCGCTCAGCCGATAATTCCAGCTGGCCCAGGCCAGTTTGCGGTCCGGCAGCAGGCGTGTGTCGGTGTGCAGCACCACGTCATTGTCGGCATAGGGCAATGCGCCGAGGATCT
+
GGAGGGGGIGAGGG.GGGGGG<AGIGGGGGGIIGGGIGIGGGIIGGGGGGGGAGGGGGGGG.A.<GG<AGGAAGAAGGIIGA<GA<G<GGGGIGGGIIIII<AGGIIIIGAGGIIIGIIGGGGGIAGGGGGGIIIGGIGGGAGAA.<GGA
@D00547:847:HYHNTBCXX:1:1101:10000:11336 BX:Z:CAGGTATCAGCGTAAG-1
TTCAGCAAGCGCAGCTTGATGGCATCCCGCAATTGATACAGGTTCAGGCCGGTGGTGTTGATCTTGAGGTCGGCCAGATCGATGATCGGTCCCAGCAGCGAGGTCTCGTCCTCGATGGCTTCGGCCA
+
A<.AGGGGGGGAGAAGGGIGGAAGGGGGGGGGGGGIGGGGGGIGGGIG.GGGAGGGGGGGG.GGGGGAGGGGG.AGA.AGG<GG<AAGGGGGGGGGAGGG.GGGIGGAG<<<AGAAGAAG<GGGGG.
@D00547:847:HYHNTBCXX:1:1101:10000:11336 BX:Z:CAGGTATCAGCGTAAG-1
GTCCTGGATGAACAAAACGGGCAGTAATCATGCGTTTGATCATCGTCAGCGGCCGCTCCGGCTCGGGTAAAAGCACCGCCCTCAACGTCCTTGAAGACAACGGCTTTTATTGCATCGACAACCTTCCCGCCGGTTTGCTGCCGGAGTTGG
+
G<..GGGGGGGGA<GGA.<GGGGGAGGG.GGAGAGAGG.AGG<<<..<...AAA<AA<G....<GAGGIGIGAGAG.....<AGG..<AAAA<..<...<.<...<<AGG...GGG<GAAAAGG..<A<AGGAG.<<.G...<7<<7GGG
@D00547:847:HYHNTBCXX:1:1101:10000:11336 BX:Z:CAGGTATCAGCGTAAG-1
GTCCTGGATGAACAAAACGGGCAGTAATCATGCGTTTGATCATCGTCAGCGGCCGCTCCGGCTCGGGTAAAAGCACCGCCCTCAACGTCCTTGAAGACAACGGCTTTTATTGCATCGACAACCTTCCCGCCGGTTTGCTGCCGGAGTTGG
+
G<..GGGGGGGGA<GGA.<GGGGGAGGG.GGAGAGAGG.AGG<<<..<...AAA<AA<G....<GAGGIGIGAGAG.....<AGG..<AAAA<..<...<.<...<<AGG...GGG<GAAAAGG..<A<AGGAG.<<.G...<7<<7GGG
```

### Performance

A full exploration of performance metrics will be available with the manuscript. For example, for metagenomics-sized dataset such as the full MOCK5 dataset (97 million reads approx., see below), deconvolution at a maximum search distance of 20,000 bp takes 7 hours and 33 minutes, 339 GB RAM, and 30 CPUs. 

## Datasets

The MOCK5 dataset used in the paper may be downloaded from [AWS](https://s3.us-east-2.amazonaws.com/readclouds/cloudspades_data.tar.gz).

## Credits

This algorithm was developed by Waris Barakzai and myself, and tested by myself with help from Dmitrii Meleshko, David Danko, and Iman Hajirasouliha.
