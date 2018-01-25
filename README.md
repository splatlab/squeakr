# squeakr
Squeakr: An Exact and Approximate k -mer Counting System

This work was published in Bioinformatics. If you use this software please cite us:
```
@article{doi:10.1093/bioinformatics/btx636,
author = {Pandey, Prashant and Bender, Michael A and Johnson, Rob and Patro, Rob},
title = {Squeakr: An Exact and Approximate k-mer Counting System},
journal = {Bioinformatics},
volume = {},
number = {},
pages = {btx636},
year = {2017},
doi = {10.1093/bioinformatics/btx636},
URL = { + http://dx.doi.org/10.1093/bioinformatics/btx636},
eprint = {/oup/backfile/content_public/journal/bioinformatics/pap/10.1093_bioinformatics_btx636/1/btx636.pdf}
}
```

Overview
--------

Squeakr is a k-mer-counting and multiset-representation system using the
recently-introduced counting quotient filter (CQF) Pandey et al. (2017), a
feature-rich approximate membership query (AMQ) data structure.

Squeakr is memory-efficient, consuming 1.5X–4.3X less memory than the
state-of-the-art. It offers competitive counting performance, in fact, it is
faster for larger k-mers, and answers queries about a particular k-mer over an
order-of- magnitude faster than other systems. The Squeakr representation of the
k-mer multiset turns out to be immediately useful for downstream processing
(e.g., De Bruijn graph traversal) because it supports fast queries and dynamic
k-mer insertion, deletion, and modification.

k-mer counts can be validated by hooking into the C++ level query API. An
example query program is also available in "kmer_query.cc".

API
--------
* 'squeakr-count': count k-mers in a read dataset.
* 'squeakr-query': query k-mers in the Squeakr representation.
* 'squeakr-inner-prod': compute inner products of two Squeakr representations.

Build
-------

Library dependencies (given version or higher):
 - libboost-dev 1.58.0.1ubuntu1
 - libssl-dev 1.0.2g-1ubuntu4.6
 - zlib1g-dev 1:1.2.8.dfsg-2ubuntu4
 - bzip2 1.0.6-8

Squeakr currently only supports fastq files. If any other file formats are
passed as input then it will throw a segmentation fault.

The CQF code uses two new instructions to implement select on machine words
introduced in intel's Haswell line of CPUs. However, there is also an alternate
implementation of select on machine words to work on CPUs older than Haswell.
To build on an older hardware (older than Haswell) use "NH=1" as a make argument.

```bash
 $ make squeakr-count
 $ ./squeakr-count -f -k 28 -s 20 -t 1 -o ./ test.fastq
```

The usage of `squeakr-count` is as follows:

```bash
SYNOPSIS
        ./squeakr-count (-f|-g|-b) -k <k-size> -s <log-slots> -t <num-threads> [-o <out-dir>] <files>... [-h]

OPTIONS
        format of the input
            -f      plain fastq
            -g      gzip compressed fastq
            -b      bzip2 compressed fastq

        <k-size>    length of k-mers to count
        <log-slots> log of number of slots in the CQF

        <num-threads>
                    number of threads to use to count

        <out-dir>   directory where output should be written (default = "./")
        <files>...  list of files to be counted
        -h, --help  show help
```

 Following are the arguments to squeakr-count:
 - file format: 0 - plain fastq, 1 - gzip compressed fastq, 2 - bzip2 compressed fastq
 - k-mer size: the size of the k-mer
 - CQF size: the log of the number of slots in the CQF
 - num of threads: number of threads to count
 - prefix: path to the directory where the output should go
 - file(s): "filename" or "dirname/*" for all the files in a directory

squeakr-count creates a files with the extension ".ser" which is the k-mer representation.

`lognumslots.sh` script can be used to estimate the `log of number of slots in the CQF` argument. The script takes as input the path to the output file of 'ntCard' (https://github.com/bcgsc/ntCard). It then calculates log of the number of slots needed by Squeakr to count k-mers.

```bash
 $ make squeakr-query
 $ ./squeakr-query -f test.fastq.ser -k 28 -n 1000 -r 0
```
The usage of `squeakr-query` is as follows:

```bash
SYNOPSIS
        ./squeakr-query -f <cqf-file> -k <k-size> -n <num-query> -r <random-queries> [-h]

OPTIONS
        <cqf-file>  input CQF file
        <k-size>    length of k-mers to query. Must be same the as the size of counted k-mers
        <num-query> number of queries

        <random-queries>
                    random queries

        -h, --help  show help
```

 Following are the arguments to squeakr-query:
 - file: dataset Squeakr representation
 - k-mer size: the size of the k-mer
 - num of queries: number of queries
 - random: 0 - query for existing k-mers, 1 - query for random k-mers

```bash
 $ make squeakr-inner-prod
 $ ./squeakr-inner-prod -a test.fastq.ser -b test.fastq.ser
```
 The usage of `squeakr-inner-product` is as follows:

```bash
SYNOPSIS
        ./squeakr-inner-prod -a <cqf-file-first> -b <cqf-file-second> [-h]

OPTIONS
        <cqf-file-first>
                    first input CQF file

        <cqf-file-second>
                    second input CQF file

        -h, --help  show help
```
 
 Following are the arguments to squeakr-inner-prod:
 - test.fastq.ser: dataset 1 Squeakr representation
 - test.fastq.ser: dataset 2 Squeakr representation

Contributing
------------
Contributions via GitHub pull requests are welcome.


Authors
-------
- Prashant Pandey <ppandey@cs.stonybrook.edu>
- Rob Patro <rob.patro@cs.stonybrook.edu>
- Rob Johnson <rob@cs.stonybrook.edu>
