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

Release notes
--------

Squeakr now has a new k-mer representation (version 2) based on the new version of the CQF and some few Squeakr specific changes. The new version of Squeakr is not compatible with the old version. We have added some new features to the 'squeakr count' command and a couple of new commands.


* Squeakr count command now supports auto-resizing. However, auto-resizing only works when the count command is run with a single thread.
* Squeakr count command can now filter out k-mers in the final representation below a certain count value.
* Squeakr count command can now exclude counts in the final representation and only keep k-mers.

* Squeakr list: to list k-mers present in a Squeakr representation. This command only works when the representation is exact.

* Squeakr info: to get the infomation about the Squeakr representation. For example, version, k-mer size, number of k-mers, CQF specific info, etc.

API
--------
* 'squeakr count': count k-mers in a read dataset.
* 'squeakr query': query k-mers in the Squeakr representation.
* 'squeakr inner-prod': compute inner products of two Squeakr representations.
* 'squeakr list': list k-mers in the Squeakr representation. Only in exact
  representation.
* 'squeakr info': get information about the Squeakr representation.

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
 $ make squeakr
 $ ./squeakr count -f -k 28 -s 20 -t 1 -o data/tmp.squeakr data/test.fastq
```

The usage of `./squeakr count` is as follows:

```bash
SYNOPSIS
        squeakr count [-e] -k <k-size> [-c <cutoff>] [-n] [-s <log-slots>] [-t <num-threads>] -o <out-file> <files>...

OPTIONS
        -e, --exact squeakr-exact (default is Squeakr approximate)
        <k-size>    length of k-mers to count
        <cutoff>    only output k-mers with count greater than or equal to cutoff (default = 1)

        -n, --no-counts
                    only output k-mers and no counts (default = false)

        <log-slots> log of number of slots in the CQF. (Size argument is only optional when numthreads is exactly 1.)

        <num-threads>
                    number of threads to use to count (default = number of hardware threads)

        <out-file>  file in which output should be written
        <files>...  list of files to be counted (supported files: fastq and compressed gzip or bzip2 fastq files)
```

squeakr-count creates a file <out-file> which is the k-mer representation.

`lognumslots.sh` script can be used to estimate the `log of number of slots in the CQF` argument. The script takes as input the path to the output file of 'ntCard' (https://github.com/bcgsc/ntCard). It then calculates log of the number of slots needed by Squeakr to count k-mers.

```bash
 $ ./squeakr query -f data/tmp.squeakr -q data/query_file -o data/query.output
```
The usage of `./squeakr query` is as follows:

```bash
SYNOPSIS
        squeakr query -f <squeakr-file> -q <query-file> -o <output-file>
OPTIONS
        <squeakr-file>
                    input squeakr file

        <query-file>
                    input query file

        <output-file>
                    output file
```

```bash
 $ make squeakr-inner-prod
 $ ./squeakr inner_prod data/tmp.squeakr data/tmp.squeakr
```
 The usage of `./squeakr inner_prod` is as follows:

```bash
SYNOPSIS
        squeakr inner_prod <first-input> <second-input>

OPTIONS
        <first-input>
                    first input squeakr file

        <second-input>
                    second input squeakr file

```

Contributing
------------
Contributions via GitHub pull requests are welcome.


Authors
-------
- Prashant Pandey <ppandey@cs.stonybrook.edu>
- Rob Patro <rob.patro@cs.stonybrook.edu>
- Rob Johnson <rob@cs.stonybrook.edu>
