# squeakr
Squeakr: An Exact and Approximate k -mer Counting System

Overview
--------

Squeakr is a k-mer-counting and multiset-representation system using the
recently-introduced counting quotient filter (CQF) Pandey et al. (2017), a
feature-rich approximate membership query (AMQ) data structure.

Squeakr is memory-efficient, consuming 1.5X–4.3X less memory than the
state-of-the-art. It offers competitive counting performance, and answers
queries about a particular k-mer over an order-of- magnitude faster than other
systems. The Squeakr representation of the k-mer multiset turns out to be
immediately useful for downstream processing (e.g., De Bruijn graph traversal)
because it supports fast queries and dynamic k-mer insertion, deletion, and
modification.

API
--------
* 'main': count k-mers in a read dataset
* 'query': query k-mers in the Squeakr representation
* 'inner-prod': compute inner products of two Squeakr representations

Build
-------

Library dependencies (given version or higher):
 - libboost-dev 1.58.0.1ubuntu1
 - libssl-dev 1.0.2g-1ubuntu4.6
 - zlib1g-dev 1:1.2.8.dfsg-2ubuntu4
 - bzip2 1.0.6-8

Squeakr currently only supports fastq files. If any other file formats are passed as input then it will 
throw a segmentation fault.

The CQF code uses two new instructions to implement select on machine words introduced in intel's Haswell line of CPUs. However, there is also an alternate implementation of select on machine words to work on CPUs older than Haswell.
To build on an older hardware (older than Haswell) use "NH=1" as a make argument.

```bash
 $ make main
 $ ./main 0 20 1 test.fastq
```

 Following are the arguments to main:
 - file format: 0 - plain fastq, 1 - gzip compressed fastq, 2 - bzip2 compressed fastq
 - CQF size: the log of the number of slots in the CQF
 - num of threads: number of threads to count
 - file(s): "filename" or "dirname/*" for all the files in a directory

The main program creates a files with the extension ".ser" which is the k-mer representation.

```bash
 $ make kmer_query
 $ ./kmer_query test.fastq.ser 10000 0
```

 Following are the arguments to query:
 - file: dataset Squeakr representation
 - num of queries: number of queries
 - random: 0 - query for existing k-mers, 1 - query for random k-mers

```bash
 $ make kmer_inner_prod
 $ ./kmer_inner_prod test.fastq.ser test.fastq.ser
```
 
 Following are the arguments to inner-prod:
 - file1: dataset 1 Squeakr representation
 - file2: dataset 2 Squeakr representation

Contributing
------------
Contributions via GitHub pull requests are welcome.


Authors
-------
- Prashant Pandey <ppandey@cs.stonybrook.edu>
- Rob Patro <rob.patro@cs.stonybrook.edu>
- Rob Johnson <rob@cs.stonybrook.edu>
