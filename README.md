# HQAlign

### Contents <a id='contents'></a>

* <a href='#intro'>Introduction</a>
* <a href='#pub'>Publication</a>
* <a href='#setup'>Setup</a>
* <a href='#use'>Usage</a>

---

### Introduction <a id='intro'></a>

Detection of structural variants (SV) from the alignment of sample DNA reads to the reference genome is an important problem in understanding human diseases. Long reads that can span repeat regions, along with an accurate alignment of these long reads play an important role in identifying novel SVs. Long read sequencers such as nanopore sequencing can address this problem by providing very long reads but with high error rates, making accurate alignment challenging. Many errors induced by nanopore sequencing have a bias because of the physics of the sequencing process and proper utilization of these error characteristics can play an important role in designing a robust aligner for SV detection problems. In this paper, we design and evaluate HQAlign, an aligner for SV detection using nanopore sequenced reads.

The key ideas of HQAlign include (i) using basecalled nanopore reads along with the nanopore physics to improve alignments for SVs (ii) incorporating SV specific changes to the alignment pipeline (iii) adapting these into existing state-of-the-art long read aligner pipeline, minimap2 (v2.24), for efficient alignments.

---

### Publication <a id='pub'></a>

If you find HQAlign is helpful, we appreciate your citation of its pre-print (https://doi.org/10.1101/2023.01.08.523172):

Dhaivat Joshi, Suhas Diggavi, Mark J.P. Chaisson, Sreeram Kannan, bioRxiv 2023.01.08.523172.

---

### Setup <a id='Installation'></a>

HQAlign has been tested under Ubuntu 16.04. Please follow the steps to setup:

###### Step 1
Download QAlign.
```
git clone https://github.com/joshidhaivat/HQAlign.git
```

### Usage <a id='use'></a>

###### Usage:
```python hqalign.py [-h] -r REF -i READS -o OUTPUT [-t THREADS] [-k KMER]```
```
python hqalign.py convert -r [/path/to/input/fasta/reference]
                          -i [/path/to/input/fasta/reads/folder]
                          -o [/path/to/output/directory]
                          -t [number of threads (default=4)]
                          -k [minimizer length for the hybrid step (default=18)]
```
