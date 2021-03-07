# LncCirSmk
 Snakemake workflow for LncRNA, Circular RNA, and novel mRNA idntification.

# Dependencies

* conda v3.8.5
* snakemake > v5.7.0
* fastp (version: 0.20.1)
* hisat2 (version: 2.2.1)
* bowtie2 (version: 2.3.5.1)
* samtools (version: 1.7)
* stringtie (version: 2.1.4)
* gffcompare (version: 0.11.2)
* gffread (version: 0.12.1)
* emboss transeq (version: 6.6.0)
* ...

Conda can be downloaded as part of the Anaconda or the Miniconda plattforms (Python 3.8). We recommend to install miniconda3. Using Linux you can get it with:

```sh
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```

Snakemake can be install with:

```sh
conda install -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Detail installation guid from: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

# Usage

```sh
snakemake -p -s LncCir.smk -j <threads>
```

* threads: The threads number for use. Less than max threads.

