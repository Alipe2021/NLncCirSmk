# NLncCirSmk
 Snakemake workflow for Novel mRNA, LncRNA, and Circular RNA identification.

Author: Peng Liu

Email: sxliulian2012@hotmail.com

# Pre-installation

An effective software management framework is recommended. 

```txt
/opt
├── biosoft
│   ├── bin
│   │   ├── fastp -> /opt/biosoft/fastp/fastp
│   │   └── fastqc -> /opt/biosoft/FastQC/fastqc
├── miniconda3
│   ├── bin
│   ├── envs
│   │   ├── biotools
│   │   ├── NLncCirSmk
│   │   │   ├── bin
│   │   │   │   ├── bowtie
│   │   │   │   ├── bowtie2
│   │   │   │   ├── snakemake
```

# Installation

**Dependencies**

* conda > v3.8.5
* snakemake > v5.7.0
* fastp (version: >= 0.20.1)
* hisat2 (version: 2.2.1)
* bowtie2 (version: 2.3.5.1)
* samtools (version: 1.9)
* stringtie (version: 2.1.4)
* gffcompare (version: 0.11.2)
* gffread (version: 0.12.1)
* emboss transeq (version: 6.6.0)
* ...

Conda can be downloaded as part of the Anaconda or the Miniconda plattforms. We recommend to install miniconda3. Using Linux you can get it with:

```sh
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```

Snakemake can be install with:

```sh
conda create -c conda-forge -c bioconda -n NLncCirSmk python=3 snakemake
```

Detail installation guide from: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

other packages can be installed by conda:

```sh
# install normal softwares
conda install -n NLncCirSmk -c bioconda -c conda-forge fastp bowtie bowtie2 hisat2 bwa
conda install -n NLncCirSmk -c bioconda -c conda-forge stringtie samtools bedtools
conda install -n NLncCirSmk -c bioconda -c conda-forge gffread gffcompare emboss

# Pfam_scan
conda create -n pfam_scan_env -c bioconda pfam_scan

# CPC2 installation
conda create -n cpc2_py3_env -c conda-forge python=3.8
conda install -n cpc2_py3_env -c bioconda biopython=1.78
wget -c https://github.com/gao-lab/CPC2_standalone/archive/refs/tags/v1.0.1.tar.gz
tar zxvf CPC2_standalone-1.0.1.tar.gz

# copy all files in CPC2_standalone-1.0.1/ to your miniconda3 directory
cp -r CPC2_standalone-1.0.1/* /opt/miniconda3/envs/cpc2_py3_env/
cd /opt/miniconda3/envs/cpc2_py3_env/
export CPC_HOME="$PWD"
cd libs/libsvm
tar zxvf libsvm-3.18.tar.gz
cd libsvm-3.18
make clean && make

# CNCI installation
conda create -n cnci_py2_env -c conda-forge python=2.7
git clone git@github.com:www-bioinfo-org/CNCI.git
cp -r CNCI/* /opt/miniconda3/envs/cnci_py2_env/bin/
cd /opt/miniconda3/envs/cnci_py2_env/bin/
unzip libsvm-3.0.zip
cd libsvm-3.0
make
cd ..
chmod 755 CNCI.py compare.py filter_novel_lincRNA.py
chmod 755 Gtf.py gtf2Bed.pl draw_class_pie.R Table.py

# FEElnc installation
conda create -n feelnc_env -c bioconda feelnc

# CPAT3 installation
conda create -n cpat_py3_env -c conda-forge python=3.9 
conda install -n cpat_py3_env -c bioconda -c conda-forge pysam numpy r-base=4.1
pip3 install CPAT

# CIRI packages
## CIRIquant
conda create -n ciri_py2_env -c bioconda -c conda-forge python=2.7
conda install -n ciri_py2_env -c bioconda -c conda-forge bwa hisat2 stringtie samtools=1.9
source activate ciri_py2_env
pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple # for Chinese user
pip install CIRIquant
## CIRI2
wget -c https://sourceforge.net/projects/ciri/files/CIRI2/CIRI_v2.0.6.zip
unzip CIRI_v2.0.6.zip
## CIRI_AS
wget -c https://pilotfiber.dl.sourceforge.net/project/ciri/CIRI-AS/CIRI_AS_v1.2.pl -O /opt/miniconda3/envs/ciri_py2_env/bin/CIRI_AS.pl
## CIRI_Vis
wget -c https://phoenixnap.dl.sourceforge.net/project/ciri/CIRI-vis/CIRI-vis_v1.4.jar -O /opt/miniconda3/envs/ciri_py2_env/bin/CIRI-vis.jar
## CIRI_Full
wget -c https://sourceforge.net/projects/ciri/files/CIRI-full/CIRI-full_v2.0.zip
wget -c https://sourceforge.net/projects/ciri/files/CIRI-full/CIRI_Full_v2.1.1.jar/download -O CIRI_Full.jar

# CIRCexplorer2
conda create -n circexplorer2_py2_env -c bioconda -c conda-forge python=2.7 circexplorer2
conda install -n circexplorer2_py2_env -c bioconda  bwa star samtools bedtools 
conda install -n circexplorer2_py2_env -c conda-forge pysam pybedtools pandas docopt scipy
conda install -n circexplorer2_py2_env -c bioconda ucsc-genepredtogtf ucsc-gtftogenepred 
conda install -n circexplorer2_py2_env -c bioconda ucsc-bedgraphtobigwig ucsc-bedtobigbed
```



# Usage

## Data preparation

```sh
# 1. download maize reference genome fasta
wget -c ftp://ftp.gramene.org/pub/gramene/release-63/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
wget -c ftp://ftp.gramene.org/pub/gramene/release-63/fasta/zea_mays/cdna/Zea_mays.AGPv4.cdna.all.fa.gz
wget -c ftp://ftp.gramene.org/pub/gramene/release-63/fasta/zea_mays/cds/Zea_mays.AGPv4.cds.all.fa.gz
wget -c ftp://ftp.gramene.org/pub/gramene/release-63/fasta/zea_mays/ncrna/Zea_mays.AGPv4.ncrna.fa.gz
wget -c ftp://ftp.gramene.org/pub/gramene/release-63/fasta/zea_mays/pep/Zea_mays.AGPv4.pep.all.fa.gz
wget -c ftp://ftp.gramene.org/pub/gramene/release-63/gtf/zea_mays/Zea_mays.B73_RefGen_v4.48.gtf.gz
# 2. uncompression
gzip -cd Zea_mays.AGPv4.dna.toplevel.fa.gz |cut -f1 > /Pub/DataBase/Species/Zea_Mays/B73v4/zma_dna_v4.fa
gzip -cd Zea_mays.AGPv4.cdna.all.fa.gz |cut -f1 > /Pub/DataBase/Species/Zea_Mays/B73v4/zma_cdna_v4.fa
gzip -cd Zea_mays.AGPv4.cds.all.fa.gz |cut -f1 > /Pub/DataBase/Species/Zea_Mays/B73v4/zma_cds_v4.fa
gzip -cd Zea_mays.AGPv4.ncrna.fa.gz |cut -f1 > /Pub/DataBase/Species/Zea_Mays/B73v4/zma_ncrna_v4.fa
gzip -cd Zea_mays.AGPv4.pep.all.fa.gz |cut -f1 > /Pub/DataBase/Species/Zea_Mays/B73v4/zma_pep_v4.fa
gzip -cd Zea_mays.B73_RefGen_v4.48.gtf.gz > /Pub/DataBase/Species/Zea_Mays/B73v4/zma.v4.48.gtf
# 3. build index
bowtie-build --threads 24 genome.fa genome.fa 		# build bowtie index for genome
bowtie2-build --threads 24 genome.fa genome.fa  	# build bowtie2 index for genome
## build hisat2 index for genome
hisat2_extract_splice_sites.py zma.v4.48.gtf > zma.ss
hisat2_extract_exons.py zma.v4.48.gtf > zma.exon
hisat2-build -p 24 --ss zma.ss --exon zma.exon genome.fa genome.fa
## build bwa index
bwa index genome.fa
## ctrate genome dict
samtools dict -o zma_dna_v4.fa.dict zma_dna_v4.fa
```

Use some tools such as `bedparse`  to convert gtf format to bed12 format.

Fetch protein coding transcripts to build `zma_mrna.v4.44.gtf` and lncRNAs to build `zma_lncrna.v4.44.gtf`.

## Configure

**Example:**

```yaml
## fasta
dna: /Pub/DataBase/Species/Zea_Mays/B73v4/zma_dna_v4.fa
cds: /Pub/DataBase/Species/Zea_Mays/B73v4/zma_cds_v4.fa
cdna: /Pub/DataBase/Species/Zea_Mays/B73v4/zma_cdna_v4.fa
ncrna: /Pub/DataBase/Species/Zea_Mays/B73v4/zma_ncrna_v4.fa
rRNA: /Pub/DataBase/RNAcentral/v16.rRNA/PlantrRNA_RNACentralv16.fa
mirna: /Pub/DataBase/miRBase/22.1/zma_mature.fa
## annotation
bed: /Pub/DataBase/Species/Zea_Mays/B73v4/zma_v4.44.bed
gtf: /Pub/DataBase/Species/Zea_Mays/B73v4/zma.v4.44.gtf
mrna_gtf: /Pub/DataBase/Species/Zea_Mays/B73v4/zma_mrna.v4.44.gtf
lncrna_gtf: /Pub/DataBase/Species/Zea_Mays/B73v4/zma_lncrna.v4.44.gtf
go: /Pub/DataBase/GO/zma.go
## index
dict: /Pub/DataBase/Species/Zea_Mays/B73v4/zma_dna_v4.fa.dict
splice_sites: /Pub/DataBase/Species/Zea_Mays/B73v4/hisat2Index/zma_v4.44.ss
rRNA_bowtie2_index: /Pub/DataBase/RNAcentral/v16.rRNA/PlantrRNA_RNACentralv16.fa
genome_bowtie2_index: /Pub/DataBase/Species/Zea_Mays/B73v4/Bowtie2Index/genome.fa
genome_hisat2_index: /Pub/DataBase/Species/Zea_Mays/B73v4/hisat2Index/genome.fa
genome_bwa_index: /Pub/DataBase/Species/Zea_Mays/B73v4/BwaIndex/genome.fa
## db
pfamdb: /Pub/DataBase/Pfam/33.1/Pfam-A/Pfam-A.hmm
pfamdt: /Pub/DataBase/Pfam/33.1/Pfam-A/Pfam-A.hmm.dat
## DEG Analysis needed
compare_pairs: /Pub/Project/NLncCirSmk/compare_pairs.txt
sample_groups: /Pub/Project/NLncCirSmk/sample_groups.txt
## CIRI2 needed
ciri_quant_cfg: /Pub/Project/NLncCirSmk/ciri_quant_cfg.yaml
### Path to a folder where intermediate files will be written.
output_dir: /Pub/Project/NLncCirSmk/output.v1/

# Path to a YAML file with samples and their corresponding FASTQ files.
sample_list: /Pub/Project/NLncCirSmk/sample.yaml
```

**Note:**

1. yaml format
2. set only you need

## Samples

**Example:**

```yaml
# sample.yaml
CK_BML1234_0h1:
- /Pub/SeqData/RawData/CK280h1_R1.fastq.gz
- /Pub/SeqData/RawData/CK280h1_R2.fastq.gz
CK_BML1234_0h2:
- /Pub/SeqData/RawData/CK280h2_R1.fastq.gz
- /Pub/SeqData/RawData/CK280h2_R2.fastq.gz
```

## compare_pairs

**Example:**

```txt
SampleA SampleB
CK_BML1234_0h   CK_BML1234_6h
CK_BML1234_0h   CK_BML1234_18h
CK_BML1234_0h   CK_BML1234_36h
CK_BML1234_0h   T_BML1234_6h
CK_BML1234_0h   T_BML1234_18h
CK_BML1234_0h   T_BML1234_36h
```

## sample_groups

```txt
Sample  Group
CK_BML1234_0h1  CK_BML1234_0h
CK_BML1234_0h2  CK_BML1234_0h
CK_BML1234_6h1  CK_BML1234_6h
CK_BML1234_6h2  CK_BML1234_6h
CK_BML1234_18h1 CK_BML1234_18h
CK_BML1234_18h2 CK_BML1234_18h
CK_BML1234_36h1 CK_BML1234_36h
CK_BML1234_36h2 CK_BML1234_36h
```

## CIRI2 Needed

```yaml
# ciri_quant_cfg.yaml
name: zma
tools:
  bwa: /opt/miniconda3/envs/LncCirSmk/bin/bwa
  hisat2: /opt/miniconda3/envs/LncCirSmk/bin/hisat2
  stringtie: /opt/miniconda3/envs/LncCirSmk/bin/stringtie
  samtools: /opt/miniconda3/envs/LncCirSmk/bin/samtools

reference:
  fasta: /Pub/DataBase/Species/Zea_Mays/B73v4/zma_dna_v4.fa
  gtf: /Pub/DataBase/Species/Zea_Mays/B73v4/zma.v4.44.gtf
  bwa_index: /Pub/DataBase/Species/Zea_Mays/B73v4/BwaIndex/genome.fa
  hisat_index: /Pub/DataBase/Species/Zea_Mays/B73v4/hisat2Index/genome.fa
```

## Run pipeline

```sh
source activate NLncCirSmk
snakemake -p -s LncCir.smk -j <threads> --latency-wait 20 
```

* threads: The threads number for use. Less than max threads.

# FAQs

**Q:** error while loading shared libraries: libreadline.so.6

**A:** link `libreadline.so.?` to `libreadline.so.6`. example: `ln -s  /opt/miniconda3/envs/feelnc_env/lib/libreadline.so.8  /opt/miniconda3/envs/feelnc_env/lib/libreadline.so.6`

# Further

More features will be added soon ...