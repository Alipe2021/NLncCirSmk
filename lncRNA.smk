#!/usr/bin/env python

import re, os, yaml
import pandas as pd
import json
import sys

''' Snakemake file for lncRNA and novel mRNA Analysis Workflow
    Author: Alipe
    E-mail: sxliulian2012@hotmail.com
      Date: 2021-01-12 22:55
'''

##========= Globals ==========
configfile: 'config.yaml'

## Set sample files
SAMPLEFILES = yaml.load(open(config['SampleList']), Loader=yaml.FullLoader)
SAMPLES = sorted(SAMPLEFILES.keys())

## set output Dir
RESULTDIR = config["ResultDir"]

## Set reference file
DNA = config["dna"]
GTF = config["gtf"]
CDS = config["cds"]
NCRNA = config["ncrna"]
TRANSCRIPT = config["transcript"]

rRNA = config["rRNA"]
PFAMDB = os.path.dirname(config["pfamdb"])

HISAT2_SPLICE_SITES = config["splice_sites"]
HISAT2_DNA_INDEX = config["genome_hisat2_index"]
BOWTIE2_rRNA_INDEX = config["rRNA_bowtie2_index"]

## Tinity Needed
TINITY_SAMPLE_LIST = config["TrinitySampleList"]

## DEGs Needs
COMPARE_PAIRS = config['compare_pairs']
SAMPLE_GROUPS = config['sample_groups']

##======== Rules ============
rule all:
    input:
    ## Step 00: Prepare data
        expand( RESULTDIR + "Step00.Prepare/{sample}.R1.fq.gz", sample=SAMPLES ),
        expand( RESULTDIR + "Step00.Prepare/{sample}.R2.fq.gz", sample=SAMPLES ),
    ## Step 01: Filter fastq
        expand( RESULTDIR + "Step01.FastqFilter/{sample}/{sample}.R1.fq.gz", sample=SAMPLES ),
        expand( RESULTDIR + "Step01.FastqFilter/{sample}/{sample}.R2.fq.gz", sample=SAMPLES ),
        expand( RESULTDIR + "Step01.FastqFilter/{sample}/{sample}.json", sample=SAMPLES ),
        expand( RESULTDIR + "Step01.FastqFilter/{sample}/{sample}.html", sample=SAMPLES ),
    ## Step 02: Filter rRNA
        # expand( RESULTDIR + "Step02.FilterrRNA/{sample}.bam", sample=SAMPLES ),
        # expand( RESULTDIR + "Step02.FilterrRNA/{sample}/un-conc-mate.1", sample=SAMPLES ),
        # expand( RESULTDIR + "Step02.FilterrRNA/{sample}/un-conc-mate.2", sample=SAMPLES ),
        # expand( RESULTDIR + "Step02.FilterrRNA/{sample}/rRNA_filter.ok", sample=SAMPLES ),
    ## Step 03: Alianment to Genome
        expand( RESULTDIR + "Step03.Hisat2Genome/{sample}.bam", sample=SAMPLES ),
        expand( RESULTDIR + "Step03.Hisat2Genome/{sample}.summary.txt", sample=SAMPLES ),
        expand( RESULTDIR + "Step03.Hisat2Genome/{sample}/un-conc-mate.ok", sample=SAMPLES ),
        expand( RESULTDIR + "Step03.Hisat2Genome/{sample}/un-conc-mate.1", sample=SAMPLES ),
        expand( RESULTDIR + "Step03.Hisat2Genome/{sample}/un-conc-mate.2", sample=SAMPLES ),
        expand( RESULTDIR + "Step03.Hisat2Genome/un-conc-mate/{sample}.R1.fq.gz", sample=SAMPLES ),
        expand( RESULTDIR + "Step03.Hisat2Genome/un-conc-mate/{sample}.R2.fq.gz", sample=SAMPLES ),
    ## Step 04-1: De novo Assembly Using Trinity
        # RESULTDIR + "Step04.TrinityAssembly/trinity_outdir/trinity_outdir.Trinity.fasta",
        # RESULTDIR + "Step04.TrinityAssembly/trinity_outdir/trinity_outdir.Trinity.ok
    ## Step 04-2: Reference Genome Guid Assembly
        expand( RESULTDIR + "Step04.StringtieAssembly/{sample}.gtf", sample=SAMPLES ),
    ## Step 05: Merge Assemblyed GTF
        RESULTDIR + "Step05.StringtieMerge/MergedList.txt",
        RESULTDIR + "Step05.StringtieMerge/StringtieMerged.gtf",
        RESULTDIR + "Step05.StringtieMerge/StringtieMerged.fa",
    ## Step 06: ReAssembly Transcripts Quantify
        expand( RESULTDIR + "Step06.ReStringtieAssembly/{sample}.gtf", sample=SAMPLES ),
        expand( RESULTDIR + "Step06.ReStringtieAssembly/{sample}.coverage.cov", sample=SAMPLES ),
        expand( RESULTDIR + "Step06.ReStringtieAssembly/{sample}.GeneAbund.txt", sample=SAMPLES ),
    ## Step 07: Gene and Transcripts Expression analysis
        RESULTDIR + "Step07.Expression/AssembledGtfList.txt",
        RESULTDIR + "Step07.Expression/gene_count_matrix.csv",
        RESULTDIR + "Step07.Expression/transcript_count_matrix.csv",
        RESULTDIR + "Step07.Expression/transcript_tpm_matrix.tsv",
    ## Step 08: DEG Analysis
        # RESULTDIR + "Step08.DEGs/transcript_count_filter_matrix.csv",
        # RESULTDIR + "Step08.DEGs/transcript_tpm_filter_matrix.csv",
        # RESULTDIR + "Step08.DGEbyEdgeR/DGEbyEdgeR.Result.ok",
    ## Step 09-1: LncRNA and novel mRNA identify by FEELnc
        RESULTDIR + "Step09.FEELncIdentify/Candidate_lncRNA_flt.gtf",
        RESULTDIR + "Step09.FEELncIdentify/Candidate_lncRNA_codpot.lncRNA.gtf",
        RESULTDIR + "Step09.FEELncIdentify/Candidate_lncRNA_classes.txt",
    ## Step 09-2-1: LncRNA and novel mRNA identify by CAPT
        RESULTDIR + "Step09.CPATIdentify/Spacies_Hexamer.tsv",
        RESULTDIR + "Step09.CPATIdentify/CpatMaize.feature.xls",
        RESULTDIR + "Step09.CPATIdentify/CpatMaize.logit.RData",
        RESULTDIR + "Step09.CPATIdentify/CpatMaize.make_logitModel.r",
        RESULTDIR + "Step09.CPATIdentify/CpatMaize.make_logitModel.ok",
        RESULTDIR + "Step09.CPATIdentify/CpatMaize.ORF_seqs.fa",
        RESULTDIR + "Step09.CPATIdentify/CpatMaize.ORF_prob.tsv",
        RESULTDIR + "Step09.CPATIdentify/CpatMaize.ORF_prob.best.tsv",
        RESULTDIR + "Step09.CPATIdentify/CpatMaize.no_ORF.txt",
        RESULTDIR + "Step09.CPATIdentify/CpatMaize.r",
        RESULTDIR + "Step09.CPATIdentify/CpatMaize.DetectORF.ok",
    ## Step 10: Compare Merged GTF to Reference GTF
        RESULTDIR + "Step10.GffCompare/GffCompared.tracking",
        RESULTDIR + "Step10.GffCompare/GffCompared.annotated.gtf",
        RESULTDIR + "Step10.GffCompare/GffCompared.stats",
        RESULTDIR + "Step10.GffCompare/GffCompared.loci",
        RESULTDIR + "Step10.GffCompare/GffCompared.ok",
    ## Step 11: LncRNA identify by CPC2,CNCI,Pfam ---- Fetch lncRNA class_code ixou, gtf, pep
        RESULTDIR + "Step11.LncRNAIdentify/Candidate_lncRNA.gtf",
        RESULTDIR + "Step11.LncRNAIdentify/Candidate_lncRNA.fa",
        RESULTDIR + "Step11.LncRNAIdentify/Candidate_lncRNA.pep.fa",
    ## Step 12: LncRNA identify by CPC2,CNCI,Pfam ---- LncRNA Protein Coding detection -- (CPC2,CNCI,PFAM)
        RESULTDIR + "Step12.CPC2LncRNAProteinPrediction/LncRNA_CPC2_Out.txt",
        RESULTDIR + "Step12.CPC2LncRNAProteinPrediction/LncRNA_CPC2_Noncoding.txt",
        RESULTDIR + "Step12.CNCILncRNAProteinPrediction/LncRNA_CNCI/CNCI.index",
        RESULTDIR + "Step12.CNCILncRNAProteinPrediction/LncRNA_CNCI_Noncoding.txt",
        RESULTDIR + "Step12.PfamScanLncRNAPrediction/LncRNA_PfamScan.Result.txt",
        RESULTDIR + "Step12.PfamScanLncRNAPrediction/LncRNA_PfamScan_Coding.txt",
    ## Step 13: Novel mRNA identify by CPC2,CNCI,Pfam ---- Fetch novel mRNA class_code uj, orf, pep
        RESULTDIR + "Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.gtf",
        RESULTDIR + "Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.fa",
        RESULTDIR + "Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.orf_cds.fa",
        RESULTDIR + "Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.orf_pep.fa",
    ## Step 14: Novel mRNA identify by CPC2,CNCI,Pfam ---- Novel mRNA Protein Coding detection (CPC2,CNCI,PFAM)
        RESULTDIR + "Step14.CPC2NovelmRNAProteinPrediction/NovelmRNA_CPC2_Out.txt",
        RESULTDIR + "Step14.CPC2NovelmRNAProteinPrediction/NovelmRNA_CPC2_Coding.txt",
        RESULTDIR + "Step14.CNCINovelmRNAProteinPrediction/NovelmRNA_CNCI/CNCI.index",
        RESULTDIR + "Step14.CNCINovelmRNAProteinPrediction/NovelmRNA_CNCI_Coding.txt",
        RESULTDIR + "Step14.PfamScanNovelmRNAPrediction/NovelmRNA_PfamScan.Result.txt",
        RESULTDIR + "Step14.PfamScanNovelmRNAPrediction/NovelmRNA_PfamScan_Coding.txt",
    ## Step 15: Cis and trans target gene
    ## Step 16: Target by known miRNA
    ## Step 17: 
    ## Step 18: 
    ## Step 19: 

    ## Report S01: Fastq Filter Result
        # RESULTDIR + "Report/S01.FastqFilter.jsonlist.txt",
        # RESULTDIR + "Report/S01.FastqFilter-State.tsv"
#########################################################################################
## ======== Step 0: Prepare Work ========
rule Analysis_00_PrepareFastQ:
    input:
        lambda wildcards:SAMPLEFILES[wildcards.sample]
    output:
       R1 =  RESULTDIR + "Step00.Prepare/{sample}.R1.fq.gz",
       R2 =  RESULTDIR + "Step00.Prepare/{sample}.R2.fq.gz"
    shell:
        "ln -sf {input[0]} {output.R1} && "
        "ln -sf {input[1]} {output.R2} "
## ======== Step 01: raw fastq filter ========
rule Analysis_01_FastqFilter:
    input:
        R1 =  RESULTDIR + "Step00.Prepare/{sample}.R1.fq.gz",
        R2 =  RESULTDIR + "Step00.Prepare/{sample}.R2.fq.gz"
    output:
        R1 = RESULTDIR + "Step01.FastqFilter/{sample}/{sample}.R1.fq.gz",
        R2 = RESULTDIR + "Step01.FastqFilter/{sample}/{sample}.R2.fq.gz",
        json = RESULTDIR + "Step01.FastqFilter/{sample}/{sample}.json",
        html = RESULTDIR + "Step01.FastqFilter/{sample}/{sample}.html"
    message:
        "Begin to filter fastq!"
    log:
        RESULTDIR + "logs/Step01.FastqFilter/{sample}.FastqFilter.log"
    params:
        "--detect_adapter_for_pe --qualified_quality_phred 15 --length_required 15 "
        "--unqualified_percent_limit 40 --n_base_limit 5 --average_qual 20 "
    threads:
        8
    shell:
        """
        fastp {params} -w {threads} -i {input.R1} -I {input.R2} -o {output.R1} \
        -O {output.R2} -j {output.json} -h {output.html} 2> {log}
        """
## ======== Step 02-1: filter rRNA using bowtie2 ========
rule Analysis_02_1_FilterrRNA:
    input:
        R1 = RESULTDIR + "Step01.FastqFilter/{sample}/{sample}.R1.fq.gz",
        R2 = RESULTDIR + "Step01.FastqFilter/{sample}/{sample}.R2.fq.gz",
        idx = BOWTIE2_rRNA_INDEX
    output:
        bam = temp(RESULTDIR + "Step02.FilterrRNA/{sample}.bam"),
        stat = RESULTDIR + "Step02.FilterrRNA/{sample}/rRNA_filter.ok",
        R1 = RESULTDIR + "Step02.FilterrRNA/{sample}/un-conc-mate.1",
        R2 = RESULTDIR + "Step02.FilterrRNA/{sample}/un-conc-mate.2"
    message:
        "Start to filter rRNA baseed on RNACentral plant rRNA database by bowtie2."
    log:
        RESULTDIR + "logs/Step02.FilterrRNA/{sample}.align.logs"
    threads:
        8
    params:
        opt = "-q --phred33 --sensitive --end-to-end --fr",
        out = RESULTDIR + "Step02.FilterrRNA/{sample}"
    shell:
        """
        bowtie2 {params.opt} -p {threads} --un-conc-gz {params.out} -x {input.idx} \
        -1 {input.R1} -2 {input.R2} 2>{log} | samtools sort -n -O Bam -@ 4 \
        -m 5G -o {output.bam} && echo Success > {output.stat}
        """
## ======== Step 02-2: Rename bowtie2 un-conc-gz ========
rule Analysis_02_2_ReNameBowtie2out:
    input:
        stat = RESULTDIR + "Step02.FilterrRNA/{sample}/rRNA_filter.ok",
        R1 = RESULTDIR + "Step02.FilterrRNA/{sample}/un-conc-mate.1",
        R2 = RESULTDIR + "Step02.FilterrRNA/{sample}/un-conc-mate.2"
    output:
        R1 = RESULTDIR + "Step02.FilterrRNA/un-conc-mate/{sample}.R1.fq.gz",
        R2 = RESULTDIR + "Step02.FilterrRNA/un-conc-mate/{sample}.R2.fq.gz"
    shell:
        """
            ln -sf {input.R1} {output.R1} && \
            ln -sf {input.R2} {output.R2}
        """
## ======== Step 03-1: Alignment to Genome ========
rule Analysis_03_1_Hisat2Genome:
    input:
        R1 = RESULTDIR + "Step02.FilterrRNA/un-conc-mate/{sample}.R1.fq.gz",
        R2 = RESULTDIR + "Step02.FilterrRNA/un-conc-mate/{sample}.R2.fq.gz",
        gtf = GTF,
        ss = HISAT2_SPLICE_SITES, # if hisat2-build with ss and exon, skip this step.
        idx = HISAT2_DNA_INDEX
    output:
        bam  = RESULTDIR + "Step03.Hisat2Genome/{sample}.bam",
        stat = RESULTDIR + "Step03.Hisat2Genome/{sample}.summary.txt",
        R1 = RESULTDIR + "Step03.Hisat2Genome/{sample}/un-conc-mate.1",
        R2 = RESULTDIR + "Step03.Hisat2Genome/{sample}/un-conc-mate.2",
        ok = RESULTDIR + "Step03.Hisat2Genome/{sample}/un-conc-mate.ok"
    message:
        "Start to map to genome with hisat2."
    log:
        RESULTDIR + "logs/Step03.Hisat2Genome/{sample}.align.logs"
    threads:
        8
    params:
        opt = "-q --phred33 --min-intronlen 20 --max-intronlen 150000 --dta",
        out = RESULTDIR + "Step03.Hisat2Genome/{sample}"
    shell:
        """
        hisat2 {params.opt} -p {threads} --summary-file {output.stat} -x {input.idx} \
         --known-splicesite-infile {input.ss} --un-conc-gz {params.out} \
         -1 {input.R1} -2 {input.R2} 2>{log} | samtools sort -O Bam \
         -@ 4 -m 5G -o {output.bam} && echo SUCCESS > {output.ok}
        """
## ======== Step 03-2: ReName un-conc-gz ========
rule Analysis_03_2_ReNameHisat2Out:
    input:
        ok = RESULTDIR + "Step03.Hisat2Genome/{sample}/un-conc-mate.ok",
        R1 = RESULTDIR + "Step03.Hisat2Genome/{sample}/un-conc-mate.1",
        R2 = RESULTDIR + "Step03.Hisat2Genome/{sample}/un-conc-mate.2"
    output:
        R1 = RESULTDIR + "Step03.Hisat2Genome/un-conc-mate/{sample}.R1.fq.gz",
        R2 = RESULTDIR + "Step03.Hisat2Genome/un-conc-mate/{sample}.R2.fq.gz"
    shell:
        """
            ln -sf {input.R1} {output.R1} && \
            ln -sf {input.R2} {output.R2}
        """
## ======== Step 04-1: Trinity Assembly ========
rule Analysis_04_1_TrinityAssembly:
    input:
        samList = TINITY_SAMPLE_LIST
    output:
        fa = RESULTDIR + "Step04.TrinityAssembly/trinity_outdir/trinity_outdir.Trinity.fasta",
        ok = RESULTDIR + "Step04.TrinityAssembly/trinity_outdir/trinity_outdir.Trinity.ok"
    log:
        RESULTDIR + "logs/Step04.TrinityAssembly/TrinityAssembly.log"
    message:
        "Start De novo Assembly useing trinity!"
    resources:
        mem_mb=100000
    threads:
        20
    params:
        opt = "--seqType fq --max_memory 100G --full_cleanup",
        outdir = RESULTDIR + "Step04.TrinityAssembly/trinity_outdir/"
    shell:
        """
        source activate tinity_env && \
        Trinity {params.opt} --samples_file {input.samList} --CPU {threads} --output {params.outdir} && \
        if [[ -s {output.fa} ]]; then echo SUCCESS ! > {output.ok}; fi && \
        conda deactivate
        """
## ======== Step 04-2: Stringtie Assembly ========
rule Analysis_04_2_StringtieAssembly:
    input:
        bam = RESULTDIR + "Step03.Hisat2Genome/{sample}.bam",
        gtf = GTF,
    output:
        gtf = RESULTDIR + "Step04.StringtieAssembly/{sample}.gtf"
    message:
        "Start to assembly use stringtie."
    log:
        RESULTDIR + "logs/Step04.StringtieAssembly/{sample}.assembly.logs"
    threads:
        8
    params:
        "-m 200 -l STRG -a 10 --conservative -g 50 -u"
    shell:
        """
        stringtie {input.bam} -G {input.gtf} -o {output.gtf} -p {threads} {params} 
        """
## ======== Step 05-1: Merge transcript ========
rule Analysis_05_1_MakeMergeList:
    input:
        expand(RESULTDIR + "Step04.StringtieAssembly/{sample}.gtf", sample = SAMPLES)
    output:
        lst = RESULTDIR + "Step05.StringtieMerge/MergedList.txt"
    run:
        with open(output.lst, 'w') as f:
            for gtf in input:
                print(gtf, file=f)
## ======== Step 05-2: StringTie Merge ========
rule Analysis_05_2_StingtieMerge:
    input:
        lst = RESULTDIR + "Step05.StringtieMerge/MergedList.txt",
        gtf = GTF,
        fasta = DNA
    output:
        gtf = RESULTDIR + "Step05.StringtieMerge/StringtieMerged.gtf"
    log:
        RESULTDIR + "logs/Step05.StringtieMerge/StringtieMerge.logs"
    threads:
        8
    params:
        "-m 200 -c 3"
    shell:
        """
        stringtie --merge {params} -p {threads} -G {input.gtf} -o {output.gtf} {input.lst} 2>{log} 
        """
## ======== Step 05-3: StringTie Merged fasta Get ========
rule Analysis_05_3_StingtieMergeFastaGet:
    input:
        refDNA = DNA,
        gtf = RESULTDIR + "Step05.StringtieMerge/StringtieMerged.gtf"
    output:
        fna = RESULTDIR + "Step05.StringtieMerge/StringtieMerged.fa"
    log:
        RESULTDIR + "logs/Step05.StringtieMerge/StringtieMergeFastaGet.logs"
    threads:
        8
    shell:
        """
        gffread -w {output.fna} -g {input.refDNA} {input.gtf} 2> {log}
        """
## ======== Step 06: Re Assembly by Stringtie ========
rule Analysis_06_ReStringtieAssemble:
    input:
        bam = RESULTDIR + "Step03.Hisat2Genome/{sample}.bam",
        gtf = RESULTDIR + "Step05.StringtieMerge/StringtieMerged.gtf",
    output:
        gtf = RESULTDIR + "Step06.ReStringtieAssembly/{sample}.gtf",
        cov = RESULTDIR + "Step06.ReStringtieAssembly/{sample}.coverage.cov",
        abu = RESULTDIR + "Step06.ReStringtieAssembly/{sample}.GeneAbund.txt"
    message:
        "ReStart to assembly use stringtie."
    log:
        RESULTDIR + "logs/Step06.ReStringtieAssembly/{sample}.assembly.logs"
    priority: 
        50
    threads:
        8
    params:
        opt = "-m 200 -l STRG -a 10 -c 1 -e",
        bal = RESULTDIR + "Step06.ReStringtieAssembly/Ballgown/{sample}"
    shell:
        """
        stringtie {input.bam} -G {input.gtf} -o {output.gtf} -p {threads} {params.opt} \
        -C {output.cov} -A {output.abu} -b {params.bal}
        """
## ======== Step 07-1: Espression, Get assembled gtf list ========
rule Analysis_07_1_GetAssembledGtfList:
    input:
        expand(RESULTDIR + "Step06.ReStringtieAssembly/{sample}.gtf", sample = SAMPLES)
    output:
        lst = RESULTDIR + "Step07.Expression/AssembledGtfList.txt"
    threads:
        1
    run:
        with open(output.lst, 'w') as f:
            for gtf in input:
                sample = re.search(r'Step06\.ReStringtieAssembly/(.*).gtf', gtf).group(1)
                outline = "\t".join((sample, gtf))
                print(outline, file=f)
## ======== Step 07-2: Espression, Get gene and transcripts count matrix ========
rule Analysis_07_2_GetCountAndTPMMatrix:
    input:
        RESULTDIR + "Step07.Expression/AssembledGtfList.txt"
    output:
        gene = RESULTDIR + "Step07.Expression/gene_count_matrix.csv",
        mrna = RESULTDIR + "Step07.Expression/transcript_count_matrix.csv",
        tran = RESULTDIR + "Step07.Expression/transcript_tpm_matrix.tsv"
    log:
        cnt = RESULTDIR + "logs/Step07.Expression/MakeCountMatrix.log",
        tpm = RESULTDIR + "logs/Step07.Expression/MakeTPMMatrix.log"
    params:
        "-l 145 -s MSTRG"
    shell:
        """
        prepDE.py -i {input} -g {output.gene} -t {output.mrna} {params} 2> {log.cnt} && \
        perl Scripts/GetTPMFromSringtieGtfList.pl {input} > {output.tran} 2> {log.tpm}
        """
## ======== Step 08: Different Expression Analysis by edgeR ========
rule Analysis_08_DGEbyDEGSeq2:
    input:
        gcm = RESULTDIR + "Step07.Expression/gene_count_matrix.csv",
        mcm = RESULTDIR + "Step07.Expression/transcript_count_matrix.csv",
        cp = COMPARE_PAIRS,
        sg = SAMPLE_GROUPS
    output:
        ok = RESULTDIR + "Step08.DGEbyDEGSeq2/DGEbyDEGSeq2.Result.ok"
    log:
        RESULTDIR + "logs/Step08.DGEbyDEGSeq2/DGEbyDEGSeq2.log"
    threads:
        8
    params:
        opt = "--foldchange 2 --fdr 0.05",
        out_dir = RESULTDIR + "Step08.DGEbyDEGSeq2/"
    shell:
        """
        source activate r_env && \
        Rscript Scripts/DGEs.R {params} -t {threads} --matrix {input.gcm} --group {input.sg} --compare {input.cp} -o {params.out_dir} && \
        Rscript Scripts/DGEs.R {params} -t {threads} -m {input.mcm} -g {input.sg} -c {input.cp} -o {params.out_dir} && \
        echo SUCCESS > {output.ok} && \
        conda deactivate
        """
## ======== Step 09-1-1: lncRNA and novel mRNA prediction by FEElnc ========
rule Analysis_09_1_1_FEELnc_filter:
    input:
        refGtf = GTF,
        newGtf = RESULTDIR + "Step05.StringtieMerge/StringtieMerged.gtf"
    output:
        flt = RESULTDIR + "Step09.FEELncIdentify/Candidate_lncRNA_flt.gtf"
    log:
        RESULTDIR + "logs/Step09.FEELncIdentify/FEELncPrediction.flt.log",
    threads:
        20
    params:
        "-s 200 -b transcript_biotype=protein_coding --monoex=0 --biex=25 -l FALSE"        
    shell:
        """
        source activate /aucluster/auhpc1_data/pub/miniconda3/envs/feelnc_install_dir && \
        FEELnc_filter.pl {params} -p {threads} -i {input.newGtf} -a {input.refGtf} -o {log} > {output.flt} && \
        conda deactivate
        """
## ======== Step 09-1-2: lncRNA and novel mRNA prediction by FEElnc ========
rule Analysis_09_1_2_FEELnc_codpot:
    input:
        refFa  = DNA,
        refGtf = GTF,
        flt = RESULTDIR + "Step09.FEELncIdentify/Candidate_lncRNA_flt.gtf"
    output:
        cod = RESULTDIR + "Step09.FEELncIdentify/Candidate_lncRNA_codpot.lncRNA.gtf"
    log:
        RESULTDIR + "logs/Step09.FEELncIdentify/FEELncPrediction.codpot.log"
    threads:
        20
    params:
        cod_opt1 = "-b transcript_biotype=protein_coding --mode=shuffle --sizeinter=0.75",
        cod_opt2 = "--learnorftype=3 --testorftype=3 --ntree 500 --seed=1234",
        cod_dir = RESULTDIR + "Step09.FEELncIdentify/",
        cod_name = "Candidate_lncRNA_codpot"
    shell:
        """
        source activate /aucluster/auhpc1_data/pub/miniconda3/envs/feelnc_install_dir && \
        FEELnc_codpot.pl {params.cod_opt1} {params.cod_opt2} -p {threads} -i {input.flt} -a {input.refGtf} -g {input.refFa} --outdir={params.cod_dir} --outname={params.cod_name} 2> {log} && \
        conda deactivate
        """
## ======== Step 09-1-3: lncRNA and novel mRNA prediction by FEElnc ========
rule Analysis_09_1_3_FEELnc_classifier:
    input:
        refFa  = DNA,
        refGtf = GTF,
        cod = RESULTDIR + "Step09.FEELncIdentify/Candidate_lncRNA_codpot.lncRNA.gtf"
    output:        
        cla = RESULTDIR + "Step09.FEELncIdentify/Candidate_lncRNA_classes.txt"
    log:
        RESULTDIR + "logs/Step09.FEELncIdentify/FEELncPrediction.class.log"
    threads:
        20
    params:
        ""
    shell:
        """
        source activate /aucluster/auhpc1_data/pub/miniconda3/envs/feelnc_install_dir && \
        FEELnc_classifier.pl {params} -i {input.cod} -a {input.refGtf} > {output.cla} 2> {log} && \
        conda deactivate
        """
## ======== Step 09-2-1: LncRNA and Novel mRNA Identify by CPAT -- Build Hexamer Table ========
rule Analysis_09_2_1_CPATBuildHexamerTable:
    input:
        cds = CDS,
        ncrna = NCRNA
    output:
        RESULTDIR + "Step09.CPATIdentify/Spacies_Hexamer.tsv"
    log:
        RESULTDIR + "logs/Step09.CPATIdentify/CPATBuildHexamerTable.log"
    threads:
        1
    params:
        ""
    shell:
        """
        source activate cpat_env && \
        make_hexamer_tab.py {params} -c {input.cds} -n {input.ncrna} > {output} 2> {log} && \
        conda deactivate
        """
## ======== Step 09-2-2: LncRNA and Novel mRNA Identify by CPAT -- Build Logit Model ========
rule Analysis_09_2_2_CPATBuildLogitModel:
    input:
        transcript = TRANSCRIPT,
        ncrna = NCRNA,
        hex = RESULTDIR + "Step09.CPATIdentify/Spacies_Hexamer.tsv"
    output:
        feature = RESULTDIR + "Step09.CPATIdentify/CpatMaize.feature.xls",
        rdata = RESULTDIR + "Step09.CPATIdentify/CpatMaize.logit.RData",
        rscipt = RESULTDIR + "Step09.CPATIdentify/CpatMaize.make_logitModel.r",
        ok = RESULTDIR + "Step09.CPATIdentify/CpatMaize.make_logitModel.ok"
    log:
        RESULTDIR + "logs/Step09.CPATIdentify/CPATBuildLogitModel.log"
    params:
        opt = "--min-orf=30",
        prefix = RESULTDIR + "Step09.CPATIdentify/CpatMaize"
    shell:
        """
        source activate cpat_env && \
        make_logitModel.py {params.opt} -x {input.hex} -c {input.transcript} -n {input.ncrna} -o {params.prefix} --log-file={log} && \
        if [[ -s {output.rdata} ]]; then echo SUCCESS! > {output.ok};  fi && \
        conda deactivate
        """
## ======== Step 09-2-3: LncRNA and Novel mRNA Identify by CPAT -- Detect ORF ========
rule Analysis_09_2_3_CPATtoDetectORF:
    input:
        hex = RESULTDIR + "Step09.CPATIdentify/Spacies_Hexamer.tsv",
        rdata = RESULTDIR + "Step09.CPATIdentify/CpatMaize.logit.RData",
        fasta = RESULTDIR + "Step05.StringtieMerge/StringtieMerged.fa"
    output:
        seqs = RESULTDIR + "Step09.CPATIdentify/CpatMaize.ORF_seqs.fa",
        peob = RESULTDIR + "Step09.CPATIdentify/CpatMaize.ORF_prob.tsv",
        best = RESULTDIR + "Step09.CPATIdentify/CpatMaize.ORF_prob.best.tsv",
        noorf = RESULTDIR + "Step09.CPATIdentify/CpatMaize.no_ORF.txt",
        rscript = RESULTDIR + "Step09.CPATIdentify/CpatMaize.r",
        ok = RESULTDIR + "Step09.CPATIdentify/CpatMaize.DetectORF.ok"
    log:
        RESULTDIR + "logs/Step09.CPATIdentify/CPATtoDetectORF.log"
    params:
        opt = "--top-orf=5 --antisense --min-orf=75 --width=100 --best-orf=p",
        prefix = RESULTDIR + "Step09.CPATIdentify/CpatMaize"
    threads:
        1
    shell:
        """
        source activate cpat_env && \
        cpat.py {params.opt} -x {input.hex} -d {input.rdata} -g {input.fasta} -o {params.prefix} --log-file={log} && \
        if [[ -s {output.best} ]]; then echo SUCCESS! > {output.ok} ; fi && \
        conda deactivate
        """
## ======== Step 10: Compare merged gtf to the reference gtf with gffcompare ========
rule Analysis_10_GffCompare2ReferenceGTF:
    input:
        ref_gtf = GTF,
        new_gtf = RESULTDIR + "Step05.StringtieMerge/StringtieMerged.gtf"
    output:
        out_gtf = RESULTDIR + "Step10.GffCompare/GffCompared.annotated.gtf",
        stats = RESULTDIR + "Step10.GffCompare/GffCompared.stats",
        track = RESULTDIR + "Step10.GffCompare/GffCompared.tracking",
        loci = RESULTDIR + "Step10.GffCompare/GffCompared.loci",
        ok = RESULTDIR + "Step10.GffCompare/GffCompared.ok"
    log:
        RESULTDIR + "logs/Step10.GffCompare/GffComapare.logs"
    conda:
        "envs/lncRNA_smk.yaml"
    threads:
        4
    params:
        opt = "-T -R",
        pfx = RESULTDIR + "Step10.GffCompare/GffCompared",
    shell:
        """
        gffcompare {params.opt} -r {input.ref_gtf} -o {params.pfx} {input.new_gtf} 2> {log} && \
        echo Success > {output.ok}
        """
## ======== Step 11-1: Fetch lncRNA gtf ========
rule Analysis_11_1_FetchCandidatelncRNAGTF:
    input:
        gtf = RESULTDIR + "Step10.GffCompare/GffCompared.annotated.gtf"
    output:
        gtf = RESULTDIR + "Step11.LncRNAIdentify/Candidate_lncRNA.gtf"
    threads:
        1
    params:
        "iuxo"
    shell:
        """
            perl Scripts/Fetch.class_code_gtf.pl {input} {params} > {output}
        """
## ======== Step 11-2: Fetch lncRNA fasta ========
rule Analysis_11_2_FetchCandidatelncRNAFasta:
    input:
        fasta = DNA,
        gtf = RESULTDIR + "Step11.LncRNAIdentify/Candidate_lncRNA.gtf"
    output:
        fna = RESULTDIR + "Step11.LncRNAIdentify/Candidate_lncRNA.fa",
        tmp = temp(RESULTDIR + "Step11.LncRNAIdentify/Candidate_lncRNA.pep"),
        pep = RESULTDIR + "Step11.LncRNAIdentify/Candidate_lncRNA.pep.fa"
    log:
        RESULTDIR + "logs/Step11.LncRNAIdentify/Candaidate_lncRNA_fa_fetch.log"
    threads:
        1
    shell:
        """
        gffread -w {output.fna} -g {input.fasta} {input.gtf} && transeq {output.fna} \
        {output.tmp} && sed 's/*//g' {output.tmp} > {output.pep} 2> {log}
        """
## ======== Step 12-1: lncRNA fasta protein prodiction CPC2 ========
rule Analysis_12_1_CPC2LncRNAProteinPrediction_1:
    input:
        fa = RESULTDIR + "Step11.LncRNAIdentify/Candidate_lncRNA.fa"
    output:
        cpc = RESULTDIR + "Step12.CPC2LncRNAProteinPrediction/LncRNA_CPC2_Out.txt"
    log:
        RESULTDIR + "logs/Step12.CPC2LncRNAProteinPrediction/CPC2_Prediction.log"
    conda:
        "envs/lncRNA_smk.yaml"
    threads:
        1
    params:
        wjj = RESULTDIR + "Step12.CPC2LncRNAProteinPrediction/",
        pfx = RESULTDIR + "Step12.CPC2LncRNAProteinPrediction/LncRNA_CPC2_Out"
    shell:
        """
        source activate lncRNA_smk && \
        if [[ ! -d {params.wjj} ]]; then mkdir -p {params.wjj}; fi && \
        CPC2.py -i {input.fa} -o {params.pfx} 2> {log} ; fi && \
        conda deactivate
        """

## ======== Step 12-1: lncRNA fasta protein prodiction CPC2 ========
rule Analysis_12_1_CPC2LncRNAProteinPrediction_2:
    input:
        RESULTDIR + "Step12.CPC2LncRNAProteinPrediction/LncRNA_CPC2_Out.txt"    
    output:
        RESULTDIR + "Step12.CPC2LncRNAProteinPrediction/LncRNA_CPC2_Noncoding.txt"
    shell:
        """
        grep -w "noncoding" {input} |cut -f1 > {output}
        """
## ======== Step 12-2: lncRNA fasta protein prediction CNCI ========
rule Analysis_12_2_CNCILncRNAProteinPrediction:
    input:
        RESULTDIR + "Step11.LncRNAIdentify/Candidate_lncRNA.fa"
    output:
        cnci = protected(RESULTDIR + "Step12.CNCILncRNAProteinPrediction/LncRNA_CNCI/CNCI.index"),
        log  =  RESULTDIR + "Step12.CNCILncRNAProteinPrediction/LncRNA_CNCI.log",
        noc  = RESULTDIR + "Step12.CNCILncRNAProteinPrediction/LncRNA_CNCI_Noncoding.txt"
    log:
        RESULTDIR + "logs/Step12.CNCILncRNAProteinPrediction/LncRNA_CNCI_Prediction.log"
    conda:
        "envs/lncRNA_smk.yaml"
    threads:
        10
    params:
        opt = "-m ve",
        out = RESULTDIR + "Step12.CNCILncRNAProteinPrediction/LncRNA_CNCI"
    shell:
        """
        CNCI.py -f {input} -o {params.out} -p {threads} {params.opt} 2> {log} && \
        grep -w "noncoding" {output.cnci} | cut -f1 > {output.noc}
        """
## ======== Step 12-3-1: lncRNA fasta protein prediction PfamScan ========
rule Analysis_12_3_PfamScanLncRNAPrediction_1:
    input:
        fa = RESULTDIR + "Step11.LncRNAIdentify/Candidate_lncRNA.pep.fa"
    output:
        res = RESULTDIR + "Step12.PfamScanLncRNAPrediction/LncRNA_PfamScan.Result.txt"
    log:
        RESULTDIR + "logs/Step12.PfamScanLncRNAPrediction/LncRNA_PfamScan.log"
    threads:
        20
    params:
        db = PFAMDB
    shell:
        """
        source activate pfam_scan_env && \
        pfam_scan.pl -cpu {threads} -fasta {input.fa} -dir {params.db} -outfile {output.res} 2> {log} && \
        conda deactivate
        """
## ======== Step 12-3-2: lncRNA fasta protein prediction PfamScan ========
rule Analysis_12_3_PfamScanLncRNAPrediction_2:
    input:
        RESULTDIR + "Step12.PfamScanLncRNAPrediction/LncRNA_PfamScan.Result.txt"
    output:
        RESULTDIR + "Step12.PfamScanLncRNAPrediction/LncRNA_PfamScan_Coding.txt"
    shell:
        """
        grep -v "#" {input} |cut -f1 > {output}
        """
## ======== Step 13-1: Fetch Candidate Novel mRNA GTF ========
rule Analysis_13_1_FetchCandidateNovelmRNAGTF:
    input:
        gtf = RESULTDIR + "Step10.GffCompare/GffCompared.annotated.gtf"
    output:
        gtf = RESULTDIR + "Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.gtf"
    log:
        RESULTDIR + "logs/Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.gtf.log"
    threads:
        1
    params:
        "ju"
    shell:
        """
            perl Scripts/Fetch.class_code_gtf.pl {input} {params} > {output} 2> {log}
        """      
## ======== Step 13-2: Fetch Candidate Novel mRNA fasta ========
rule Analysis_13_2_FetchCandidateNovelmRNAFasta:
    input:
        fasta = DNA,
        gtf = RESULTDIR + "Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.gtf"
    output:
        fna = RESULTDIR + "Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.fa"
    log:
        RESULTDIR + "logs/Step13.NovelmRNAIdentify/Candidate_Novel_mRNA_fetch.log"
    threads:
        1
    shell:
        """
        gffread -w {output.fna} -g {input.fasta} {input.gtf} 2> {log}
        """
## ======== Step 13-3: Fetch Candidate Novel mRNA ORF ========
rule Analysis_13_3_TransDecoderLongOrfs:
    input:
        fna = RESULTDIR + "Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.fa"
    output:
        fna = protected(RESULTDIR + "Step13.NovelmRNAIdentify/TransdecoderORF/longest_orfs.cds"),
        pep = protected(RESULTDIR + "Step13.NovelmRNAIdentify/TransdecoderORF/longest_orfs.pep"),
        gff = protected(RESULTDIR + "Step13.NovelmRNAIdentify/TransdecoderORF/longest_orfs.gff3")
    log:
        RESULTDIR + "logs/Step13.NovelmRNAIdentify/TransdecoderORF_Fetch.log"
    threads:
        1
    params:
        opt = "-m 100",
        odr = RESULTDIR + "Step13.NovelmRNAIdentify/TransdecoderORF/"
    shell:
        """
        TransDecoder.LongOrfs -t {input.fna} {params.opt} --output_dir {params.odr} 2> {log}
        """
## ======== Step 13-4: Make ln to Candidate Novel mRNA ORF ========
rule Analysis_13_4_MakeLink2ORF:
    input:
        fna = RESULTDIR + "Step13.NovelmRNAIdentify/TransdecoderORF/longest_orfs.cds",
        pep = RESULTDIR + "Step13.NovelmRNAIdentify/TransdecoderORF/longest_orfs.pep"
    output:
        cds = RESULTDIR + "Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.orf_cds.fa",
        orf = RESULTDIR + "Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.orf_pep.fa"
    shell:
        """
        ln -sf {input.fna} {output.cds} && \
        ln -sf {input.pep} {output.orf}        
        """
## ======== Step 14-1: Candidate Novel mRNA fasta protein prodict CPC2 ========
rule Analysis_14_1_CPC2NovelmRNAProteinPrediction_1:
    input:
        cds = RESULTDIR + "Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.orf_cds.fa"
    output:
        cpc = RESULTDIR + "Step14.CPC2NovelmRNAProteinPrediction/NovelmRNA_CPC2_Out.txt"
    log:
        RESULTDIR + "logs/Step14.CPC2NovelmRNAProteinPrediction/CPC2_Prediction.log"
    conda:
        "envs/lncRNA_smk.yaml"
    threads:
        1
    params:
        wjj = RESULTDIR + "Step12.CPC2NovelmRNAProteinPrediction/",
        pfx = RESULTDIR + "Step14.CPC2NovelmRNAProteinPrediction/NovelmRNA_CPC2_Out"
    shell:
        """
        source activate lncRNA_smk && \
        if [[ !-d {params.wjj} ]]; then mkdir -p {params.wjj}; fi && \
        CPC2.py -i {input.cds} -o {params.pfx} 2> {log} && \
        conda deactivate
        """
## ======== Step 14-1: Candidate Novel mRNA fasta protein prodict CPC2 ========
rule Analysis_14_1_CPC2NovelmRNAProteinPrediction_2:
    input:
        RESULTDIR + "Step14.CPC2NovelmRNAProteinPrediction/NovelmRNA_CPC2_Out.txt"
    output:
        RESULTDIR + "Step14.CPC2NovelmRNAProteinPrediction/NovelmRNA_CPC2_Coding.txt"
    shell:
        """
        grep -w "coding" {input} |cut -f1 > {output}
        """
## ======== Step 14-2: Candidate Novel mRNA fasta protein prediction CNCI ========
rule Analysis_14_2_CNCINovelmRNAProteinPrediction_1:
    input:
        RESULTDIR + "Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.fa"
    output:
        cnci = protected(RESULTDIR + "Step14.CNCINovelmRNAProteinPrediction/NovelmRNA_CNCI/CNCI.index"),
        log =  RESULTDIR + "Step14.CNCINovelmRNAProteinPrediction/NovelmRNA_CNCI.log"
    log:
        RESULTDIR + "logs/Step14.CNCINovelmRNAProteinPrediction/NovelmRNA_CNCI_Prediction.log"
    conda:
        "envs/lncRNA_smk.yaml"
    threads:
        10
    params:
        opt = "-m ve",
        out = RESULTDIR + "Step14.CNCINovelmRNAProteinPrediction/NovelmRNA_CNCI"
    shell:
        """
        source activate lncRNA_smk && \
        CNCI.py -f {input} -o {params.out} -p {threads} {params.opt} 2> {log} && \
        conda deactivate
        """
## ======== Step 14-2: Candidate Novel mRNA fasta protein prediction CNCI ========
rule Analysis_14_2_CNCINovelmRNAProteinPrediction_2:
    input:
        RESULTDIR + "Step14.CNCINovelmRNAProteinPrediction/NovelmRNA_CNCI/CNCI.index"
    output:
        RESULTDIR + "Step14.CNCINovelmRNAProteinPrediction/NovelmRNA_CNCI_Coding.txt"
    shell:
        """
        grep -w "coding" {output} | cut -f1 > {output}
        """
## ======== Step 14-3: Candidate Novel mRNA fasta protein prediction Pfam ========
rule Analysis_14_3_PfamScanNovelmRNAPrediction_1:
    input:
        RESULTDIR + "Step13.NovelmRNAIdentify/Candidate_Novel_mRNA.orf_pep.fa"
    output:
        RESULTDIR + "Step14.PfamScanNovelmRNAPrediction/NovelmRNA_PfamScan.Result.txt"
    log:
        RESULTDIR + "logs/Step14.PfamScanNovelmRNAPrediction/NovelmRNA_PfamScan.log"
    conda:
        "envs/pfam_scan_env.yaml"
    threads:
        20
    params:
        db = PFAMDB
    shell:
        """
        source activate pfam_scan_env && \
        pfam_scan.pl -cpu {threads} -fasta {input} -dir {params.db} -outfile {output} 2> {log} && \
        conda deactivate
        """
## ======== Step 14-3: Candidate Novel mRNA fasta protein prediction Pfam ========
rule Analysis_14_3_PfamScanNovelmRNAPrediction_2:
    input:
        RESULTDIR + "Step14.PfamScanNovelmRNAPrediction/NovelmRNA_PfamScan.Result.txt"
    output:
        RESULTDIR + "Step14.PfamScanNovelmRNAPrediction/NovelmRNA_PfamScan_Coding.txt"
    shell:
        """
        grep -v "#" {input} |cut -f1 > {output}
        """

## ======== Step 15: LncRNA target gene predection ========




#########################################################################################
## -------- Report 01: Fastq Filter Result --------
rule Report_01_FastqFilterState:
    input:
        expand( RESULTDIR + "Step01.FastqFilter/{sample}/{sample}.json", sample=SAMPLES )
    output:
        lst = RESULTDIR + "Report/S01.FastqFilter.jsonlist.txt"
    run:
        with open(output.lst, 'w') as jo:
            for filename in input:
                print(filename, file=jo)

## -------- XXXXXXXXXXXXXXXXXXXXXXXXXXX --------------------
rule Report_01_FastqFilterState2:
    input:
        RESULTDIR + "Report/S01.FastqFilter.jsonlist.txt"
    output:
        tsv = RESULTDIR + "Report/S01.FastqFilter-State.tsv"
    shell:
        """
        source activate data_env && \
        python Scripts/GetFastqFilterState.py {input} {output}
        """
## -------- Report 02: rRNA Filter Result --------
## -------- Report 03: Genome alianment Result --------
## -------- Report 04: Genome guid assembly Result --------
## -------- Report 05: Genome guid assembly Result --------
## -------- Report 06: Genome guid assembly Result --------
## -------- Report 07: Genome guid assembly Result --------
## -------- Report 08: Genome guid assembly Result --------
## -------- Report 09: Genome guid assembly Result --------
## -------- Report 10: Genome guid assembly Result --------
## -------- Report 00: Output html report --------
