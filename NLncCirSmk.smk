#!/usr/bin/env python

"""
    Title:  Snakemake File for Novel mRNA, LncRNA, and Circular RNA Analysis
   Author:  Liu Peng
   E-mail:  sxliulian2012@hotmail.com
  Version:  1.0
     Date:  2021-05-08 16:31
   Update:  2021-09-14 09:50
"""
#
import re, os, yaml, json, sys
from typing_extensions import ParamSpecArgs
#
# 1. Raw Reads --> Clean Reads --[filter rRNA]--> rRNA-free Reads
# 2. rRNA-free Reads --[mapp to genome]--> Assembly --> Merged Info -
#    -> Re-assembly --> LncRNA/Novel mRNA Identification --> DELs/DEMs -
#    -> Target / TFs Prediction / WGCNA / miRNA Sponge ...
# 3. rRNA-free Reads --[mapp to genome]--> CircRNA Identification - 
#    -> CircRNA Quantification --> DECs --> miRNA Sponge ...
# 4. CircRNA Identification
#   a. CIRI2
#   b. find_circRNA
#   c. Psirc
#   ...
## ================ Load Global Configurations ================= ##
configfile: 'config.yaml'             # load config file
#
OUTPUTDIR = config["output_dir"]      # output directory
SAMPLEFILES = yaml.load(open(config['sample_list']), Loader=yaml.FullLoader)  # sample files
SAMPLES = sorted(SAMPLEFILES.keys())  # sample names
#
# Need for Prediction of lncRNA, novel mRNA, and circRNA
DNA = config["dna"]     # reference genome fasta
GTF = config["gtf"]     # reference gene annotation
#
# Needed for rRNA filtering
rRNA = config["rRNA"]   # rRNAs from RNACentral or Specise-Spacial rRNA
BOWTIE2_rRNA_INDEX = config["rRNA_bowtie2_index"]   # rRNA bowtie2 index
#
# Specifically Needed for LncRNA and Novel mRNA Identification
HISAT2_SPLICE_SITES = config["splice_sites"]        # load hisat2 splice site information
HISAT2_DNA_INDEX    = config["genome_hisat2_index"]   # load hisat2 index for genome maaping
#
BED = config["bed"]     # refernece bed files
CDS = config["cds"]     # cds fasta
CDNA = config["cdna"]   # cdna fasta
NCRNA = config["ncrna"] # ncRNA fasta
MIRNA = config["mirna"] # mature miRNA fasta
MRNA_GTF = config["mrna_gtf"]       # known mRNA gtf, need for FEELnc
LNCRNA_GTF = config["lncrna_gtf"]   # known lncRNA gtf, need for FEElnc
PFAMDB = os.path.dirname(config["pfamdb"]) # PFAM database, need for Protein Prediction
#
# Specifically Needed for CircRNA Identification and Quantification
BWA_DNA_INDEX   = config["genome_bwa_index"]    # bwa index for genome mapping, need for CIRI2
CIRI_QUANT_CFG  = config["ciri_quant_cfg"]      # yaml format config file for CIRIquant
BOWTIE2_DNA_INDEX = config["genome_bowtie2_index"]
#
# DEGs Analysis Needed
COMPARE_PAIRS = config['compare_pairs']         # comparation paris
SAMPLE_GROUPS = config['sample_groups']         # sample group information
#
## ============================================================ ##
#                    Start Steps for each Module                 #
#                        2021-08-19 09:30                        #
#                          Edit by Alipe                         #    
##################################################################
#                                                                #
#                                                                #
#                             _oo0oo_                            #
#                            o8888888o                           #
#                            88" . "88                           #
#                            (| -_- |)                           #
#                            0\  =  /0                           #
#                          ___/`---'\___                         #
#                        .' \\|     |# '.                        #
#                       / \\|||  :  |||# \                       #
#                      / _||||| -:- |||||- \                     #
#                     |   | \\\  -  #/ |   |                     #
#                     | \_|  ''\---/''  |_/ |                    #
#                     \  .-\__  '-'  ___/-. /                    #
#                   ___'. .'  /--.--\  `. .'___                  #
#                ."" '<  `.___\_<|>_/___.' >' "".                #
#               | | :  `- \`.;`\ _ /`;.`/ - ` : | |              #
#               \  \ `_.   \_ __\ /__ _/   .-` /  /              #
#           =====`-.____`.___ \_____/___.-`___.-'=====           #
#                             `=---='                            #
#                                                                #
#         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          #
#                                                                #
#                     佛祖保佑         永无BUG                    #
#                                                                #
##################################################################
#                                                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rule all:
    input:
  ## -------------- Part 01 Data Preprocessing ---------------- ##
    # Step 00: Prepare data
        expand( OUTPUTDIR + "./Preprocess/00.DataPrepare/{sample}.R1.fq.gz", sample=SAMPLES ),
        expand( OUTPUTDIR + "./Preprocess/00.DataPrepare/{sample}.R2.fq.gz", sample=SAMPLES ),
    # Step 01: Quality Control
        expand( OUTPUTDIR + "./Preprocess/01.FastqFilter/{sample}/{sample}.R1.fq.gz", sample=SAMPLES ),
        expand( OUTPUTDIR + "./Preprocess/01.FastqFilter/{sample}/{sample}.R2.fq.gz", sample=SAMPLES ),
        expand( OUTPUTDIR + "./Preprocess/01.FastqFilter/{sample}/{sample}.json", sample=SAMPLES ),
        expand( OUTPUTDIR + "./Preprocess/01.FastqFilter/{sample}/{sample}.html", sample=SAMPLES ),
    # Step 02: Filter rRNA
        expand( OUTPUTDIR + "./Preprocess/02.FilterrRNA/{sample}.bam", sample=SAMPLES ),
        expand( OUTPUTDIR + "./Preprocess/02.FilterrRNA/{sample}/un-conc-mate.1", sample=SAMPLES ),
        expand( OUTPUTDIR + "./Preprocess/02.FilterrRNA/{sample}/un-conc-mate.2", sample=SAMPLES ),
        expand( OUTPUTDIR + "./Preprocess/02.FilterrRNA/{sample}/rRNA_filter.ok", sample=SAMPLES ),
    # Step 03: Rename bowtie2 output un-conc-mate
        expand( OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R1.fq.gz", sample=SAMPLES ),
        expand( OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R2.fq.gz", sample=SAMPLES ),
    # Step 04: FastQC
        expand( OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/{sample}.R1_fastqc.zip", sample=SAMPLES ),
        expand( OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/{sample}.R1_fastqc.html", sample=SAMPLES ),
        expand( OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/{sample}.R2_fastqc.zip", sample=SAMPLES ),
        expand( OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/{sample}.R2_fastqc.html", sample=SAMPLES ),
        expand( OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/{sample}.fastqc.ok", sample=SAMPLES ),
        expand( OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/{sample}.fastqc.script", sample=SAMPLES ),
    # Step 05: Prepare-MultiQC-List
        OUTPUTDIR + "./Preprocess/05.MultiQC/FastQC_Out_List/fastqc_output_list.txt",
        OUTPUTDIR + "./Preprocess/05.MultiQC/FastQC_Out_List/make_fastqc_output_list.ok",
    # Step06： MultiQC
        OUTPUTDIR + "./Preprocess/05.MultiQC/Report/multiqc_report.html",
        OUTPUTDIR + "./Preprocess/05.MultiQC/Report/multiqc_report.ok",
  ##
  ## -------------- Part 02 Mapping and Assembly -------------- ##
    # Step 01: Align to reference genome
        expand( OUTPUTDIR + "./MappingAndAssembly/01.Hisat2Genome/{sample}.bam", sample=SAMPLES ),
        expand( OUTPUTDIR + "./MappingAndAssembly/01.Hisat2Genome/{sample}.summary.txt", sample=SAMPLES ),
        expand( OUTPUTDIR + "./MappingAndAssembly/01.Hisat2Genome/{sample}/un-conc-mate.1", sample=SAMPLES ),
        expand( OUTPUTDIR + "./MappingAndAssembly/01.Hisat2Genome/{sample}/un-conc-mate.2", sample=SAMPLES ),
    # Step 02: ReName un-conc-gz
        expand( OUTPUTDIR + "./MappingAndAssembly/02.Hisat2Rename/{sample}.R1.fq.gz", sample=SAMPLES ),
        expand( OUTPUTDIR + "./MappingAndAssembly/02.Hisat2Rename/{sample}.R2.fq.gz", sample=SAMPLES ),
    # Step 03: Stringtie Assembly
        expand( OUTPUTDIR + "./MappingAndAssembly/03.StringtieAssembly/{sample}.gtf", sample=SAMPLES ),
    # Step 04: Make gtf List
        OUTPUTDIR  + "MappingAndAssembly/04.MakeGtfMergeList/MergedList.txt",
    # Step 05: Merge transcript
        OUTPUTDIR + "./MappingAndAssembly/05.StringtieMerge/StringtieMerged.gtf",
    # Step 06: Compare to reference annotation
        OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.annotated.gtf",
        OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.stats",
        OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.tracking",
        OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.loci",
        OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.ok",
  ##
  ## -------------- Part 03 LncRNA Identification ------------- ##
    # Step 01: Fetch class code "i", "o", "x", and "u"
        OUTPUTDIR + "./LncRNA/01.CandidateLncRNAGtf/GffCompared.ioux.gtf",
    # Step 02: Fetch Candidate lncRNA fasta
        OUTPUTDIR + "./LncRNA/02.CandidatelncRNAFas/GffCompared.ioux.fa",
        OUTPUTDIR + "./LncRNA/02.CandidatelncRNAFas/GffCompared.ioux.pep.fa",
    # Step 03: lncRNA protein coding potential prodiction with CPC2
        OUTPUTDIR + "./LncRNA/03.CPC2_Predict/CPC2PredictOut.txt",
        OUTPUTDIR + "./LncRNA/03.CPC2_Predict/CPC2_Noncoding.txt",
        OUTPUTDIR + "./LncRNA/03.CPC2_Predict/CPC2_Noncoding.ok",
    # Step 04: LncRNA protein coding potential prediction with CNCI
        OUTPUTDIR + "./LncRNA/04.CNCI_Predict/CNCI_Predict/CNCI.index",
        OUTPUTDIR + "./LncRNA/04.CNCI_Predict/CNCI_Noncoding.txt",
        OUTPUTDIR + "./LncRNA/04.CNCI_Predict/CNCI_Noncoding.ok",
    # Step 05: LncRNA protein coding potential prediction with PfamScan
        OUTPUTDIR + "./LncRNA/05.Pfam_Predict/PfamPredictOut.txt",
        OUTPUTDIR + "./LncRNA/05.Pfam_Predict/Pfam_Coding.txt",
    # Step 06: LncRNA Identification by FEElnc, S1: filter
        OUTPUTDIR + "./LncRNA/06.FEELnc_filter/Candidate_lncRNA_flt.gtf",
    # Step 07: LncRNA Identification by FEElnc, S2: codpot predict
        OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/Candidate_lncRNA_codpot.lncRNA.gtf",
        OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/Candidate_lncRNA_codpot.mRNA.gtf",
        OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/Candidate_lncRNA_codpot.noORF.gtf",
        OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/Candidate_lncRNA_codpot_RF_learningData.txt",
        OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/Candidate_lncRNA_codpot_RF_statsLearn_CrossValidation.txt",
        OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/Candidate_lncRNA_codpot_RF.txt",
    # Step 08: LncRNA Identification by FEElnc, S3: classifier
        OUTPUTDIR + "./LncRNA/08.FEELnc_classifier/Candidate_lncRNA_classes.txt",
  ##
  ## -------------- Part 04 Novel mRNA Identification --------- ##
    # Step 01: Fetch Candidate Novel mRNA GTF
        OUTPUTDIR + "./Novel_mRNA/01.FetchCandidateNovelmRNAGtf/GffCompared.ju.gtf",
    # Step 02: Fetch Candidate Novel mRNA fasta
        OUTPUTDIR + "./Novel_mRNA/02.FetchCandidateNovelmRNAFas/GffCompared.ju.fa",
    # Step 03: Fetch Gene<\t>Transcript table
        OUTPUTDIR + "./Novel_mRNA/03.Gene2Tanscript/GffCompared_ju.g2t.txt",
    # Step 04: Fetch Candidate Novel mRNA ORF
        OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.cds",
        OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.rm_NN.cds",
        OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.pep",
        OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.ok",
    # Step 05: Novel mRNA protein coding potential prodict by CPC2
        OUTPUTDIR + "./Novel_mRNA/05.CodingPotential_CPC2/CPC2PredictOut.txt",
        OUTPUTDIR + "./Novel_mRNA/05.CodingPotential_CPC2/CPC2_Coding.txt",
        OUTPUTDIR + "./Novel_mRNA/05.CodingPotential_CPC2/CPC2_Coding.ok",
    # Step 06: Novel mRNA protein coding potential prodict by CNCI
        OUTPUTDIR + "./Novel_mRNA/06.CodingPotential_CNCI/CNCI_Predict/CNCI.index",
        OUTPUTDIR + "./Novel_mRNA/06.CodingPotential_CNCI/CNCI_Coding.txt",
        OUTPUTDIR + "./Novel_mRNA/06.CodingPotential_CNCI/CNCI_Coding.ok",
    # Step 07: Novel mRNA protein coding potential prodict by Pfam
        OUTPUTDIR + "./Novel_mRNA/07.CodingPotential_Pfam/PfamPredictOut.txt",
        OUTPUTDIR + "./Novel_mRNA/07.CodingPotential_Pfam/Pfam_Coding.txt",
    # Step 08: Random fetch fas for CPAT
        OUTPUTDIR + "./Novel_mRNA/08.RandomFetchFas/Traning_CDS.fa",
        OUTPUTDIR + "./Novel_mRNA/08.RandomFetchFas/Traning_Noc.fa",
    # Step 09: Novel mRNA protein coding potential prodict by CPAT-BuildHexamerTable
        OUTPUTDIR + "./Novel_mRNA/09.CPATBuildHexamerTable/Maize_Hexamer.tsv",
    # Step 10: Novel mRNA protein coding potential prodict by CPAT-Build Logit Model
        OUTPUTDIR + "./Novel_mRNA/10.CPATBuildLogitModel/CpatMaize.feature.xls",
        OUTPUTDIR + "./Novel_mRNA/10.CPATBuildLogitModel/CpatMaize.logit.RData",
        OUTPUTDIR + "./Novel_mRNA/10.CPATBuildLogitModel/CpatMaize.make_logitModel.r",
    # Step 11: Novel mRNA protein coding potential prodict by CPAT-Detect ORF
        OUTPUTDIR + "./Novel_mRNA/11.CPATtoDetectORF/CPAT_Predicet_Out.ORF_seqs.fa",
        OUTPUTDIR + "./Novel_mRNA/11.CPATtoDetectORF/CPAT_Predicet_Out.ORF_prob.tsv",
        OUTPUTDIR + "./Novel_mRNA/11.CPATtoDetectORF/CPAT_Predicet_Out.ORF_prob.best.tsv",
        OUTPUTDIR + "./Novel_mRNA/11.CPATtoDetectORF/CPAT_Predicet_Out.no_ORF.txt",
        OUTPUTDIR + "./Novel_mRNA/11.CPATtoDetectORF/CPAT_Predicet_Out.r",
  ##
  # -------------- Part 05 Expression Analysis --------------- ##
    # Step 01: Re Assembly by Stringtie
        expand( OUTPUTDIR + "./Expression_Analysis/01.ReStringtieAssemble/{sample}.gtf", sample = SAMPLES),
        expand( OUTPUTDIR + "./Expression_Analysis/01.ReStringtieAssemble/{sample}.coverage.cov", sample = SAMPLES),
        expand( OUTPUTDIR + "./Expression_Analysis/01.ReStringtieAssemble/{sample}.GeneAbund.txt", sample = SAMPLES),
    # Step 02: Get assembled gtf list
        OUTPUTDIR + "./Expression_Analysis/02.GetAssembledGtfList/AssembledGtfList.txt",
    # Step 03: Get gene and transcripts count matrix
        OUTPUTDIR + "./Expression_Analysis/03.GetCountAndTPMMatrix/gene_count_matrix.csv",
        OUTPUTDIR + "./Expression_Analysis/03.GetCountAndTPMMatrix/transcript_count_matrix.csv",
        OUTPUTDIR + "./Expression_Analysis/03.GetCountAndTPMMatrix/transcript_tpm_matrix.tsv",
    # Step 04: Different Expression Analysis by edgeR
        # OUTPUTDIR + "./Expression_Analysis/04.DGEbyEdgeR2/DGEbyEdgeR2.Result.ok",    
  #
  ## -------------- Part 06 CircRNA Identification ------------ ##
    # Step 01: Trim fastq to equal reads
        expand( OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.R1P.fq.gz", sample=SAMPLES),
        expand( OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.R2P.fq.gz", sample=SAMPLES),
        expand( OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.R1U.fq.gz", sample=SAMPLES),
        expand( OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.R2U.fq.gz", sample=SAMPLES),
        expand( OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.summary.txt", sample=SAMPLES),
        expand( OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.summary.ok", sample=SAMPLES),
    # Step 01: Map to Genome with BWA-MEM
        expand( OUTPUTDIR + "./CircRNA/01.BWA2Genome/{sample}.sam", sample=SAMPLES ),
    # Step 02: CircRNA identification with CIRI2
        expand( OUTPUTDIR + "./CircRNA/02.CIRI2_Prediction/{sample}.ciri", sample=SAMPLES ),
#     # Step 03: CircRNA quantitation with CIRIquant
#         expand( OUTPUTDIR + "./CircRNA/03.CircRNA_Quantitation/{sample}/{sample}.gtf", sample=SAMPLES ),
#     # Step 04: Identify circRNA by find_circ -- 1.mapping
#         expand( OUTPUTDIR + "./CircRNA/04.Bowtie2ToGenome/{sample}.sorted.bam", sample=SAMPLES),
#     # Step 05: Identify circRNA by find_circ -- 2.Fetch unmapped read with bowtie2
#         expand( OUTPUTDIR + "./CircRNA/05.UnmappedBam/{sample}.unmapped.bam", sample=SAMPLES),
#     # Step 06: Identify circRNA by find_circ -- 3.Convert bam to qfa
#         expand( OUTPUTDIR + "./CircRNA/06.Bam2Anchors/{sample}/unmapped_anchors.fq.gz", sample=SAMPLES),
#     # Step 07: Identify circRNA by find_circ -- 4.Find circRNA
#         expand( OUTPUTDIR + "./CircRNA/07.FindCircRNA/{sample}/spliced_reads.fa", sample = SAMPLES),
#         expand( OUTPUTDIR + "./CircRNA/07.FindCircRNA/{sample}/spliced_reads.bed", sample = SAMPLES),
#         expand( OUTPUTDIR + "./CircRNA/07.FindCircRNA/{sample}/stat.txt", sample = SAMPLES),
#     # Step 08: Identify circRNA by find_circ -- 5.Merge all samples bed
#         OUTPUTDIR + "./CircRNA/08.MergeAllSamplesBed/merged_spliced_reads.bed",
#         OUTPUTDIR + "./CircRNA/08.MergeAllSamplesBed/merged_stat.txt",
#     # Step 09: Identify circRNA by find_circ -- 6.Fetch Good circRNA
#         OUTPUTDIR + "./CircRNA/09.FinalCircRNA/circ_candidates.bed", 
#   ##
#   ## -------------- Part 07 Final Analysis Report  ------------ ##
#     # Report S01: Fastq Filter Result
#         OUTPUTDIR + "./Report/S01.FastqFilter.jsonlist.txt",
#         OUTPUTDIR + "./Report/S01.FastqFilter-State.tsv",
#     # Report S03: Read distribution by RseQC 3.0
#         expand( OUTPUTDIR + "./Report/S03.ReadDistribution/{sample}.state.tsv", sample=SAMPLES ),
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
##################################################################
## ================ Part 01 Data Preprocessing ================ ##
##################################################################
#
# Step 00: Prepare data
rule Part01_Preprocess_00_PrepareFastqFile:
    input:
        lambda wildcards:SAMPLEFILES[wildcards.sample]
    output:
       R1 =  OUTPUTDIR + "./Preprocess/00.DataPrepare/{sample}.R1.fq.gz",
       R2 =  OUTPUTDIR + "./Preprocess/00.DataPrepare/{sample}.R2.fq.gz"
    threads:
        1
    shell:
        """
        ln -sf {input[0]} {output.R1} && ln -sf {input[1]} {output.R2}
        """        
#
# Step 01: Quality Control
rule Part01_Preprocess_01_FastqFilter:
    input:
        R1 =  OUTPUTDIR + "./Preprocess/00.DataPrepare/{sample}.R1.fq.gz",
        R2 =  OUTPUTDIR + "./Preprocess/00.DataPrepare/{sample}.R2.fq.gz"
    output:
        R1   = OUTPUTDIR + "./Preprocess/01.FastqFilter/{sample}/{sample}.R1.fq.gz",
        R2   = OUTPUTDIR + "./Preprocess/01.FastqFilter/{sample}/{sample}.R2.fq.gz",
        json = OUTPUTDIR + "./Preprocess/01.FastqFilter/{sample}/{sample}.json",
        html = OUTPUTDIR + "./Preprocess/01.FastqFilter/{sample}/{sample}.html"
    message:
        "Begin to filter fastq!"
    log:
        OUTPUTDIR + "./AllLogs/Preprocess/01.FastqFilter/{sample}.FastqFilter.log"
    params:
        "--detect_adapter_for_pe --average_qual 15 --length_required 50"
    threads:
        8
    shell:
        """
        source activate qc_env && \
        fastp {params} -w {threads} -i {input.R1} -I {input.R2} -o {output.R1} \
        -O {output.R2} -j {output.json} -h {output.html} 2> {log}
        """
#
# Step 02: rRNA filter
rule Part01_Preprocess_02_FilterrRNA:
    input:
        R1 = OUTPUTDIR + "./Preprocess/01.FastqFilter/{sample}/{sample}.R1.fq.gz",
        R2 = OUTPUTDIR + "./Preprocess/01.FastqFilter/{sample}/{sample}.R2.fq.gz",
        idx = BOWTIE2_rRNA_INDEX
    output:
        bam  = OUTPUTDIR + "./Preprocess/02.FilterrRNA/{sample}.bam",
        stat = OUTPUTDIR + "./Preprocess/02.FilterrRNA/{sample}/rRNA_filter.ok",
        R1   = OUTPUTDIR + "./Preprocess/02.FilterrRNA/{sample}/un-conc-mate.1",
        R2   = OUTPUTDIR + "./Preprocess/02.FilterrRNA/{sample}/un-conc-mate.2"
    message:
        "Start to filter rRNA baseed on RNACentral plant rRNA database by bowtie2."
    log:
        OUTPUTDIR + "./AllLogs/Preprocess/02.FilterrRNA/{sample}.align.logs"
    threads:
        8
    params:
        opt = "-q --phred33 --sensitive --end-to-end --fr",
        out = OUTPUTDIR + "./Preprocess/02.FilterrRNA/{sample}"
    shell:
        """
        bowtie2 {params.opt} -p {threads} --un-conc-gz {params.out} -x {input.idx} \
        -1 {input.R1} -2 {input.R2} 2>{log} | samtools sort -n -O Bam -@ 4 \
        -m 5G -o {output.bam} && echo Success > {output.stat}
        """
#
# Step 03: Rename bowtie2 un-conc-gz
rule Part01_Preprocess_03_RenameBowtieOut:
    input:
        stat = OUTPUTDIR + "./Preprocess/02.FilterrRNA/{sample}/rRNA_filter.ok",
        R1   = OUTPUTDIR + "./Preprocess/02.FilterrRNA/{sample}/un-conc-mate.1",
        R2   = OUTPUTDIR + "./Preprocess/02.FilterrRNA/{sample}/un-conc-mate.2"
    output:
        R1 = OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R1.fq.gz",
        R2 = OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R2.fq.gz"
    threads:
        1
    shell:
        """
        ln -sf {input.R1} {output.R1} && ln -sf {input.R2} {output.R2}
        """
# 
# Step 04: FastQC
rule Part01_Preprocess_04_FastQC:
    input: 
        R1 = OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R1.fq.gz",
        R2 = OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R2.fq.gz",
    output:
        R1 = OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/{sample}.R1_fastqc.zip",
        H1 = OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/{sample}.R1_fastqc.html",
        R2 = OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/{sample}.R2_fastqc.zip",
        H2 = OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/{sample}.R2_fastqc.html",
        ok = OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/{sample}.fastqc.ok",
        src = OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/{sample}.fastqc.script"
    threads:
        4
    params:
        dir = OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/"
    run:
        import os
        import subprocess

        if not os.path.exists(params.dir):
            os.makedirs(params.dir)

        cmd = """
        source activate qc_env && \
        fastqc -Xmx5g -Xms1g -o {o} -t {t} -f fastq {r1} {r2}
        """.format(o=params.dir, t=threads, r1=input.R1, r2=input.R2)

        with open(output.src, "w") as f:
            print(cmd, file=f)

        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.H1) and os.path.getsize(output.H2):
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)
# 
# Step 05: MultiQC-Prepare
rule Part01_Preprocess_05_Make_MultiQC_List:
    input:
        expand(OUTPUTDIR + "./Preprocess/04.FastQC/{sample}/{sample}.fastqc.ok", sample = SAMPLES)
    output:
        lst = OUTPUTDIR + "./Preprocess/05.MultiQC/FastQC_Out_List/fastqc_output_list.txt",
        ok  = OUTPUTDIR + "./Preprocess/05.MultiQC/FastQC_Out_List/make_fastqc_output_list.ok",
    threads:
        1
    params:
        dir = OUTPUTDIR + "./Preprocess/05.MultiQC/FastQC_Out_List/"
    run:
        import os
        import subprocess

        if not os.path.exists(params.dir):
            os.makedirs(params.dir)

        with open(output.lst, 'w') as f:
            for qc_ok in input:
                qc_dir = os.path.dirname(qc_ok)
                print(qc_dir, file=f)

        if os.path.getsize(output.lst) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)  
# 
# # Step 06: MultiQC
rule Part01_Preprocess_06_MultiQC:
    input:
        lst = OUTPUTDIR + "./Preprocess/05.MultiQC/FastQC_Out_List/fastqc_output_list.txt",
    output:
        rpt = OUTPUTDIR + "./Preprocess/05.MultiQC/Report/multiqc_report.html",
        ok = OUTPUTDIR + "./Preprocess/05.MultiQC/Report/multiqc_report.ok",
    threads:
        1
    params:
        dir = OUTPUTDIR + "./Preprocess/05.MultiQC/Report/"
    run:
        import os
        import subprocess

        if not os.path.exists(params.dir):
            os.makedirs(params.dir)

        # output_dir = os.path.dirname(output.rpt)
        cmd = """
        source activate qc_env && \
        multiqc --quiet --file-list {lst} -n {rpt}
        """.format(lst=input.lst, rpt=output.rpt)

        print(cmd)

        subprocess.call(cmd, shell=True)
        if os.path.getsize(output.rpt) > 0:
            subprocess.call("echo SUCCESS > {ok}".format(ok=output.ok), shell=True)
# 
##################################################################
## ================ Part 02 Mapping and Assembly ============== ##
##################################################################
#
# Step 01: Align to reference genome
rule Part02_MappingAndAssembly_01_Hisat2Genome:
    input:
        R1 = OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R1.fq.gz",
        R2 = OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R2.fq.gz",
        SS = HISAT2_SPLICE_SITES, # if hisat2-build with ss and exon, skip this step.
        GTF = GTF,
        IDX = HISAT2_DNA_INDEX
    output:
        BAM = OUTPUTDIR + "./MappingAndAssembly/01.Hisat2Genome/{sample}.bam",
        SUM = OUTPUTDIR + "./MappingAndAssembly/01.Hisat2Genome/{sample}.summary.txt",
        R1  = OUTPUTDIR + "./MappingAndAssembly/01.Hisat2Genome/{sample}/un-conc-mate.1",
        R2  = OUTPUTDIR + "./MappingAndAssembly/01.Hisat2Genome/{sample}/un-conc-mate.2"
    message:
        "Start to map  genome with hisat2."
    log:
        OUTPUTDIR + "./AllLogs/MappingAndAssembly/01.Hisat2Genome/{sample}.align.logs"
    threads:
        8
    params:
        OPT = "-q --phred33 --min-intronlen 20 --max-intronlen 150000 --dta",
        DIR = OUTPUTDIR + "./MappingAndAssembly/01.Hisat2Genome/{sample}"
    shell:
        """
        hisat2 {params.OPT} -p {threads} --summary-file {output.SUM} -x {input.IDX} \
            --known-splicesite-infile {input.SS} --un-conc-gz {params.DIR} \
            -1 {input.R1} -2 {input.R2} 2>{log} | samtools sort -O Bam \
            -@ {threads} -m 5G -o {output.BAM}
        """
#
# Step 02: ReName un-conc-gz
rule Part02_MappingAndAssembly_02_ReNameHisat2Out:
    input:
        R1 = OUTPUTDIR + "./MappingAndAssembly/01.Hisat2Genome/{sample}/un-conc-mate.1",
        R2 = OUTPUTDIR + "./MappingAndAssembly/01.Hisat2Genome/{sample}/un-conc-mate.2",
    output:
        R1 = OUTPUTDIR + "./MappingAndAssembly/02.Hisat2Rename/{sample}.R1.fq.gz",
        R2 = OUTPUTDIR + "./MappingAndAssembly/02.Hisat2Rename/{sample}.R2.fq.gz"
    threads:
        1
    shell:
        """
        ln -sf {input.R1} {output.R1} && ln -sf {input.R2} {output.R2}
        """
#
# Step 03: Stringtie Assembly
rule Part02_MappingAndAssembly_04_StringtieAssembly:
    input:
        bam = OUTPUTDIR + "./MappingAndAssembly/01.Hisat2Genome/{sample}.bam",
        gtf = GTF,
    output:
        gtf = OUTPUTDIR + "./MappingAndAssembly/03.StringtieAssembly/{sample}.gtf"
    message:
        "Start to assembly use stringtie."
    log:
        OUTPUTDIR + "./AllLogs/MappingAndAssembly/03.StringtieAssembly/{sample}.assembly.logs"
    threads:
        8
    params:
        "-m 200 -l STRG -a 10 --conservative -g 50 -u"
    shell:
        """
        stringtie {input.bam} -G {input.gtf} -o {output.gtf} -p {threads} {params}
        """
#
# Step 04: Make gtf List
rule Part02_MappingAndAssembly_05_MakeMergeList:
    input:
        expand(OUTPUTDIR + "./MappingAndAssembly/03.StringtieAssembly/{sample}.gtf", sample = SAMPLES)
    output:
        lst = OUTPUTDIR  + "MappingAndAssembly/04.MakeGtfMergeList/MergedList.txt",
        ok = OUTPUTDIR  + "MappingAndAssembly/04.MakeGtfMergeList/MergedList.ok",
    threads:
        1
    params:
        dir = OUTPUTDIR  + "MappingAndAssembly/04.MakeGtfMergeList/"
    run:
        import os
        import subprocess

        if not os.path.exists(params.dir):
            os.makedirs(params.dir)

        with open(output.lst, 'w') as f:
            for gtf in input:
                print(gtf, file=f)

        if os.path.getsize(output.lst) > 0:
            subprocess.call("echo Success > {ok}".format(ok=output.ok), shell=True)        
#
# Step 05: Merge transcript
rule Part02_MappingAndAssembly_06_StringtieMerge:
    input:
        ok = OUTPUTDIR  + "MappingAndAssembly/04.MakeGtfMergeList/MergedList.ok",
        lst = OUTPUTDIR + "./MappingAndAssembly/04.MakeGtfMergeList/MergedList.txt",
        gtf = GTF
    output:
        gtf = OUTPUTDIR + "./MappingAndAssembly/05.StringtieMerge/StringtieMerged.gtf"
    log:
        OUTPUTDIR + "./AllLogs/MappingAndAssembly/05.StringtieMerge/StringtieMerge.logs"
    threads:
        8
    params:
        "-m 200 -c 3"
    shell:
        """
        stringtie --merge {params} -p {threads} -G {input.gtf} -o {output.gtf} {input.lst} 2> {log}
        """
#
# Step 06: Compare to reference annotation
rule Part02_MappingAndAssembly_07_Compare2Ref:
    input:
        RefGtf = GTF,
        ComGtf = OUTPUTDIR + "./MappingAndAssembly/05.StringtieMerge/StringtieMerged.gtf"
    output:
        gtf   = OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.annotated.gtf",
        stats = OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.stats",
        track = OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.tracking",
        loci  = OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.loci",
        ojbk  = OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.ok"
    log:
        OUTPUTDIR + "./AllLogs/MappingAndAssembly/06.Compare2Ref/GffCompared.logs"
    threads:
        1
    params:
        opt = "-T -R",
        pfx = OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared",
    shell:
        """
        gffcompare {params.opt} -r {input.RefGtf} -o {params.pfx} {input.ComGtf} 2> {log} && \
        echo SUCCESS > {output.ojbk}
        """
#
##################################################################
## ================ Part 03 LncRNA Identification ============= ##
##################################################################
#
# Step 01: Fetch class code "i", "o", "x", and "u"
rule Part03_LncRNA_Identification_01_FetchCandidateLncRNAGtf:
    input:
        gtf  = OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.annotated.gtf",
        ojbk = OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.ok"
    output:
        gtf = OUTPUTDIR + "./LncRNA/01.CandidateLncRNAGtf/GffCompared.ioux.gtf"
    threads:
        1
    params:
        "iuxo"
    shell:
        """
        perl Scripts/Fetch.class_code_gtf.pl {input} {params} > {output}
        """
#
# Step 02: Fetch Candidate lncRNA fasta
rule Part03_LncRNA_Identification_02_FetchCandidatelncRNAFas:
    input:
        dna = DNA,
        gtf = OUTPUTDIR + "./LncRNA/01.CandidateLncRNAGtf/GffCompared.ioux.gtf"
    output:
        fna = OUTPUTDIR + "./LncRNA/02.CandidatelncRNAFas/GffCompared.ioux.fa",
        tmp = OUTPUTDIR + "./LncRNA/02.CandidatelncRNAFas/GffCompared.ioux.tmp",
        pep = OUTPUTDIR + "./LncRNA/02.CandidatelncRNAFas/GffCompared.ioux.pep.fa"
    log:
        OUTPUTDIR + "./AllLogs/LncRNA/02.CandidatelncRNAFas/Candaidate_lncRNA_fa_fetch.log"
    threads:
        1
    shell:
        """
        gffread -w {output.fna} -g {input.dna} {input.gtf} && transeq {output.fna} \
        {output.tmp} && sed 's/*//g' {output.tmp} > {output.pep} 2> {log}
        """
#
# Step 03: lncRNA protein coding potential prodiction with CPC2
rule Part03_LncRNA_Identification_03_CPC2_Predict:
    input:
        fas = OUTPUTDIR + "./LncRNA/02.CandidatelncRNAFas/GffCompared.ioux.fa",
    output:
        cpc = protected(OUTPUTDIR + "./LncRNA/03.CPC2_Predict/CPC2PredictOut.txt"),
        noc = protected(OUTPUTDIR + "./LncRNA/03.CPC2_Predict/CPC2_Noncoding.txt"),
        ok = OUTPUTDIR + "./LncRNA/03.CPC2_Predict/CPC2_Noncoding.ok",
    log:
        OUTPUTDIR + "./AllLogs/LncRNA/03.CPC2_Predict/CPC2Predict.log"
    threads:
        1
    params:
        dir = OUTPUTDIR + "./LncRNA/03.CPC2_Predict/",
        pfx = OUTPUTDIR + "./LncRNA/03.CPC2_Predict/CPC2PredictOut",
    run:
        import os
        import subprocess
        if not os.path.exists(params.dir):
            os.makedirs(params.dir)

        cmd = """
        source activate cpc2_py3_env && \
        CPC2.py -i {fas} -o {pfx} 2> {lo} && grep -w "noncoding" {cpc} | cut -f1 > {noc}
        """.format(fas=input.fas, lo=log, pfx=params.pfx, cpc=output.cpc, noc=output.noc)
        
        print(cmd)
        subprocess.call(cmd, shell = True)

        if os.path.getsize(output.noc) > 0:
            subprocess.call("echo SUCCESS > {ok}".format(ok=output.ok), shell = True)
#
# Step 04: LncRNA protein coding potential prediction with CNCI
rule Part03_LncRNA_Identification_04_CNCI_Predict:
    input:
        fas = OUTPUTDIR + "./LncRNA/02.CandidatelncRNAFas/GffCompared.ioux.fa",
    output:
        cnci = protected(OUTPUTDIR + "./LncRNA/04.CNCI_Predict/CNCI_Predict/CNCI.index"),
        noc  = protected(OUTPUTDIR + "./LncRNA/04.CNCI_Predict/CNCI_Noncoding.txt"),
        ok   = OUTPUTDIR + "./LncRNA/04.CNCI_Predict/CNCI_Noncoding.ok",
    log:
        OUTPUTDIR + "./AllLogs/LncRNA/04.CNCI_Predict/CNCI_Predict.log"
    threads:
        10
    params:
        opt = "-m pl",
        dir = OUTPUTDIR + "./LncRNA/04.CNCI_Predict/CNCI_Predict",
    run:
        import os
        import subprocess

        if not os.path.exists(params.dir):
            os.makedirs(params.dir)
       
        cmd = """source activate cnci_py2_env && CNCI.py -f {f} -o {out} -p {p} \
        {opt} 2> {lo} && grep -w "noncoding" {cnci} | cut -f1 > {noc}
        """.format(f=input.fas,out=params.dir, p=threads, opt=params.opt, lo=log, cnci=output.cnci, noc=output.noc)

        print(cmd)
        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.noc) > 0:
            subprocess.call("echo SUCESS > {ok}".format(ok=output.ok), shell=True)
#
# Step 05: LncRNA protein coding potential prediction with PfamScan
rule Part03_LncRNA_Identification_05_Pfam_Predict:
    input:
        OUTPUTDIR + "./LncRNA/02.CandidatelncRNAFas/GffCompared.ioux.pep.fa",
    output:
        res = protected(OUTPUTDIR + "./LncRNA/05.Pfam_Predict/PfamPredictOut.txt"),
        cod = protected(OUTPUTDIR + "./LncRNA/05.Pfam_Predict/Pfam_Coding.txt"),
    log:
        OUTPUTDIR + "./AllLogs/LncRNA/05.Pfam_Predict/PfamPredict.log"
    threads:
        20
    params:
        db = PFAMDB
    shell:
        """
        source activate pfam_scan_env && \
        pfam_scan.pl -cpu {threads} -fasta {input} -dir {params.db} -outfile {output.res} 2> {log} && \
        grep -v "#" {output.res} |cut -d' ' -f1 > {output.cod}
        """
#
# Step 06: LncRNA Identification by FEElnc, S1: filter
rule Part03_LncRNA_Identification_06_FEELnc_filter:
    input:
        RefGtf = GTF,
        ComGtf = OUTPUTDIR + "./LncRNA/01.CandidateLncRNAGtf/GffCompared.ioux.gtf"
    output:
        flt = OUTPUTDIR + "./LncRNA/06.FEELnc_filter/Candidate_lncRNA_flt.gtf"
    log:
        OUTPUTDIR + "./AllLogs/LncRNA/06.FEELnc_filter/FEELncPrediction.flt.log",
    threads:
        20
    params:
        "-s 200 -b transcript_biotype=protein_coding --monoex=0 --biex=25 -l FALSE"        
    shell:
        """
        source activate feelnc_env && \
        FEELnc_filter.pl {params} -p {threads} -i {input.ComGtf} -a {input.RefGtf} -o {log} > {output.flt}
        """
#
# Step 07: LncRNA Identification by FEElnc, S2: codpot predict
rule Part03_LncRNA_Identification_07_FEELnc_codpot:
    input:
        RefFas = DNA,
        RefGtf = GTF,
        mRNAGtf   = MRNA_GTF,
        lncRNAGtf = LNCRNA_GTF,
        flt = OUTPUTDIR + "./LncRNA/06.FEELnc_filter/Candidate_lncRNA_flt.gtf"
    output:
        cod = OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/Candidate_lncRNA_codpot.lncRNA.gtf",
        rna = OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/Candidate_lncRNA_codpot.mRNA.gtf",
        orf = OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/Candidate_lncRNA_codpot.noORF.gtf",
        rld = OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/Candidate_lncRNA_codpot_RF_learningData.txt",
        rsc = OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/Candidate_lncRNA_codpot_RF_statsLearn_CrossValidation.txt",
        rft = OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/Candidate_lncRNA_codpot_RF.txt"
    log:
        OUTPUTDIR + "./AllLogs/LncRNA/07.FEELnc_codpot/FEELncPrediction.codpot.log"
    threads:
        20
    params:
        cod_opt1 = "-b transcript_biotype=protein_coding --mode=shuffle --sizeinter=0.75",
        cod_opt2 = "--learnorftype=3 --testorftype=3 --ntree 500 --seed=1234",
        cod_dir = OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/",
        cod_name = "Candidate_lncRNA_codpot"
    shell:
        """
        source activate feelnc_env && \
        FEELnc_codpot.pl {params.cod_opt1} {params.cod_opt2} -p {threads} -i {input.flt} -a {input.RefGtf} \
            -g {input.RefFas} -l {input.lncRNAGtf} --outdir={params.cod_dir} --outname={params.cod_name} 2> {log}
        """
#
# Step 08: LncRNA Identification by FEElnc, S3: classifier
rule Part03_LncRNA_Identification_08_FEELnc_classifier:
    input:
        RefGtf = GTF,
        CodGtf = OUTPUTDIR + "./LncRNA/07.FEELnc_codpot/Candidate_lncRNA_codpot.lncRNA.gtf"
    output:        
        OUTPUTDIR + "./LncRNA/08.FEELnc_classifier/Candidate_lncRNA_classes.txt"
    log:
        OUTPUTDIR + "./AllLogs/LncRNA/08.FEELnc_classifier/FEELncPrediction.class.log"
    threads:
        1
    shell:
        """
        source activate feelnc_env && \
        FEELnc_classifier.pl -i {input.CodGtf} -a {input.RefGtf} > {output} -l {log}
        """
#
# Step 09: Fetch Final identified lncRNAs
rule Part03_LncRNA_Identification_09_FinalCantidatelncRNAs:
    input:
        cpc2 = OUTPUTDIR + "./LncRNA/03.CPC2_Predict/CPC2_Noncoding.txt",
        cnci = OUTPUTDIR + "./LncRNA/04.CNCI_Predict/CNCI_Noncoding.txt",
        pfam = OUTPUTDIR + "./LncRNA/05.Pfam_Predict/Pfam_Coding.txt",
        feel = OUTPUTDIR + "./LncRNA/08.FEELnc_classifier/Feelnc_Noncoding.txt"
    output:
    threads:
        1
    run:
        import os
        import subprocess
        cmd = """
        cat {c} {n} {p} {f} | sort |uniq -c |awk "$1==4 {print $2}" > {o}
        """.format(c=input.cpc2, n=input.cnci, p=input.pfam, f=input.feel, o=output.lst)

        print(cmd)
        subprocess.call(cmd, shell=True)
#
##################################################################
## ================ Part 04 Novel mRNA Identification ========= ##
##################################################################
#
# Step 01: Fetch Candidate Novel mRNA GTF
rule Part04_NovelmRNA_Identification_01_FetchCandidateNovelmRNAGtf:
    input:
        gtf  = OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.annotated.gtf",
        ojbk = OUTPUTDIR + "./MappingAndAssembly/06.Compare2Ref/GffCompared.ok"
    output:
        OUTPUTDIR + "./Novel_mRNA/01.FetchCandidateNovelmRNAGtf/GffCompared.ju.gtf"
    log:
        OUTPUTDIR + "./AllLogs/Novel_mRNA/01.FetchCandidateNovelmRNAGtf/FetchCandidateNovelmRNAGtf.log"
    threads:
        1
    params:
        "ju"
    shell:
        """
        perl Scripts/Fetch.class_code_gtf.pl {input.gtf} {params} > {output} 2> {log}
        """      
#
# Step 02: Fetch Candidate Novel mRNA fasta
rule Part04_NovelmRNA_Identification_02_FetchCandidateNovelmRNAFas:
    input:
        dna = DNA,
        gtf = OUTPUTDIR + "./Novel_mRNA/01.FetchCandidateNovelmRNAGtf/GffCompared.ju.gtf"
    output:
        OUTPUTDIR + "./Novel_mRNA/02.FetchCandidateNovelmRNAFas/GffCompared.ju.fa"
    log:
        OUTPUTDIR + "./AllLogs/Novel_mRNA/02.FetchCandidateNovelmRNAFas/FetchCandidateNovelmRNAFas.log"
    threads:
        1
    shell:
        """
        gffread -w {output} -g {input.dna} {input.gtf} 2> {log}
        """
#
# Step 03: Fetch Gene<\t>Transcript table
rule Part04_NovelmRNA_Identification_03_FetchGene2Tanscript:
    input:
        OUTPUTDIR + "./Novel_mRNA/01.FetchCandidateNovelmRNAGtf/GffCompared.ju.gtf"
    output:
        OUTPUTDIR + "./Novel_mRNA/03.Gene2Tanscript/GffCompared_ju.g2t.txt"
    threads:
        1    
    shell:
        """
        perl Scripts/Fetch.gene2transcript.pl {input} > {output}
        """
#
# Step 04: Fetch Candidate Novel mRNA ORF
rule Part04_NovelmRNA_Identification_04_TransDecoderLongOrfs:
    input:
        fna = OUTPUTDIR + "./Novel_mRNA/02.FetchCandidateNovelmRNAFas/GffCompared.ju.fa",
        g2t = OUTPUTDIR + "./Novel_mRNA/03.Gene2Tanscript/GffCompared_ju.g2t.txt"
    output:
        cds = protected(OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.cds"),
        fds = protected(OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.rm_NN.cds"),
        pep = protected(OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.pep"),
        gff = protected(OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.gff3"),
        ok = OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.ok",
    log:
        OUTPUTDIR + "./AllLogs/Novel_mRNA/04.TransDecoderLongOrfs/TransdecoderORF_Fetch.log"
    threads:
        1
    params:
        opt = "-m 100 -S",
        dir = OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/"
    run:
        import os,re
        import subprocess
        from Bio import SeqIO

        cmd = """
        TransDecoder.LongOrfs {opt} -t {fna} --gene_trans_map {g2t} --output_dir {dir} 2> {log}
        """.format(opt=params.opt, fna=input.fna, g2t=input.g2t, dir=params.dir, log=log)

        print(cmd)
        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.cds) > 0:
            records = SeqIO.parse(input.cds, "fasta")
            filtered = (rec for rec in records if not re.search('N', str(rec.seq), re.IGNORECASE))
            # SeqIO.write(filtered, output.fds, 'fasta')
            # same to cut -d' ' -f1 
            with open(output.fds, "w") as f:
                for rec in filtered:
                    out = ">" + str(rec.id) + "\n" + str(rec.seq)
                    print(out, file=f)

        if os.path.getsize(output.fds) > 0:
            subprocess.call("echo SUCCESS >{ok}".format(ok=output.ok), shell=True)
#     
# Step 05: Novel mRNA protein coding potential prodict by CPC2
rule Part04_NovelmRNA_Identification_05_CodingPotential_CPC2:
    input:
        cds = OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.rm_NN.cds",
        ok  = OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.ok",
    output:
        cpc = protected(OUTPUTDIR + "./Novel_mRNA/05.CodingPotential_CPC2/CPC2PredictOut.txt"),
        cod = protected(OUTPUTDIR + "./Novel_mRNA/05.CodingPotential_CPC2/CPC2_Coding.txt"),
        ok = OUTPUTDIR + "./Novel_mRNA/05.CodingPotential_CPC2/CPC2_Coding.ok"
    log:
        OUTPUTDIR + "./AllLogs/Novel_mRNA/05.CodingPotential_CPC2/CPC2Prediction.log"
    threads:
        1
    params:
        dir = OUTPUTDIR + "./Novel_mRNA/05.CodingPotential_CPC2/",
        pfx = OUTPUTDIR + "./Novel_mRNA/05.CodingPotential_CPC2/CPC2PredictOut"
    run:
        import os
        import subprocess
        if not os.path.exists(params.dir):
            os.makedirs(params.dir)
        
        cmd = """
        source activate cpc2_py3_env && \
        CPC2.py -i {cds} -o {pfx} 2> {lo} && grep -w "coding" {cpc} | cut -f1 > {cod}
        """.format(cds=input.cds, pfx=params.pfx, lo=log, cpc=output.cpc, cod=output.cod)

        print(cmd)
        subprocess.call(cmd, shell = True)

        sz = os.path.getsize(output.cod)
        if sz > 0:
            subprocess.call("echo SUCCESS > {ok}".format(ok=output.ok), shell = True)
#
# Step 06: Novel mRNA protein coding potential prodict by CNCI
rule Part04_NovelmRNA_Identification_06_CodingPotential_CNCI:
    input:
        cds = OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.rm_NN.cds",
        ok  = OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.ok",
    output:
        cnci = protected(OUTPUTDIR + "./Novel_mRNA/06.CodingPotential_CNCI/CNCI_Predict/CNCI.index"),
        cod = protected(OUTPUTDIR + "./Novel_mRNA/06.CodingPotential_CNCI/CNCI_Coding.txt"),
        ok = OUTPUTDIR + "./Novel_mRNA/06.CodingPotential_CNCI/CNCI_Coding.ok",
    log:
        OUTPUTDIR + "./AllLogs/Novel_mRNA/06.CodingPotential_CNCI/CNCI_Prediction.log"
    threads:
        10
    params:
        opt = "-m pl",
        dir = OUTPUTDIR + "./Novel_mRNA/06.CodingPotential_CNCI/CNCI_Predict"
    run:
        import os
        import subprocess

        if not os.path.exists(params.dir):
            os.makedirs(params.dir)

        cmd = """source activate cnci_py2_env && CNCI.py -f {cds} -o {dir} \
        -p {th} {opt} 2> {log} && grep -w "coding" {cnci} | cut -f1 > {cod}
        """.format(cds=input.cds, dir=params.dir, th=threads, opt=params.opt, 
        log=log, cnci=output.cnci, cod=output.cod)

        print(cmd)
        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.cod) > 0:
            subprocess.call("echo SUCCESS > {ok}".format(ok=output.ok), shell=True)
#
# Step 07: Novel mRNA protein coding potential prodict by Pfam
rule Part04_NovelmRNA_Identification_07_CodingPotential_Pfam:
    input:
        OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.pep"
    output:
        res = OUTPUTDIR + "./Novel_mRNA/07.CodingPotential_Pfam/PfamPredictOut.txt",
        cod = OUTPUTDIR + "./Novel_mRNA/07.CodingPotential_Pfam/Pfam_Coding.txt",
    log:
        OUTPUTDIR + "./AllLogs/Novel_mRNA/07.CodingPotential_Pfam/PfamPredictOut.log"
    threads:
        20
    params:
        db = PFAMDB
    shell:
        """
        source activate pfam_scan_env && \
        pfam_scan.pl -cpu {threads} -fasta {input} -dir {params.db} -outfile {output.res} 2> {log} && \
        grep -v "#" {output.res} |cut -d' ' -f1 > {output.cod}
        """
#
# Step 08: Random fetch fas for CPAT
rule Part04_NovelmRNA_Identification_08_RandomFetchFas:
    input:
        cds = CDS,
        noc = NCRNA,
    output:
        cds = OUTPUTDIR + "./Novel_mRNA/08.RandomFetchFas/Traning_CDS.fa",
        noc = OUTPUTDIR + "./Novel_mRNA/08.RandomFetchFas/Traning_Noc.fa",
    threads:
        1
    shell:
        """
        perl Scripts/RandomFetchFas.pl {input.cds} > {output.cds} && \
        perl Scripts/RandomFetchFas.pl {input.noc} > {output.noc}
        """    
#
# Step 09: Novel mRNA protein coding potential prodict by CPAT-BuildHexamerTable
rule Part04_NovelmRNA_Identification_09_CPATBuildHexamerTable:
    input:
        cod = OUTPUTDIR + "./Novel_mRNA/08.RandomFetchFas/Traning_CDS.fa",
        noc = OUTPUTDIR + "./Novel_mRNA/08.RandomFetchFas/Traning_Noc.fa",
    output:
        tsv = OUTPUTDIR + "./Novel_mRNA/09.CPATBuildHexamerTable/Maize_Hexamer.tsv"
    log:
        OUTPUTDIR + "./AllLogs/Novel_mRNA/09.CPATBuildHexamerTable/CPATBuildHexamerTable.log"
    threads:
        1
    run:
        import os
        import subprocess

        if not os.path.exists(os.path.dirname(output.tsv)):
            os.makedirs(os.path.dirname(output.tsv))

        cmd = """
        source activate cpat_py3_env && \
        make_hexamer_tab.py -c {cod} -n {noc} > {tsv} 2> {log}
        """.format(cod=input.cod, noc=input.noc, tsv=output.tsv, log=log)

        print(cmd)
        subprocess.call(cmd, shell=True)
#
# Step 10: Novel mRNA protein coding potential prodict by CPAT-Build Logit Model
rule Part04_NovelmRNA_Identification_10_CPATBuildLogitModel:
    input:
        cod = OUTPUTDIR + "./Novel_mRNA/08.RandomFetchFas/Traning_CDS.fa",
        noc = OUTPUTDIR + "./Novel_mRNA/08.RandomFetchFas/Traning_Noc.fa",
        hex = OUTPUTDIR + "./Novel_mRNA/09.CPATBuildHexamerTable/Maize_Hexamer.tsv",
    output:
        feature = OUTPUTDIR + "./Novel_mRNA/10.CPATBuildLogitModel/CpatMaize.feature.xls",
        rdata = OUTPUTDIR + "./Novel_mRNA/10.CPATBuildLogitModel/CpatMaize.logit.RData",
        rscipt = OUTPUTDIR + "./Novel_mRNA/10.CPATBuildLogitModel/CpatMaize.make_logitModel.r"
    log:
        OUTPUTDIR + "./AllLogs/Novel_mRNA/10.CPATBuildLogitModel/CPATBuildLogitModel.log"
    params:
        pfx = OUTPUTDIR + "./Novel_mRNA/10.CPATBuildLogitModel/CpatMaize"
    run:
        import os
        import subprocess

        if not os.path.exists(os.path.dirname(output.feature)):
            os.makedirs(os.path.dirname(output.feature))

        cmd = """
        source activate cpat_py3_env && \
        make_logitModel.py -x {hex} -c {cod} -n {noc} -o {pfx} --log-file={log}
        """.format(hex=input.hex, cod=input.cod, noc=input.noc, pfx=params.pfx, log=log)

        print(cmd)
        subprocess.call(cmd, shell=True)
#
# Step 11: Novel mRNA protein coding potential prodict by CPAT -- Detect ORF
rule Part04_NovelmRNA_Identification_11_CPATtoDetectORF:
    input:
        fas = OUTPUTDIR + "./Novel_mRNA/02.FetchCandidateNovelmRNAFas/GffCompared.ju.fa",
        # OUTPUTDIR + "./Novel_mRNA/04.TransDecoderLongOrfs/longest_orfs.rm_NN.cds",
        hex = OUTPUTDIR + "./Novel_mRNA/09.CPATBuildHexamerTable/Maize_Hexamer.tsv",
        rdata = OUTPUTDIR + "./Novel_mRNA/10.CPATBuildLogitModel/CpatMaize.logit.RData",
    output:
        seqs = protected(OUTPUTDIR + "./Novel_mRNA/11.CPATtoDetectORF/CPAT_Predicet_Out.ORF_seqs.fa"),
        peob = protected(OUTPUTDIR + "./Novel_mRNA/11.CPATtoDetectORF/CPAT_Predicet_Out.ORF_prob.tsv"),
        best = protected(OUTPUTDIR + "./Novel_mRNA/11.CPATtoDetectORF/CPAT_Predicet_Out.ORF_prob.best.tsv"),
        norf = protected(OUTPUTDIR + "./Novel_mRNA/11.CPATtoDetectORF/CPAT_Predicet_Out.no_ORF.txt"),
        rscript = protected(OUTPUTDIR + "./Novel_mRNA/11.CPATtoDetectORF/CPAT_Predicet_Out.r")
    log:
        OUTPUTDIR + "./AllLogs/Novel_mRNA/11.CPATtoDetectORF/CPATtoDetectORF.log"
    params:
        opt = "--top-orf=5 --min-orf=75 --width=100 --best-orf=p",
        pfx = OUTPUTDIR + "./Novel_mRNA/11.CPATtoDetectORF/CPAT_Predicet_Out"
    threads:
        1
    run:
        import os
        import subprocess

        cmd = """
        source activate cpat_py3_env && cpat.py -g {fas} -o {pfx} -x {hex} -d {r} --log-file={log}
        """.format(hex=input.hex, r=input.rdata, fas=input.fas, pfx=params.pfx, log=log)

        print(cmd)
        subprocess.call(cmd, shell=True)
#
# Step 12: Get Final Novel mRNAs
# rule Part04_NovelmRNA_Identification_12_FinalCandidatemRNAs:
#     input:
#         cpc2 = OUTPUTDIR + "./Novel_mRNA/05.CodingPotential_CPC2/CPC2_Coding.txt",
#         cnci = OUTPUTDIR + "./Novel_mRNA/06.CodingPotential_CNCI/CNCI_Coding.txt",
#         pfam = OUTPUTDIR + "./Novel_mRNA/07.CodingPotential_Pfam/Pfam_Coding.txt",
#         # cpat = ,
#     output:
#         lst = protected(OUTPUTDIR + "./Novel_mRNA/12.FinalCandidatemRNAs/Coding_mRNAs_list.txt"),
#     threads:
#         1
#     run:
#         import os
#         import subprocess

#         # cmd = """
#         # cat {c} {n} {p} {a}| sort |uniq -c | awk "$1==4 {print $2}" > {o}
#         # """.format(c=input.cpc2, n=input.cnci, p=input.pfam, a=input.cpat, o=output.lst)

#         # print(cmd)
#         # subprocess.call(cmd, shell=True)
#
##################################################################
## ================ Part 05 Expression Analysis =============== ##
##################################################################
#
# Step 01: Re Assembly by Stringtie
rule Part05_Expression_Analysis_01_ReStringtieAssemble:
    input:
        bam = OUTPUTDIR + "./MappingAndAssembly/01.Hisat2Genome/{sample}.bam",
        gtf = OUTPUTDIR + "./MappingAndAssembly/05.StringtieMerge/StringtieMerged.gtf"
    output:
        gtf = OUTPUTDIR + "./Expression_Analysis/01.ReStringtieAssemble/{sample}.gtf",
        cov = OUTPUTDIR + "./Expression_Analysis/01.ReStringtieAssemble/{sample}.coverage.cov",
        abu = OUTPUTDIR + "./Expression_Analysis/01.ReStringtieAssemble/{sample}.GeneAbund.txt"
    message:
        "Start to assembly use stringtie."
    log:
        OUTPUTDIR + "./AllLogs/Expression_Analysis/01.ReStringtieAssemble/{sample}.assembly.logs"
    priority:
        50
    threads:
        8
    params:
        opt = "-m 200 -l STRG -a 10 -c 1 -e",
        bal = OUTPUTDIR + "./Expression_Analysis/01.ReStringtieAssemble/Ballgown/{sample}"
    shell:
        """
        stringtie {input.bam} -G {input.gtf} -o {output.gtf} -p {threads} {params.opt} \
        -C {output.cov} -A {output.abu} -b {params.bal}
        """
#
# Step 02: Get assembled gtf list
rule Part05_Expression_Analysis_02_GetAssembledGtfList:
    input:
        expand(OUTPUTDIR + "./Expression_Analysis/01.ReStringtieAssemble/{sample}.gtf", sample = SAMPLES)
    output:
        lst = OUTPUTDIR + "./Expression_Analysis/02.GetAssembledGtfList/AssembledGtfList.txt"
    threads:
        1
    run:
        with open(output.lst, 'w') as f:
            for gtf in input:
                sample = re.search(r'01\.ReStringtieAssemble/(.*).gtf', gtf).group(1)
                outline = "\t".join((sample, gtf))
                print(outline, file=f)
#
# Step 03: Get gene and transcripts count matrix
rule Part05_Expression_Analysis_03_GetCountAndTPMMatrix:
    input:
        OUTPUTDIR + "./Expression_Analysis/02.GetAssembledGtfList/AssembledGtfList.txt"    
    output:
        gene = OUTPUTDIR + "./Expression_Analysis/03.GetCountAndTPMMatrix/gene_count_matrix.csv",
        mrna = OUTPUTDIR + "./Expression_Analysis/03.GetCountAndTPMMatrix/transcript_count_matrix.csv",
        tran = OUTPUTDIR + "./Expression_Analysis/03.GetCountAndTPMMatrix/transcript_tpm_matrix.tsv"
    log:
        cnt = OUTPUTDIR + "./AllLogs/Expression_Analysis/03.GetCountAndTPMMatrix/MakeCountMatrix.log",
        tpm = OUTPUTDIR + "./AllLogs/Expression_Analysis/03.GetCountAndTPMMatrix/MakeTPMMatrix.log"
    params:
        "-l 145 -s MSTRG"
    shell:
        """
        prepDE.py -i {input} -g {output.gene} -t {output.mrna} {params} 2> {log.cnt} && \
        perl Scripts/GetTPMFromSringtieGtfList.pl {input} > {output.tran} 2> {log.tpm}
        """
# #
# # Step 04: Different Expression Analysis by edgeR
# rule Part05_Expression_Analysis_04_DGEbyEdgeR2:
#     input:
#         gcm = OUTPUTDIR + "./Expression_Analysis/03.GetCountAndTPMMatrix/gene_count_matrix.csv",
#         mcm = OUTPUTDIR + "./Expression_Analysis/03.GetCountAndTPMMatrix/transcript_count_matrix.csv",
#         cp = COMPARE_PAIRS,
#         sg = SAMPLE_GROUPS
#     output:
#         ok = OUTPUTDIR + "./Expression_Analysis/04.DGEbyEdgeR2/DGEbyEdgeR2.Result.ok"
#     log:
#         OUTPUTDIR + "./AllLogs/Expression_Analysis/04.DGEbyEdgeR2/DGEbyEdgeR2.log"
#     threads:
#         8
#     params:
#         opt = "--foldchange 2 --fdr 0.05",
#         out_dir = OUTPUTDIR + "./Expression_Analysis/04.DGEbyEdgeR2/"
#     shell:
#         """
#         source activate r_edger_env && \
#         # Rscript Scripts/DGEs.R {params} -t {threads} --matrix {input.gcm} --group {input.sg} --compare {input.cp} -o {params.out_dir} && \
#         # Rscript Scripts/DGEs.R {params} -t {threads} -m {input.mcm} -g {input.sg} -c {input.cp} -o {params.out_dir} && \
#         echo SUCCESS > {output.ok} 
#         """
#
##################################################################
## ================ Part 06 circRNA Identification ============ ##
##################################################################
#
# Step 01: Trim the read to a specified length [CIRI-Full needed]
rule Part06_CircRNA_Analysis_01_Trim2EqualReads:
    input:
        R1 = OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R1.fq.gz",
        R2 = OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R2.fq.gz",
    output:
        R1P = OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.R1P.fq.gz",
        R2P = OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.R2P.fq.gz",
        R1U = OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.R1U.fq.gz",
        R2U = OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.R2U.fq.gz",
        sum = OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.summary.txt",
        ok = OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.summary.ok",
    log:
        OUTPUTDIR + "./AllLogs/CircRNA/01.Trim2EqualReads/TrimReads_{sample}.log"
    threads:
        8
    params:
        opt = "-Xmx10g -Xms2g PE -phred33",
        len = "CROP:145 MINLEN:145"
    run:
        import os
        import subprocess

        cmd = """
        source activate qc_env && \
        trimmomatic {opt} -threads {t} -summary {sum} \
        {R1} {R2} {R1P} {R1U} {R2P} {R2U} {len} 2> {log}
        """.format(opt=params.opt, t=threads, sum=output.sum, 
        R1=input.R1, R2=input.R2, R1P=output.R1P, R1U=output.R1U, 
        R2P=output.R2P, R2U=output.R2U, len=params.len, log=log)

        print(cmd)
        subprocess.call(cmd, shell=True)

        if os.path.getsize(output.sum) > 0:
            subprocess.call("echo SUCCESS > {ok}".format(ok=output.ok), shell=True)
#
# Step 01: Map rRNA free reads to genome with bwa-mem
rule Part06_CircRNA_Analysis_02_BWA2Genome:
    input:
        idx = BWA_DNA_INDEX,
        R1P = OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.R1P.fq.gz",
        R2P = OUTPUTDIR + "./CircRNA/01.Trim2EqualReads/{sample}.R2P.fq.gz",
    output:
        sam = OUTPUTDIR + "./CircRNA/02.BWA2Genome/{sample}.sam"
    log:
        OUTPUTDIR + "./AllLogs/CircRNA/02.BWA2Genome/{sample}.bwa_mem.log"
    message:
        "Start to filter rRNA baseed on RNACentral plant rRNA database by bowtie2."
    params:
        "-T 19"    
    threads:
        8    
    shell:
        """
        bwa mem {params} -o {output.sam} -t {threads} {input.idx} {input.R1P} {input.R2P} 2> {log}
        """
#
# Step 02: CircRNA identification with CIRI2
rule Part06_CircRNA_Analysis_02_CICR2_Prediction:
    input:
        ref = DNA,
        gtf = GTF,
        sam = OUTPUTDIR + "./CircRNA/02.BWA2Genome/{sample}.sam"
    output:
        ciri = OUTPUTDIR + "./CircRNA/03.CIRI2_Prediction/{sample}.ciri"
    log:
        OUTPUTDIR + "./AllLogs/CircRNA/03.CIRI2_Prediction/{sample}.ciri.log"
    params:
        "-M Mt -Q"
    threads:
        8    
    shell:
        """
        source activate ciri_py2_env && \
        perl /opt/biosoft/CIRI-full_v2.0/bin/CIRI_v2.0.6/CIRI2.pl {params} -T {threads} -I {input.sam} \
        -O {output.ciri} -F {input.ref} -A {input.gtf} --log {log} 
        """        
#
# Step 03: CircRNA Quantitation
rule Part06_CircRNA_Analysis_03_CircRNA_Quantitation:
    input:
        conf = CIRI_QUANT_CFG,
        ciri = OUTPUTDIR + "./CircRNA/02.CIRI2_Prediction/{sample}.ciri",
        R1 = OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R1.fq.gz",
        R2 = OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R2.fq.gz"
    output:
        gtf = OUTPUTDIR + "./CircRNA/03.CircRNA_Quantitation/{sample}/{sample}.gtf"
    log:
        OUTPUTDIR + "./AllLogs/CircRNA/03.CircRNA_Quantitation/{sample}.CIRIquant.log"
    params:
        opt = "--tool CIRI2",
        dir = OUTPUTDIR + "./CircRNA/03.CircRNA_Quantitation/{sample}",
        perfix = "{sample}",
    threads:
        8
    shell:
        """
        source activate ciri_py2_env && \
        CIRIquant --config {input.conf} {params.opt} -1 {input.R1} -2 {input.R2} -o {params.dir} \
            -p {params.perfix} -t {threads} --circ {input.ciri} -e {log}
        """
#
# Step 04: Identify circRNA by find_circ -- 1.mapping
rule Part06_CircRNA_Analysis_04_Bowtie2ToGenome:
    input:
        idx = BOWTIE2_DNA_INDEX,
        R1 = OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R1.fq.gz",
        R2 = OUTPUTDIR + "./Preprocess/03.rRNAFreeFastq/{sample}.R2.fq.gz"
    output:
        bam  = OUTPUTDIR + "./CircRNA/04.Bowtie2ToGenome/{sample}.sorted.bam",
    log:
        OUTPUTDIR + "./AllLogs/CircRNA/04.Bowtie2ToGenome/{sample}.align.logs"
    threads:
        8
    params:
        bwt = "-q --very-sensitive --score-min=C,-15,0 --mm",
        sam = "-O bam -@ 4"
    shell:
        """
        bowtie2 {params.bwt} -p {threads} -x {input.idx} -1 {input.R1} -2 {input.R2} 2> {log} \
            | samtools view -hbuS - | samtools sort {params.sam} -o {output.bam}
        """
#
# Step 05: Identify circRNA by find_circ -- 2.Fetch unmapped read with bowtie2
rule Part06_CircRNA_Analysis_05_FetchUnmappedBam:
    input:
        bam = OUTPUTDIR + "./CircRNA/04.Bowtie2ToGenome/{sample}.sorted.bam"
    output:
        bam = OUTPUTDIR + "./CircRNA/05.UnmappedBam/{sample}.unmapped.bam",
        bai = OUTPUTDIR + "./CircRNA/05.UnmappedBam/{sample}.unmapped.bai",
    threads:
        4        
    shell:
        """
        samtools view -hf4 -@ {threads} {input.bam} | samtools view -@ {threads} -Sb - > {output.bam} && \
        samtools index -b -@ {threads} {output.bam} {output.bai}
        """
#
# Step 06: Identify circRNA by find_circ -- 3.Convert bam to qfa
rule Part06_CircRNA_Analysis_06_Bam2Anchors:
    input:
        bam = OUTPUTDIR + "./CircRNA/05.UnmappedBam/{sample}.unmapped.bam",
        bai = OUTPUTDIR + "./CircRNA/05.UnmappedBam/{sample}.unmapped.bai"
    output:
        fq = OUTPUTDIR + "./CircRNA/06.Bam2Anchors/{sample}/unmapped_anchors.fq.gz"
    params:
        DIR = OUTPUTDIR + "./CircRNA/06.Bam2Anchors/{sample}"
    run:
        import os
        import subprocess
        if not os.path.exists(params.DIR):
            os.makedirs(params.DIR)

        cmd = ("source activate find_circ_env && unmapped2anchors.py {bam} | "
        " gzip > {fq} && conda deactivate").format(bam=input.bam, fq=output.fq)

        print(cmd)
        subprocess.call(cmd, shell=True)
#
# Step 07: Identify circRNA by find_circ -- 4.Find circRNA
rule Part06_CircRNA_Analysis_07_FindCircRNA:
    input:
        fq = OUTPUTDIR + "./CircRNA/06.Bam2Anchors/{sample}/unmapped_anchors.fq.gz",
        dna = DNA,
        idx = BOWTIE2_DNA_INDEX
    output:
        fa = OUTPUTDIR + "./CircRNA/07.FindCircRNA/{sample}/spliced_reads.fa",
        bed = OUTPUTDIR + "./CircRNA/07.FindCircRNA/{sample}/spliced_reads.bed",
        stat = OUTPUTDIR + "./CircRNA/07.FindCircRNA/{sample}/stat.txt"
    log:
        OUTPUTDIR + "./AllLogs/CircRNA/07.FindCircRNA/{sample}.align.logs"
    threads:
        8
    params:
        bwt = "--score-min=C,-15,0 --reorder --mm -q",
        fic = "--prefix=zma_{sample}_  --name={sample}",
    run:
        import subprocess

        cmd = ("source activate find_circ_env && bowtie2 {bwt} -p {t} -U {fq} -x {idx} 2> {log} | "
        "find_circ.py --genome={dna} {fic} --stats={stat} --reads={fa} > {bed} && conda deactivate").format(bwt=params.bwt, 
        t=threads, fq=input.fq, idx=input.idx, log=log, dna=input.dna, fic=params.fic, stat=output.stat, fa=output.fa, bed=output.bed)

        print(cmd)
        subprocess.call(cmd, shell=True)
#
# Step 08: Identify circRNA by find_circ -- 5.Merge all samples bed
rule Part06_CircRNA_Analysis_08_MergeAllSamplesBed:
    input:
        expand(OUTPUTDIR + "./CircRNA/07.FindCircRNA/{sample}/spliced_reads.bed", sample = SAMPLES)
    output:
        bed = OUTPUTDIR + "./CircRNA/08.MergeAllSamplesBed/merged_spliced_reads.bed",
        stat = OUTPUTDIR + "./CircRNA/08.MergeAllSamplesBed/merged_stat.txt"
    run:
        import subprocess

        cmd = ("source activate find_circ_env && merge_bed.py {i} -s {s} > {bed}").format(i=input, 
        s=output.stat, bed=output.bed)
        
        print(cmd)
        subprocess.call(cmd, shell=True)
#
# Step 09: Identify circRNA by find_circ -- 6.Fetch Good circRNA
rule Part06_CircRNA_Analysis_09_FetchFinalCircRNA:
    input:
        OUTPUTDIR + "./CircRNA/08.MergeAllSamplesBed/merged_spliced_reads.bed"
    output:
        OUTPUTDIR + "./CircRNA/09.FinalCircRNA/circ_candidates.bed"
    threads:
        1
    params:
        DIR = OUTPUTDIR + "./CircRNA/09.FinalCircRNA/"
    run:
        import os
        import subprocess

        if not os.path.exists(params.DIR):
            os.makedirs(params.DIR)

        cmd = ("source activate find_circ_env && grep CIRCULAR {i} | awk '$5>=2' | "
                "grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE | maxlength.py 100000 > {o}"
                ).format(i=input, o=output)
        
        print(cmd)
        subprocess.call(cmd, shell=True)
#
##################################################################
## ================ Part 07 Final Analysis Report ============= ##
##################################################################
#
# Report 01: Fastq Filter Result
    # rule Report_01_FastqFilterState1:
    #     input:
    #         expand( OUTPUTDIR + "./Step01.FastqFilter/{sample}/{sample}.json", sample=SAMPLES )
    #     output:
    #         lst = OUTPUTDIR + "./Report/S01.FastqFilter.jsonlist.txt"
    #     run:
    #         with open(output.lst, 'w') as jo:
    #             for filename in input:
    #                 print(filename, file=jo)

    # ## -------- XXXXXXXXXXXXXXXXXXXXXXXXXXX --------------------
    # rule Report_01_FastqFilterState2:
    #     input:
    #         OUTPUTDIR + "./Report/S01.FastqFilter.jsonlist.txt"
    #     output:
    #         tsv = OUTPUTDIR + "./Report/S01.FastqFilter-State.tsv"
    #     shell:
    #         """
    #         source activate data_env && \
    #         python Scripts/GetFastqFilterState.py {input} {output}
    #         """
    # ## -------- Report 02: rRNA Filter Result --------
    # ## -------- Report 03: Genome alianment Result --------
    # rule Report_03_ReadDistribution:
    #     input:
    #         bed = BED,
    #         bam = OUTPUTDIR + "./Step03.Hisat2Genome/{sample}.bam"        
    #     output:
    #         OUTPUTDIR + "./Report/S03.ReadDistribution/{sample}.state.tsv"
    #     log:
    #         OUTPUTDIR + "./logs/Report/S03.GenomeAlianmentStatebyRseQC.{sample}.log"
    #     params:
    #         prx = OUTPUTDIR + "./Report/S03.ReadDistribution/",
    #         log_dir = OUTPUTDIR + "./logs/Report/"
    #     shell:
    #         """
    #         source activate py3.7 && \
    #         if [ ! -d {params.prx} ]; then mkdir -p {params.prx}; fi && \
    #         if [ ! -d {params.log_dir} ]; then mkdir -p {params.log_dir}; fi && \
    #         read_distribution.py -i {input.bam} -r {input.bed} > {output} 2> {log} && \
    #         conda deactivate
    #         """


    # ## -------- Report 04: Genome guid assembly Result --------
    # ## -------- Report 05: Genome guid assembly Result --------
    # ## -------- Report 06: Genome guid assembly Result --------
    # ## -------- Report 07: Genome guid assembly Result --------
    # ## -------- Report 08: Genome guid assembly Result --------
    # ## -------- Report 09: Genome guid assembly Result --------
    # ## -------- Report 10: Genome guid assembly Result --------
    # ## -------- Report 00: Output html report --------
#

