import subprocess
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
import matplotlib.pyplot as plt
import random
from pysam import VariantFile, AlignmentFile, TabixFile
from config.py import bam_path, vcf_path, newbam_name


def add_unintersecting_records_to_bam(bed_path, bam, new_bam):
    ref_bam = AlignmentFile(bam)
    unintersect_df = pd.read_csv(bed_path, sep='\t', names=["chromosome", "start", "end", "id", "quality", "strand"], usecols = range(6))

    for row_number, row in unintersect_df.iterrows():
        records = ref_bam.fetch(str(row["chromosome"]), row["start"], row["end"])

        for record in records:
            new_bam.write(record)

    return new_bam


def extract_genotype(snp_info):
    return snp_info.split(":")[0]

def insert_snp(seq, pos, ref, alts, genotype):

    if genotype == "0|0":
        new_seq = seq[:pos] + ref + seq[pos+1:]

    elif genotype == "1|1":
        new_seq = seq[:pos] + alts + seq[pos+1:]

    elif genotype == "0|1":
        alleles = genotype.split("|")
        rand = random.choice(alleles)

        if rand == "0":
            new_seq = seq[:pos] + ref + seq[pos+1:]

        elif rand == "1":
            new_seq = seq[:pos] + alts + seq[pos+1:]

    return new_seq


def change_snps_and_save_to_bam(bed_path, bam, new_bam, vcf_tabix):
    ref_bam = AlignmentFile(bam)
    ref_tabix = TabixFile(vcf_tabix)
    intersect_df = pd.read_csv(bed_path, sep='\t', names=["chromosome", "start", "end", "id", "quality", "strand"], usecols = range(6))

    for row_number, row in intersect_df.iterrows():
        alignments = ref_bam.fetch(str(row["chromosome"]), row["start"], row["end"])

        for alignment in alignments:
            snps = ref_tabix.fetch(str(row["chromosome"]), row["start"], row["end"], parser=pysam.asTuple())
            seq = alignment.seq

            for snp in snps:

                pos = snp[1]
                ref = snp[3]
                alts = snp[4]
                snp_info = snp[9]

                if (len(ref) == 1 and len(alts) == 1):
                    change_pos = int(pos) - int(alignment.pos)
                    genotype = extract_genotype(snp_info)

                    alignment.seq = insert_snp(seq, change_pos, ref, alts, genotype)
                    new_bam.write(alignment)
    return new_bam

def make_bam(unintersect_bed, intersect_bed, bam, vcf_tabix, newbam_name):
    ref_bam = AlignmentFile(bam)
    ref_tabix = TabixFile(vcf_tabix)
    new_bam = AlignmentFile(newbam_name, "wb", template = ref_bam)
    add_unintersecting_records_to_bam(unintersect_bed, ref_bam, new_bam)
    change_snps_and_save_to_bam(intersect_bed, ref_bam, new_bam, ref_tabix)
    new_bam.close()
    return new_bam


def add_to_vcf(bed_path, vcf, newvcf_name):
    ref_vcf = VariantFile(vcf)
    new_vcf = VariantFile(newvcf_name, "wb", header=ref_vcf.header)
    intersect_df = pd.read_csv(bed_path, sep='\t', names=["chromosome", "start", "end", "id", "quality", "strand"], usecols = range(6))

    for row_number, row in unintersect_df.iterrows():
        records = ref_vcf.fetch(str(row["chromosome"]), row["start"], row["end"])

        for record in records:
            new_vcf.write(record)

    return new_vcf


def pipeline(input_bam, input_vcf, newbam_name):
    print("Starting pipeline...")
    sort_command = f"samtools sort -o sort.bam {input_bam}"
    subprocess.run(sort_command, shell=True)

    print("Bam sorted....")
    subprocess.run("samtools index sort.bam", shell=True)

    print("Bam indexed....")
    vcf_intersect_command = f"bedtools intersect -a sort.bam -sorted -b {input_vcf} -sorted -wa -bed > vcf_intersect.bed"
    subprocess.run(vcf_intersect_command, shell=True)

    print("Ready to reduce the vcf...")
    
    smaller_vcf = add_to_vcf("vcf_intersect.bed", input_vcf, "smaller_vcf.vcf")
    subprocess.run("bcftools sort -o sorted_smaller_vcf.vcf smaller_vcf.vcf", shell=True)
    subprocess.run("bgzip sorted_smaller_vcf.vcf, shell=True)
    

    print("Vcf reduced...")
    intersect_command = f"bedtools intersect -a sort.bam -sorted -b sorted_smaller_vcf.vcf.gz -sorted -wa -bed > intersect.bed"
    subprocess.run(intersect_command, shell=True)

    print("Bam & vcf intersection created....")
    unintersect_command = f"bedtools intersect -a sort.bam -sorted -b sorted_smaller_vcf.vcf.gz -sorted -v -wa -bed > unintersect.bed"
    subprocess.run(unintersect_command, shell=True)

    print("Bam & vcf unintersection created....")
    vcf_index_command = f"bcftools index -f -t sorted_smaller_vcf.vcf.gz"
    subprocess.run(vcf_index_command, shell=True)

    print("Everything is ready, bam creation has started...")
    newbam = make_bam("unintersect.bed", "intersect.bed", input_bam, "sorted_smaller_vcf.vcf.gz.tbi", newbam_name)
    # remove temporary files

pipeline(bam_path, vcf_path, newbam_name)
