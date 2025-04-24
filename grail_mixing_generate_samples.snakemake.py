import csv

SAMPLES_INFO = "files/grail_mixing_metadata.txt"
# SAMPLES_INFO = "files/grail_mixing_metadata_test.txt"
EXONS = "files/UCSC_hg19_canonical_cds_grail_panel.bed"

# REPLACE WITH YOUR OWN PATH
# GENOME = "/mnt/Data01/genome/hg19.fa"
GENOME = "/path/to/genome/hg19.fa"

ATAC = "files/TCGA_ATAC_peak.all.probeMap.simple.hg19.bed"
TFBS = "files/Homo_sapiens_meta_clusters_hg19_midpoint_top5k_sorted.bed.gz"

def samples_to_list(sample_file):
	with open(sample_file) as f:
		samples = f.read().splitlines()
	return samples

def samples_to_dict(sample_info):
        info_dict = {}
        with open(sample_info, mode = "r") as f:
            reader = csv.reader(f, delimiter = "\t")
            next(reader)
            for row in reader:
                key = row[0]
                values = list(row[1:])
                info_dict[key] = values
        return info_dict

# order of list is cancer_sample[0], normal_sample[1], cancer_ratio[2], normal_ratio[3], cancer_reads[4], normal_reads[5]
SAMPLES_DICT = samples_to_dict(SAMPLES_INFO)

SAMPLES = list(SAMPLES_DICT.keys())

rule all:
    input:
        expand("data/grail/mixed_samples/{sample}.bed.gz", sample = SAMPLES),
        "output/mixing_data/feature_tables/se.rds",
        "output/mixing_data/feature_tables/depth.rds",
        "output/mixing_data/feature_tables/frag_bins.rds",
        "output/mixing_data/feature_tables/small_frags.rds",
        "output/mixing_data/feature_tables/full_gene_depth.rds",
        "output/mixing_data/feature_tables/mds.rds",
        "output/mixing_data/feature_tables/TFBS_entropy.rds",
        "output/mixing_data/feature_tables/ATAC_entropy.rds"

rule downsample_cancer:
    input:
        sample_bed = lambda wcs: "data/grail/beds/" + SAMPLES_DICT[wcs.sample][0] + ".bed.gz",
        exons = EXONS
    output:
        temp("data/grail/mixed_samples/{sample}_cancer_reads.bed.gz")
    params:
        cancer_reads = lambda wcs: SAMPLES_DICT[wcs.sample][4]
    shell:
        """
        zcat {input.sample_bed} \
        | shuf -n {params.cancer_reads} \
        | bedtools intersect -a stdin -b {input.exons} -wa -wb \
        | gzip > {output}
        """

rule downsample_normal:
    input:
        sample_bed = lambda wcs: "data/grail/beds/" + SAMPLES_DICT[wcs.sample][1] + ".bed.gz",
        exons = EXONS
    output:
        temp("data/grail/mixed_samples/{sample}_normal_reads.bed.gz")
    params:
        normal_reads = lambda wcs: SAMPLES_DICT[wcs.sample][5]
    shell:
        """
        zcat {input.sample_bed} \
        | shuf -n {params.normal_reads} \
        | bedtools intersect -a stdin -b {input.exons} -wa -wb \
        | gzip > {output}
        """

rule cat_and_overlap_exons:
    input:
        exons = EXONS,
        cancer_bed = "data/grail/mixed_samples/{sample}_cancer_reads.bed.gz",
        normal_bed = "data/grail/mixed_samples/{sample}_normal_reads.bed.gz"
    output:
        "data/grail/mixed_samples/{sample}.bed.gz"
    shell:
        """
        zcat {input.cancer_bed} {input.normal_bed} \
        | cut -f 1,2,3,8,9,10,11,12 \
        | sort -k1,1 -k2,2n \
        | gzip > {output}
        """


rule get_SE_fragstats:
    input:
        "data/grail/mixed_samples/{sample}.bed.gz"
    output:
        "output/mixing_data/fragstats/SE_files/{sample}.txt.gz"
    shell:
        """
        zcat {input} \
        | awk '{{print $6, "\t", $7, "\t", ($3-$2)}}' \
        | sort -k1,1 -k2,2n -k3,3n \
        | uniq -c \
        | awk '{{print $2, "\t", $3, "\t", $4, "\t", $1}}' \
        | gzip > {output}
        """


rule get_depth_fragstats:
    input:
        "data/grail/mixed_samples/{sample}.bed.gz"
    output:
        "output/mixing_data/fragstats/depth_files/{sample}.txt.gz"
    shell:
        """
        zcat {input} \
        | cut -f 6,7 \
        | sort -k1,1 -k2,2n \
        | uniq -c \
        | awk '{{print $2, "\t", $3, "\t", $1 }}' > {output}
        """
        
rule calculate_SE:
    input:
       "output/mixing_data/fragstats/SE_files/{sample}.txt.gz"
    output:
       "output/mixing_data/metrics/se/{sample}.SE.tsv"
    shell:
       """
       Rscript scripts/calculate_SE.R {input} {output}
       """

rule calculate_normalized_depth:
    input:
       "output/mixing_data/fragstats/depth_files/{sample}.txt.gz"
    output:
       "output/mixing_data/metrics/depth/{sample}.depth.tsv"
    shell:
       """
       Rscript scripts/calculate_depth.R {input} {output}
       """

rule calculate_frag_bins:
    input:
        "output/mixing_data/fragstats/SE_files/{sample}.txt.gz"
    output:
        "output/mixing_data/metrics/frag_bins/{sample}.fragbins.tsv"
    shell:
        """
        Rscript scripts/calculate_frag_bins.R {input} {output}
        """

rule calculate_small_frags:
    input:
        "output/mixing_data/fragstats/SE_files/{sample}.txt.gz"
    output:
        "output/mixing_data/metrics/small_frags/{sample}.smallfrag.tsv"
    shell:
        """
        Rscript scripts/calculate_small_frags.R {input} {output}
        """

rule calculate_full_gene_depth:
    input:
        "output/mixing_data/fragstats/depth_files/{sample}.txt.gz"
    output:
        "output/mixing_data/metrics/full_gene_depth/{sample}.fullgenedepth.tsv"
    shell:
        """
        Rscript scripts/calculate_full_gene_depth.R {input} {output}
        """

rule get_left_4mer:
    input:
        sample_reads = "data/grail/mixed_samples/{sample}.bed.gz",
        genome = GENOME
    output:
        "output/mixing_data/metrics/mds/{sample}_left_4mers.txt.gz"
    shell:
        """
        zcat {input.sample_reads} \
        | awk -F "\t" '{{OFS=FS}}; {{print $1, $2, $2+4, $4, $5, $6, $7, $8}}' \
        | bedtools getfasta -fi {input.genome} -bed stdin -bedOut \
        | awk -F "\t" '{{OFS=FS}}; {{print $6, $7, toupper($9)}}' \
        | sort -k1,1 -k2,2n -k3,3 \
        | uniq -c \
        | gzip > {output}
        """

rule get_right_4mer:
    input:
        sample_reads = "data/grail/mixed_samples/{sample}.bed.gz",
        genome = GENOME
    output:
        "output/mixing_data/metrics/mds/{sample}_right_4mers.txt.gz"
    shell:
        """
        zcat {input.sample_reads} \
        | awk -F "\t" '{{OFS=FS}}; {{print $1, $3-4, $3, $4, $5, "-", $6, $7, $8}}' \
        | bedtools getfasta -fi {input.genome} -bed stdin -bedOut -s \
        | awk -F "\t" '{{OFS=FS}}; {{print $7, $8, toupper($10)}}' \
        | sort -k1,1 -k2,2n -k3,3 \
        | uniq -c \
        | gzip > {output}
        """

rule calculate_MDS:
    input:
        left = "output/mixing_data/metrics/mds/{sample}_left_4mers.txt.gz",
        right = "output/mixing_data/metrics/mds/{sample}_right_4mers.txt.gz",
        mds_script = "scripts/calculate_mds.R"
    output:
        "output/mixing_data/metrics/mds/{sample}_mds.txt"
    shell:
        """
        Rscript {input.mds_script} {input.left} {input.right} {output}
        """

# top 5k sites
rule overlap_TFBS:
    input:
        sample_reads = "data/grail/mixed_samples/{sample}.bed.gz",
        tfbs = TFBS
    output:
        "output/mixing_data/metrics/TFBS_entropy/{sample}_TFBS_frag_count.txt.gz"
    shell:
        """
        zcat {input.sample_reads} \
        | bedtools intersect -a stdin -b {input.tfbs} -wa -wb \
        | awk -F"\t" 'BEGIN {{OFS=FS}} {{print $3-$2, $12}}' \
        | sort -k2,2 -k1,1n \
        | uniq -c \
        | gzip  > {output}
        """

rule calculate_TFBS_entropy:
    input:
        "output/mixing_data/metrics/TFBS_entropy/{sample}_TFBS_frag_count.txt.gz"
    output:
        "output/mixing_data/metrics/TFBS_entropy/{sample}_TFBS_entropy.txt"
    shell:
        """
        Rscript scripts/calculate_TFBS_entropy.R {input} {output}
        """

rule overlap_ATAC:
    input:
        sample_reads = "data/grail/mixed_samples/{sample}.bed.gz",
        atac = ATAC
    output:
        "output/mixing_data/metrics/ATAC_entropy/{sample}_ATAC_frag_count.txt.gz"
    shell:
        """
        zcat {input.sample_reads} \
        | bedtools intersect -a stdin -b {input.atac} -wa -wb \
        | awk -F"\t" '{{OFS=FS}}; {{split($12, arr, "_"); print arr[1], $3-$2}}' \
        | sort -k1,1 -k2,2n \
        | uniq -c \
        | gzip > {output}
        """

rule calculate_ATAC_entropy:
    input:
        "output/mixing_data/metrics/ATAC_entropy/{sample}_ATAC_frag_count.txt.gz"
    output:
        "output/mixing_data/metrics/ATAC_entropy/{sample}_ATAC_entropy.txt"
    shell:
        """
        Rscript scripts/calculate_ATAC_entropy.R {input} {output}
        """

rule build_feature_tables:
    input:
        # all of the .SE.tsv files for each 
        expand("output/mixing_data/metrics/se/{sample}.SE.tsv", sample = SAMPLES),
        expand("output/mixing_data/metrics/depth/{sample}.depth.tsv", sample = SAMPLES),
        expand("output/mixing_data/metrics/frag_bins/{sample}.fragbins.tsv", sample = SAMPLES),
        expand("output/mixing_data/metrics/small_frags/{sample}.smallfrag.tsv", sample = SAMPLES),
        expand("output/mixing_data/metrics/full_gene_depth/{sample}.fullgenedepth.tsv", sample = SAMPLES),
        expand("output/mixing_data/metrics/mds/{sample}_mds.txt", sample = SAMPLES),
        expand("output/mixing_data/metrics/TFBS_entropy/{sample}_TFBS_entropy.txt", sample = SAMPLES),
        expand("output/mixing_data/metrics/ATAC_entropy/{sample}_ATAC_entropy.txt", sample = SAMPLES)
    output:
        "output/mixing_data/feature_tables/se.rds",
        "output/mixing_data/feature_tables/depth.rds",
        "output/mixing_data/feature_tables/frag_bins.rds",
        "output/mixing_data/feature_tables/small_frags.rds",
        "output/mixing_data/feature_tables/full_gene_depth.rds",
        "output/mixing_data/feature_tables/mds.rds",
        "output/mixing_data/feature_tables/TFBS_entropy.rds",
        "output/mixing_data/feature_tables/ATAC_entropy.rds"
    shell:
        """
        Rscript scripts/build_feature_tables_mixing_data.R
        """
