# static paths
SAMPLES_PATH = "files/grail_samples_test.txt"
EXONS = "files/UCSC_hg19_canonical_cds_grail_panel.bed"

ATAC = "files/TCGA_ATAC_peak.all.probeMap.simple.hg19.bed"
TFBS = "files/Homo_sapiens_meta_clusters_hg19_midpoint_top5k_sorted.bed.gz"

# REPLACE WITH YOUR OWN PATH
GENOME = "/path/to/genome/hg19.fa"
# GENOME = "/mnt/Data01/genome/hg19.fa"


def samples_to_list(sample_file):
	with open(sample_file) as f:
		samples = f.read().splitlines()
	return samples

# sample names without the ".bed.gz" at the end of the file name
SAMPLES = samples_to_list(SAMPLES_PATH)

LEVELS = ["1M", "5M", "10M", "25M", "50M", "100M"]

DOWNSAMPLE_DICT = {"1M":1000000, "5M":5000000, "10M":10000000, "25M":25000000, "50M":50000000, "100M":100000000}

rule all:
    input:
        expand("output/downsample/cds_overlap/{level}/{sample}.bed.gz", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/fragstats/SE_files/{level}/{sample}.txt.gz", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/fragstats/depth_files/{level}/{sample}.txt.gz", sample = SAMPLES, level = LEVELS),
        # generating metrics
        expand("output/downsample/metrics/se/{level}/{sample}.SE.tsv", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/depth/{level}/{sample}.depth.tsv", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/frag_bins/{level}/{sample}.fragbins.tsv", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/small_frags/{level}/{sample}.smallfrag.tsv", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/full_gene_depth/{level}/{sample}.fullgenedepth.tsv", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/mds/{level}/{sample}_mds.txt", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/TFBS_entropy/{level}/{sample}_TFBS_entropy.txt", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/ATAC_entropy/{level}/{sample}_ATAC_entropy.txt", sample = SAMPLES, level = LEVELS),

        # building feature tables
        expand("output/downsample/feature_tables/grail_se_{level}.rds", level = LEVELS),
        expand("output/downsample/feature_tables/grail_depth_{level}.rds", level = LEVELS),
        expand("output/downsample/feature_tables/grail_frag_bins_{level}.rds", level = LEVELS),
        expand("output/downsample/feature_tables/grail_small_frags_{level}.rds", level = LEVELS),
        expand("output/downsample/feature_tables/grail_full_gene_depth_{level}.rds", level = LEVELS),
        expand("output/downsample/feature_tables/grail_mds_{level}.rds", level = LEVELS),
        expand("output/downsample/feature_tables/grail_TFBS_entropy_{level}.rds", level = LEVELS),
        expand("output/downsample/feature_tables/grail_ATAC_entropy_{level}.rds", level = LEVELS),
        
        


rule downsample:
    input:
        sample_bed = "data/grail/beds/{sample}.bed.gz"
    output:
        "output/downsample/all_seqs/{level}/{sample}.bed.gz"
    params:
        downsample_level = lambda wcs: DOWNSAMPLE_DICT[wcs.level]
    shell:
        """
        zcat {input} | shuf -n {params.downsample_level} | sort -k1,1 -k2,2n | gzip > {output}
        """
        
rule overlap_exons:
    input:
        exons = EXONS,
        sample_bed = "output/downsample/all_seqs/{level}/{sample}.bed.gz"
    output:
        "output/downsample/cds_overlap/{level}/{sample}.bed.gz"
    shell:
        """
        bedtools intersect -a {input.sample_bed} -b {input.exons} -wa -wb \
        | cut -f 1,2,3,8,9,10,11,12 \
        | sort -k1,1 -k2,2n \
        | gzip > {output}
        """
        
        
rule get_SE_fragstats:
    input:
        "output/downsample/cds_overlap/{level}/{sample}.bed.gz"
    output:
        "output/downsample/fragstats/SE_files/{level}/{sample}.txt.gz"
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
        "output/downsample/cds_overlap/{level}/{sample}.bed.gz"
    output:
        "output/downsample/fragstats/depth_files/{level}/{sample}.txt.gz"
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
       "output/downsample/fragstats/SE_files/{level}/{sample}.txt.gz"
    output:
       "output/downsample/metrics/se/{level}/{sample}.SE.tsv"
    shell:
       """
       Rscript scripts/calculate_SE.R {input} {output}
       """

rule calculate_normalized_depth:
    input:
       "output/downsample/fragstats/depth_files/{level}/{sample}.txt.gz"
    output:
       "output/downsample/metrics/depth/{level}/{sample}.depth.tsv"
    shell:
       """
       Rscript scripts/calculate_depth.R {input} {output}
       """

rule calculate_frag_bins:
    input:
        "output/downsample/fragstats/SE_files/{level}/{sample}.txt.gz"
    output:
        "output/downsample/metrics/frag_bins/{level}/{sample}.fragbins.tsv"
    shell:
        """
        Rscript scripts/calculate_frag_bins.R {input} {output}
        """

rule calculate_small_frags:
    input:
        "output/downsample/fragstats/SE_files/{level}/{sample}.txt.gz"
    output:
        "output/downsample/metrics/small_frags/{level}/{sample}.smallfrag.tsv"
    shell:
        """
        Rscript scripts/calculate_small_frags.R {input} {output}
        """

rule calculate_full_gene_depth:
    input:
        "output/downsample/fragstats/depth_files/{level}/{sample}.txt.gz"
    output:
        "output/downsample/metrics/full_gene_depth/{level}/{sample}.fullgenedepth.tsv"
    shell:
        """
        Rscript scripts/calculate_full_gene_depth.R {input} {output}
        """

rule get_left_4mer:
    input:
        sample_reads = "output/downsample/cds_overlap/{level}/{sample}.bed.gz",
        genome = GENOME
    output:
        "output/downsample/metrics/mds/{level}/{sample}_left_4mers.txt.gz"
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
        sample_reads = "output/downsample/cds_overlap/{level}/{sample}.bed.gz",
        genome = GENOME
    output:
        "output/downsample/metrics/mds/{level}/{sample}_right_4mers.txt.gz"
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
        left = "output/downsample/metrics/mds/{level}/{sample}_left_4mers.txt.gz",
        right = "output/downsample/metrics/mds/{level}/{sample}_right_4mers.txt.gz",
        mds_script = "scripts/calculate_mds.R"
    output:
        "output/downsample/metrics/mds/{level}/{sample}_mds.txt"
    shell:
        """
        Rscript {input.mds_script} {input.left} {input.right} {output}
        """

# top 5k sites
rule overlap_TFBS:
    input:
        sample_reads = "output/downsample/cds_overlap/{level}/{sample}.bed.gz",
        #bedtools = BEDTOOLS,
        tfbs = TFBS
    output:
        "output/downsample/metrics/TFBS_entropy/{level}/{sample}_TFBS_frag_count.txt.gz"
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
        "output/downsample/metrics/TFBS_entropy/{level}/{sample}_TFBS_frag_count.txt.gz"
    output:
        "output/downsample/metrics/TFBS_entropy/{level}/{sample}_TFBS_entropy.txt"
    shell:
        """
        Rscript scripts/calculate_TFBS_entropy.R {input} {output}
        """

rule overlap_ATAC:
    input:
        sample_reads = "output/downsample/cds_overlap/{level}/{sample}.bed.gz",
        atac = ATAC
    output:
        "output/downsample/metrics/ATAC_entropy/{level}/{sample}_ATAC_frag_count.txt.gz"
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
        "output/downsample/metrics/ATAC_entropy/{level}/{sample}_ATAC_frag_count.txt.gz"
    output:
        "output/downsample/metrics/ATAC_entropy/{level}/{sample}_ATAC_entropy.txt"
    shell:
        """
        Rscript scripts/calculate_ATAC_entropy.R {input} {output}
        """

rule build_feature_tables:
    input:
        # all of the .SE.tsv files for each 
        expand("output/downsample/metrics/se/{level}/{sample}.SE.tsv", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/depth/{level}/{sample}.depth.tsv", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/frag_bins/{level}/{sample}.fragbins.tsv", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/small_frags/{level}/{sample}.smallfrag.tsv", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/full_gene_depth/{level}/{sample}.fullgenedepth.tsv", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/mds/{level}/{sample}_mds.txt", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/TFBS_entropy/{level}/{sample}_TFBS_entropy.txt", sample = SAMPLES, level = LEVELS),
        expand("output/downsample/metrics/ATAC_entropy/{level}/{sample}_ATAC_entropy.txt", sample = SAMPLES, level = LEVELS)
    output:
        "output/downsample/feature_tables/grail_se_{level}.rds",
        "output/downsample/feature_tables/grail_depth_{level}.rds",
        "output/downsample/feature_tables/grail_frag_bins_{level}.rds",
        "output/downsample/feature_tables/grail_small_frags_{level}.rds",
        "output/downsample/feature_tables/grail_full_gene_depth_{level}.rds",
        "output/downsample/feature_tables/grail_mds_{level}.rds",
        "output/downsample/feature_tables/grail_TFBS_entropy_{level}.rds",
        "output/downsample/feature_tables/grail_ATAC_entropy_{level}.rds"
    shell:
        """
        Rscript scripts/build_feature_tables_downsample.R
        """

