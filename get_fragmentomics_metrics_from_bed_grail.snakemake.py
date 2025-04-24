# static paths
SAMPLES_PATH = "files/grail_samples_test.txt"
ATAC = "files/TCGA_ATAC_peak.all.probeMap.simple.hg19.bed"
TFBS = "files/Homo_sapiens_meta_clusters_hg19_midpoint_top5k_sorted.bed.gz"

# REPLACE WITH YOUR OWN PATH
GENOME = "/mnt/Data01/genome/hg19.fa"
# GENOME = "path/to/genome.fa"

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
    

# sample names without the ".bed.gz" at the end of the file name
SAMPLES = samples_to_list(SAMPLES_PATH)

COMM_PANELS = ["default", "tempus", "guardant", "foundationOne"]

rule all:
    input:
        # Filtered Commercial Panel Data
        expand("data/{comm_panel}/{sample}.bed.gz", sample = SAMPLES, comm_panel = COMM_PANELS),
        # SE
        expand("output/fragstats/SE_files/{comm_panel}/{sample}.txt.gz", sample = SAMPLES, comm_panel = COMM_PANELS),
        expand("output/metrics/se/{comm_panel}/{sample}.SE.tsv", sample = SAMPLES, comm_panel = COMM_PANELS),
        # Depth
        expand("output/fragstats/depth_files/{comm_panel}/{sample}.txt", sample = SAMPLES, comm_panel = COMM_PANELS),
        expand("output/metrics/depth/{comm_panel}/{sample}.depth.tsv", sample = SAMPLES, comm_panel = COMM_PANELS),
        # frag bins
        expand("output/metrics/frag_bins/{comm_panel}/{sample}.fragbins.tsv", sample = SAMPLES, comm_panel = COMM_PANELS),
        # small frags
        expand("output/metrics/small_frags/{comm_panel}/{sample}.smallfrag.tsv", sample = SAMPLES, comm_panel = COMM_PANELS),
        # full gene depth
        expand("output/metrics/full_gene_depth/{comm_panel}/{sample}.fullgenedepth.tsv", sample = SAMPLES, comm_panel = COMM_PANELS),
        # MDS
        expand("output/metrics/mds/{comm_panel}/{sample}_mds.txt", sample = SAMPLES, comm_panel = COMM_PANELS),
        # TFBS
        expand("output/metrics/TFBS_entropy/{comm_panel}/{sample}_TFBS_frag_count.txt.gz", sample = SAMPLES, comm_panel = COMM_PANELS),
        expand("output/metrics/TFBS_entropy/{comm_panel}/{sample}_TFBS_entropy.txt", sample = SAMPLES, comm_panel = COMM_PANELS),
        # ATAC
        expand("output/metrics/ATAC_entropy/{comm_panel}/{sample}_ATAC_frag_count.txt.gz", sample = SAMPLES, comm_panel = COMM_PANELS),
        expand("output/metrics/ATAC_entropy/{comm_panel}/{sample}_ATAC_entropy.txt", sample = SAMPLES, comm_panel = COMM_PANELS),
        # Feature Tables
        expand("output/feature_tables/grail/{comm_panel}/se.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/depth.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/frag_bins.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/small_frags.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/full_gene_depth.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/mds.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/TFBS_entropy.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/ATAC_entropy.rds", comm_panel = COMM_PANELS),

        expand("output/feature_tables/grail/{comm_panel}/se_E1.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/depth_E1.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/small_frags_E1.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/mds_E1.rds", comm_panel = COMM_PANELS),

        expand("output/feature_tables/grail/{comm_panel}/all_combined.rds", comm_panel = COMM_PANELS)



rule filter_comm_panel_genes:
    input:
        "data/{sample}.bed.gz"
    output:
        # "data/{comm_panel}/{sample}.bed",
        "data/default/{sample}.bed",
        "data/tempus/{sample}.bed",
        "data/foundationOne/{sample}.bed",
        "data/guardant/{sample}.bed"
    shell:
        """
        Rscript scripts/filter_comm_panels.R {input}
        """

rule gzip_comm_panel_genes:
    input:
        "data/{comm_panel}/{sample}.bed"
    output:
        "data/{comm_panel}/{sample}.bed.gz"
    shell:
        """
        gzip {input}
        """


rule get_SE_fragstats:
    input:
        "data/{comm_panel}/{sample}.bed.gz"
    output:
        "output/fragstats/SE_files/{comm_panel}/{sample}.txt.gz"
    priority: 20
    shell:
        """
        zcat {input} \
        | awk '{{print $6, "\t", $7, "\t", ($3-$2)}}' \
        | sort -k1,1 -k2,2n -k3,3n \
        | uniq -c \
        | awk '{{print $2, "\t", $3, "\t", $4, "\t", $1}}' \
        | gzip > {output}
        """

rule calculate_SE:
    input:
        "output/fragstats/SE_files/{comm_panel}/{sample}.txt.gz"
    output:
        "output/metrics/se/{comm_panel}/{sample}.SE.tsv"
    shell:
        """
        Rscript scripts/calculate_SE.R {input} {output}
        """

rule get_depth_fragstats:
    input:
        "data/{comm_panel}/{sample}.bed.gz"
    output:
        "output/fragstats/depth_files/{comm_panel}/{sample}.txt"
    priority: 20
    shell:
        """
        zcat {input} \
        | cut -f 6,7 \
        | sort -k1,1 -k2,2n \
        | uniq -c \
        | awk '{{print $2, "\t", $3, "\t", $1 }}' > {output}
        """

rule calculate_normalized_depth:
    input:
        "output/fragstats/depth_files/{comm_panel}/{sample}.txt"
    output:
        "output/metrics/depth/{comm_panel}/{sample}.depth.tsv"
    shell:
        """
        Rscript scripts/calculate_depth.R {input} {output}
        """

rule calculate_frag_bins:
    input:
        "output/fragstats/SE_files/{comm_panel}/{sample}.txt.gz"
    output:
        "output/metrics/frag_bins/{comm_panel}/{sample}.fragbins.tsv"
    shell:
        """
        Rscript scripts/calculate_frag_bins.R {input} {output}
        """

rule calculate_small_frags:
    input:
        "output/fragstats/SE_files/{comm_panel}/{sample}.txt.gz"  
    output:
        "output/metrics/small_frags/{comm_panel}/{sample}.smallfrag.tsv"
    shell:
        """
        Rscript scripts/calculate_small_frags.R {input} {output}
        """

rule calculate_full_gene_depth:
    input:
        "output/fragstats/depth_files/{comm_panel}/{sample}.txt"
    output:
        "output/metrics/full_gene_depth/{comm_panel}/{sample}.fullgenedepth.tsv"
    shell:
        """
        Rscript scripts/calculate_full_gene_depth.R {input} {output}
       """

rule get_left_4mer:
    input:
        sample_reads = "data/{comm_panel}/{sample}.bed.gz",
        genome = GENOME
    output:
        "output/metrics/mds/{comm_panel}/{sample}_left_4mers.txt.gz"
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
        sample_reads = "data/{comm_panel}/{sample}.bed.gz",
        genome = GENOME
    output:
        "output/metrics/mds/{comm_panel}/{sample}_right_4mers.txt.gz"
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
        left = "output/metrics/mds/{comm_panel}/{sample}_left_4mers.txt.gz",
        right = "output/metrics/mds/{comm_panel}/{sample}_right_4mers.txt.gz",
        mds_script = "scripts/calculate_mds.R"
    output:
        "output/metrics/mds/{comm_panel}/{sample}_mds.txt"
    shell:
        """
        Rscript {input.mds_script} {input.left} {input.right} {output}
        """

# top 5k sites
rule overlap_TFBS:
    input:
        sample_reads = "data/{comm_panel}/{sample}.bed.gz",
        #bedtools = BEDTOOLS,
        tfbs = TFBS
    output:
        "output/metrics/TFBS_entropy/{comm_panel}/{sample}_TFBS_frag_count.txt.gz"
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
        "output/metrics/TFBS_entropy/{comm_panel}/{sample}_TFBS_frag_count.txt.gz"
    output:
        "output/metrics/TFBS_entropy/{comm_panel}/{sample}_TFBS_entropy.txt"
    shell:
        """
        Rscript scripts/calculate_TFBS_entropy.R {input} {output}
        """

rule overlap_ATAC:
    input:
        sample_reads = "data/{comm_panel}/{sample}.bed.gz",
        atac = ATAC
    output:
        "output/metrics/ATAC_entropy/{comm_panel}/{sample}_ATAC_frag_count.txt.gz"
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
        "output/metrics/ATAC_entropy/{comm_panel}/{sample}_ATAC_frag_count.txt.gz"
    output:
        "output/metrics/ATAC_entropy/{comm_panel}/{sample}_ATAC_entropy.txt"
    shell:
        """
        Rscript scripts/calculate_ATAC_entropy.R {input} {output}
        """

rule build_feature_tables:
    input:
        # all of the .SE.tsv files for each 
        expand("output/metrics/se/{comm_panel}/{sample}.SE.tsv", sample = SAMPLES, comm_panel = COMM_PANELS),
        expand("output/metrics/depth/{comm_panel}/{sample}.depth.tsv", sample = SAMPLES, comm_panel = COMM_PANELS),
        expand("output/metrics/frag_bins/{comm_panel}/{sample}.fragbins.tsv", sample = SAMPLES, comm_panel = COMM_PANELS),
        expand("output/metrics/small_frags/{comm_panel}/{sample}.smallfrag.tsv", sample = SAMPLES, comm_panel = COMM_PANELS),
        expand("output/metrics/full_gene_depth/{comm_panel}/{sample}.fullgenedepth.tsv", sample = SAMPLES, comm_panel = COMM_PANELS),
        expand("output/metrics/mds/{comm_panel}/{sample}_mds.txt", sample = SAMPLES, comm_panel = COMM_PANELS),
        expand("output/metrics/TFBS_entropy/{comm_panel}/{sample}_TFBS_entropy.txt", sample = SAMPLES, comm_panel = COMM_PANELS),
        expand("output/metrics/ATAC_entropy/{comm_panel}/{sample}_ATAC_entropy.txt", sample = SAMPLES, comm_panel = COMM_PANELS)
    output:
        "output/feature_tables/{comm_panel}/se.rds",
        "output/feature_tables/{comm_panel}/depth.rds",
        "output/feature_tables/{comm_panel}/frag_bins.rds",
        "output/feature_tables/{comm_panel}/small_frags.rds",
        "output/feature_tables/{comm_panel}/full_gene_depth.rds",
        "output/feature_tables/{comm_panel}/mds.rds",
        "output/feature_tables/{comm_panel}/TFBS_entropy.rds",
        "output/feature_tables/{comm_panel}/ATAC_entropy.rds"
    shell:
        """
        Rscript scripts/build_feature_tables.R
        """


rule extract_first_exon:
    input:
         expand("output/feature_tables/{comm_panel}/se.rds", comm_panel = COMM_PANELS),
         expand("output/feature_tables/{comm_panel}/depth.rds", comm_panel = COMM_PANELS),
         expand("output/feature_tables/{comm_panel}/frag_bins.rds", comm_panel = COMM_PANELS),
         expand("output/feature_tables/{comm_panel}/small_frags.rds", comm_panel = COMM_PANELS),
         expand("output/feature_tables/{comm_panel}/full_gene_depth.rds", comm_panel = COMM_PANELS),
         expand("output/feature_tables/{comm_panel}/mds.rds", comm_panel = COMM_PANELS),
         expand("output/feature_tables/{comm_panel}/TFBS_entropy.rds", comm_panel = COMM_PANELS),
         expand("output/feature_tables/{comm_panel}/ATAC_entropy.rds", comm_panel = COMM_PANELS)
    output:
        "output/feature_tables/{comm_panel}/se_E1.rds",
        "output/feature_tables/{comm_panel}/depth_E1.rds",
        "output/feature_tables/{comm_panel}/small_frags_E1.rds",
        "output/feature_tables/{comm_panel}/mds_E1.rds"
    shell:
        """
        Rscript scripts/filter_first_exons.R
        """

rule build_combined_ft:
    input:
        expand("output/feature_tables/{comm_panel}/se.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/{comm_panel}/depth.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/{comm_panel}/frag_bins.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/{comm_panel}/small_frags.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/{comm_panel}/full_gene_depth.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/{comm_panel}/mds.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/{comm_panel}/TFBS_entropy.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/{comm_panel}/ATAC_entropy.rds", comm_panel = COMM_PANELS)
    output:
        "output/feature_tables/{comm_panel}/all_combined.rds"
    shell:
        """
        Rscript scripts/generate_combined_ft.R
        """

rule move_comm_panels:
    input:
        expand("output/feature_tables/{comm_panel}/se_E1.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/{comm_panel}/depth_E1.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/{comm_panel}/small_frags_E1.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/{comm_panel}/mds_E1.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/{comm_panel}/all_combined.rds", comm_panel = COMM_PANELS)
    output:
        expand("output/feature_tables/grail/{comm_panel}/se.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/depth.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/frag_bins.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/small_frags.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/full_gene_depth.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/mds.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/TFBS_entropy.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/ATAC_entropy.rds", comm_panel = COMM_PANELS),

        expand("output/feature_tables/grail/{comm_panel}/se_E1.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/depth_E1.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/small_frags_E1.rds", comm_panel = COMM_PANELS),
        expand("output/feature_tables/grail/{comm_panel}/mds_E1.rds", comm_panel = COMM_PANELS),

        expand("output/feature_tables/grail/{comm_panel}/all_combined.rds", comm_panel = COMM_PANELS)
    shell:
        """
        mv output/feature_tables/default/* output/feature_tables/grail/default/;
        mv output/feature_tables/tempus/* output/feature_tables/grail/tempus/;
        mv output/feature_tables/guardant/* output/feature_tables/grail/guardant/;
        mv output/feature_tables/foundationOne/* output/feature_tables/grail/foundationOne/
        """