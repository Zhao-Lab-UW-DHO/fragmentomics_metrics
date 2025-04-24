
SAMPLES_PATH = "files/grail_samples.txt"

def samples_to_list(sample_file):
	with open(sample_file) as f:
		samples = f.read().splitlines()
	return samples

# sample names without the "raw.bam" at the end of the file name
SAMPLES = samples_to_list(SAMPLES_PATH)

rule all:
	input:
		expand("data/{sample}.bed.gz", sample=SAMPLES)

rule bam2bed:
	input:
		bam = "data/grail_bams/{sample}raw.bam"
	output:
		"data/{sample}.bed.gz"
	threads:
		8
	shell:
		"""
		{input.sambamba} sort -n -t {threads} -m 20G -o /dev/stdout {input.bam} \
		| {input.samtools} view -f3 -F2308 -b - \
		| {input.bedtools} bamtobed -bedpe -i stdin \
		| cut -f 1,2,6- \
		| awk -F "\t" 'BEGIN {{ OFS=FS }}; {{ if($6 != $7) {{ print $1,$2,$3,$4=substr($4, length($4)-12) }} }}' \
		| {input.bedops_sort} --max-mem 20G --unique - \
		| uniq \
		| gzip \
		> {output}
		"""
