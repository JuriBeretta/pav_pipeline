	configfile:"config.yaml"
genome=config['genome']
threads_tot=config['threads_config']
samples=[]
with open('sample_list.txt') as f:
        samples = f.read().splitlines()

rule all:
	input:
		expand("{sample}_1.fastq", sample=sample),
        expand("{sample}_2.fastq", sample=sample),
        expand("{sample}_1.fastq.gz", sample=sample),
        expand("{sample}_2.fastq.gz", sample=sample)
		expand("trimmed/{sample}_{read}.trimmed.fastq.gz", sample=samples, read=['1', '2']),
		expand("coverage/{sample}.sorted.mapping.depth", sample=samples),
		#expand("coverage/{sample}.sorted.mapping.depth.d",)		
		expand("coverage/{sample}.sorted.mapping.depth.d/{sample}.sorted.mapping.depth.coverage.png", sample=samples),
		expand("coverage/{sample}.sorted.mapping.depth.d/{sample}.sorted.mapping.depthgenes.coverage.csv", sample=samples),
		expand("coverage/{sample}.sorted.mapping.depth.d/{sample}.sorted.mapping.depth.absent.list", sample=samples)


rule download_reads:
    output:
        reads_1="{sample}_1.fastq",
        reads_2="{sample}_2.fastq"


    
        


rule trimming:
	input:
                read1="{sample}_1.fastq.gz",
                read2="{sample}_2.fastq.gz"
	output:
		trimmed1="trimmed/{sample}_1.trimmed.fastq.gz",
		trimmed2="trimmed/{sample}_2.trimmed.fastq.gz"
	group: "pipeline"
	threads: threads_tot
	conda:
		"fastp"

	shell: "/home/thanatos/miniforge3/envs/fastp/bin/fastp -i {input.read1} -I {input.read2} -o {output.trimmed1} -O {output.trimmed2} --detect_adapter_for_pe -V -w {threads} -x -g -n 2 -5 -3 -p -l 75 -M 24"

rule mapping:
	input:
		genome=expand("{genome}", genome = genome),
		trimmed1="trimmed/{sample}_1.trimmed.fastq.gz",
		trimmed2="trimmed/{sample}_2.trimmed.fastq.gz"
	output:
		mapping="mapped/{sample}.mapping.bam"
	group: "pipeline"
	threads: threads_tot
	conda:
		"bwa"
	shell:
		"/home/thanatos/miniforge3/envs/bwa/bin/bwa mem -t {threads} -M {input.genome} {input.trimmed1} {input.trimmed2} | /home/thanatos/miniforge3/envs/bwa/bin/samtools view -bS - > {output.mapping}"

rule sorting:
	input:
		mapping="mapped/{sample}.mapping.bam"
	output:
		sorted="mapped/{sample}.sorted.mapping.bam"
	threads: threads_tot
	conda:
		"bwa"
	shell:
		"/home/thanatos/miniforge3/envs/bwa/bin/samtools sort -@ {threads} -O bam -o {output.sorted} {input.mapping}"
rule depth:
	input:
		sorted="mapped/{sample}.sorted.mapping.bam"

	output:
		coverage="coverage/{sample}.sorted.mapping.depth"

	group: "pipeline"

	conda:
		"bwa"
	shell:
		"/home/thanatos/miniforge3/envs/bwa/bin/samtools depth -aa {input.sorted} > {output.coverage}"

rule analysis:
	input:
		coverage="coverage/{sample}.sorted.mapping.depth"
	output:
		results_dir = directory("coverage/{sample}.sorted.mapping.depth.d"),
		results_fig= "coverage/{sample}.sorted.mapping.depth.d/{sample}.sorted.mapping.depth.coverage.png",
		results_tab = "coverage/{sample}.sorted.mapping.depth.d/{sample}.sorted.mapping.depthgenes.coverage.csv",
		resultas_lista = "coverage/{sample}.sorted.mapping.depth.d/{sample}.sorted.mapping.depth.absent.list"
		#tempdir = "{sample}.sorted.mapping.depth.d"

	group: "pipeline"

	shell:
		"python coverage.py ./exons_coordinates_transcriptome.tsv ./complete_busco_coordinates_transcriptome.tsv {input.coverage}"
