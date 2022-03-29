# GDS pipe
# Jules GILET <jules.gilet@inserm.fr>

configfile: "config.yaml"
print("GDS is : ", config["gse"])

import yaml
with open('config.yaml') as file:
	try:
		CONFIG = yaml.safe_load(file)
	except yaml.YAMLError as mess:
		print(mess)


import os
if os.path.exists("SRRlist.txt"):
	print("SRRlist.txt file is present")
else:
	print("Generating SRR file list...")
	os.system("scripts/download.sh "+CONFIG["gse"]+" SRRlist.txt")


from numpy import loadtxt
def get_srr_list(wildcards):
        lst=loadtxt("SRRlist.txt", dtype=str).tolist()
        return lst
SRRs=get_srr_list("srr")
print("SRR ids are: ", SRRs)


rule all:
	input:
                expand("fastq/{srr}_S1_L001_R1_001.fastq.gz", srr=SRRs),
                expand("fastq/{srr}_S1_L001_R2_001.fastq.gz", srr=SRRs),
		expand("outs/seurat/{srr}.seurat", srr=SRRs),
		"outs/seurat/integrated/UMAPs_"+config["gse"]+".pdf"


rule dump_fq:
	input:
		file="SRRlist.txt"
	output:
		R1="fastq/{srr}_1.fastq",
		R2="fastq/{srr}_2.fastq",
		R3="fastq/{srr}_3.fastq"
	threads: 1
	shell:
		"""
		fastq-dump --split-files {wildcards.srr} --outdir fastq/
		"""


rule compress_fq:
	input:
                R1="fastq/{srr}_1.fastq",
                R2="fastq/{srr}_2.fastq",
                R3="fastq/{srr}_3.fastq"

	output:
                CR1="fastq/{srr}_2.fastq.gz",
                CR2="fastq/{srr}_3.fastq.gz"
	threads: 1
	message: "Compressing fastqs, that might take a while..."
	shell:
		"""
		rm {input.R1}
		gzip {input.R2}
		gzip {input.R3}
		"""


rule namefix_fq:
	input:
		R2="fastq/{srr}_2.fastq.gz",
                R3="fastq/{srr}_3.fastq.gz"
	output:
		CRFQ1=protected("fastq/{srr}_S1_L001_R1_001.fastq.gz"),
		CRFD2=protected("fastq/{srr}_S1_L001_R2_001.fastq.gz")
	threads: 1
	shell:
		"""
		FILE1={input.R2}
		FILE2={input.R3}
		mv ${{FILE1}} ${{FILE1//_2/_S1_L001_R1_001}}
		mv ${{FILE2}} ${{FILE2//_3/_S1_L001_R2_001}}
		"""


FASTQ=glob_wildcards("fastq/{smp}_S1_L001_R2_001.fastq.gz").smp
print("Available fastq are: ", FASTQ)
rule check_fq:
	input:
		R1="fastq/{srr}_S1_L001_R1_001.fastq.gz",
		R2="fastq/{srr}_S1_L001_R2_001.fastq.gz"

	output:
		"fastq/{srr}.fqcheck"
	threads: 1
	shell:
		"""
		if [ -f {input.R1} ] && [ -f {input.R2} ]
		then
			if [ -s {input.R1} ] && [ -s {input.R2} ]
			then
				touch {output}
			fi
		fi
		"""


CHECK=glob_wildcards("fastq/{chk}.fqcheck").chk
print("Checked samples are: ", CHECK)
rule cr_mapping:
	input:
		"fastq/{srr}.fqcheck"
	output:
		protected("outs/cellranger/cr_{srr}/outs/filtered_feature_bc_matrix.h5")

	log:
		"outs/cellranger/crmap_{srr}.log"
	params:
		path = config["cr_path"],
		ref = config["ref"],
	threads: 16
	shell:
		"""
		cd outs/cellranger
		rm -fr cr_{wildcards.srr}
		{params.path} count --id=cr_{wildcards.srr} --transcriptome={params.ref} --fastqs=../../fastq --sample={wildcards.srr} &> {log}
		cd ../..
		"""


MAPPED=glob_wildcards("outs/cellranger/cr_{map}/outs/filtered_feature_bc_matrix.h5").map
print("Mapped samples are: ", MAPPED)
rule preprocess:
	input:
		"outs/cellranger/cr_{srr}/outs/filtered_feature_bc_matrix.h5"
	output:
		file="outs/seurat/{srr}.seurat",
		plot="outs/seurat/UMAP_{srr}.pdf",
		log="outs/seurat/ppseurat_{srr}.log"
	params:
		meta=config["meta_select"]
	threads: 16
	script:
		"scripts/postprocess.R"


PPED=glob_wildcards("outs/seurat/{pp}.seurat").pp
print("Preprocessed samples are: ", PPED)
print("GSE is :", config["gse"])
rule integrate:
	input:
		list = expand("outs/seurat/{pp}.seurat", pp=FASTQ)
	output:
		file = "outs/seurat/integrated/seu_int_"+config["gse"]+".seurat",
		plot = "outs/seurat/integrated/UMAPs_"+config["gse"]+".pdf",
		log = "outs/seurat/integrated/integration_"+config["gse"]+".log",
		matrix = protected("outs/seurat/integrated/readCountMatrix_"+config["gse"]+".csv")
	threads: 32
	params:
		gse = config.get("gse")
	script:
		"scripts/merging.R"


onsuccess:
	print("Succesfully downloaded and prepared fastq files.")
	print("Succefully preprocessed and integrated datasets.")

onerror:
	print("WARNING: an error occured somewhere...")
	shell("{log} > ERROR.log")

