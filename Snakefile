# GDS pipe
# Jules GILET <jules.gilet@inserm.fr>

configfile: "config.yaml"
print("GDS is : ", config["gse"])


from numpy import loadtxt
def get_srr_list(wildcards):
        lst=loadtxt("SRRlist.txt", dtype=str).tolist()
        return lst
SRRs=get_srr_list("srr")
print("SRR ids are: ", SRRs)


rule all:
	input:
		"SRRlist.txt",
                expand("data/{srr}_S1_L001_R1_001.fastq.gz", srr=SRRs),
                expand("data/{srr}_S1_L001_R2_001.fastq.gz", srr=SRRs),
		expand("outs/seurat/{srr}.seurat", srr=SRRs),
		"outs/seurat/integrated/UMAPs_"+config["gse"]+".pdf"

rule get_SRRs:
	input: "config.yaml"
	params: config["gse"]
	output: "SRRlist.txt"
	threads: 1
	shell:
		"""
		echo "Config file is:"
		ls {input}
		echo "Fetching SRRs..."
		SRX=$(esearch -db gds -query {params} | efetch -format native | grep SRX | cut -d "=" -f 2)
		SRR=$(for A in ${{SRX[@]}}; do esearch -db sra -query ${{A}} | efetch -format runinfo | grep -e "SRR" | cut -d "," -f 1; done)
		echo ${{SRR}} > {output}
		"""

rule download_SRAs:
	input: 
		file="SRRlist.txt",
		id=get_srr_list
	output: temp("data/NCBI/{srr}/{srr}.sra")
	threads: 1
	shell:
		"""
		prefetch --max-size 100GB {input.id} && vdb-validate {input.id}
		"""

rule dump_fq:
	input: "data/NCBI/{srr}/{srr}.sra"
	output:
		R1=temp("data/{srr}_1.fastq"),
		R2="data/{srr}_2.fastq",
		R3="data/{srr}_3.fastq"
	threads: 1
	shell:
		"""
		cd data
		fastq-dump --split-files */*.sra
		cd ..
		"""

rule compress_fq:
	input:
                R1="data/{srr}_1.fastq",
                R2="data/{srr}_2.fastq",
                R3="data/{srr}_3.fastq"

	output:
                CR1="data/{srr}_2.fastq.gz",
                CR2="data/{srr}_3.fastq.gz"
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
		R2="data/{srr}_2.fastq.gz",
                R3="data/{srr}_3.fastq.gz"
	output:
		CRFQ1=protected("data/{srr}_S1_L001_R1_001.fastq.gz"),
		CRFD2=protected("data/{srr}_S1_L001_R2_001.fastq.gz")
	threads: 1
	shell:
		"""
		FILE1={input.R2}
		FILE2={input.R3}
		mv ${{FILE1}} ${{FILE1//_2/_S1_L001_R1_001}}
		mv ${{FILE2}} ${{FILE2//_3/_S1_L001_R2_001}}
		"""

FASTQ=glob_wildcards("data/{smp}_S1_L001_R2_001.fastq.gz").smp
print("Available fastq are: ", FASTQ)
rule check_fq:
	input:
		R1="data/{srr}_S1_L001_R1_001.fastq.gz",
		R2="data/{srr}_S1_L001_R2_001.fastq.gz"

	output:
		"data/{srr}.fqcheck"
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

CHECK=glob_wildcards("data/{chk}.fqcheck").chk
print("Checked samples are: ", CHECK)
rule cr_mapping:
	input:
		"data/{srr}.fqcheck"
	output:
		protected("data/cr_{srr}/outs/filtered_feature_bc_matrix.h5")

	log:
		"crmap_{srr}.log"
	params:
		path = config["cr_path"],
		ref = config["ref"],
	threads: 16
	shell:
		"""
		cd data/
		rm -fr cr_{wildcards.srr}
		{params.path} count --id=cr_{wildcards.srr} --transcriptome={params.ref} --fastqs=./ --sample={wildcards.srr} &> {log}
		cd ..
		"""

MAPPED=glob_wildcards("data/cr_{map}/outs/filtered_feature_bc_matrix.h5").map
print("Mapped samples are: ", MAPPED)
rule preprocess:
	input:
		"data/cr_{srr}/outs/filtered_feature_bc_matrix.h5"
	output:
		file="outs/seurat/{srr}.seurat",
		plot="outs/seurat/UMAP_{srr}.pdf",
		log="outs/seurat/ppseurat_{srr}.log"
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
		log = "outs/seurat/integrated/integration_"+config["gse"]+".log"
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

