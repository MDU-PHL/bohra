
import pathlib
configfile: 'config.yaml'
localrules: all,combine_seqdata, qc_snippy, index_reference, calculate_iqtree_command_core,combine_assembly_metrics,assembly_statistics,collate_report,write_html_report

def get_collation_input(pipeline):
	
	s = "{sample}/seqdata.toml"
	k = "{sample}/kraken.toml"
	output = [s,k]
	a = [
		"{sample}/mlst.toml",
		"{sample}/assembly.toml",
		"{sample}/resistome.toml",
		"{sample}/prokka.toml"
		]
	s = [
		"{sample}/snippy.toml",
		"{sample}/snippy_qc.toml", 
		]
		
	if pipeline == 'sa' or pipeline == 'all':
		for x in s:
			output.append(x)
		for y in a:
			output.append(y)
	elif pipeline == 'a':
		for y in a:
			output.append(y)
	elif pipeline == 's':
		for x in s:
			output.append(x)
	
	return output


def get_report_tomls(pipeline):
	output = ["seqdata.toml",
		"kraken.toml"]
	a = ["mlst.toml", 
		"assembly.toml", 
		"resistome.toml"]
	s = ['gubbins.toml', 'snippy_core.toml', 'iqtree.toml', 'distances.toml']
	r = ["roary.toml",
		'pan_genome.toml']
	
	if pipeline == 'sa':
		for x in s:
			output.append(x)
		for y in a:
			output.append(y)
	elif pipeline == 'all':
		for x in s:
			output.append(x)
		for y in a:
			output.append(y)
		for z in r:
			output.append(z)
	elif pipeline == 'preview':
		output = ['report.toml', 'preview.toml']
	
	return output

def final_output(tomls):

	output = [f"{pathlib.Path('report', t)}" for t in tomls]
	output.append(f"{pathlib.Path('report', 'index.html')}")
	# if 'report.html' not in output:
	# 	output.append(f"{pathlib.Path('report', 'report.toml')}")
	return output

PREFILLPATH = config['prefill_path']
SAMPLE = config['isolates'].split()
SNIPPY_SINGULARITY = config['snippy_singularity']
# ASSEMBLER_SINGULARITY = config['assembler_singularity']
ABRITAMR_SINGULARITY = config['abritamr_singularity']
MIN_ALN = int(config['min_perc'])
REFERENCE = config['reference']
GUBBINS = config['gubbins']
PIPELINE = config['pipeline']
ALL_TOMLS = get_collation_input(pipeline = PIPELINE)
print(ALL_TOMLS)
REPORT_TOMLS = get_report_tomls(pipeline = PIPELINE)
print(REPORT_TOMLS)
FINAL_OUTPUT = final_output(tomls = REPORT_TOMLS)
print(FINAL_OUTPUT)
WORKDIR = config['workdir']
JOB_ID = config['job_id']
ASSEMBLER = config['assembler']
TEMPLATE_PATH = config['template_path']
SCRIPT_PATH = config['script_path']
MASK_STRING = config['mask_string'] if config['mask_string'] != '' else 'nomask'
PREVIEW = config['preview']	
KRAKEN_DB=config['kraken_db']
KRAKEN_THREADS = int(config['kraken_threads'])
MIN_COV = config['min_cov']
print(f"Kraken threads {KRAKEN_THREADS}")
rule all:
	input:
		FINAL_OUTPUT,
		expand("{sample}/final.toml", sample = SAMPLE)

rule estimate_coverage:
	input:
		r1="{sample}/R1.fq.gz",
		r2="{sample}/R2.fq.gz"
	output:
		"{sample}/mash.toml"
	params:
		script_path = SCRIPT_PATH
	
	script:
		"mash.py"

rule seqdata:
	input:
		r1 = '{sample}/R1.fq.gz',
		r2 = '{sample}/R2.fq.gz',
		# mash = '{sample}/mash.toml'
	output:
		"{sample}/seqdata.toml"
	params:
		script_path=SCRIPT_PATH,
		mincov = MIN_COV,
		reference = REFERENCE
	
	script:
		"seqdata.py"

rule combine_seqdata:
	input:
		expand("{sample}/seqdata.toml", sample = SAMPLE)
	output:
		"seqdata.toml"
	params:
		script_path=SCRIPT_PATH
	script:
		"combine_seqdata.py"

if PREVIEW:
	rule preview:
		input:
			expand("{sample}/mash.toml", sample = SAMPLE)
		output:
			"preview.toml"
		params:
			reference = REFERENCE,
			script_path = SCRIPT_PATH
		script:
			"preview.py"

	rule combine_preview_tomls:
		input:
			"{sample}/mash.toml"
		output:
			"{sample}/final.toml"
		shell:
			"""
			cp {input} {output}
			"""
	rule compile:
		input:
			"preview.toml"
		output: #this is where I stopped
			"report.toml", "report.html"
		params:
			pipeline = 'preview',
			script_path = SCRIPT_PATH,
			job_id = JOB_ID,
			template_path = TEMPLATE_PATH,
			assembler = ASSEMBLER
		script:
			"compile.py"
else:
	rule kraken:
		input: 
			r1='{sample}/R1.fq.gz',
			r2='{sample}/R2.fq.gz'
		output: 
			'{sample}/kraken.toml'
		params:
			prefill_path = PREFILLPATH,
			kraken_db = KRAKEN_DB,
			script_path = SCRIPT_PATH,
		threads:
			KRAKEN_THREADS
		script:	
			"kraken.py"

	rule combine_kraken:
		input:
			expand("{sample}/kraken.toml", sample = SAMPLE)
		output:
			"kraken.toml"
		params:
			script_path = SCRIPT_PATH
		script:
			"combine_kraken.py"
	
	rule snippy:
		input:
			'{sample}/seqdata.toml'	
		output:
			'{sample}/snippy.toml'
		threads:
			8
		singularity: SNIPPY_SINGULARITY
		params:
			script_path=SCRIPT_PATH,
			reference = REFERENCE
		script:
			"snippy.py"		

	rule qc_snippy: 
		input:
			'{sample}/snippy.toml'
			
		output:
			'{sample}/snippy_qc.toml'
		params:
			script_path = SCRIPT_PATH,
			minaln = MIN_ALN
		script:
			"snippy_qc.py"

	rule run_snippy_core:
		input:
			expand("{sample}/snippy_qc.toml", sample = SAMPLE)
		output:
			'snippy_core.toml'
		singularity: SNIPPY_SINGULARITY
		params:
			mask_string = MASK_STRING,
			script_path = SCRIPT_PATH,
			reference = REFERENCE,
			
		script:
			"snippy_core.py"

	rule run_gubbins:
		input:
			'snippy_core.toml'
		output:
			'gubbins.toml'
		params:
			script_path = SCRIPT_PATH,	
			gubbins = GUBBINS
		script:
			"gubbins.py"

	rule run_snpdists:
		input:
			'gubbins.toml'
		output:
			'distances.toml' 
		params:
			script_path = SCRIPT_PATH
		singularity: SNIPPY_SINGULARITY
		script:
			"snp_dists.py"
		

	# rule index_reference:
	# 	input:
	# 		REFERENCE
	# 	output:
	# 		"ref.fa",
	# 		"ref.fa.fai"
	# 	run:
	# 		from Bio import SeqIO
	# 		import pathlib, subprocess
	# 		ref = f"{output[0]}"
	# 		idx = f"{output[1]}"
	# 		if '.fa' not in REFERENCE:
	# 			# print(f"converting {REFERENCE}")
	# 			SeqIO.convert(f"{input[0]}", 'genbank', ref	, 'fasta')
	# 			# print(f"converted {REFERENCE}")
	# 		else:
	# 			subprocess.run(f"ln -sf {REFERENCE} {ref}", shell = True)
	# 		subprocess.run(f"samtools faidx {ref}", shell =True)

	rule run_iqtree_core:
		input:
			gubbins = 'gubbins.toml', 
			# ref = 'ref.fa', 
			# idx = 'ref.fa.fai'
		
		output:
			'iqtree.toml',
		params:
			script_path = SCRIPT_PATH,
			ref = "ref.fa",
			idx = "ref.fa.fai"
		script:
			"run_iqtree.py"
	
	rule assemble:
		input:
			'{sample}/seqdata.toml'
		output:
			'{sample}/assembly.toml'
		threads:
			16
		params:
			prefill_path = PREFILLPATH,
			assembler = ASSEMBLER, 
			script_path = SCRIPT_PATH
		script:
			"assemble.py"

	rule assembly_statistics:
		input:
			"{sample}/assembly.toml"
		output:
			"{sample}/assembly_stats.toml"
		params:
			script_path = SCRIPT_PATH,
			minsize = 500,
			prefill_path = PREFILLPATH
		script:
			"assembly_stat.py"
		
	rule run_prokka:
		input:
			assembly = "{sample}/assembly_stats.toml",
			seqdata = "{sample}/seqdata.toml"
		output:
			"{sample}/prokka.toml"
		params:
			script_path = SCRIPT_PATH
		threads: 4
		script:
			"prokka.py"
	

	rule combine_assembly_metrics:
		input:
			prokka = expand("{sample}/prokka.toml",sample = SAMPLE), 
			# assembly = expand("{sample}/assembly_stats.toml",sample = SAMPLE)
		output:
			"assembly.toml"
		params:
			script_path= SCRIPT_PATH
		script:
			"assembly_combine.py"


	rule resistome:
		input:
			assembly = "{sample}/assembly.toml",
			seqdata = '{sample}/seqdata.toml'
		output:
			'{sample}/resistome.toml'
		params:
			script_path= SCRIPT_PATH,
			work_dir = WORKDIR,
			job_id = JOB_ID
		singularity:ABRITAMR_SINGULARITY
		script:
			"resistome.py"

	rule combine_resistome:
		input:
			expand("{sample}/resistome.toml", sample = SAMPLE)
		output:
			'resistome.toml'
		params:
			script_path=SCRIPT_PATH
		script:
			"combine_resistome.py"

	rule mlst:
		input:
			assembly = "{sample}/assembly.toml",
			seqdata = "{sample}/seqdata.toml"
		output:
			'{sample}/mlst.toml'
		params:
			script_path=SCRIPT_PATH
		script:
			"mlst.py"

	rule combine_mlst:
		input:
			expand('{sample}/mlst.toml', sample = SAMPLE)
		output:
			"mlst.toml"
		params:
			script_path = SCRIPT_PATH
		script:
			"combine_mlst.py"	
	

	if PIPELINE == 'all':
	
		rule run_roary:
			input:
				prokka = expand("{sample}/prokka.toml",sample = SAMPLE)
			output:
				"roary.toml"
			params:
				script_path = SCRIPT_PATH
			script:
				"roary.py"

		rule pan_figure:
			input:
				"roary.toml"
			output:
				"pan_genome.toml"
			params:
				script_path = SCRIPT_PATH
			script:
				"pan_figure.py"

	rule collate_isolate_tomls:
		input:
			ALL_TOMLS
		output:
			"{sample}/final.toml"
		params:
			script_path = SCRIPT_PATH
		script:	
			"collate_tomls.py"
		
	rule compile_report_toml:
		input:
			REPORT_TOMLS
		output:
			"report.html", 
			"report.toml"
		params:
			pipeline = PIPELINE,
			template_path = TEMPLATE_PATH,
			assembler = ASSEMBLER, 
			job_id = JOB_ID
		script:
			"compile.py"


rule move_outputs:
	input:
		"report.html",
		REPORT_TOMLS
	output:
		FINAL_OUTPUT
	params:
		script_path= SCRIPT_PATH
	script:
		"move_outputs.py"