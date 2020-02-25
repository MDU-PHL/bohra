workdir:{% raw %}{% endraw %} '{{workdir}}'{% raw %}
configfile: 'config.yaml'
localrules: all, generate_yield, combine_seqdata, qc_snippy, index_reference, calculate_iqtree_command_core,combine_assembly_metrics,assembly_statistics,collate_report,write_html_report

def get_final_output(pipeline):
	output = [
		# expand("{sample}/seqdata.tab", sample = sample_list),
		"report/seqdata.tab",
		"report/index.html",
		"excluded_isolates.txt",
		"core_isolates.txt",
		"species_identification.tab",
		"report/species_identification.tab"
	]
		# expand("{sample}/kraken.tab",sample = sample_list)]
	a = [
		# expand("{sample}/{sample}.fa", sample = sample_list),
		# expand("{sample}/resistome.tab", sample = sample_list),
		# expand("{sample}/{sample}.gff", sample = sample_list),
		# expand("{sample}/{sample}.txt", sample = sample_list),
		"mlst.tab", 
		"denovo.tab", 
		"assembly.tab", 
		"summary_matches.csv",
		"report/assembly.tab",
		"report/mlst.tab", 
		"report/summary_matches.tab"]
	s = ["ref.fa",
		"ref.fa.fai",
		# expand("{sample}/snps.vcf", sample = sample_list),
		# expand("{sample}/snps.aligned.fa", sample = sample_list),
		"core.vcf", 
		"distances.tab",
		"core.treefile", 
		"report/core_genome.tab", 
		"report/core.treefile", 
		"report/distances.tab",
		"report/core.tab"]
	r = ["pan_genome.svg", 
		"report/pan_genome.svg"]
	
	if pipeline == 'sa':
		output.extend(a)
		output.extend(s)
	elif pipeline == 'a':
		output.extend(a)
	elif pipeline == 's':
		output.extend(s)
	elif pipeline == 'all':
		output.extend(s)
		output.extend(a)
		output.extend(r)
	
	return output

def get_report_input(pipeline):
	output = ["seqdata.tab"
		"species_identification.tab"]
	a = ["mlst.tab", 
		"denovo.tab", 
		"assembly.tab", 
		"summary_matches.csv"]
	s = ['core.txt', 
		'core.treefile', 
		'core.tab', 
		'distances.tab']
	r = ["pan_genome.svg",
		'roary/summary_statistics.txt']
	
	if pipeline == 'sa':
		output.extend(a)
		output.extend(s)
	elif pipeline == 'a':
		output.extend(a)
	elif pipeline == 's':
		output.extend(s)
	elif pipeline == 'all':
		output.extend(s)
		output.extend(a)
		output.extend(r)
	
	return output


def get_report_output(pipeline):
	output = ["report/seqdata.tab"
		"report/species_identification.tab"]
	a = ["report/mlst.tab", 
		"report/assembly.tab", 
		"report/summary_matches.tab"]
	s = ['report/core_genome.tab', 
		'report/core.treefile', 
		'report/core.tab', 
		'report/distances.tab']
	r = ["report/pan_genome.svg",
		'report/summary_statistics.txt']
	
	if pipeline == 'sa':
		output.extend(a)
		output.extend(s)
	elif pipeline == 'a':
		output.extend(a)
	elif pipeline == 's':
		output.extend(s)
	elif pipeline == 'all':
		output.extend(s)
		output.extend(a)
		output.extend(r)
	
	return output



PREFILLPATH = config['prefill_path']
SAMPLE = config['isolates'].split()
print(SAMPLE)
MIN_ALN = int(config['min_perc'])
REFERENCE = config['reference']
GUBBINS = config['gubbins']
PIPELINE = config['pipeline']
REPORT_INPUT = get_report_input(pipeline = PIPELINE)
REPORT_OUTPUT= get_report_output(pipeline = PIPELINE)
FINAL_OUTPUT = get_final_output(pipeline=PIPELINE)
print(FINAL_OUTPUT)
WORKDIR = config['workdir']
JOB_ID = config['job_id']
ASSEMBLER = config['assembler']
TEMPLATE_PATH = config['template_path']
SCRIPT_PATH = config['script_path']
MASK_STRING = config['mask_string']
		
if GUBBINS:
	CORE_OUTPUT = 'gubbins.aln'
else:
	CORE_OUTPUT = 'core.aln'
	
rule all:
	input:
		# expand('{sample}/seqdata.tab', sample = SAMPLE ),
		# FINAL_OUTPUT
		expand("{sample}/seqdata.tab", sample = SAMPLE),
		"report/seqdata.tab",
		"report/index.html",
		"excluded_isolates.txt",
		"core_isolates.txt",
		"species_identification.tab",
		"report/species_identification.tab",
		expand("{sample}/{sample}.fa", sample = sample_list),
		expand("{sample}/resistome.tab", sample = sample_list),
		expand("{sample}/{sample}.gff", sample = sample_list),
		expand("{sample}/{sample}.txt", sample = sample_list),
		"mlst.tab", 
		"denovo.tab", 
		"assembly.tab", 
		"summary_matches.csv",
		"report/assembly.tab",
		"report/mlst.tab", 
		"report/summary_matches.tab",
		"ref.fa",
		"ref.fa.fai",
		expand("{sample}/snps.vcf", sample = sample_list),
		expand("{sample}/snps.aligned.fa", sample = sample_list),
		"core.vcf", 
		"distances.tab",
		"core.treefile", 
		"report/core_genome.tab", 
		"report/core.treefile", 
		"report/distances.tab",
		"report/core.tab"
		

rule run_kraken:
	input:
		'{sample}/R1.fq.gz',
		'{sample}/R2.fq.gz'
	output:
		"{sample}/kraken.tab"
	params:
		prefill_path = PREFILLPATH
	shell:
		"""
		KRAKENPATH={params.prefill_path}/{wildcards.sample}/kraken2.tab
		if [ -f $KRAKENPATH ]; then
			cp $KRAKENPATH {output}
		else
			kraken2 --paired {input[0]} {input[1]} --minimum-base-quality 13 --report {output} --memory-mapping
		fi
		"""

rule seqdata:
	input:
		'{sample}/R1.fq.gz',
		'{sample}/R2.fq.gz'
	output:
		"{sample}/seqdata.tab"
	singularity:{% endraw %}"{{singularity_dir}}/seqtk"{% raw %}
	shell:
		"""
		seqtk fqchk {input[0]} {input[1]} > {output}
		"""


rule estimate_coverage:
	input:
		"{sample}/R1.fq.gz",
		"{sample}/R2.fq.gz"
	output:
		"{sample}/mash.txt"
	singularity:{% endraw %}"{{singularity_dir}}/mash_kmc"{% raw %}
	shell:
		"""
		mash sketch -r {input[0]} {input[1]} -m 3 -k 31 -o mash  &> {output}
		"""

rule generate_yield:
	input:
		"{sample}/mash.txt",
		"{sample}/seqdata.tab"
	output:
		"{sample}/yield.tab"
	params:
		script_path = SCRIPT_PATH
	shell:
		"""
		python3 {params.script_path}/generate_yield.py {input[1]} {input[0]} {output}
		"""

rule combine_seqdata:
	input:
		expand("{sample}/yield.tab", sample = SAMPLE)
	output:
		"seqdata.tab"
	run:
		import pathlib, pandas, numpy
		sdfiles = f"{input}".split()
		seqdata = pandas.DataFrame()
		for sd in sdfiles:
			p = pathlib.Path(sd)
			df = pandas.read_csv(sd, sep = "\t")
			print(df)
			df['Isolate'] = f"{p.parts[0]}"
			
			if seqdata.empty:
				seqdata = df
			else:
				seqdata = seqdata.append(df, sort = True)
		seqdata['Quality'] = numpy.where(seqdata['Estimated depth'] >= 40, 'PASS','FAIL')
		seqdata = seqdata[['Isolate','Reads','Yield','GC content','Min len','Avg len','Max len','Avg Qual','Estimated depth', 'Quality']]
		seqdata.to_csv(f"{output}", sep = '\t', index = False)

if PIPELINE != 'a':
	rule snippy:
		input:
			'{sample}/R1.fq.gz',
			'{sample}/R2.fq.gz'
		output:
			'{sample}/snps.vcf',
			'{sample}/snps.aligned.fa'
		threads:
			8
		singularity:{% endraw %}"{{singularity_dir}}/snippy"{% raw %}
		shell:
			"""
			snippy --outdir {wildcards.sample} --ref {REFERENCE} --R1 {input[0]} --R2 {input[1]} --force --cpus {threads}
			"""
		

	rule qc_snippy: 
		input:
			expand('{sample}/snps.aligned.fa', sample = SAMPLE)
			
		output:
			'core_isolates.txt'
			
		run:
			from Bio import SeqIO
			import pathlib
			import pandas
			# create an output
			isolate_list = []
			excluded_list = []
			outfile = pathlib.Path(f"{output[0]}")
			# get input file list
			input_list = f"{input}".split()
			# set the log path
			logpath = pathlib.Path('isolates.log')
			for i in input_list: # for each input file
				# get the isolate name
				p = pathlib.Path(f"{i}")
				isolate = p.parts[-2]
				if p.exists(): # if the file exists open it
					fasta = p.open()
					for i in SeqIO.parse(fasta,'fasta'): # use BioPython to determine percent alignment
						length = len(i.seq)
						nocov = i.seq.count('-')
						lowcov = i.seq.count('N')
						het = i.seq.count('n')
						unaln = nocov + lowcov + het
						perc_aln = 100*(length - unaln) / length
						# if the percent alignement is greater than the min alignment
						if perc_aln > min_aln:
							isolate_list.append(f"{isolate}")
						else:
							excluded_list.append(isolate)
							print(f"{isolate} has been excluded from the analysis due to poor alignement with reference")
							
			isolate_list = list(set(isolate_list))
			with open(outfile, 'w') as f:
				f.write('\n'.join(isolate_list))
			# get log if the excluded list has any isolates in
			if excluded_list != []:
				if logpath.exists():
					lf = pandas.read_csv(logpath, sep = '	', index_col = False)
					for e in excluded_list:
						lf.loc[lf['Isolate'] == e.strip('#'), 'Status'] = f"(FAILED ALIGNMENT (<{min_aln}% ALIGNMENT))"
						lf.loc[lf['Isolate'] == e.strip('#'), 'Date'] = f"{config['day']}"
						lf.to_csv(logpath, sep = '	', index=False)

		

	rule run_snippy_core:
		input:
			'core_isolates.txt'
		output:
			'core.vcf',
			'core.txt',
			'core.aln', 
			'core.full.aln',
			'core.tab'
		singularity:{% endraw %}"{{singularity_dir}}/snippy"{% raw %}
		params:
			mask_string = MASK_STRING
		shell:
			"""
			snippy-core {params.mask_string} --ref {REFERENCE}  $(cat core_isolates.txt)
			"""

	if GUBBINS:

		rule run_gubbins:
			input:
				'core.full.aln'
			output:
				'clean.full.aln',
				'gubbins.aln'
				
			
			shell:
				"""
				snippy-clean_full_aln {input} > {output[0]}
				run_gubbins.py -c 36  --prefix core {output[0]}
				snp-sites -c core.filtered_polymorphic_sites.fasta > {output[1]}
				"""

	rule run_snpdists:
		input:
			CORE_OUTPUT
		output:
			'distances.tab' 
		singularity:{% endraw %}"{{singularity_dir}}/snippy"{% raw %}
		shell:
			"""
			snp-dists {input} > {output}
			"""
		

	rule index_reference:
		input:
			REFERENCE
		output:
			"ref.fa",
			"ref.fa.fai"
		run:
			from Bio import SeqIO
			import pathlib, subprocess
			ref = f"{output[0]}"
			idx = f"{output[1]}"
			print(type(ref))
			print(type(idx))
			if '.fa' not in REFERENCE:
				print(f"converting {REFERENCE}")
				SeqIO.convert(f"{input[0]}", 'genbank', ref	, 'fasta')
				print(f"converted {REFERENCE}")
			else:
				subprocess.run(f"ln -sf {REFERENCE} {ref}", shell = True)
			subprocess.run(f"samtools faidx {ref}", shell =True)


	rule calculate_iqtree_command_core:
		input:
			CORE_OUTPUT,
			"ref.fa"
		output:
			'run_iqtree_core.sh'
		params:
			script_path = SCRIPT_PATH
		shell:
			"bash {params.script_path}/iqtree_generator.sh {input[1]} {input[0]} core 20 > {output}"

	rule run_iqtree_core:
		input:
			'run_iqtree_core.sh'
		
		output:
			'core.iqtree',
			'core.treefile',
			
		singularity:{% endraw %}"{{singularity_dir}}/iqtree"{% raw %}
		shell:
			"""	
			bash run_iqtree_core.sh
			rm *.ckp.gz *.contree *.bionj
			"""
			
	
if PIPELINE != 's':
	rule assemble:
		input:
			'{sample}/R1.fq.gz',
			'{sample}/R2.fq.gz'
		output:
			'{sample}/{sample}.fa'
		threads:
			16
		params:
			prefill_path = PREFILLPATH
		singularity:{% endraw %}"{{singularity_dir}}/assemblers"{% raw %}
		shell:
			"""
			ASSEMBLEPATH={params.prefill_path}/{wildcards.sample}
			if [ -f $ASSEMBLEPATH/contigs.fa ]; then
				cp $ASSEMBLEPATH/contigs.fa {output}

			else
				echo No assembly found. Assembling {wildcards.sample} with shovill
				
				shovill --outdir {wildcards.sample} --R1 {input[0]} --R2 {input[1]} --force --minlen 500 --cpus {threads}
				mv {wildcards.sample}/contigs.fa {output}

			fi		
			"""
	
	rule collate_assemblies:
		input:
			expand('{sample}/{sample}.fa', sample = SAMPLE)
		output:
			'isolates_abritamr.tab'
		params:
			work_dir=WORKDIR
		run:
			import pathlib
			lines = []
			for i in {input}:
				line = f"{wildcards.sample}\t{params.work_dir}/{i}"
				lines.append(line)
			with open(f"{output}", 'w') as f:
				f.write('\n'.join(lines))

	rule resistome:
		input:
			'isolates_abritamr.tab'
		output:
			'summary_matches.csv'
		# singularity:{% endraw %}"{{singularity_dir}}/abricate"{% raw %}
		shadow:
			"full"
		shell:
			"""
			abriTAMR -c {input}
			"""

	rule mlst:
		input:
			expand('{sample}/{sample}.fa', sample = SAMPLE)
		output:
			'mlst.tab'
		singularity:{% endraw %}"{{singularity_dir}}/mlst"{% raw %}
		shell:
			"""
			mlst --nopath {input} | sed 's/\.fa//g' | sed '1iIsolate\tScheme\tST\tAlleles' > {output}
			"""		

	rule assembly_statistics:
		input:
			expand("{sample}/{sample}.fa", sample = SAMPLE)
		output:
			"denovo.tab"
		params:
			script_path = SCRIPT_PATH
		shell:
			"""
			python3 {params.script_path}/assembly_stat.py {input} -m 500 > {output}
			"""
		
	rule run_prokka:
		input:
			"{sample}/{sample}.fa"
		output:
			"{sample}/{sample}.gff","{sample}/{sample}.txt"
		singularity:{% endraw %}"{{singularity_dir}}/prokka"{% raw %}
		shell:
			"""
			prokka --outdir {wildcards.sample} --prefix {wildcards.sample} --mincontiglen 500 --notrna --fast --force {input} --cpus 1
			rm {wildcards.sample}/*.err {wildcards.sample}/*.faa {wildcards.sample}/*.ffn {wildcards.sample}/*.fsa {wildcards.sample}/*.sqn {wildcards.sample}/*.tbl {wildcards.sample}/*.tsv
			"""
	

	rule combine_assembly_metrics:
		input:
			prokka = expand("{sample}/{sample}.txt",sample = SAMPLE), 
			assembly = "denovo.tab"
		output:
			"assembly.tab"
		run:
			import pandas, pathlib

			prokka = f"{input.prokka}".split()
			gff = pandas.DataFrame()
			
			for p in prokka:
				g = pathlib.Path(p)
				df = pandas.read_csv(g, sep = ':', header = None, names = ['cond', f"{g.parts[1]}"])
				
				if gff.empty:
						gff = df
				else:
						gff = gff.merge(df, how = 'outer')
			gff = gff[gff['cond'].isin(['CDS', 'rRNA'])]
			gff = gff.T
			gff.columns = gff.iloc[0]
			gff = gff.iloc[1:]
		
			d = pathlib.Path(f"{input.assembly}")
			df = pandas.read_csv(d, sep = '\t')
		
			assembly = df.merge(gff, left_on = ['Name'], right_on= gff.index)
			assembly = assembly.rename(columns={'Name':'Isolate'})
			assembly.to_csv(f"{output}", sep = '\t', index = False)


if PIPELINE == 'all':
 
	rule run_roary:
		input:
			expand("{sample}/{sample}.gff", sample = SAMPLE)
		output:
			"roary/gene_presence_absence.csv", "roary/summary_statistics.txt"
		threads:
			36
		singularity:
			"{% endraw %}{{singularity_dir}}/roary"{% raw %}
{% endraw %}

		shell:
			"""
			roary -p {threads} -f roary {input}
			mv roary_*/* roary
			rm -r roary_*
			"""

	rule pan_figure:
		input:
			"roary/gene_presence_absence.csv"
		output:
			"pan_genome.svg"
		params:
			script_path = SCRIPT_PATH
		shell:
			"""
			perl {params.script_path}/roary2svg.pl {input} > {output}
			"""

rule collate_report:
	input:
		REPORT_INPUT
	output:
		REPORT_OUTPUT
	params:
		pipeline = PIPELINE
	run:	
		import pandas, pathlib, subprocess, numpy
		if {params.pipeline} != 'a':
			df = pandas.read_csv(pathlib.Path(f"core.txt"), sep = '\t')
			df['% USED'] = 100 * (df['LENGTH'] - df['UNALIGNED'])/ df['LENGTH']
			df['% USED'] = df['% USED'].round(2)
			df = df.rename(columns={'ID':'Isolate'})
			df.to_csv(f"report/core_genome.tab", sep='\t', index = False)


		if {params.pipeline} != 's':
			# calculate mean + 2SD and use as cutoff for quality of contigs and fix column names
			dfass = pandas.read_csv(pathlib.Path(f"assembly.tab"), sep = '\t')
			cut = dfass['# Contigs'].mean() + (2* dfass['# Contigs'].std())
			dfass['Quality'] = numpy.where(dfass['# Contigs'] <= cut, 'PASS','FAIL')
			dfass = dfass.rename(columns={'rRNA':'# rRNA', 'CDS':'# CDS'})
			dfass.to_csv(f"report/assembly.tab", sep= '\t', index=False)

			abritamr_summary = pandas.read_csv("summary_matches.csv")
			abritamr_summary.to_csv(f"report/summary_matches.tab", sep= '\t', index=False)
		# cmd = []
		cmd = [f"""
cp seqdata.tab report/seqdata.tab
cp species_identification.tab report/sequence_identification.tab
"""]
		s = ["""
cp core.treefile report/core.treefile
cp distances.tab report/distances.tab
cp core.tab report/core.tab
"""]
		a = ["""
cp mlst.tab report/mlst.tab
"""]
		r = ["""
cp pan_genome.svg report/pan_genome.svg
cp roary/summary_statistics.txt report/summary_statistics.txt
"""]
		if {params.pipeline} == 's':
			cmd.extend(s)
		elif {params.pipeline} == 'a':
			cmd.extend(a)
		elif {params.pipeline} == 'sa':
			cmd.extend(s)
			cmd.extend(a)
		elif {params.pipeline} == 'all':
			cmd.extend(s)
			cmd.extend(a)
			cmd.extend(r)
		subprocess.run(cmd, shell = True)

rule write_html_report:
	input:
		REPORT_OUTPUT
	output:
		'report/report.html'
	params:
		work_dir = WORKDIR,
		jobid = JOB_ID,
		assembler = ASSEMBLER,
		template_path = TEMPLATE_PATH,
		script_path = SCRIPT_PATH,
		pipeline = PIPELINE
	shell:
		"""
		python3 {params.script_path}/write_report.py {params.work_dir} {params.template_path} {params.pipeline} {params.job_id} {params.assembler}
		"""