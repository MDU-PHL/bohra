
class MakeWorkflow():
	'''
	A Class to write out the snakemkae file for each specific job

	'''

	# TODO add in logic for species and exp for validation and prefill path
	def write_config_params(self, rerun = False):
		
		return(f"""
job = config['name']
min_aln = int(config['min_perc'])
now = config['now']
REFERENCE = config['reference']


""")

	
	def write_all(self, pipeline = 'sa'):
		base = f"""	expand(\"{{sample}}/seqdata.tab\", sample = SAMPLE),
		\"report/seqdata.tab\""""

		sstring = f"""
		expand(\"{{sample}}/snps.vcf\", sample = SAMPLE),
 		expand(\"{{sample}}/snps.aligned.fa\", sample = SAMPLE),
		\"core.vcf\", 
		\"distances.tab\",
		\"core.treefile\", 
		\"report/core_genome.tab\", 
		\"report/core.treefile\", 
		\"report/distances.tab\",
		\"report/core.tab\""""

		astring = f"""
		expand(\"{{sample}}/{{sample}}.fa\", sample = SAMPLE),
		expand(\"{{sample}}/kraken.tab\", sample = SAMPLE),
		expand(\"{{sample}}/resistome.tab\", sample = SAMPLE),
		expand(\"prokka/{{sample}}/{{sample}}.gff\", sample = SAMPLE),
		expand(\"prokka/{{sample}}/{{sample}}.txt\", sample = SAMPLE),
		\"mlst.tab\", 
		\"species_identification.tab\",
		\"denovo.tab\", 
		\"assembly.tab\", 
		\"resistome.tab\",
		\"report/assembly.tab\",
		\"report/mlst.tab\", 
		\"report/species_identification.tab\", 
		\"report/resistome.tab\""""

		rstring= f"""
		\"roary/gene_presence_absence.csv\", 
		\"pan_genome.svg\", 
		\"report/pan_genome.svg\",
		\"report/summary_statistics.txt\""""

		report_string = f"""
		\"report/report.html\""""
		i_dict = {
			's': f"""{base},{sstring},{report_string}
""",
			'a': f"""{base},{astring},{report_string}
""",
			'sa': f"""{base},{astring},{sstring},{report_string}
""",
			'all': f"""{base},{astring},{sstring},{rstring},{report_string}

"""
		}

		return(f"""
rule all:
	input:
	{i_dict[pipeline]}
""")


	def write_seqdata(self, prefillpath = False):
		return(f"""
rule seqdata:
	input:
		'READS/{{sample}}/R1.fq.gz',
		'READS/{{sample}}/R2.fq.gz'
	output:
		"{{sample}}/seqdata.tab"
	shell:
		\"""
		seqtk fqchk {{input[0]}} {{input[1]}} > {{output}}
		\"""
""")

	def write_kraken(self, prefillpath = ''):
		return(f"""
rule kraken:
	input:
		'READS/{{sample}}/R1.fq.gz',
		'READS/{{sample}}/R2.fq.gz'
	output:
		"{{sample}}/kraken.tab"

	shell:
		\"""
		KRAKENPATH={prefillpath}/{{wildcards.sample}}/kraken2.tab
		if [ -f $KRAKENPATH ]; then
			cp $KRAKENPATH {{output}}
		else
			kraken2 --paired {{input[0]}} {{input[1]}} --memory-mapping --minimum-base-quality 13 --report {{output}} 
		fi
		\"""
		
""")

	def write_combine_kraken(self):
		return(f"""
rule combine_kraken:
	input: 
		expand("{{sample}}/kraken.tab", sample = SAMPLE)
	output:
		"species_identification.tab"
	run:
		import pandas, pathlib, subprocess
		kfiles = f"{{input}}".split()
		id_table = pandas.DataFrame()
		for k in kfiles:
			kraken = pathlib.Path(k)
			
			df = pandas.read_csv(kraken, sep = "\t", header =None, names =  ['percentage', 'frag1', 'frag2','code','taxon','name'])
			df['percentage'] = df['percentage'].apply(lambda x:float(x.strip('%')) if isinstance(x, str) == True else float(x)) #remove % from columns
			df = df.sort_values(by = ['percentage'], ascending = False)
			df = df[df['code'].isin(['U','S'])]     
			df = df.reset_index(drop = True) 
			tempdf = pandas.DataFrame()
			d = {{'Isolate': f"{{kraken.parts[0]}}",    
					'#1 Match': df.ix[0,'name'].strip(), '%1': df.ix[0,'percentage'],
					'#2 Match': df.ix[1,'name'].strip(), '%2': df.ix[1,'percentage'],       
					'#3 Match': df.ix[2,'name'].strip(), '%3': df.ix[2,'percentage'] ,
					'#4 Match': df.ix[3,'name'].strip(), '%4': df.ix[3,'percentage']
					}}
			tempdf = pandas.DataFrame(data = d, index= [0])
			if id_table.empty:
					id_table = tempdf
			else:
					id_table = id_table.append(tempdf)
		id_table.to_csv(f"{{output}}", sep = "\t", index = False)
		subprocess.run(f"sed -i 's/%[0-9]/%/g' {{output}}", shell=True)
""")
	def write_estimate_coverage(self):
		return(f"""
rule estimate_coverage:
	input:
		\"READS/{{sample}}/R1.fq.gz\",
		\"READS/{{sample}}/R2.fq.gz\"
	output:
		\"{{sample}}/mash.txt\"
	shell:
		\"""
		mash sketch -r {{input[0]}} {{input[1]}} -m 3 -k 31 -o mash  &> {{output}}
		\"""
""")
	def write_generate_yield(self, script_path):
		return(f"""
rule generate_yield:
	input:
		\"{{sample}}/mash.txt\",
		\"{{sample}}/seqdata.tab\"
	output:
		\"{{sample}}/yield.tab"
	shell:
		\"""
		python3 {script_path}/generate_yield.py {{input[1]}} {{input[0]}} {{output}}
		\"""

""")

	def write_combine_seqdata(self):
		
		return(f"""
rule combine_seqdata:
	input:
		expand("{{sample}}/yield.tab", sample = SAMPLE)
	output:
		"seqdata.tab"
	run:
		import pathlib, pandas, numpy
		sdfiles = f"{{input}}".split()
		seqdata = pandas.DataFrame()
		for sd in sdfiles:
			p = pathlib.Path(sd)
			df = pandas.read_csv(sd, sep = "\t")
			print(df)
			df['Isolate'] = f\"{{p.parts[0]}}\"
			
			if seqdata.empty:
				seqdata = df
			else:
				seqdata = seqdata.append(df)
		seqdata['Quality'] = numpy.where(seqdata['Estimated depth'] >= 40, 'PASS','FAIL')
		seqdata = seqdata[['Isolate','Reads','Yield','GC content','Min len','Avg len','Max len','Avg Qual','Estimated depth', 'Quality']]
		seqdata.to_csv(f"{{output}}", sep = '\t', index = False)
		

""")

	def write_snippy(self):
		
		return(f"""
rule snippy:
	input:
		'READS/{{sample}}/R1.fq.gz',
		'READS/{{sample}}/R2.fq.gz'
	output:
		'{{sample}}/snps.vcf',
		'{{sample}}/snps.aligned.fa'
	threads:
		8
	shell:
		\"""
		snippy --outdir {{wildcards.sample}} --ref {{REFERENCE}} --R1 {{input[0]}} --R2 {{input[1]}} --force --cpus {{threads}}
		\"""
	""")


	def write_qc_snippy_initial(self):
		return(f"""
rule qc_snippy: 
	input:
		expand('{{sample}}/snps.aligned.fa', sample = SAMPLE)
		
	output:
		'core_isolates.txt'
		
	run:
		from Bio import SeqIO
		import pathlib
		import pandas
		# create an output
		isolate_list = []
		excluded_list = []
		outfile = pathlib.Path(f"{{output[0]}}")
		# get input file list
		input_list = f"{{input}}".split()
		# set the log path
		logpath = pathlib.Path('isolates.log')
		for i in input_list: # for each input file
			# get the isolate name
			p = pathlib.Path(f"{{i}}")
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
						isolate_list.append(f"{{isolate}}")
					else:
						excluded_list.append(isolate)
						print(f"{{isolate}} has been excluded from the analysis due to poor alignement with reference")
						
		isolate_list = list(set(isolate_list))
		with open(outfile, 'w') as f:
			f.write('\\n'.join(isolate_list))
		# get log if the excluded list has any isolates in
		if excluded_list != []:
			if logpath.exists():
				lf = pandas.read_csv(logpath, sep = '\t', index_col = False)
				for e in excluded_list:
					lf.loc[lf['Isolate'] == e.strip('#'), 'Status'] = f"(FAILED ALIGNMENT (<{{min_aln}}% ALIGNMENT))"
					lf.loc[lf['Isolate'] == e.strip('#'), 'Date'] = f"{{config['day']}}"
					lf.to_csv(logpath, sep = '\t', index=False)

	""")


	def write_snippy_core(self, mask):
			return(f"""
rule run_snippy_core:
	input:
		'core_isolates.txt'
	output:
		'core.vcf',
		'core.txt',
		'core.aln', 
		'core.full.aln',
		'core.tab'
	
	shell:
		\"""
		snippy-core --ref {{REFERENCE}} {mask} $(cat core_isolates.txt)
		
		\"""
	""")

	def write_snp_dists(self):
		return(f"""
rule run_snpdists:
	input:
		'core.aln' 
	output:
		'distances.tab' 
	
	shell:
		\"""
		snp-dists {{input}} > {{output}}
		\"""
	""")


	def write_tree(self, script_path, alntype):
		return(f"""
rule calculate_iqtree_command_{alntype}:
	input:
		'{alntype}.aln',
		REFERENCE
	output:
		'run_iqtree_{alntype}.sh'
	run:
		from Bio import SeqIO
		import pathlib
		ref = pathlib.Path(REFERENCE)
		name = ref.stem
		if '.fa' not in REFERENCE:
			print(f"converting {{REFERENCE}}")
			SeqIO.convert(f"{{input[1]}}", 'genbank', f"{{name}}.fa", 'fasta')
			ref = f"{{name}}.fa"
		else:
			ref = f"{{ref}}"
		
		shell("bash {script_path}/iqtree_generator.sh {{ref}} {{input[0]}} {alntype} 20 > {{output}}")

rule run_iqtree_{alntype}:
	input:
		'run_iqtree_{alntype}.sh'
	
	output:
		'{alntype}.iqtree',
		'{alntype}.treefile',
		
	
	shell:
		\"""	
		bash run_iqtree_{alntype}.sh
		
		rm *.ckp.gz *.contree *.bionj
		\"""
		
	""")



	def write_assemblies(self, prefillpath = '', assembler = 'skesa'):
		
		if assembler == 'skesa':
			assemble = f"skesa --fastq {{input[0]}},{{input[1]}} --vector_percent 1 --use_paired_ends --cores {{threads}} > {{output}}"
		elif assembler == 'shovill':
			assemble = f""""
			shovill --outdir {{wildcards.sample}} --R1 {{input[0]}} --R2 {{input[1]}} --force --minlen 500 --cpus {{threads}}
			mv {{wildcards.sample}}/contigs.fa {{output}}
"""
		elif assembler == 'spades':
			assemble = f"""
			spades.py -o {{wildcards.sample}} -1 {{input[0]}} -2 {{input[1]}} --threads {{threads}}
			mv {{wildcards.sample}}/contigs.fasta {{output}}
"""
		
		return(f"""
rule assemble:
	input:
		'READS/{{sample}}/R1.fq.gz',
		'READS/{{sample}}/R2.fq.gz'
	output:
		'{{sample}}/{{sample}}.fa'
	threads:
		16
	shell:
		\"""
		ASSEMBLEPATH={prefillpath}/{{wildcards.sample}}
		if [ -f $ASSEMBLEPATH/contigs.fa ]; then
			cp $ASSEMBLEPATH/contigs.fa {{output}}

		else
			echo No assembly found. Assembling {{wildcards.sample}} with {assembler}
			{assemble}
		fi		
		\"""
	""") # TODO Will need to add in assembly options here

	def write_mlst(self):
		return(f"""
rule mlst:
	input:
		expand('{{sample}}/{{sample}}.fa', sample = SAMPLE)
	output:
		'mlst.tab'
	
	shell:
		\"""
		mlst --nopath {{input}} | sed 's/\.fa//g' | sed '1iIsolate\tScheme\tST\tAlleles' > {{output}}
		\"""

	""")

	def write_resistome(self):
		return(f"""
rule resistome:
	input:
		'{{sample}}/{{sample}}.fa'
	output:
		'{{sample}}/resistome.tab'
	
	shell:
		\"""
		abricate --nopath {{input}} > {{output}}
		\"""

	""")


	def write_combine(self):
		return(f"""
rule combine_results:
	input:
		expand('{{sample}}/resistome.tab', sample = SAMPLE)
		
	output:
		'resistome.tab'
		
	
	shell:
		\"""
		abricate --summary {{input}} | sed 's/\/resistome.tab//g' | sed 's/\#FILE/Isolate/g' > {{output[0]}}
		
		\"""
	""")

	def write_prokka(self):
		return(f"""
rule run_prokka:
	input:
		"{{sample}}/{{sample}}.fa"
	output:
		"prokka/{{sample}}/{{sample}}.gff","prokka/{{sample}}/{{sample}}.txt"
	shell:
		\"""
		prokka --outdir prokka/{{wildcards.sample}} --prefix {{wildcards.sample}} --mincontiglen 500 --notrna --fast --force {{input}}
		\"""
""")

	def write_roary(self):
		return(f"""
rule run_roary:
	input:
		expand("prokka/{{sample}}/{{sample}}.gff", sample = SAMPLE)
	output:
		"roary/gene_presence_absence.csv", "roary/summary_statistics.txt"
	threads:
		36
	shell:
		\"""
		roary -p {{threads}} -f roary {{input}}
		mv roary_*/* roary
		rm -r roary_*
		\"""
""")

	def write_pan_graph(self, script_path):
		return(f"""
rule pan_figure:
	input:
		"roary/gene_presence_absence.csv"
	output:
		"pan_genome.svg"
	shell:
		\"""
		perl {script_path} roary2svg.pl {{input}} > {{output}}
		\"""
""")

	def write_assembly_stats(self, script_path):
		return(f"""
rule assembly_statistics:
	input:
		expand("{{sample}}/{{sample}}.fa", sample = SAMPLE)
	output:
		"denovo.tab"
	shell:
		\"""
		 python3 {script_path}/assembly_stat.py {{input}} -m 500 > {{output}}
		\"""
	
""")

	def write_gff_summary(self):
		return(f"""
rule combine_assembly_metrics:
	input:
		prokka = expand("prokka/{{sample}}/{{sample}}.txt",sample = SAMPLE), 
		assembly = "denovo.tab"
	output:
		"assembly.tab"
	run:
		import pandas, pathlib

		prokka = f"{{input.prokka}}".split()
		gff = pandas.DataFrame()
		
		for p in prokka:
			g = pathlib.Path(p)
			df = pandas.read_csv(g, sep = ':', header = None, names = ['cond', f"{{g.parts[1]}}"])
			
			if gff.empty:
					gff = df
			else:
					gff = gff.merge(df, how = 'outer')
		gff = gff[gff['cond'].isin(['CDS', 'rRNA'])]
		gff = gff.T
		gff.columns = gff.iloc[0]
		gff = gff.iloc[1:]
	
		d = pathlib.Path(f"{{input.assembly}}")
		df = pandas.read_csv(d, sep = '\\t')
	
		assembly = df.merge(gff, left_on = ['Name'], right_on= gff.index)
		assembly = assembly.rename(columns={{'Name':'Isolate'}})
		assembly.to_csv(f"{{output}}", sep = '\\t', index = False)


""")
	def io_collation(self):
		i = f"'seqdata.tab'"
		o = f"'report/seqdata.tab'"
		# snps input string
		si = f"'core.txt', 'core.treefile', 'core.tab', 'distances.tab', 'core.tab'"
		# snps output string
		so = f"'report/core_genome.tab', 'report/core.treefile','report/distances.tab','report/core.tab'"
		# assembly input string
		ai = f"'assembly.tab', 'mlst.tab', 'species_identification.tab', 'resistome.tab'"
		# assembly output string
		ao = f"'report/assembly.tab', 'report/mlst.tab', 'report/species_identification.tab', 'report/resistome.tab'"
		# roary input
		ri = f"'pan_genome.svg', 'roary/summary_statistics.txt'"
		# roary output 
		ro = f"'report/pan_genome.svg', 'report/summary_statistics.txt'"
		# input string based on pipeline value
		return(i,o,si,so,ai,ao,ri,ro)
	def run_setup(self):
		
		run_string = f"""
		import pandas, pathlib, subprocess, numpy
"""
		coretxt = f"""
		
		# for core.txt
		df = pandas.read_csv(pathlib.Path(f"core.txt"), sep = '\\t')
		df['% USED'] = 100 * (df['LENGTH'] - df['UNALIGNED'])/ df['LENGTH']
		df['% USED'] = df['% USED'].round(2)
		df = df.rename(columns={{'ID':'Isolate'}})
		df.to_csv(f"report/core_genome.tab", sep='\\t', index = False)

"""
		assemblystatstring = f"""

		# calculate mean + 2SD and use as cutoff for quality of contigs and fix column names
		dfass = pandas.read_csv(pathlib.Path(f"assembly.tab"), sep = '\\t')
		cut = dfass['# Contigs'].mean() + (2* dfass['# Contigs'].std())
		dfass['Quality'] = numpy.where(dfass['# Contigs'] <= cut, 'PASS','FAIL')
		dfass = dfass.rename(columns={{'rRNA':'# rRNA', 'CDS':'# CDS'}})
		dfass.to_csv(f"report/assembly.tab", sep= '\\t', index=False)
"""
		cmdstring = f"""cp seqdata.tab report/seqdata.tab
"""
		scmd = f"""cp core.treefile report/core.treefile
cp distances.tab report/distances.tab
cp core.tab report/core.tab
"""
		acmd = f"""cp mlst.tab report/mlst.tab
cp species_identification.tab report/species_identification.tab
cp resistome.tab report/resistome.tab
"""
		rcmd = f"""cp pan_genome.svg report/pan_genome.svg
cp roary/summary_statistics.txt report/summary_statistics.txt
"""
		cmd = f"""
"""
		r_dict = {
			's' : f"""
		{run_string}
		{coretxt}
		cmd = f\"""
{cmdstring}{scmd}
		\"""
		subprocess.run(cmd, shell = True)
""",
			'a' : f"""
		{run_string}
		{assemblystatstring}
		cmd = f\"""
{cmdstring}{acmd}
		\"""
		subprocess.run(cmd, shell = True)
""",
			'sa' : f"""
		{run_string}
		{coretxt}
		{assemblystatstring}
		cmd = f\"""
{cmdstring}{scmd}{acmd}
		\"""
		subprocess.run(cmd, shell = True)
""",
			'all' : f"""
		{run_string}
		{coretxt}
		{assemblystatstring}
		cmd = f\"""
{cmdstring}{scmd}{acmd}{rcmd}
		\"""
		subprocess.run(cmd, shell = True)
"""}
		return(r_dict)


	def write_report_collation(self, pipeline = 'sa'):
		# input string always starts with seqdata and output string always starts with seqdata
		i,o,si,so,ai,ao,ri,ro = self.io_collation()
		
		io_dict = {
			's' : {
				'input':f"{i}, {si}",
				'output':f"{o}, {so}"
			},
			'a' : {
				'input':f"{i}, {ai}",
				'output':f"{o}, {ao}"
			},
			'sa' : {
				'input':f"{i}, {ai}, {si}",
				'output':f"{o}, {ao}, {so}"
			},
			'all' : {
				'input':f"{i}, {ai}, {si}, {ri}",
				'output':f"{o}, {ao}, {so}, {ro}"
			}
		}
		r_dict = self.run_setup()
		input_string = io_dict[pipeline]['input']
		output_string = io_dict[pipeline]['output']
		run_string = r_dict[pipeline]
		return(f"""
rule collate_report:
	input:
		{input_string}
	output:
		{output_string}
	run:
		{run_string}
""")

	def write_html(self, workdir, resources, job_id, script_path,assembler, pipeline = 'sa'):
		seq = f"'report/seqdata.tab'"
		s = f"'report/core_genome.tab', 'report/core.treefile', 'report/distances.tab'"
		a = f"'report/assembly.tab','report/mlst.tab', 'report/species_identification.tab','report/resistome.tab'"
		r = f"'report/pan_genome.svg', 'report/summary_statistics.txt'"
		if pipeline == 's':
			input_string = f"{seq},{s}"
		elif pipeline == 'a':
			input_string = f"{seq}, {a}"
		elif pipeline == 'sa':
			input_string = f"{seq}, {a}, {s}"
		elif pipeline == 'all':
			input_string = f"{seq}, {a}, {s}, {r}"
		wd_string = f"{workdir / job_id}"
		return(f"""

rule write_html_report:
	input:
		{input_string}
	output:
		'report/report.html'
	
	shell:
		\"""
		python3 {script_path}/write_report.py {wd_string} {resources} {pipeline} {job_id} {assembler}
		\"""
""")
	
	def write_gubbins(self, gubbins,job):
		if gubbins == True:
			return("""
rule run_gubbins:
		input:
			'core.full.aln',
			'core.treefile'
		
		output:
			'clean.full.aln', 
			'gubbins.aln'
			
		# TODO fix gubbins!!
		shell:
			\"""
			snippy-clean_full_aln {{input[0]}} > {{output[0]}}
			run_gubbins.py -c 36 --starting_tree {{input[1]}} {{output[0]}}
			snp-sites -c clean.full.filtered_polymorphic_sites.fasta > {{output[1]}}
			\"""
""")
		else:
			return('')
	
