import toml, pathlib, subprocess, sys, re
from Bio import SeqIO
from snakemake import shell

def generate_seqdata_cmd(r1, r2, isolate):
    
    cmd = f"seqtk fqchk {r1} {r2} > {isolate}/seqdata.tab"
    
    return cmd

def run_cmd(cmd):
    
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding = 'utf-8')
    
    return p.returncode


def get_data(cmd_output):
#    print(seqtkdata[0]
    output = pathlib.Path(f"{cmd_output}").open().readlines()
    summary_pat = re.compile(r'(\w+):\s(\d+\.?\d?\d?);')
    distinct_values = re.compile(r';\s(\d+)\s\w+')
    # print(distinct_values)
    record = re.compile(r'(\d+\.?\d?\d?)')
    summary_dict = dict((k, float(v)) for k, v in summary_pat.findall(output[0]))
    # print(distinct_values.findall(output[0]))
    summary_dict['distinct_error_codes'] = int(distinct_values.findall(output[0])[0])
    header = output[1].strip().replace('%', '').replace('#', '').replace("POS\t", "").split("\t")
    line_one = [float(v) for v in record.findall(output[2])]
    summary_dict.update(dict(zip(header, line_one)))
    header.insert(0, 'position')
    summary_dict['Reads'] = dict(zip(header, [float(v) for v in record.findall(output[3])]))['bases']
    summary_dict['GC content'] = round((summary_dict['C'] + summary_dict['G']),2)
    median_position = (summary_dict['Reads'] + 1) / 2
    for line in output[4:]:
        parsed_record = dict(zip(header, [float(v) for v in record.findall(line)]))
        if parsed_record['bases'] > median_position:
            continue
        else:
            summary_dict['med_len'] = parsed_record['position'] - 1
            break
    return(summary_dict)

def get_length(fasta):

    length = 0
    for i in SeqIO.parse(fasta,'fasta'): # use BioPython to determine percent alignment
        l = len(i.seq)
        print(l)
        length = length + l
        print(length)
    print(length)
    return length

def get_coverage(reference, isolate, yld):

    if '.fa' in reference:
        length = get_length(reference)
    elif pathlib.Path('ref.fa').exists():
        length = get_length('ref.fa')
    elif '.gbk' in reference:
        SeqIO.convert(reference, 'genbank', 'ref.fa', 'fasta')
        length = get_length('ref.fa')
        
    
    cov = float(yld)/length

    return cov


def open_toml(tml):

    data = toml.load(tml)

    return data

def write_toml(data, output):
    
    with open(output, 'wt') as f:
        toml.dump(data, f)
    
def main(r1, r2, isolate, mincov, reference):
    
    # msh = open_toml(mash)
    # print(msh)
    data = {}
    data[isolate] = {}
    data[isolate]['seqdata'] = {}
    data[isolate]['seqdata']['R1'] = r1
    data[isolate]['seqdata']['R2'] = r2
    data[isolate]['seqdata']['file'] = f"{pathlib.Path(isolate,'seqdata.tab')}"


    cmd = generate_seqdata_cmd(r1 = r1, r2 = r2, isolate = isolate)
    p = run_cmd(cmd)
    
    if p == 0:
        # print(cov)
        data[isolate]['seqdata']['data'] = get_data(cmd_output = data[isolate]['seqdata']['file'])
        cov = get_coverage(reference = reference, isolate = isolate, yld = data[isolate]['seqdata']['data']['bases'])
        data[isolate]['seqdata']['data']['Estimated_coverage'] = round(float(cov),2)
        data[isolate]['seqdata']['data']['Quality'] = 'PASS' if data[isolate]['seqdata']['data']['Estimated_coverage'] > int(mincov) else 'FAIL - will be removed from further analysis'
        write_toml(data = data, output = f"{isolate}/seqdata.toml")   


r1 = snakemake.input.r1
r2 = snakemake.input.r2
# mash = snakemake.input.mash
isolate = snakemake.wildcards.sample
mincov = snakemake.params.mincov
output = snakemake.output
reference = snakemake.params.reference

main(r1 = r1, r2 = r2, isolate = isolate, mincov = mincov, reference = reference)



# mash triangle -C *.msh

# mash sketch -m 5 -s 10000 -r -o 2019-12803-6/sketch -I 2019-12803-6 -C 2019-12803-6/R1.fq.gz 2019-12803-6/R1.fq.gz