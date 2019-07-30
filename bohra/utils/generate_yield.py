import pandas, re, pathlib, sys




def get_data(output):
#    print(seqtkdata[0])
    summary_pat = re.compile(r'(\w+):\s(\d+\.?\d?\d?);')
    distinct_values = re.compile(r';\s(\d+)\s\w+')
    record = re.compile(r'(\d+\.?\d?\d?)')
    summary_dict = dict((k, float(v)) for k, v in summary_pat.findall(output[0]))
    summary_dict['distinct_error_codes'] = int(distinct_values.findall(output[0])[0])
    header = output[1].replace('%', '').replace('#', '').replace("POS\t", "").split("\t")
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

def get_coverage(mashdata):

    return(float(mashdata[1].strip().split(':')[1].strip()))


def main(pathtoseqtkdata, pathtomashdata, outputpath):
    output = pathlib.Path(f"{pathtoseqtkdata}").open().readlines()
    seqtkdata = get_data(output)
    df = pandas.DataFrame(data = seqtkdata, index = [0])
    df = df[['Reads','bases','GC content','min_len', 'avg_len', 'max_len','avgQ']]
    mash = pathlib.Path(f"{pathtomashdata}").open().readlines()
    df['Estimated depth'] = int(get_coverage(mash))
    df = df.rename(columns={'bases':'Yield', 'min_len': 'Min len', 'avg_len': 'Avg len', 'max_len':'Max len', 'avgQ': 'Avg Qual'})
    print(df.head())
    df = df[['Reads','Yield','GC content','Min len','Avg len','Max len','Avg Qual','Estimated depth']]
    df.to_csv(pathlib.Path(f"{outputpath}"), index = False, sep = '\t')

  

if __name__ == '__main__':
    
    # print(sys.argv[1], sys.argv[2], sys.argv[3])
    main(f"{sys.argv[1]}", f"{sys.argv[2]}", f"{sys.argv[3]}")