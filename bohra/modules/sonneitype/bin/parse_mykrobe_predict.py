#!/usr/bin/env python3
import json
import sys
import pandas as pd
from argparse import ArgumentParser
"""
This script was cloned from https://github.com/katholt/sonneityping 10th May 2021
"""
def get_arguments():
    parser = ArgumentParser(description='Parse mykrobe predict JSON files')

    # job submission options
    parser.add_argument('--jsons', required=True, nargs='+', help='JSON files output from mykrobe predict')
    parser.add_argument('--alleles', required=True, help='alleles.txt file that defines each genotype, its marker, and the human-readable name')
    parser.add_argument('--prefix', required=True, help='prefix for output files')

    return parser.parse_args()

def extract_qrdr_info(genome_data, genome_name, spp_call):

    # set up empty dict for final output
    qrdr_out_dict = {}
    # make a list of all possible mutations
    qrdr_possible = ['parC_S80I', 'gyrA_S83L', 'gyrA_S83A', 'gyrA_D87G', 'gyrA_D87N', 'gyrA_D87Y']

    # can only extract qrdr info if sample is sonnei, otherwise set 'unknown' for all calls
    if spp_call == "Shigella_sonnei":
        try:
            qrdr_data = genome_data["susceptibility"]["quinolones"]
        # if nothiing has been called, just set everything to NA
        except KeyError:
            qrdr_out_dict['genome'] = [genome_name]
            for allele in qrdr_possible:
                qrdr_out_dict[allele] = ['NA']
            qrdr_out_dict['num QRDR'] = ['NA']
            # make table
            qrdr_out_df = pd.DataFrame(qrdr_out_dict, columns=['genome'] + qrdr_possible + ['num QRDR'])
            return qrdr_out_df
    else:
        qrdr_out_dict['genome'] = [genome_name]
        for allele in qrdr_possible:
            qrdr_out_dict[allele] = ['NA']
        qrdr_out_dict['num QRDR'] = ['NA']
        # make table
        qrdr_out_df = pd.DataFrame(qrdr_out_dict, columns=['genome'] + qrdr_possible + ['num QRDR'])
        return qrdr_out_df

    # set up the list of calls in this genome
    qrdr_calls = []
    # only need to parse if the predict is R - if it's S all values will be 0
    if qrdr_data["predict"] == "R":
        calls = qrdr_data["called_by"]
        # loop through each mutation
        for mutation in calls:
            # get the mutation name and add it to our list of calls
            qrdr_calls.append(mutation.split('-')[0])
    # now create our pandas dataframe for this genome
    qrdr_out_dict['genome'] = [genome_name]
    num_qrdr_calls = 0
    for allele in qrdr_possible:
        if allele in qrdr_calls:
            qrdr_out_dict[allele] = [1]
            num_qrdr_calls += 1
        else:
            qrdr_out_dict[allele] = [0]
    # add column with total number of qrdr calls
    qrdr_out_dict['num QRDR'] = [num_qrdr_calls]
    # make table
    qrdr_out_df = pd.DataFrame(qrdr_out_dict, columns=['genome'] + qrdr_possible + ['num QRDR'])
    return qrdr_out_df

def inspect_calls(full_lineage_data):
    # for all the genotypes that mykrobe is calling, inspect the values listed next to the actual call heirarhcy
    genotype_details = full_lineage_data['calls_summary']
    genotype_list = list(genotype_details.keys())
    # intialise empty values for best score and corresponding genotype
    best_score = 0
    best_genotype = None
    # node support is the number of 'good nodes' called by mykrobe
    # needs to be X/Y, where Y is the total number of levels at that genotype
    node_support = None
    # loop through each genotype
    for genotype in genotype_list:
        # the maximum score we can get is the total depth of the heirarchy
        max_score = genotype_details[genotype]['tree_depth']
        # the actual score is a sum of the values within each levels call
        actual_score = sum(list(genotype_details[genotype]['genotypes'].values()))
        # if actual score is equal to the max score, then this is the best genotype
        # BUT WE NEED TO DEAL WITH INSTANCES WHERE WE MIGHT HAVE A GREAT CALL 3.7.25 say AND A LESS GREAT CALL 3.6.1.1.2 say BUT BECACUSE THE HEIRARHCY IS BIGGER FOR 3.6.1.1.2 WE INADVERTANTLY CALL THAT
        if actual_score == max_score:
            best_score = actual_score
            best_genotype = genotype
            #node_support = genotype_details[best_genotype]['good_nodes']
        # if actual score is < max score, but is greater than the current best score,
        # then this is our best genotype for the moment
        elif actual_score < max_score and actual_score > best_score:
            best_score = actual_score
            best_genotype = genotype
            #node_support = genotype_details[best_genotype]['good_nodes']

    # determine any quality issues and split them amongst our various columns
    # if the node support for the best genotype is not '1' at all positions, we need to report that
    best_calls = genotype_details[best_genotype]['genotypes']
    # get a list of all the calls for calculating confidence later
    best_calls_vals = list(best_calls.values())
    poorly_supported_markers = []
    # dict that keeps percentage supports, key=percent support; value=level
    lowest_within_genotype_percents = {}
    # list that gives supports for all markers in best call
    final_markers = []
    for level in best_calls.keys():
        # regardless of the call, get info
        call_details = full_lineage_data['calls'][best_genotype][level]
        # check that there is something there
        if call_details:
            # need to do this weird thing to grab the info without knowing the key name
            call_details = call_details[list(call_details.keys())[0]]
            ref = call_details['info']['coverage']['reference']['median_depth']
            alt = call_details['info']['coverage']['alternate']['median_depth']
            # calculate percent support for marker
            try:
                percent_support = alt / (alt + ref)
            except ZeroDivisionError:
                percent_support = 0
            marker_string = level + ' (' + str(best_calls[level]) + '; ' + str(alt) + '/' + str(ref) + ')'
            final_markers.append(marker_string)
        # if the value is null, just report 0 (indicates that no SNV detected, either ref or alt?)
        else:
            # note we do not have a markers for lineage5.1 so don't report this as 0
            if level != 'lineage5.1':
                lowest_within_genotype_percents[0] = level
                marker_string = level + ' (0)'
                # its returned null so is therefore by definition poorly supported
                poorly_supported_markers.append(marker_string)
                # add to the final markers list too
                final_markers.append(marker_string)
        # if call is 1 then that is fine, count to determine confidence later
        # if call is 0.5, then get info
        # if call is 0, there will be no info in the calls section, so just report 0s everywhere
        if best_calls[level] < 1 and level != 'lineage5.1':
            # then it must be a 0 or a 0.5
            # report the value (0/0.5), and also the depth compared to the reference
            lowest_within_genotype_percents[percent_support] = level
            # add this to the list of poorly supported markers
            poorly_supported_markers.append(marker_string)

    # determining final confidence is based ONLY on the actual genotype, not incongruent genotype calls
    # if there is a lineage5.1 in the call, then we need to remove one of the 0s in best_calls_vals
    if 'lineage5.1' in best_calls.keys():
        best_calls_vals.remove(0)
    # strong = quality 1 for all calls
    if best_calls_vals.count(0) == 0 and best_calls_vals.count(0.5) == 0:
        confidence = 'strong'
        lowest_support_val = ''
    # moderate = quality 1 for all calls bar 1 (which must be quality 0.5 with percent support > 0.5)
    elif best_calls_vals.count(0) == 0 and best_calls_vals.count(0.5) == 1:
        if min(lowest_within_genotype_percents.keys()) > 0.5:
            confidence = 'moderate'
        else:
            confidence = 'weak'
        lowest_support_val = round(min(lowest_within_genotype_percents.keys()), 3)
    # weak = more than one quality < 1, or the single 0.5 call is < 0.5% support
    else:
        confidence = 'weak'
        lowest_support_val = round(min(lowest_within_genotype_percents.keys()), 3)

    # make a list of all possible quality issues (incongruent markers, or not confident calls within the best geno)
    non_matching_markers = []
    non_matching_supports = []
    # we now want to report any additional markers that aren't congruent with our best genotype
    #ie if 3.6.1 is the best genotype, but we also have a 3.7.29 call, we need to report the 3.7 and 3.29 markers as incongruent
    if len(genotype_list) > 1:
        # remove the best genotype from the list
        genotype_list.remove(best_genotype)
        # loop through each incongruent phenotype
        for genotype in genotype_list:
            # extract the calls for that genotype
            other_calls = genotype_details[genotype]['genotypes']
            # for every call in our BEST calls, we're only interested
            # in calls that are incongruent
            for call in other_calls.keys():
                if call not in best_calls.keys():
                    call_info = full_lineage_data['calls'][genotype][call]
                    # check that there is something there (if the value is null, don't report it for incongruent calls)
                    if call_info:
                        # need to do this weird thing to grab the info without knowing the key name
                        call_info = call_info[list(call_info.keys())[0]]
                        ref_depth = call_info['info']['coverage']['reference']['median_depth']
                        alt_depth = call_info['info']['coverage']['alternate']['median_depth']
                        # only keep the call if the alternate has a depth of > 1
                        # this is because mykrobe fills in intermediate levels of the heirarchy with 0s
                        # if a lower level SNV marker is detected
                        if alt_depth >= 1:
                            percent_support = alt_depth / (alt_depth + ref_depth)
                            non_matching_supports.append(percent_support)
                            marker_string = call + ' (' + str(other_calls[call]) + '; ' + str(alt_depth) + '/' + str(ref_depth) + ')'
                            non_matching_markers.append(marker_string)

    # get max value for non matching supports, if empty, return ''
    if len(non_matching_supports) > 0:
        max_non_matching = round(max(non_matching_supports), 3)
    else:
        max_non_matching = ''

    return best_genotype, confidence, lowest_support_val, poorly_supported_markers, max_non_matching, non_matching_markers, final_markers

def extract_lineage_info(lineage_data, genome_name, lineage_name_dict):

    # first, extract phylo spp information - if the spp is not sonnei then don't proceed
    spp_data = lineage_data['species']
    spp_call = list(spp_data.keys())[0]
    # if spp is unknown, then this is not sonnei, exit this function
    if spp_call == "Unknown":
        out_dict = {'genome':[genome_name], 'species':['not S. sonnei'], 'name':['NA'],'final genotype':['NA'],'confidence':['NA'],'lowest support for genotype marker':[''], 'poorly supported markers':[''], 'node support':[''], 'max support for additional markers':[''], 'additional markers':['']}
        out_df = pd.DataFrame(out_dict, columns=['genome', 'final genotype', 'name', 'confidence', 'lowest support for genotype marker', 'poorly supported markers', 'node support', 'max support for additional markers', 'additional markers'])
        return out_df, spp_call
    else:
        # if it is sonnei, then get the percentage
        spp_percentage = spp_data["Shigella_sonnei"]["percent_coverage"]
        # if the percentage is <90, then exit this function as it's likely not sonnei
        if spp_percentage < 90:
            out_dict = {'genome':[genome_name], 'species':['not S. sonnei'], 'name':['NA'],'final genotype':['NA'],'confidence':['NA'],'lowest support for genotype marker':[''], 'poorly supported markers':[''], 'node support':[''], 'max support for additional markers':[''], 'additional markers':['']}
            out_df = pd.DataFrame(out_dict, columns=['genome', 'final genotype', 'name', 'confidence', 'lowest supprot for genotype marker', 'poorly supported markers', 'node support', 'max support for additional markers', 'additional markers'])
            return out_df, "Unknown"

    # if we are sonnei, then continue
    # set up dictionary for final table output
    lineage_out_dict = {'genome': [genome_name]}

    # this try/except statement deals with instances where for some reason there is no lineage output
    # in the json file
    try:
        genotype_calls = lineage_data['lineage']['lineage']
    except KeyError:
        out_dict = {'genome':[genome_name], 'species':['S. sonnei'], 'name':['NA'],'final genotype':['uncalled'],'confidence':['NA'],'lowest support for genotype marker':[''], 'poorly supported markers':[''], 'max support for additional markers':[''], 'additional markers':[''], 'node support':['']}
        out_df = pd.DataFrame(out_dict, columns=['genome', 'final genotype', 'name', 'confidence', 'lowest support for genotype marker', 'poorly supported markers', 'node support', 'max support for additional markers', 'additional markers'])
        return out_df, spp_call
    # if there are no calls, populate with none
    if len(genotype_calls) == 0:
        lineage_out_dict['final_genotype'] = ['uncalled']
        lineage_out_dict['name']  = ['NA']
        lineage_out_dict['confidence'] = ['NA']
        lineage_out_dict['lowest support for genotype marker']  = ['']
        lineage_out_dict['poorly supported markers']  = ['']
        lineage_out_dict['max support for additional markers']  = ['']
        lineage_out_dict['additional markers']  = ['']
        lineage_out_dict['node support'] = ['']
    # if there is just one call, list that - but CHECK that all values in heirarchy are good (>0.5)
    # only report up to level in heirarchy where we have a good call
    if len(genotype_calls) == 1:
        # inspect calls for genotype
        best_genotype, confidence, lowest_support_val, poorly_supported_markers, non_matching_support, non_matching_markers, final_markers = inspect_calls(lineage_data['lineage'])
        # extract genotype from lineage thing, add that to the table
        genotype = best_genotype.split('lineage')[1]
        lineage_out_dict['final genotype'] = [genotype]
        lineage_out_dict['name'] = [lineage_name_dict[best_genotype]]

        # fill in all details on confidence/support
        lineage_out_dict['confidence'] = [confidence]
        lineage_out_dict['lowest support for genotype marker'] = [lowest_support_val]
        lineage_out_dict['poorly supported markers'] = ['; '.join(poorly_supported_markers)]
        lineage_out_dict['max support for additional markers'] = [non_matching_support]
        lineage_out_dict['additional markers'] = ['; '.join(non_matching_markers)]
        lineage_out_dict['node support'] = ['; '.join(final_markers)]
        #lineage_out_dict['all_genotype_calls'] = genotype_calls
    # if there is more than one call, we want to report the best, but also other calls
    elif len(genotype_calls) > 1:
        # get the info for each call
        # work out which lineage has the best call
        best_genotype, confidence, lowest_support_val, poorly_supported_markers, non_matching_support, non_matching_markers, final_markers = inspect_calls(lineage_data['lineage'])
        # now write out info
        lineage_out_dict['final genotype'] = [best_genotype.split('lineage')[1]]
        lineage_out_dict['name'] = [lineage_name_dict[best_genotype]]
        lineage_out_dict['confidence'] = [confidence]
        lineage_out_dict['lowest support for genotype marker'] = [lowest_support_val]
        lineage_out_dict['poorly supported markers'] = ['; '.join(poorly_supported_markers)]
        lineage_out_dict['max support for additional markers'] = [non_matching_support]
        lineage_out_dict['additional markers'] = ['; '.join(non_matching_markers)]
        lineage_out_dict['node support'] = ['; '.join(final_markers)]
    # add species info
    lineage_out_dict['species']=['S. sonnei']
    lineage_out_df = pd.DataFrame(lineage_out_dict, columns=['genome', 'species', 'final genotype', 'name', 'confidence', 'lowest support for genotype marker', 'poorly supported markers', 'node support', 'max support for additional markers', 'additional markers'])

    return lineage_out_df, spp_call

def main():
    print("doing something")
    args = get_arguments()
    # list of tables for each result type
    results_tables = []

    # want a single table of outputs
    # create lineage_name_dict, key=mykrobe lineage name, value=human readable name
    lineage_name_dict = {}
    with open(args.alleles, 'r') as lineage_names:
        for line in lineage_names:
            fields = line.strip().split('\t')
            lineage_name_dict[fields[4]] = fields[3]

    # read in json files
    for json_file in args.jsons:
        with open(json_file) as f:
            #with open("ERR017671.mykrobe.json") as f:
            myk_result = json.load(f)
        # get genome name (should first and ONLY key at start of json file)
        if len(list(myk_result.keys())) > 1:
            print("More than one result in mykrobe output file " + json_file + ", quitting")
            sys.exit()
        genome_name = list(myk_result.keys())[0]
        # extract all the data for this genome
        genome_data = myk_result[genome_name]
        # extract the genotyping information
        lineage_data = genome_data["phylogenetics"]
        lineage_table, spp_call = extract_lineage_info(lineage_data, genome_name, lineage_name_dict)
        #all_lineage_tables.append(lineage_table)
        # extract the qrdr information, only if sonnei
        genome_qrdr_table = extract_qrdr_info(genome_data, genome_name, spp_call)
        #all_qrdr_tables.append(genome_qrdr_table)
        # combine the two tables together
        all_info_table = pd.merge(lineage_table, genome_qrdr_table, on="genome", how="inner")
        results_tables.append(all_info_table)


    # concatenate, re-order columns, and write out
    final_results = pd.concat(results_tables, sort=True)
    final_results.to_csv(args.prefix + "_predictResults.tsv", index=False, sep="\t", columns=["genome", "species", "final genotype", "name", "confidence", "num QRDR", "parC_S80I", "gyrA_S83L", "gyrA_S83A", "gyrA_D87G", "gyrA_D87N", "gyrA_D87Y", "lowest support for genotype marker", "poorly supported markers", "max support for additional markers", "additional markers", "node support"])

if __name__ == '__main__':
    main()
