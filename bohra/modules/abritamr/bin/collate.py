#!/usr/bin/env python3
import pathlib, pandas, math, sys,  re


class Collate:

    """
    a base class for collation of amrfinder results - to be used when not doing MDU QC
    """
    

    REFGENES = pathlib.Path(__file__).parent / "db" / "refgenes_latest.csv"
    MATCH = ["ALLELEX", "BLASTX", "EXACTX", "POINTX"]
    NONRTM = [
        "Amikacin/Gentamicin/Kanamycin/Tobramycin",
        "Amikacin/Kanamycin",
        "Amikacin/Kanamycin/Tobramycin",
        "Amikacin/Quinolone",
        "Amikacin/Tobramycin",
        "Aminoglycosides",
        "Gentamicin",
        "Gentamicin/Kanamycin/Tobramycin",
        "Gentamicin/Tobramcyin",
        "Kanamycin",
        "Kanamycin/Tobramycin",
        "Spectinomycin",
        "Streptogramin",
        "Streptomycin",
        "Streptomycin/Spectinomycin",
        "Tobramycin",
    ]
    MACROLIDES = [
        "Erythromycin",
        "Erythromycin/Telithromycin",
        "Lincosamides",
        "Lincosamides/Streptogramin",
        "Macrolides",
        "Streptogramin",
    ]
    RTM = "Other aminoglycoside resistance (non-RMT)"
    MAC = "Macrolide, lincosamide and/or streptogramin resistance"
    
    
    def __init__(self, isolate, amr_output):
        self.workdir = pathlib.Path.cwd()
        self.isolate = isolate
        self.amr_output = amr_output

    def joins(self, dict_for_joining):
        """
        make them a comma separated list
        """
        for i in dict_for_joining:
            if i != "Isolate":
                dict_for_joining[i] = list(set(dict_for_joining[i]))
                dict_for_joining[i] = ",".join(dict_for_joining[i])

        return dict_for_joining

    def get_drugclass(self, reftab, row, colname):

        """
        if the enhanced subclass is in either NONRTM or MACROLIDES then then use the groups specified by Norelle. If it is empty (-) then fall back on the AMRFinder subclass, else report the extended subclass
        """
        # print(row)
        gene_id_col = "Gene symbol" if colname != "refseq_protein_accession" else "Accession of closest sequence"
        # print(gene_id_col)
        # print(reftab[reftab[colname] == row[1][gene_id_col]])
        d = reftab[reftab[colname] == row[1][gene_id_col]]['enhanced_subclass'].values[0]
        # print(d)
        if d in self.NONRTM:
            return self.RTM
        elif d in self.MACROLIDES:
            return self.MAC
        elif d == "-":
            return reftab[reftab[colname] == row[1]["Gene symbol"]][
                "subclass"
            ].unique()[0]
        else:
            return d

    def extract_bifunctional_name(self, protein, reftab):
        """
        extract the joint name of bifunctional genes
        """
        return reftab[reftab["refseq_protein_accession"] == protein]["gene_family"].values[0]

    def extract_gene_name(self, protein, reftab):

        if reftab[reftab["refseq_protein_accession"] == protein]["allele"].values[0] != '-':
            return reftab[reftab["refseq_protein_accession"] == protein]["allele"].values[0]
        else:
            return reftab[reftab["refseq_protein_accession"] == protein]["gene_family"].values[0]
            
    def setup_dict(self, drugclass_dict, reftab, row):
        """
        return the dictionary for collation
        """
        # drugname = 'x'
        # print(row[1])
        if row[1]["Gene symbol"] in list(reftab["allele"]):
            # print('gene symbol is an allele')
            drugclass = self.get_drugclass(
                    reftab=reftab, row=row, colname="allele"
                    )
            drugname = self.extract_gene_name(protein = row[1]["Accession of closest sequence"], reftab = reftab)

        elif row[1]["Gene symbol"] in list(reftab["gene_family"]):
            drugclass = self.get_drugclass(
                reftab=reftab, row=row, colname="refseq_protein_accession"
            )
            drugname = f"{self.extract_gene_name(protein = row[1]['Accession of closest sequence'], reftab = reftab)}*" if not row[1]["Method"] in ["EXACTX", "ALLELEX"] else f"{self.extract_gene_name(protein = row[1]['Accession of closest sequence'], reftab = reftab)}"
            
        elif row[1]["Accession of closest sequence"] in list(reftab["refseq_protein_accession"]):
            drugclass = self.get_drugclass(
                reftab = reftab, row = row, colname = "refseq_protein_accession"
            )
            drugname = self.extract_bifunctional_name(protein = row[1]['Accession of closest sequence'], reftab = reftab)
        else:
            drugname = row[1]["Gene symbol"]
            drugclass = "Unknown"

        if drugclass in drugclass_dict:
            drugclass_dict[drugclass].append(drugname)
        elif drugclass not in drugclass_dict:
            drugclass_dict[drugclass] = [drugname]

        return drugclass_dict

    def get_per_isolate(self, reftab, df, isolate):
        """
        make two dictionaries for each isolate that contain the drug class assignments for each match that is one of ALLELEX, EXACTX or BLASTX and then another dictionary which lists all partial mathces
        """
        drugclass_dict = {"Isolate": isolate}
        partials = {"Isolate": isolate}
        
        for row in df.iterrows():
            # print(row[1]["Gene symbol"])
            
            if row[1]['Element subtype'] == 'POINT':
                print(row[1])
            # if the match is good then generate a drugclass dict
            if row[1]["Gene symbol"] == "aac(6')-Ib-cr" and row[1]["Method"] in ["EXACTX", "ALLELEX"]: # restrict the calling of quinolone resistance - if not an exact match then is should be considered a partial match
                partials = self.setup_dict(drugclass_dict = partials, reftab = reftab, row = row)
            elif row[1]["Method"] in self.MATCH and row[1]["Scope"] == "core" and row[1]["Element type"] == "AMR":
                drugclass_dict = self.setup_dict(drugclass_dict = drugclass_dict, reftab = reftab, row = row)
            elif row[1]["Method"] not in self.MATCH and row[1]["Scope"] == "core" and row[1]["Element type"] == "AMR":
                partials = self.setup_dict(drugclass_dict = partials, reftab = reftab, row = row)
        
        drugclass_dict = self.joins(dict_for_joining=drugclass_dict)
        partials = self.joins(dict_for_joining=partials)
        
        return drugclass_dict, partials

    def get_amr_output(self, path):
        """
        check that AMR finder output is present in the correct format.
        """
        amr_output = sorted(path.glob("*.out"))
        
        if amr_output != []:
            return amr_output
        else:
            print(
                "You do not have any AMR finder output. Please check your settings and try again."
            )
            raise SystemExit

    def save_files(self, path, tosave):
        # print(tosave)
        tosave.to_csv(path, sep = '\t')

    def collate(self):
        """
        if the refgenes.csv is present then proceed to collate data and save the csv files.
        """
        
    
    
        if self.REFGENES.exists():
            reftab = pandas.read_csv(self.REFGENES)
            reftab = reftab.fillna("-")
            header = list(reftab['enhanced_subclass'].unique())
            # print(reftab.head())
            # p = pathlib.Path.cwd()
            # print(p)
            # amr_output = self.get_amr_output(path=p)
            # print(amr_output)
            summary_drugs = pandas.DataFrame(columns = header)

            summary_partial = pandas.DataFrame(columns = header)
            # print(summary_drugs)
            # for a in amr_output:
            
            # print(a)
            # print(isolate)
            try:
                df = pandas.read_csv(self.amr_output, sep="\t")
                # print(df)
                drug, partial = self.get_per_isolate(
                    reftab=reftab, df=df, isolate=self.isolate
                )
                temp_drug = pandas.DataFrame(drug, index=[0])
                # print(temp_drug)
                temp_partial = pandas.DataFrame(partial, index=[0])
                # if summary_drugs.empty:
                #     summary_drugs = temp_drug
                # else:
                summary_drugs = summary_drugs.append(temp_drug)

                # if summary_partial.empty:
                #     summary_partial = temp_partial
                # else:
                summary_partial = summary_partial.append(temp_partial)
                summary_drugs = summary_drugs.set_index("Isolate")
                summary_partial = summary_partial.set_index("Isolate")
                
            
                return summary_drugs, summary_partial
            except:
                print(f"Something has gone wrong with {self.isolate} there appears to be no genes or mutations.")        
        else:
            print(f"The refgenes DB ({self.REFGENES}) seems to be missing.")
            raise SystemExit

    def run(self):

        summary_drugs, summary_partial = self.collate( )

        self.save_files(path="summary_matches.txt", tosave=summary_drugs)
        self.save_files(path="summary_partials.txt", tosave=summary_partial)


class MduCollate(Collate):
    def __init__(self, qc, db, runid, matches, partials):
        self.workdir = pathlib.Path.cwd()
        self.mduqc = self.workdir / f"{qc}"
        self.db = db
        self.check_for_mduqc()
        self.mduqctab = self.mdu_qc_tab()
        self.passed_partials = partials
        self.passed_match = matches
        self.runid = runid

        self.NONE_CODES = {
            "Salmonella":"CPase_ESBL_AmpC_16S_NEG",
            "Shigella":"CPase_ESBL_AmpC_16S_NEG",
            "Staphylococcus":"Mec_VanAB_Linez_NEG",
            "Enterococcus":"Van_Linez_NEG",
            "Other":"Cpase_16S_mcr_NEG"
        }

    def mdu_qc_tab(self):
        
        tab = pandas.read_csv(self.mduqc)
        
        pos = pandas.DataFrame(data = {"ISOLATE": "9999-99888", "TEST_QC" : "PASS", "SPECIES_EXP":"Staphylococcus aureus", "SPECIES_OBS":"Staphylococcus aureus" }, index = [0])
        
        return tab.append(pos)

    def check_for_mduqc(self):

        if self.mduqc.exists():
            return True
        else:
            print(
                f"It appears you are running mdu-amr in the context of MDU QC, the {self.mduqc} file is not present in {self.workdir}. Please check your settings and try again."
            )
            raise SystemExit

    def strip_bla(self, gene):
        '''
        strip bla from front of genes except
        '''
        if gene.startswith("bla") and len(gene) >6 and gene.endswith("*"):
            gene = gene.replace("bla", "")
        elif gene.startswith("bla") and len(gene) >5 and not gene.endswith("*"):
            gene = gene.replace("bla", "")
        return gene

    def get_passed_isolates(self, qc_name):
        """
        generate lists of isolates that passed QC and need AMR, failed QC and should have AMR and all other isolates
        """
        
        # print(self.mdu_qc_tab)
        
        failed = list(
            self.mduqctab[self.mduqctab["TEST_QC"] == False]["ISOLATE"]
        )  # isolates that failed qc and should have amr
        passed = list(
            self.mduqctab[self.mduqctab["TEST_QC"] == True]["ISOLATE"]
        )  # isolates that failed qc and should have amr
        # print(passed)
        return (passed, failed)


    # Carbapenemase to be reported for all cases
    # Carbapenemase (MBL) all in all HOWEVER if blaL1 should  not be reported in Stenotrophomonas maltophilia
    # Carbapenemase (OXA-51 family) REPORTED IN ALL except in Acinetobacter baumannii,Acinetobacter calcoaceticus,Acinetobacter nosocomialis,Acinetobacter pittii,Acinetobacter baumannii complex,
    # All ESBL in Salmonella and Shigella
    # ESBL (Amp C type) in Salmonella and Shigella
    # Aminoglycosides (Ribosomal methyltransferases) ALL
    # Colistin ALL
    # Oxazolidinone & phenicol resistance if Genus = Enterococcus or Staphylococcus aureus and Staphylococcus argenteus
    # report vanA*, B*, C*, vanD*, vanE*, vanG*, vanL*, vanM*, vanN* 
    # Methicillin ALL

    def get_all_genes(self, row):
        all_genes = []
        for r in row[1]:
            if isinstance(r, str):
                if len(r.split(",")) > 1:
                    for a in r.split(","):
                        all_genes.append(a)
                else:
                    all_genes.append(r)
        return all_genes

    def none_replacement_code(self, genus):

        if genus in self.NONE_CODES:
            return self.NONE_CODES[genus]
        else:
            return "GENRL_AMR_NEG1"


    def reporting_logic(self, row, species, neg_code = True):
        # get all genes found
        all_genes = self.get_all_genes(row)
        # print(row)
        isodict = row[1].to_dict()
        # determine the genus EXPECTED
        genus = species.split()[0]
         # A list of reportable genes - TODO move to a class variable
        reportable = [
            "Carbapenemase",
            "Carbapenemase (MBL)",
            "Carbapenemase (OXA-51 family)",
            "ESBL",
            "ESBL (AmpC type)",
            "Aminoglycosides (Ribosomal methyltransferase)",
            "Colistin",
            "Oxazolidinone & phenicol resistance",
            "Vancomycin",
            "Methicillin"
        ]
        # print(reportable)
        non_caveat_reportable = [
            "Carbapenemase",
            "Aminoglycosides (Ribosomal methyltransferase)",
            "Colistin"
        ]

        abacter_excluded = [
            "Acinetobacter baumannii",
            "Acinetobacter calcoaceticus",
            "Acinetobacter nosocomialis",
            "Acinetobacter pittii",
            "Acinetobacter baumannii complex"
        ]

        
        van_match = re.compile("van[A,B,C,D,E,G,L,M,N][\S]*")
        mec_match = re.compile("mec[^IR]")
        
        genes_reported = []  # genes for reporting
        genes_not_reported = []  # genes found but not reportable
        for i in isodict:
            # print(i)
                        
            genes = []
            if not isinstance(isodict[i], float):
                genes = isodict[i].split(',')
                # print(genes)
            # print(isodict[i])
            if genes != []: # for each bin we do things to genes
                
                if i in reportable:
                    
                    # print(isodict[i])
                    if i in non_caveat_reportable:
                        genes_reported.extend(genes)
                    elif i == "Carbapenemase (MBL)" and species != "Stenotrophomonas maltophilia":
                        genes_reported.extend(genes)
                    elif i == "Carbapenemase (MBL)" and species == "Stenotrophomonas maltophilia":
                         # if species is "Stenotrophomonas maltophilia" don't report blaL1
                        genes_reported.extend([g for g in genes if not g.startswith("blaL1")])
                        genes_not_reported.extend([g for g in genes if g.startswith("blaL1")])
                    elif i == "Carbapenemase (OXA-51 family)" and species not in abacter_excluded:
                        genes_reported.extend(genes)
                    elif i in ["ESBL","ESBL (AmpC type)"] and genus in ["Salmonella"]: 
                        genes_reported.extend(genes)
                    elif i in ["ESBL","ESBL (AmpC type)"] and genus in ["Shigella"]: 
                        genes_reported.extend([g for g in genes if "blaEC" not in g])
                        genes_not_reported.extend([g for g in genes if "blaEC" in g]) # don't report blaEC for shigella
                    elif i == "Oxazolidinone & phenicol resistance":
                        if species in ["Staphylococcus aureus","Staphylococcus argenteus"] or genus == "Enterococcus":
                            genes_reported.extend(genes)
                        else:
                            genes_not_reported.extend(genes)
                    elif i == "Vancomycin":
                        genes_reported.extend([g for g in genes if van_match.match(g)])
                        genes_not_reported.extend([g for g in genes if not van_match.match(g)])
                    elif i == "Methicillin":
                        genes_reported.extend([g for g in genes if mec_match.match(g)])
                        genes_not_reported.extend([g for g in genes if not mec_match.match(g)])
                    else:
                        genes_not_reported.extend(genes)

                else:
                    genes_not_reported.extend(genes)
            # print(genes_reported)
        if genes_reported == []:
            genes_reported = [self.none_replacement_code(genus= genus)] if neg_code else ''
        if genes_not_reported == []:
            genes_not_reported = ["No non-reportable genes found."] if neg_code else ''
            # break
    
        return genes_reported, genes_not_reported


    def assign_itemcode(self,mduid, reg):
        m = reg.match(mduid)
        try:
            itemcode = m.group('itemcode') if m.group('itemcode') else ''
        except AttributeError:
            itemcode = ''
        return itemcode

    def assign_mduid(self, mduid, reg):
        m = reg.match(mduid)
        try:
            mduid = m.group('id')
        except AttributeError:
            mduid = mduid.split('/')[-1]
        return mduid

    def mdu_reporting(self, match, neg_code = True):

        mduidreg = re.compile(r'(?P<id>[0-9]{4}-[0-9]{5,6})-?(?P<itemcode>.{1,2})?')
        reporting_df = pandas.DataFrame()
        qc = pandas.read_csv(self.mduqc, sep=None, engine="python")
        qc = qc.rename(columns={qc.columns[0]: "ISOLATE"})
        match_df = pandas.read_csv(match, sep = '\t')
        for row in match_df.iterrows():
            isolate = row[1]['Isolate']
            item_code = self.assign_itemcode(isolate, mduidreg)
            md = self.assign_mduid(isolate, mduidreg)
            d = {"MDU sample ID": md, "Item code" : item_code}
            qcdf = self.mduqctab[self.mduqctab["ISOLATE"].str.contains(isolate)]
            exp_species = qcdf["SPECIES_EXP"].values[0]
            obs_species = qcdf["SPECIES_OBS"].values[0]
           
            species = obs_species if obs_species == exp_species else exp_species
            genes_reported, genes_not_reported = self.reporting_logic(
                row=row, species=species, neg_code=neg_code
            )
            # strip bla
            genes_not_reported = [self.strip_bla(g) for g in genes_not_reported]
            genes_reported = [self.strip_bla(g) for g in genes_reported]
            d["Resistance genes (alleles) detected"] = ",".join(genes_reported)
            d["Resistance genes (alleles) det (non-rpt)"] = ",".join(genes_not_reported)
            if qcdf["TEST_QC"].values[0] == 'FAIL': # species not needed for MDU LIMS upload
                d["Species_exp"] = exp_species
            d["Species_obs"] = obs_species
            d["Species_exp"] = exp_species
            d['db_version'] = self.db
            tempdf = pandas.DataFrame(d, index=[0])
            tempdf = tempdf.set_index("MDU sample ID")
            # print(tempdf)
            if reporting_df.empty:
                reporting_df = tempdf
            else:
                reporting_df = reporting_df.append(tempdf)
        
        return reporting_df.reindex(labels = ['Item code','Resistance genes (alleles) detected','Resistance genes (alleles) det (non-rpt)','Species_obs', 'Species_exp', 'db_version'], axis = 'columns')

    def mdu_partials(self, partials):
        '''
        split the partial matches into reportable/nonreportable
        input is a dataframe with partial type as colnames 
        '''
        partial_names = ['PARTIALX', 'PARTIAL_CONTIG_ENDX']
        partial_df = pandas.DataFrame()
        for row in partials.iterrows:
            qc = pandas.read_csv(self.mduqc, sep=None, engine="python")
            qc = qc.rename(columns={qc.columns[0]: "ISOLATE"})
            exp_species = self.mduqctab["SPECIES_EXP"].values[0]
            obs_species = self.mduqctab["SPECIES_OBS"].values[0]
            species = obs_species if obs_species == exp_species else exp_species
            d = {"MDU sample ID": row[0]}
            
    def save_spreadsheet(
        self,
        passed_match,
        passed_partials,
        
    ):
        print(self.runid)
        writer = pandas.ExcelWriter(f"{self.runid}_MMS118.xlsx", engine="xlsxwriter")

        
        passed_match.to_excel(writer, sheet_name="MMS118")
        passed_partials.to_excel(writer, sheet_name="Passed QC partial")
        

        writer.close()


    def get_single_amr(self, amr_file):
        '''
        get a single amr output
        '''
        if len(amr_file) > 0:
            data = pandas.read_csv(amr_file[0], sep = "\t")
            return data
        
   

    def run(self):

        passed_match_df = self.mdu_reporting(match=self.passed_match)
        passed_partials_df = self.mdu_reporting(match = self.passed_partials, neg_code=False)
        # failed_match_df = self.mdu_reporting(match = failed_match)
        # failed_partials_df = self.mdu_reporting(match = failed_partials, neg_code=False)
        self.save_spreadsheet(
            passed_match_df,
            passed_partials_df,

        )

if __name__ == "__main__":
# $params.qc $params.db $params.runid $abritamr_matches $abritamr_partials
    if len(sys.argv) > 4:
        c = MduCollate(qc = f"{sys.argv[1]}", db = f"{sys.argv[2]}", runid = f"{sys.argv[3]}", matches = f"{sys.argv[4]}", partials = f"{sys.argv[5]}")
    else:
        c = Collate(isolate = f"{sys.argv[1]}", amr_output = f"{sys.argv[2]}")
    c.run()

