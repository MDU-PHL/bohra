# Running `bohra`

`bohra` has 6 pipelines to choose from

1. `basic`
2. `assemble`
3. `amr_typing`
4. `comparative`
5. `full`
6. `tb`

You can see available pipelines
```
bohra run  --help
```
and help for each pipeline

```
bohra run <pipeline> --help
```
For any pipeline you can elect to run speciation using `--speciation sylph` or `--speciation kraken2` or you can turn off speciation by using `--speciation none`.


## Genome characterisation pipelines

These first three pipelines will generate single sample results, such as sequence metrics, typing and recovery of AMR mechanisms. These pipelines can be run on a variety of species in the same run, ie for quality control purposes or where only sample level results are required.

### basic

This pipeline is a basic sequence assessment pipeline. It can be run as a standalone to assess your sequence data, but is also run as part of all other pipelines.

```
bohra run basic -i input_file.tsv -j my_basic_pipeline
```
where
- `-i/--input_file` is a tab-delimited file formatted as described [here](../usage/overview.md)
- `-j/--job_id` is the name of your run. This value will appear on your report.

### assemble

This pipeline should be used if you would like to generate _de novo_ assemblies from paired-end reads. You can choose from `spades`, `skesa`, `shovill + skesa` or `shovill + spades` (default)

```
bohra run assemble -i input_file.tsv -j my_assembly_pipeline -a shovill_skesa
```
where
- `-i/--input_file` is a tab-delimited file formatted as described [here](../usage/overview.md)
- `-j/--job_id` is the name of your run. This value will appear on your report.
- `-a/--assembler` is the assembler to use (ONT coming soon)

### amr_typing

This pipeline will use `abritamr` for AMR mechanism detection and undertake species appropriate serotyping. Please note that if you do not supply a species and set `--speciation none` there will be no species specific AMR or serotyping done.
You can provide paired-end fastq and/or assemblies as input for this pipeline. If only paired-end fastq supplied, assembly will be run to generata appropriate inputs for `abritamr` and serotyping. If you would like to use a an assembly combination different from the default you will need to specify

1. where inputs are only (or mostly paired-end fastq) and you want to use a different assembler to the default of `shovill + spades`.
```
bohra run amr_typing -i input_file.tsv -j my_typing_pipeline -a shovill_skesa
```
2. where inputs are assemblies (or you are happy to stick with `shovill + spades` )
```
bohra run amr_typing -i input_file.tsv -j my_typing_pipeline
```

## Comparative pipelines

These are pipelines where dataset wide comparasions are made. For example species wide phylogenetic analysis, outbreak investigations or other investigations of the relationships between sequences in a dataset. Where you have included mixed species, please be aware you may get unexpected and frankly wrong interpretations.

## comparative

This pipeline invloves undertaking comparisons within a dataset. There are three tools available in `bohra` that can be used for this purpose
1. `snippy` (default) - a reference based approach, detecting SNPS compared to a reference genome. This option can ONLY take paired-end fastq files as input and requires the provision of a reference genome.
2. `ska2` - this is a reference free approach, inferring SNPS based on splits in kmers. You can use paired-end fastq and/or _de novo_ assemblies as input.
3. `mash` - this is a quick reference free approach which allows for approximation of genetic distances between sequences. Like `ska2` you can provide paired-end fastq and/or _de novo_ assemblies as input.