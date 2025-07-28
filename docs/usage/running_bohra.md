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

All comparative pipelines have the option to cluster the distances using a heirarchical approach. You can specifiy the algorithm with `--cluster_method` (choose from `single`, `average`, `complete`, `centroid`, `median`, `ward`, `weighted`), with `single` being the default setting. You can also choose the thresholds that you want to use as a comma separate list with `--cluster_threshold` (not that the algorithms use a `<` so if you provide a theshold of 10 - the cluster will actually defined at `<= 9` ).

### comparative

This pipeline invloves undertaking comparisons within a dataset. Please take note - it will not run any assembly based tools, like MLST or AMR. If these are required use the `full` pipeline. There are three tools available in `bohra` that can be used for this purpose
1. `snippy` (default) - a reference based approach, detecting SNPS compared to a reference genome. This option can ONLY take paired-end fastq files as input and requires the provision of a reference genome.
2. `ska2` - this is a reference free approach, inferring SNPS based on splits in kmers. You can use paired-end fastq and/or _de novo_ assemblies as input.
3. `mash` - this is a quick reference free approach which allows for approximation of genetic distances between sequences. Like `ska2` you can provide paired-end fastq and/or _de novo_ assemblies as input.

Further customisation includes selection of the data which is used to generate the tree (for `snippy` and `ska2`)

1. distances - you can generate a distance based tree, where the branch lengths are simply the SNP distances between sequences
2. alignment - this will generate a maximum liklihood phylogentic tree, using the GTR model of phylogenetic inference.

Additionally - you can select the tree builder to use. `VeryFastTree` is the default tree builder - as it is very quick. But if required you can also use `IQtree`.

1. `IQtree`
2. `VeryFastTree` (default)


For example

**snippy (alignment and iqtree) and cluster at a threshold of <= 5 and <= 25**
```
bohra run comparative -i input_file.tsv -j my_snippy_pipeline -ref <path_to_reference.fa(gbk)> --tree_builder iqtree --cluster_threshold 6,26
```
**ska2 (distance and veryfastree) and cluster at a threshold of <= 5 and <= 25, using complete linkage**
```
bohra run comparative -i input_file.tsv -j my_snippy_pipeline --comparative_tool ska2 --cluster_method complete --cluster_threshold 6,26
```

### full

Like the `comparative` pipeline, you can select the comparative tool, tree builder, cluster methods, and thresholds. In addition, this pipeline will undertake the assembly, amr and typing described above and additional pangenome analysis.

For example

```
bohra run full -i input_file.tsv -j my_snippy_pipeline -ref <path_to_reference.fa(gbk)> --tree_builder iqtree --cluster_threshold 6,26
```
