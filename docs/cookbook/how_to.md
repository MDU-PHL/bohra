# How do I

## Generate an input file?

```
bohra generate-input --reads /path/to/reads --contigs /path/to/contigs --outname my_data.txt
```
(you can rename your input as you please - the default input file if no name is supplied is `bohra_input.tsv`)

- You will need a folder with all of your sequences in it.
- `bohra` will match read 1 and read 2 and use the string before the first `_`  as the isolate identifier.
- If your input sequences are assemblies, `bohra` will use the string before the first `.` as the isolate identifier. 
- Note if these do not match (ie)

## Run `mash` to assess dataset?

```
bohra run preview -i input_file.tsv -j my_basic_pipeline --cpus N --report_outdir your_report
```

## Run the full pipeline?

```
bohra run full -i input_file.tsv --cpus N -ref your_reference.fa
```

## Change the cluster threshold?
```
bohra run full -i input_file.tsv --cpus N -ref your_reference.fa -ct 6,11,21
```

## Use ska instead of `snippy`?
```
bohra run full --comparative_tool ska -i input_file.tsv --cpus N 
```

## Just run amr and serotype?
```
bohra run amr_typing -i input_file.tsv --cpus N
```

## Remove recombination?

```
bohra run full -i input_file.tsv --cpus N -ref your_reference.fa --gubbins
```