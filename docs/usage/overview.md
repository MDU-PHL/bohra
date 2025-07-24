# Usage overview

## Input files

`bohra` requires a single tab-delimited file as input, you can find examples [here] TO ADD.

| Column name | Description | Required? | 
| :---: | :---: | :---: | 
| Isolate | This is the name of the sequence or sample and will appear throughout `bohra` outputs. It must be unique.| Yes|
| r1 | The path to read 1 | If an assembly file is not supplied you must supply reads | 
| r2 | The path to read 2 | If an assembly file is not supplied you must supply reads |
| assembly | The path to the assembly for the isolate | If reads are not supplied you must supply an assembly file |
| species | The expected species of the sample or 'control'.  | No |

### Species column

If you do not require speciation as part of the pipeline and already know the species, you can provide it here. Please note if no speciation is undertaken, `bohra` will use this value to undertake typing and AMR mechanisms/inferrence. If the species in this column is NOT accurate - unexpected results may occur. Furthermore if your analysis includes control sequences you can provide that information here (`control`) and the sequence will not be included in any comparative analysis.

### Annotation

Where you are undertaking a comparative analysis (`snippy`, `ska2`, `mash`) you may also provide additional columns of relevant metadata in your input file. `bohra` will do data validation on these columns - that is up to the user. But any additional metadata provided will be visible on the tree provided in the report file.

 
## Pipelines

`bohra`

### Basic

### Assemble

### Comparative

### AMR and typing

### TB

### Default

### Full

### C. auris

COMING SOON