# Pipeline overview

| Pipeline name | Speciation | Sequence assessment | Assemble | MLST | Serotyping | AMR | Plasmid detection | Tree | Distances/Comparative analysis | Pangenome|
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |:---: |
| [basic](../pipelines/basic.md) | Optional | Yes | No | No | No | No | No | No | No | No |
| [preview](../pipelines/preview.md) | Optional | Yes | No | No | No | No | No | Yes (default - mash) | Yes (default - mash) |No |
| [assemble](../pipelines/assemble.md) | Optional | Yes | Yes | No | No | No | No | No | No | No |
| [amr_typing](../pipelines/amr_typing.md) | Optional | Yes | Yes (if required) | Yes | Yes (if speciation is on) | Yes | Yes | No | No | No |
| [comparative](../pipelines/comparative.md) | Optional | Yes | No | No | No | No | No | Yes (default - iqtree) | Yes (default - snippy)| No |
| [full](../pipelines/full.md) | Optional | Yes | Yes (if required) | Yes | Yes (if speciation is on)| Yes | Yes | Yes (default - iqtree) | Yes (default - snippy)| Yes |
| [tb](../pipelines/tb.md) | Optional | Yes | No | No | No | Yes | No | Yes (default - iqtree) | Yes (default - snippy)| Yes |