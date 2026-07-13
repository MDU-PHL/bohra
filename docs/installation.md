# Installation



## Recommended (conda or mamba)

**1. Install conda (skip this step if you already have conda installed)**

If you do not already have `conda` installed, you can check out the documentation [here](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html). We recommend you install [`miniforge`](https://github.com/conda-forge/miniforge)


**2. Create and install bohra**

```
mamba (or conda) create -n bohra -c bioconda bohra
```

**3. Install dependencies**

```
conda activate bohra
bohra deps install
```

**4. Testing your installation**

If you would like to make sure that all the dependencies are installed properly and that the pipeline will run on your system you can run


```
bohra test --cpus X
```

(where X is the number of cpus you want to use)

This will download some fastq files and run the full pipeline.


## Install from source


**Setup the `bohra` environment and install the pipeline co-ordinator**

1. Clone the `bohra` repository

```
git clone https://github.com/MDU-PHL/bohra.git
```


. Create the environment and install `bohra`


```
cd bohra
conda env create -f environment.yml -n bohra
conda activate bohra
pip3 install .
```

**Install dependencies and setup environment variables**


This step will setup conda environments under the `path/to/conda/envs/bohra`(depending on how you have configured your `conda` installation). 


1. Activate your environment from step 2 above
```
conda activate bohra
```

2. Run the bohra dependency installation

```
bohra  deps install
```

## Updating `bohra`

Updates to the `bohra` pipeline  will often include 2 steps/

1. **Update `bohra` env**

```
conda update -n bohra --all
```

This will update the bohra pipeline and running scripts. But not any dependencies

2. **Update any changed dependencies**

```
conda activate bohra
bohra deps update
```

This will update any dependencies that have pinned versions or any new dependencies and should be run each time you update `bohra`.

### Other update options

You may also update a specific tool or dependency environment. This should be done with caution as updates that have not been adeqautely incorporated into the pipeline may cause unexpected behaviour or breakages. 

```
conda activate bohra
bohra deps update --tool <env_name> --force
```

If you would like another tool or different versions please leave an issue and we will endeavour to facillitate it.

If you need to re-install your dependencies for some reason you can run

```
conda activate bohra
bohra deps update --force
```
This will overwrite all existing dependencies and may take some time.