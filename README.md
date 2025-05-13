# project-template

This is a template for R analysis projects for the Raab lab, forked from the Phantsiel lab.

## How to use this template

1.  Create a repository from this template according to these instructions, make sure to select 'raab-lab' as owner so everyone in the lab can have access: <https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-repository-from-a-template>

2.  Then clone your new repository locally:

    ``` bash
    git clone https://github.com/{YourUserName}/{YourRepoName}.git
    ```

3.  Rename the `project-template.Rproj` file to `{YourRepoName}.Rproj`


## What to put in `README.md`

The README contains information about your project. Here you can describe your analyses, processing steps, or approaches. Most importantly, the README should contain information about your raw data. Describe 1) how the raw data was generated, 2) where it was obtained, and 3) how it can be accessed for others. This is the perfect place to describe steps that were carried out on the command line or on a remote cluster.

## Snakefile

The `Snakefile` is like a recipe book that describes the desired output files, the input files used to make them and the instructions for transforming inputs into outputs. Everytime you create an output file of any type, make sure to add it to the `Snakefile`. This ensures that changes to your scripts will produce the most up-to-date outputs after running `snakemake`.

Snakemake comes packaged with the longleaf Anaconda module. Run `module load anaconda` to access it.

## Directory organization

The project template is laid out with the following example structure.
Feel free to remove or add as you feel necessary for you project

```
project-template/
├── data
│   ├── derived_data		<- 'Final' analysis outputs, like differential binding/expression tables.
│   │   ├── atac
│   │   ├── cnr
│   │   ├── integrate
│   │   └── rna
│   ├── external		<- Location for potential public or non-pipeline processed data
│   ├── processed_data		<- Intermediate datasets, like summarizedExperiment RDS files or consensus peak calls
│   │   ├── atac
│   │   ├── cnr
│   │   └── rna
│   └── source_data		<- Here is where you should link in data from pipeline runs
│       ├── atac
│       ├── cnr
│       └── rna
├── figures			<- Any graphical outputs -- boxplots, scatterplots, etc
│   ├── atac
│   ├── cnr
│   ├── integrate
│   └── rna
├── project-template.Rproj	<- Rename this to your project name
├── README.md
├── renv
├── renv.lock			<- For use with the `renv` package. Helps everyone stay in the same R environment
├── Snakefile			<- See above
└── src				<- All code necessary to generate figures and derived data
    ├── atac
    │   └── clean_atac.R
    ├── cnr
    │   └── clean_cnr.R
    ├── create_canonical_sheet.R
    ├── integrate
    │   └── integrate.R
    └── rna
        └── clean_rna.R
```
