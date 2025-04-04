# RNASeq Quantification Workflow
This workflow is meant for use originally for isolates with STAR using a transcriptome; however additions have made it now possible to use Kallisto. If running on a MAG sample you may wish to take an alignmnet approach instead of using this method (especially if the annotation set is limited). In the case of an limited annotation set, you may wish to annotate using [anvi'o](https://anvio.org/) which provides a much more extensive annotation set then other tools tested (as of April, 2025) prior to using this workflow. 

### Installation and Environment Set Up
This workflow expects that you have [conda](https://docs.conda.io/en/latest/miniconda.html) installed prior to starting. Conda is very easy to install in general and will allow you to easily install the other dependencies needed in this workflow. then you will need to download this repository via git clone.

1. Activate your conda environment so that you are in the "base" environment. This can be achieved by running a source /path/to/conda/install/folder/bin/activate
2. Make sure you have your channels set up to allow conda-forge by running `conda config --add channels conda-forge`
3. Install snakemake into a new environment by running `conda create -n rnaquant -c bioconda snakemake=7.25.4`. Then activate it by running `conda activate rnaquant`.
4. Move into the directory `cd rnaseq-quantification-pipeline` and if desired generate the test data by running `cd test && snakemake --use-conda --cores 4 -j`. The md5sum may say the files are not downloaded on the first go; however, if you run it again the files should all be OK. Move back to the rnaseq-quantification-pipeline folder after this is done, then run `snakemake --use-conda --cores 2 -n` to test if everything is installed correctly. If not, then you'll need to trouble shoot your environment. The

That's it! If you have downloaded snakemake, the rest of the dependencies will be downloaded automatically for you via the workflow.

### Prepare the input files
In this pipeline we are using PEP which allows for easier portability between projects. You need to alter some of the files and generate two `.tsv`.

1. Edit the config/pepconfig.yml file's `raw_data:` string to be where your files are located at. Make sure to use the full path to avoid any errors.
2. We need to generate the samples and subsamples files. subsamples is just where the paired end files are input, samples is an overview file. Each row contains a different sample and it's related information. First create your sample sheet and make sure it ends with `.csv`. 
`sample_name`: The name of the sample without the ending. If it can't match it then you'll get an error.
`alternate_id`: Really doesn't matter, but if there is a different ID that the sample is under add it here. Otherwise, just set this as 1,2,3,4 etc.
`project`: The name of the project you are working on
`organism`: The organism / strain
`sample_type`: Put `rna` here.
3. Create a subsample sheet fill out the following.
`sample_name`: The name of the sample as supplied in the sample file above.
`subsample`: The same name of the sample_name. This is needed if it is paired end as there are technically two files to process.
`protocol`: Put `rna` here.
`seq_method`: Put `paired_end` here
4. Edit the config/pepconfig.yml to contain the full path to the subsample.csv and the sample.csv on lines `sample_table:` and `subsample_table:`
5. You can select to download your genome and gff3 via a url or to download it manually before hand and set the path. Make sure it is the genome and not the transcriptome or the gffread step will fail. Also make sure the gff3 is the annotation associated with the genome you have provided. You'll need to then edit the config/config.json file's lines:
`is_local`: to be True or False depending on if you wish for it to download from a provided url in the config or from the provided path.
`genome`: To be the path to where your genome.fa is.
`genome_name`: To be the exact name of your genome including the fasta.
`gff3`: To be the path where the gff3 file is location.
`run_gunc`: Put either "true" of flase, to choose to run GUNC on the reference genome being used.
`quantification_tool`: Put either "kallisto" or "star" to choose which quantification method you  with to use.

### Run The Workflow
You should be good to go. Resources are automatically downloaded for tools that need them in the `resouces/` folder. This may take an hour or so after running. At the ened you will get two combined matrices of counts and TPM from the quantification.
1. To run the workflow enter `snakemake --use-conda -j --cores 4` in the `rnaseq-quantification-pipeline` directory. Make sure you are *not* trying to run the quantification pipeline in the `test` or `workflow` directories. If the workflow is failing at a certain point, you can still generate some of the data *likely by running it with the `-k` flag as well.
NOTE:: You can also run the workflow on slurm by using `snakemake --slurm --default-resources slurm_account=osc_account_number -j --use-conda`.

TIP: If you want to run multiple samples, move them all to the same directory (or symlink them to save space), and design your input sheet to have different `project` ids. 


We are running the following programs, and are welcome to PR and input on other programs we should be running. 
Trimming: Trim-galore
Quanitifcation: Kallisto
Quality Control: 
    Fastq-screen
    Fastqc
    Multiqc

### This is an example DAG for an analysis run with 6 samples.

![Workflow DAG](/images/dag.svg)

Solid lines indicate the rules that have not been executed yet, whereas dashed lines depict completed jobs at time of the dag generation.

# Common FAQ
- `list index out of range` Error: If you recieve an error that says `list index out of range` and the traceback points to `get_r1_fastq` or `get_r2_fastq` as the reason. More often than not the problem is in how the sample name is being given in the sample and subsample.tsv files. A common error is that the R1 or R2 extension is not recognized or it is being duplicated. When inputing the sample names, you do not need to include R1 or R2 in the name. Only include the file name up until the last . or _ before the paired read name. Recognized file endings are `_R1`, `.R1`, `.r1`, `_r1`, `_1` with the ending `fastq.gz` or `fq.gz`. If the error persists please feel free to submit an issue.
- When making the input csv files, there shouldn't be any duplicates in the `sample` column of the sample.csv or else an error will be thrown by snakemake.
