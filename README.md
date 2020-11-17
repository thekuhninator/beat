# beat
BEAT - (Batch Effect Assessment Tool) is designed for researchers to assess the severity of the batch effect present in their RNA-Seq data. The tool is written in R and runs from the command line. It takes in as input gene counts, metadata, a dataset name, an output directory, and a flag denoting whether the dataset is uncorrected and outputs an html report with a pca plot, t-sne plot, comparative boxplot, as well as results form kBET. It also ouptuts a beat log file for use by multi_beat. After generating several beat reports and log files, a user can then run multi_beat on a parent directory to generate an aggregate report to compare the batch effect across the datasets using the metrics provided by BEAT. This can be used to help researchers choose which batch correction method works best for their particular dataset. 

## Installation 

The Batch Effect Assessment Tool (BEAT) package for R. Used to assess the magnitude of the batch effect present in RNA-Seq data.

To install this tool run the following commands `library(devtools) install_github('thekuhninator/beat')`

## Running beat

BEAT is a package within the R programming language. Below is an example of how to run beat:

<pre><code>input_counts <- "./dataset_1/dataset_unfiltered_gene_counts.csv"
input_annot  <- "./dataset_1/dataset_unfiltered_metadata.csv"
output_dir <- "./beat_output"
dataset_name <- "dataset_1"
original <- TRUE
beat::beat(input_counts, input_annot, output_dir, dataset_name, original)
</code></pre> 

The following is a detailed description of how the arguments used for `beat`.

### Arguments

<b> input_counts \<string\> </b>

The path to the gene counts file in comma seperated value format (.csv). The gene counts should be normalized prior to running beat. The gene counts should have the rows be gene counts and the columns be the name of the samples. See the dataset below for an example.

Example:

  genes     | sample1 | sample2 | sample3
-----|-------|--------|-------------
gene1 | val  | val     | val
gene2 | val  | val     | val
gene3 | val  | val     | val

<b> input_annot \<string\> </b>

The path to the metadata file in comma seperated value format (.csv). The metadata should have the rows be the samples and the columns be the features of the metadata. See the dataset below for an example. **It is crucial that the sample names be the first column in the metadata. Additionally, these names must match all the sample names in the gene counts file.** Otherwise, beat will not be able to run properly. **Additionally, there must be a columnn named "batch" (case sensitive) for beat to run properly.**

Example:

sample_id  | feature1 | feature2 | feature3
-----------|----------|----------|---------
sample1    |   val    | val     |   val
sample2    |   val    | val     |   val  
sample3    |   val    | val     |   val  

<b> output_dir \<string\> </b>

A directory to output the pca plot,  t-sne plot, k-bet plot, boxplot, beat log file, and html beat report.

<b> dataset_name \<string\> </b>

The name of the dataset to be used on the plots and name of the files.

<b> original </b>

A flag that should be used if the dataset is uncorrected. This will be used later when multi_beat is run. There should only be one original dataset when running multi_beat later on.

## beat Output Files

Beat produces 6 output files and outputs them in the output directory specified. The path it outputs to uses the output directory and dataset name in the form of: <output_dir>/\<dataset_name/.

#### beat Report: <output_dir>/<dataset_name>_batch_correction_report.html

The html file is the most valuable output from beat. It contains several plots that provide an overview as to how severe the batch effect is in the dataset. It contains a PCA plot, k-bet plot, t-sne plot, and comparative boxplot. Learn more about how each of the plots inform the user about the severity of the batch effect in their data below.

#### pca_plot: <output_dir>/<dataset_name>_pca_plot.png

Principal Component Analysis (PCA) is a dimensionality reduction technique that emphasizes the variation in the data and allows us to see patterns in the data. The X axis represents the first principal component and its contributor rate. The Y axis represents the second component and its contributor rate. Points represent each sample. Sample colors and shapes are according to a group the sample belongs to. If the plot shows many samples of the same color (same batch) clustering together, this means there is a strong batch effect presense in the data. If the plot shows colors well mixed the batch effect is not severe in the data.

#### kbet_plot: <output_dir>/<dataset_name>_kbet_plot.png

The K-Nearest Neighbor Batch Effect Test (kBET) is a test metric used for assessing the severity of a batch effect in the data. The algorithm creates k-nearest neighbour matrix and choses 10% of the samples to check the batch label distribution in its neighbourhood. If the local batch label distribution is sufficiently similar to the global batch label distribution, the $\chi^2$-test does not reject the null hypothesis (that is "all batches are well-mixed"). The neighbourhood size k is fixed for all tests. Next, the test returns a binary result for each of the tested samples. Finally, the result of kBET is the average test rejection rate. The lower the test result, the less bias is introduced by the batch effect. kBET is very sensitive to any kind of bias. By comparing the distribution of the expected and observed plots, one can see how severe the batch effect is in the data. The closer the two boxes are, the less severe, the batch effect. For more information about kbet, see their [their paper](https://www.nature.com/articles/s41592-018-0254-1) or [github](https://github.com/theislab/kBET/blob/master/README.md).

#### t-sne plot: <output_dir>/<dataset_name>_t_sne_plot.png

T-distributed Stochastic Neighbor Embedding (t-sne) is a machine learning algorithm for visualization. It is also a dimensionality reduction technique like PCA and is also useful in determing the severity of the batch effect by examining how strongly the colors (batches) are clustering together.

#### comparative boxplot: <output_dir>/<dataset_name>_boxplot.png

The comparative boxplot is a useful way of visualizing how the batches vary in the distribution of each gene's mean expression. Each gene's mean expression value across all samples within a batch are used as data points in constructing the comparative boxplot. If the boxes appear to be similar in their distribution the batch effect is not as severe for the dataset.

#### beat log file: <output_dir>/<dataset_name>_beat_log.beat

The beat log files contain important information that are used when mutli_beat is run. They are stored in the R data format.
___

# multi_beat

multi_beat is a function used to combine several beat reports. multi_beat takes in as input a parent directory, output_directory, and output_name. It generates as output a report of the combined pca plots, combined t-sne plots, combined comparative boxplot, as well as a scatterplot of highly variable genes retained vs kBET acceptance rate. This is best used when all the beat files used for the report belong to the same dataset, as the legend will be shared between them. This is useful in determining which of several batch correction methods works best on a dataset. 

## Running multi_beat

multi_beat comes with the beat R package. Below is an example on how to run multi_beat:

<pre><code>parent_dir <- "./"
output_dir <- "multi_beat_output"
output_name <- "dataset_1"
beat::multi_beat(parent_dir, output_dir, output_name)
</code></pre> 

Run `multi_beat` from the command line as follows.
 
The following is a detailed description of how the arguments used for `multi_beat`.

### Arguments

<b> parent_dir \<string\> </b> (required)

The path to the parent directory holding other folders containing .beat log files. multi_beat checks all subfolders recursively.

<b> output_dir \<string\> </b> (required)

The path to the output file where multi_beat will output the aggregate report.

<b> output_name \<string\> </b> (required)

The name of the dataset being examined, used for titles and file names.


### multi_beat Output Files

multi_beat produces 5 output files. The files and their names are listed below.

#### beat Report: <output_dir>/<output_name>_batch_correction_report.html

The html file is the most valuable output from multi_beat. It contains several combined plots that provide an overview as to how severe the batch effect is in each of the correction methods for the dataset. It contains a combined PCA plot, k-bet vs highly variable genes plot, combined t-sne plot, and combined comparative boxplot. Learn more about how each of the plots inform the user about the severity of the batch effect in their data below.

#### kBET vs Highly Variable Genes Plot: <output_dir>/<output_name>_kbet_hvg_scatterplot.png

**The k-BET vs Highly Variable Genes plot is perhaps the most useful plot from multi_beat**. The highly variable genes are defined as the top 10% of genes with the highest variance. multi_beat checks which genes are retained between each dataset and the original uncorrected dataset which was specified by the user when beat was first run. The percentage of highly variable genes retained after correction serves as a metric for biological preservation. k-BET serves as a metric for the severity of the batch effect. By plotting these two in a scatterplot, one can ascertain which correction method best suits their dataset.

#### Combined PCA Plot: <output_dir>/<output_name>_pca_plot.png

Principal Component Analysis (PCA) is a dimensionality reduction technique that emphasizes the variation in the data and allows us to see patterns in the data. The X axis represents the first principal component and its contributor rate. The Y axis represents the second component and its contributor rate. Points represent each sample. Sample colors and shapes are according to a group the sample belongs to. If the plot shows many samples of the same color (same batch) clustering together, this means there is a strong batch effect presense in the data. If the plot shows colors well mixed the batch effect is not severe in the data.

#### Combined t-sne plot: <output_dir>/<output_name>_t_sne_plot.png

T-distributed Stochastic Neighbor Embedding (t-sne) is a machine learning algorithm for visualization. It is also a dimensionality reduction technique like PCA and is also useful in determing the severity of the batch effect by examining how strongly the colors (batches) are clustering together.

#### Combined comparative boxplot: <output_dir>/<output_name>_boxplot.png

The combined comparative boxplot is a useful way of visualizing how the batches vary in the distribution of each gene's mean expression. Each gene's mean expression value across all samples within a batch are used as data points in constructing the comparative boxplot. If the boxes appear to be similar in their distribution the batch effect is not as severe for the dataset.

## Contact the Author

If you are having issues with `beat`, feel free to reach out to me, at roman (dot) kuhn1 (at) gmail (dot) com.

## Cite beat 

If you find `beat` useful in your research please cite the related publication:

[Assessment of Batch Correction Methods for Whole Blood RNA-Seq Data](http://google.com)
<h1> BEAT Workflow Pipeline/Use Diagram </h1>

![BEAT Pipeline](https://github.com/thekuhninator/beat/blob/master/beat_pipeline.png?raw=true)
