load_data <- function(input_counts, input_metadata, dataset_name, output_dir, factor_of_interest)
{
  # create the direcotry if it does not already exist
  dir.create(file.path(output_dir), showWarnings = FALSE)

  # check if file paths exist
  if (!file.exists(input_counts))
    stop("Cannot find gene counts file. Are you sure the path is correct?")
  if (!file.exists(input_metadata))
    stop("Cannot find metadata file. Are you sure the path is correct?")

  print(paste('Running beat for ', dataset_name, sep=""))

  # Read in data file
  gene_counts <- t(read.csv(input_counts, header = TRUE, row.names = 1, check.names=FALSE))
  #Read in annotation file
  annot <- read.csv(input_metadata, header = TRUE, row.names = 1, check.names = FALSE)

  # check if factor of interest was inputted and exists in the metadataw
  if(!is.null(factor_of_interest) && !factor_of_interest %in% names(annot))
    stop("Factor of interest not in the columns of the metadata. Are you sure it exists and is spelled correctly?")
}

generate_report <- function(dataset_name, output_dir, kbet_src, pca_src, tsne_src, boxplot_src)
{

  #dataset_name = "combined_dataset1"
  #pca_path = 'C:/Users/Roman/Documents/Work/Depression_and_Immunology/beast/dataset_1_stuff_picture.png'
  #pca_path = 'C:/Users/Roman/Documents/Work/Depression_and_Immunology/beast/dataset_1_stuff_picture.png'
  #pca_path = 'C:/Users/Roman/Documents/Work/Depression_and_Immunology/beast/dataset_1_stuff_picture.png'
  html_string =
    paste('
          <html>
          <head>
          <link rel="stylesheet" href="https://cdn.jsdelivr.net/gh/kognise/water.css@latest/dist/light.min.css">

          <style>body{ margin:0 100; background:whitesmoke; }</style>
          </head>
          <body>
          <h1>Batch Correction Report for ', dataset_name, '</h1>


          <!-- *** Section 3 *** --->
          <h2>Principal Component Analysis</h2>
          <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
          src="', pca_src, '" ></iframe> <p>Principal Component Analysis (PCA) is a dimensionality reduction technique that emphasizes the variation in the data and allows us to see patterns in the data. The X axis represents the first principal component and its contributor rate. The Y axis represents the second component and its contributor rate. Points represent each sample. Sample colors and shapes are according to a group the sample belongs to. If the plot shows many samples of the same color (same batch) clustering together, this means there is a strong batch effect presense in the data. If the plot shows colors well mixed the batch effect is not severe in the data.

          </p>

          <h2>kBET - K-Nearest Neighbour Batch Effect test</h2>
          <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
          src="', kbet_src, '"></iframe> <p>The K-Nearest Neighbor Batch Effect Test (kBET) is a test metric used for assessing the severity of a batch effect in the data. The algorithm creates k-nearest neighbour matrix and choses 10% of the samples to check the batch label distribution in its neighbourhood. If the local batch label distribution is sufficiently similar to the global batch label distribution, the chi-squared test does not reject the null hypothesis (that is "all batches are well-mixed"). The neighbourhood size k is fixed for all tests. Next, the test returns a binary result for each of the tested samples. Finally, the result of kBET is the average test rejection rate. The lower the test result, the less bias is introduced by the batch effect. kBET is very sensitive to any kind of bias. By comparing the distribution of the expected and observed plots, one can see how severe the batch effect is in the data. The closer the two boxes are, the less severe, the batch effect. For more information about kbet, see their their paper or github.

          </p>

          <!-- *** Section 1 *** --->
          <h2>T-Stochastic Neighbor Embedding (T-SNE)</h2>
          <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
          src="' , tsne_src, '"></iframe>
          <p>T-distributed Stochastic Neighbor Embedding (t-sne) is a machine learning algorithm for visualization. It is also a dimensionality reduction technique like PCA and is also useful in determing the severity of the batch effect by examining how strongly the colors (batches) are clustering together.</p>

          <!-- *** Section 2 *** --->
          <h2>Comparative BoxPlot</h2>
          <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
          src="',  boxplot_src, '"></iframe>
          <p>The comparative boxplot is a useful way of visualizing how the batches vary in the distribution of each gene mean expression value. The mene gene expression value across all samples within a batch are used as data points in constructing the comparative boxplot. If the boxes appear to be similar in their distribution the batch effect is not as severe for the dataset.

          .</p>



          </body>
          </html>', sep="")

  file_name <- paste(dataset_name, '_batch_correction_report.html', sep="")
  output_file_path <- file.path(output_dir, file_name)
  html_report <- file(output_file_path)
  writeLines(c(html_string), html_report)
  close(html_report)
}

generateLogFile <- function(dataset_name, output_dir, original, kbet_results, tsne_plot, pca_plot, boxplot_results, hvgs)
{

  file_name <- paste(dataset_name, '_beat_log.beat',sep="")
  output_file_path <- file.path(output_dir, file_name)
  save(dataset_name, original, kbet_results, tsne_plot, boxplot_results, pca_plot, hvgs, file=output_file_path)
}

toBase64 <- function(image_file) {
  uri=knitr::image_uri(image_file)
}
