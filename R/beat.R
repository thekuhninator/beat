#' Run Beat
#'
#' Takes in gene counts, metadata, and output paths, as well as a dataset name and outputs an html report.
#' @param gene_counts the path to the gene counts with rows as genes and columns as samples
#' @export
beat <- function(input_counts, input_metadata, output_dir, dataset_name, original) {

  factor_of_interest = NULL
  #
  # Let's load in the data
  #
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


  #
  # get all the results
  #
  print('Running kbet...')
  kbet_results <- kbet(gene_counts, annot, output_dir, dataset_name)
  print('Creating t-sne plot...')
  tsne_results <- tsne_batch(gene_counts, annot, output_dir, dataset_name)
  print('Creating pca plot...')
  pca_results <- pca_m3c(t(gene_counts), annot, output_dir, dataset_name, factor_of_interest)
  print('Creating boxplot...')
  boxplot_results <- grouped_boxplot(gene_counts, annot, output_dir, dataset_name);
  print('Getting HVGs...')
  hvgs <- getHvgs(gene_counts, annot)

  tsne_path <- tsne_results$path
  tsne_plot <- tsne_results$plot
  pca_path  <- pca_results$path
  pca_plot  <- pca_results$plot
  kbet_path <- kbet_results$path
  kbet_data <- kbet_results$results
  boxplot_path <- boxplot_results$path
  boxplot_data <- boxplot_results$data

  kbet_base64     <- toBase64(kbet_path)
  pca_base64      <- toBase64(pca_path)
  tsne_base64     <- toBase64(tsne_path)
  boxplot_base64  <- toBase64(boxplot_path)

  # get all the base64 data of images
  #tsne_base64 <- toBase64(tsne_path)
  #pca_base64  <- toBase64(pca_path)
  #box_base64  <- toBase64(boxplot_path)

  # write them to an html report
  generate_report(dataset_name, output_dir, kbet_base64, pca_base64, tsne_base64, boxplot_base64)
  # write them to the log file
  generateLogFile(dataset_name, output_dir, original, kbet_data, tsne_plot, pca_plot, boxplot_results, hvgs)

  print('The report and log file have succesfully been generated!')

}


#' Run MultiBeat
#'
#' Takes in gene counts, metadata, and output paths, as well as a dataset name and outputs an html report.
#' @param parent_dir the path to the gene counts with rows as genes and columns as samples
#' @param output_dir the path to the directory in which the output shall be.
#' @param output_name the name of the report that will be generated.
#' @export
multi_beat <- function(parent_dir, output_dir, output_name) {
  
  # create the direcotry if it does not already exist
  dir.create(file.path(output_dir), showWarnings = FALSE)
  
  print(paste('Running multi_beat for ', output_name, sep=""))
  
  # check if the parent_dir exists...
  if (!dir.exists(parent_dir))
    stop("Cannot find parent directory. Are you sure the path is correct?")
  
  # load in necessary data from folders
  beat_files <- list.files(path=parent_dir, recursive=TRUE, pattern="\\.beat$")
  
  original_dataset_name <- NULL
  hvgs_hash        <- hash()
  hvgs_retention   <- hash()
  kbet_accept_hash <- hash()
  pca_plots        <- hash()
  tsne_plots       <- hash()
  
  retained_hvgs <- NULL
  comparative_boxplot_data <- NULL
  datasets <- c()
  
  #print(beat_files)
  #beat_file
  for (beat_file in beat_files)
  {
    load(file=file.path(parent_dir, beat_file))
    datasets <- c(datasets, dataset_name)
    
    if(original)
      original_dataset_name <- dataset_name
    
    # save the hvgs
    hvgs_hash[[dataset_name]] <- hvgs
    kbet_accept_hash[[dataset_name]] <- kbet_results['kBET.observed'][1,]
    boxplot_data['dataset'] <- dataset_name
    pca_plots[[dataset_name]] <- pca_plot
    tsne_plots[[dataset_name]] <- tsne_plot
    
    if(is.null(comparative_boxplot_data))
    {
      comparative_boxplot_data <- boxplot_data
    }
    else
    {
      comparative_boxplot_data <- rbind(comparative_boxplot_data, boxplot_data)
    }
    
  }
  # create the output
  kbet_plot_data <- get_kbet_plot_data(hvgs_hash, original_dataset_name, kbet_accept_hash)  
  print(kbet_plot_data)
  kbet_hvg_path  <- kbet_hvg_scatterplot(kbet_plot_data, output_dir, output_name)
  boxplot_path   <- grouped_boxplot(comparative_boxplot_data, output_dir, output_name)
  pca_tile_path  <- tile_plots(pca_plots, datasets, 'pca', output_dir, output_name, 'PCA Combined Plots')
  tsne_tile_path <- tile_plots(tsne_plots,  datasets,'tsne',   output_dir, output_name, 'T-SNE Combined Plots')
  
  kbet_hvg_base64  <- toBase64(kbet_hvg_path)
  boxplot_base64   <- toBase64(boxplot_path)
  pca_tile_base64  <- toBase64(pca_tile_path)
  tsne_tile_base64 <- toBase64(tsne_tile_path)
  
  generate_aggregate_report(kbet_hvg_base64, boxplot_base64, pca_tile_base64, tsne_tile_base64, output_dir, output_name)
  
  
  print('The multibeat report and log file have succesfully been generated!')
  
}

toBase64 <- function(image_file) {
  uri=image_uri(image_file)
}
