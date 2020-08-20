# kBET function
# Input: Takes in gene conts data table, metadata data table, a foi, output_directory, and output name
# Ouptut: returns a dataframe containing data regarding the kBET test

kbet <- function(gene_counts, annot, output_dir, output_name)
{
  print('Made it to the beggninning of this function')
  # get the factor of interest (not currently used)
  #foi = annot[factor_of_interest]

  # additionall quantitative metrics
  #avedistVal = avedist(new_data, as.factor(batch) )#, as.factor(diagnosis))
  #pvcamVal = pvcam(new_data, batch, foi)
  #skewdivVal = skewdiv(new_data, as.factor(batch))
  #kldistVal = kldist(new_data, as.factor(batch))
  #sepscoreVal = sepscore(new_data, as.factor(batch) )
  # diffexprmVal =
  # cobraVal =

  data = gene_counts

  batch = annot$batch
  print('did we make it pas this far?')
  # Capitalize Function

  capitalize <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }

  # Generate the Plot title
  # plotTitle = paste(capitalize(unlist(strsplit(name, "_"))), collapse= " ")
  plotTitle = paste(output_name, " kBET Plot", sep="")

  file_path <- paste( file.path(output_dir, output_name), "_kbet_plot.png", sep="")
  png(file=file_path)
  print('running kbet')
  print(dim(data))
  print(dim(batch))
  batch.estimate <- kBET::kBET(data, batch, plot=FALSE)
  print('finmihsed running kbet')
  plot.data <- data.frame(class=rep(c('observed', 'expected'),
                                   each=length(batch.estimate$stats$kBET.observed)),
                          data =  c(batch.estimate$stats$kBET.observed,
                                    batch.estimate$stats$kBET.expected))
  print('finished plotting the data')
  print(ggplot2::ggplot(plot.data, aes(class, data)) + geom_boxplot() +
          labs(x='Test', y='Rejection rate',title=paste(plotTitle,'kBET test results',sep=' ')) +
          theme_bw() +
          scale_y_continuous(limits=c(0,1)))

  dev.off()

  results <- batch.estimate$summary

  # add extra quantitative measures to results
  #results['avedist'] = avedistVal
  #results['pvcam'] = pvcamVal
  #results['skewdiv'] = skewdivVal
  #results['kldist'] = kldistVal

  return( list("path" = file_path, "results" = results))

}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

# Capitalize Function
capitalize <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

tsne_batch <- function(gene_counts, annot, output_dir, output_name)
{
  # output PCA
  file_path <- paste( file.path(output_dir, output_name), "_tsne_plot.png", sep="")
  png(file=file_path)

  #name = paste(sapply(paste(unlist(strsplit(output_name, " ")), sep=" "), simpleCap), collapse=" ")
  #plotTitle = paste(capitalize(unlist(strsplit(name, "_"))), collapse= " ")
  plotTitle <- paste( output_name, "T-SNE", sep=" ")
  invisible(capture.output(
    g <- M3C::tsne(t(gene_counts), labels=as.factor(annot$batch), legendtitle ="Batch", dotsize = 2) +
      ggtitle(plotTitle) +
      theme(plot.title = element_text(size=20, hjust = .5))
  ))
  #print(typeof(g))
  print(g)
  dev.off()

  return (list("path" = file_path, "plot" = g))
}

pca_m3c <- function(gene_counts, annot, output_dir, output_name, foi=NULL)
{

  file_name <- paste(file.path(output_dir, output_name), '.png', sep="")
  #png(file_name)

  # output PCA
  file_path <- paste( file.path(output_dir, output_name), "_pca_picture.png", sep="")
  png(file=file_path)

  name = paste(sapply(paste(unlist(strsplit(output_name, "[_]")), sep=" "), simpleCap), collapse=" ")
  plotTitle = paste(capitalize(unlist(strsplit(name, "_"))), collapse= " ")
  plotTitle <- paste( plotTitle, " PCA", sep="")

  g <- M3C::pca(gene_counts, labels=as.factor(annot$batch), legendtitle ="Batch", dotsize = 2) +
    ggtitle(plotTitle) +
    theme(plot.title = element_text(size=20, hjust = .5)) +
    scale_size_continuous(range = c(1, 2))
  print(g)
  dev.off()

  return (list("path" = file_path, "plot" = g))


}

grouped_boxplot <- function(gene_counts, annot, output_dir, dataset_name)
{
  # get list of samples for all batchs
  # for each batch variable
  # get list of samples for this batch
  # pull out their gene expression values
  # get the means and put it in a list
  # add it to a list
  # add batch x to names

  boxplot_data = NULL


  for (x in unique(annot$batch))
  {
    batchX = annot[annot$batch == x,]
    #print(batchX)

    rownames(batchX)

    gene_counts_subset <- t(gene_counts[rownames(batchX), ])

    #print(dim(gene_counts_subset))

    means = rowMeans(gene_counts_subset)

    mean_values <- as.vector(data.frame(means = rowMeans(gene_counts_subset))$means)

    batch_name = paste(c("batch_", x), collapse="", sep="")

    if(is.null(boxplot_data))
    {
      #print('first time')
      boxplot_data = data.frame("mean" = mean_values, "batch" = rep(x, each=length(mean_values)))
    }
    else
    {
      #print('not first time')
      new_data = data.frame("mean" = mean_values, "batch" = rep(x, each=length(mean_values)))
      #print(dim(new_data))

      boxplot_data = rbind(boxplot_data, new_data)
    }
    #print(dim(boxplot_data))

  }


  #print(boxplot)data
  file_name <- paste(dataset_name, '_comparative_boxplot.png', sep ="")
  output_file_name <- file.path(output_dir, file_name)
  png(output_file_name)


  g <- ggplot2::ggplot(boxplot_data, aes(x = factor(batch), y = mean, fill=factor(batch))) +
    geom_boxplot() +
    labs(title = "Comparative Boxplot", x = "Batch", y = "Mean Gene Expression", fill = "Batch") +
    theme(plot.title   = element_text(size=19, hjust = .5),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          legend.title = element_text(size=14),
          axis.text.x =  element_text(size=16),
          axis.text.y =  element_text(size=16))


  print(g)

  dev.off()

  return (list( "path" = output_file_name, "data" = boxplot_data))
}


getHvgs <- function(gene_counts, annot)
{
  hvgs <- NULL

  # let's iterate through all the batches
  for (batch_i in unique(annot$batch))
  {
    batch_counts <- t(gene_counts[rownames(annot[annot$batch == batch_i,]),]);
    batch_counts <- as.data.frame(cbind(batch_counts, 'var'=matrixStats::rowVars(batch_counts)))
    top_percent <- .10
    n_genes <- dim(batch_counts)[1] * top_percent
    top_varying <- rownames(head(batch_counts[order(batch_counts$var,decreasing=T),],n_genes))

    if(is.null(hvgs))
    {
      hvgs <- top_varying
    }
    else
    {
      hvgs <- Reduce(intersect, list(hvgs, top_varying))
    }

  }
  return (hvgs)
}
