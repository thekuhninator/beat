
#
# Get the kbet plot data
#
get_kbet_plot_data <- function(hvgs_hash, original_dataset_name, kbet_accept_hash)
{
  original_hvg_list <- hvgs_hash[[original_dataset_name]]
  #print('length of orignal hvgs list')
  num_original_hvgs <- length(original_hvg_list)
  #print(num_original_hvgs)
  
  # get the hvgs percent retained
  for (dataset_name in datasets)
  {
    #print('dataset name')
    #print(dataset_name)
    hvgs_retained <- Reduce(intersect, list(orig = original_hvg_list, other = hvgs_hash[[dataset_name]] ))
    #print('length of hvgs retained')
    #print(length(hvgs_retained))
    percent_retained <- length(hvgs_retained) / num_original_hvgs
    hvgs_retention[[dataset_name]] <- percent_retained
    #print('percent retained')
    #print(percent_retained)
  }
  
  hvgs_retention[["dataset1_combat"]]
  
  kbet_plot_data <- data.frame("kbet_acceptance" = numeric(), "dataset_name" = character(), "hvgs_retained" = numeric(), stringsAsFactors = FALSE)
  kbet_plot_data
  
  # get the kbet value for each dataset
  for (dataset_name in datasets)
  {
    kbet_val <- 1 - kbet_accept_hash[[dataset_name]]
    retained <- hvgs_retention[[dataset_name]]
    new_data <- list(kbet_val, dataset_name, retained)
    #print('new data')
    #print(new_data)
    #print(dataset_name)
    kbet_plot_data[nrow(kbet_plot_data) + 1,] = new_data
    #kbet_plot_data[nrow(kbet_plot_data) + 1, 1] = kbet_val
    #kbet_plot_data[nrow(kbet_plot_data) , 2] = dataset_name
    #kbet_plot_data[nrow(kbet_plot_data), 3] = retained
  }
  print(kbet_plot_data)
  return(kbet_plot_data)
}


kbet_hvg_scatterplot <- function(kbet_plot_data, output_dir, output_name)
{
  # create the output file
  file_name <- paste(output_name, '_kbet_hvg_scatterplot.png')
  output_file_path <- file.path(output_dir, file_name)
  png(output_file_path)
  
  print(kbet_plot_data)
  g <- ggplot(kbet_plot_data, aes(x=as.numeric(hvgs_retained), y = as.numeric(kbet_acceptance), color = dataset_name)) + 
    geom_point(size=4) + 
    labs(title = "KBET vs HVGs Plot", x = "Percent of HVGS Retained",
         y = "kBET Acceptance Rate", color="Dataset") + 
    scale_x_continuous(limits=c(0,1)) + 
    scale_y_continuous(limits=c(0,1))
  
  print(g)
  
  dev.off()
  return(output_file_path)
  
}

#
# Genereate grouped boxplot
#
grouped_boxplot <- function(boxplot_data, output_dir, output_name)
{
  # create the output file
  file_name <- paste(output_name, '_comparative_boxplot.png')
  output_file_path <- file.path(output_dir, file_name)
  png(output_file_path)
  
  # generate the plot
  g <- ggplot(comparative_boxplot_data, aes(x=dataset,
                                            y=mean,
                                            fill=as.factor(batch))) +
    geom_boxplot() +
    labs(title = "Comparative Grouped Boxplot", x ="Dataset", y = "Gene Mean Expression", fill = "Batch") +
    theme(axis.text.x = element_text(angle=30))
  
  
  # save the output file
  print(g)
  dev.off()
  
  return(output_file_path)
}


#
# Tile plots
#
tile_plots <- function(plots, dataset_names, plot_name, output_dir, output_name, combined_title)
{
  # create the output file
  file_name <- paste(output_name, plot_name, 'plot.png', sep = "_")
  output_file_path <- file.path(output_dir, file_name)
  png(output_file_path)
  plot_list <- list()
  for(name in dataset_names)
  {
    plot <- plots[[name]]
    plot <- plot + theme_bw()
    plot_list[[name]] <-plot
  }
  
  p_no_legend <- lapply(plot_list, function(x) x + theme(legend.position = "none"))
  legend <- cowplot::get_legend(plot_list[[1]] + theme(legend.position = "bottom"))
  title <- cowplot::ggdraw() + cowplot::draw_label(combined_title, fontface = "bold")
  p_grid <- cowplot::plot_grid(plotlist = p_no_legend, ncol = 2)
  g <- cowplot::plot_grid(title, p_grid, legend, ncol = 1, rel_heights = c(0.1, 1, 0.2))
  
  print(g)
  dev.off()
  return (output_file_path)
}

generate_aggregate_report <- function(kbet_hvg_src, boxplot_src, pca_tile_src, tsne_tile_src, output_dir, output_name)
{
  html_string = 
    paste('
          <html>
          <head>
          <link rel="stylesheet" href="https://cdn.jsdelivr.net/gh/kognise/water.css@latest/dist/light.min.css">
          
          <style>body{ margin:0 100; background:whitesmoke; }</style>
          </head>
          <body>
          <h1>Batch Correction Aggregate Report for ', output_name, '</h1>
          
          
          <!-- *** Section 3 *** --->
          <h2>Principal Component Analysis</h2>
          <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
          src="', pca_tile_src, '" ></iframe> <p>Principal Component Analysis (PCA) is a dimensionality reduction technique that emphasizes the variation in the data and allows us to see patterns in the data. The X axis represents the first principal component and its contributor rate. The Y axis represents the second component and its contributor rate. Points represent each sample. Sample colors and shapes are according to a group the sample belongs to. If the plot shows many samples of the same color (same batch) clustering together, this means there is a strong batch effect presense in the data. If the plot shows colors well mixed the batch effect is not severe in the data.</p>
          
          <h2>kBET - K-Nearest Neighbour Batch Effect Test Acceptance Rate vs HVGs Retention Rate</h2>
          <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
          src="', kbet_hvg_src, '"></iframe> <p>The highly variable genes are defined as the top 10% of genes with the highest variance. multi_beat checks which genes are retained between each dataset and the original uncorrected dataset which was specified by the user when beat was first run. The percentage of highly variable genes retained after correction serves as a metric for biological preservation. k-BET serves as a metric for the severity of the batch effect. By plotting these two in a scatterplot, one can ascertain which correction method best suits their dataset.</p>
          
          <!-- *** Section 1 *** --->
          <h2>T-Stochastic Neighbor Embedding (T-SNE) Tiled Plots</h2>
          <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
          src="' , tsne_tile_src, '"></iframe>
          <p>T-distributed Stochastic Neighbor Embedding (t-sne) is a machine learning algorithm for visualization. It is also a dimensionality reduction technique like PCA and is also useful in determing the severity of the batch effect by examining how strongly the colors (batches) are clustering together.</p>
          
          <!-- *** Section 2 *** --->
          <h2>Comparative BoxPlot</h2>
          <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
          src="',  boxplot_src, '"></iframe>
          <p>The combined comparative boxplot is a useful way of visualizing how the batches vary in the distribution of each gene\'s mean expression. Each gene\'s mean expression value across all samples within a batch are used as data points in constructing the comparative boxplot. If the boxes appear to be similar in their distribution the batch effect is not as severe for the dataset.</p>
          
          
          
          </body>
          </html>', sep="")
  
  file_name <- paste(output_name, '_batch_correction_report.html', sep="")
  output_file_path <- file.path(output_dir, file_name)
  html_report <- file(output_file_path)
  writeLines(c(html_string), html_report)
  close(html_report)
}
