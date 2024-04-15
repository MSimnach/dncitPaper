#' Title: Plotting functions for rejection rates of CITs
#'
#' @param folder_path path to folder with rejection rates of all settings and all CITs
#' @param dncits vector of all DNCITs
#' @param settings data frame with all settings, columns should consist of settings$confounder, fct_relation, embedding_map, noise_x
#' @param sample_sizes vector of all sample sizes (either sample_sizes or betas has to be specified)
#' @param betas vector of all betas (either betas or sample_sizes has to be specified)
#'
#' @return returns a plot with all rejection rates of all settings and all CITs
#' @export
create_test_plot <- function(folder_path, dncits, settings, sample_sizes=NULL, betas=NULL){
  # List all files in the folder
  all_files <- list.files(folder_path, full.names = TRUE)

  ################# Extract unique parts of file names -> MANUALLY ADJUSTED ##############################
  settings_files <- unique(substr(basename(all_files), 1, 39))
  # exclude additional noise for Y
  settings_files <- settings_files[!grepl('_50_0_fas', settings_files)]
  # exclude noise 10 on embeddings, since tests not sensitive for noise on features
  settings_files <- settings_files[!grepl('10_1_0_fas', settings_files)]
  #Reorder matching with tree diagram in y axis
  settings_files <- settings_files[c(23,21,17,19,15,13,9,11,7,5,1,3,24,22,18,20,16,14,10,12,8,6,2,4)]

  # empty result tab for each DNCIT
  list_result_tabs <- list()
  for(dncit in dncits){
    dncit_char <- paste0(dncit)
    list_result_tabs[[dncit_char]] <- create_empty_result_tab(settings_files,sample_sizes=sample_sizes, betas=betas)
  }

  list_result_tabs_full <- create_reject_tabs(dncits, settings_files, all_files, list_result_tabs)

  results_plots <- list()
  for(dncit in dncits){
    if(which(dncits == dncit) == 1){
      results_plots[[dncit]] <- plot_left_dncit_col(dncit, list_result_tabs_full[[dncit]])
    }else if(which(dncits == dncit) == length(dncits)){
      results_plots[[dncit]] <- plot_right_dncit_col(dncit, list_result_tabs_full[[dncit]])
    }else{
      results_plots[[dncit]] <- plot_mid_dncit_col(dncit, list_result_tabs_full[[dncit]])
    }
  }

  y_axis <- create_y_axis(settings)

  # Combine plots
  p1 <- cowplot::plot_grid(y_axis, results_plots,
                           ncol = length(dncits)+1,  rel_widths = c(2, rep(0.5,length(dncits)-1), 1), axis= "bt")

  #y.grob <- grid::textGrob("Setting",
  #                   gp=ggplot2::gpar(fontface="bold", col="black", fontsize=15), rot=90)
  if(is.null(betas)){
    x.grob <- grid::textGrob(expression(paste("Sample size n=30,100,300,1000,3000,10000")),
                             gp=grid::gpar(fontface="bold", col="black", fontsize=15))
  }else if(is.null(sample_sizes)){
    x.grob <- grid::textGrob(expression(paste("Effect size ", beta, " in ", 10 %*% 10^{-k}, " for k=0,1,2,3")),
                       gp=grid::gpar(fontface="bold", col="black", fontsize=15))
  }

  x.grob <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p1,# left = y.grob,
                                      bottom = x.grob))
  return(x.grob)
}

create_empty_result_tab <- function(settings_files,sample_sizes=NULL, betas=NULL){
  if(is.null(betas)){
    result_tab <- data.frame(matrix(nrow=length(settings_files), ncol = length(sample_sizes)+1))
    colnames(result_tab) <- c('setting', sample_sizes)
  }else if(is.null(sample_sizes)){
    result_tab <- data.frame(matrix(nrow=length(settings_files), ncol = length(betas)))
    colnames(result_tab) <- c('setting',betas)
  }
  return(result_tab)
}

create_reject_tabs <- function(dncits, settings_files, all_files, list_result_tabs){
  count2 <- 0

  for (part in settings_files) {
    count2 = count2+1
    matching_files <- grep(part, all_files, value = TRUE)

    if (length(matching_files) > 1) {
      count=1
      for (file_path in matching_files) {

        data <- utils::read.csv(file_path)
        if(count==1){
          merged_part <- data.frame('X' = data$X)
        }
        count=count+1
        col_name <- gsub(part, "", basename(file_path))

        # Modify col_name to always start with first letter after '_'
        first_underscore <- regexpr("_", col_name)
        if (first_underscore != -1) {
          if(startsWith(col_name, 'kpc_graph')){

          }else{
            col_name <- substring(col_name, first_underscore + 1)
          }
        }

        if (ncol(data) > 1) {
          colnames(data)[-1] <- col_name
        }

        merged_part <- dplyr::full_join(merged_part, data, by = "X")
      }

      #extract rejection rates for each DNCIT
      for(dncit in dncits){
        dncit_csv <- paste0(dncit, '.csv')
        if (dncit_csv %in% colnames(merged_part)){
          print(dncit)
          list_result_tabs[[dncit]][count2,] <- c(part, merged_part[[dncit_csv]])
        }
      }
    }
  }
  return(list_result_tabs)
}

create_long_tab <- function(dncit, result_tab_dncit){
  dncit_long <- reshape2::melt(result_tab_dncit, id.vars = "setting")
  dncit_long$setting <- factor(dncit_long$setting, levels=rev(unique(dncit_long$setting)))
  colnames(dncit_long)[3] <- 'rejection_rate'
  dncit_long$rejection_rate <- as.numeric(dncit_long$rejection_rate)
  return(dncit_long)
}

plot_left_dncit_col <- function (dncit, result_tab_dncit){
  dncit_tab_long <- create_long_tab(dncit, result_tab_dncit)
  plot_dncit <- ggplot2::ggplot(dncit_tab_long, ggplot2::aes(x = variable, y = setting, fill = rejection_rate)) +
    ggplot2::geom_tile()+
    ggplot2::scale_fill_gradientn(colors = c('green', 'darkgreen', 'yellow1','orange', 'red'), values=c(0,0.05,0.1,0.15,1), limits=c(0,1)) +  # Adjust the color gradient as needed
    ggplot2::labs(x = "Sample size", y = "", title = dncit) +  # No need to repeat y-axis label
    ggplot2::theme(axis.text= ggplot2::element_blank(),  # Remove y-axis text
          axis.title= ggplot2::element_blank())  +
    ggplot2::theme(legend.position = "none")
  return(plot_dncit)
}

plot_mid_dncit_col <- function(dncit, result_tab_dncit){
  dncit_tab_long <- create_long_tab(dncit, result_tab_dncit)
  plot_dncit <- ggplot2::ggplot(dncit_tab_long, ggplot2::aes(x = variable, y = setting, fill = rejection_rate)) +
    ggplot2::geom_tile()+
    ggplot2::scale_fill_gradientn(colors = c('green', 'darkgreen', 'yellow1','orange', 'red'), values=c(0,0.05,0.1,0.15,1), limits=c(0,1)) +  # Adjust the color gradient as needed
    ggplot2::labs(x = "Sample size", y = "", title = dncit) +  # No need to repeat y-axis label
    ggplot2::theme(axis.text= ggplot2::element_blank(),  # Remove y-axis text
          axis.title= ggplot2::element_blank()) +
    ggplot2::theme(legend.position = "none")
  return(plot_dncit)
}

plot_right_dncit_col <- function(dncit, result_tab_dncit){
  dncit_tab_long <- create_long_tab(dncit, result_tab_dncit)
  plot_dncit <- ggplot2::ggplot(dncit_tab_long, ggplot2::aes(x = variable, y = setting, fill = rejection_rate)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = c('green', 'darkgreen', 'yellow1','orange', 'red'), values=c(0,0.05,0.1,0.15,1), limits=c(0,1)) +  # Adjust the color gradient as needed
    ggplot2::labs(x = "Sample size", y = "Setting", title = dncit) +
    ggplot2::theme(axis.text= ggplot2::element_blank(),  # Remove y-axis text
          axis.title= ggplot2::element_blank())
  return(plot_dncit)
}

## y-axis as tree diagram
create_y_axis <- function(settings){
  settings$pathString <- paste(" ",
                               settings$confounder,
                               settings$fct_relation,
                               settings$embedding_map,
                               settings$noise_x,
                               sep = "/")
  settings.tree <- data.tree::as.Node(settings)
  settings.phylo <- ape::as.phylo(settings.tree)
  #reorder matching rows of plots
  #ggtree::ggtree(settings.phylo) +
  #  geom_text(aes(label=node), hjust=-.3)
  settings.phylo <- ape::rotate(settings.phylo,node=25, polytom=c(26,39))
  y_axis <-ggtree::ggtree(settings.phylo, layout='rectangular') +
    #geom_nodepoint()+
    #geom_tippoint() +
    ggplot2::geom_label(ggplot2::aes(x=branch,label=label)) +
    ggplot2::labs(title=' ')
  #y_axis_b <- ggplot2::ggplot_build(y_axis)
  return(y_axis)
}
