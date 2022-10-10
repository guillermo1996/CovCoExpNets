#' Plots the RMSE
#'
#' For the given grouped results (see \link[CovCoExpNets]{groupMetrics}),
#' it plots the RMSE for the conditions and its 95% confidence
#' interval.
#'
#' @param df.group data.frame obtained from the function \link[CovCoExpNets]{groupMetrics}
#' @param output.file path the output file to save the plot. If not provided, the plot will
#'   not be stored
#' @param title string, title of the plot
#' @param ylab string, y-axis title
#' @param legend.title string, title of the legend
#'
#' @return ggplot object containing the RMSE comparison of the two groups
#' @export
plotRMSE <- function(df.group, output.file = NA, title= "", ylab = "", legend.title = "Algorithm", first.arg = "GLMNET", second.arg = "CovCoExpNets"){
  require(ggplot2)
  #df.rmse.combined = rbind(cbind(df.group.A, Group = first.arg), cbind(df.group.B, Group = second.arg))
  plot = ggplot(df.group, aes(x = reorder(Condition, rmse.test.mean), y = rmse.test.mean)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", show.legend = T, fill = "steelblue") +
    geom_errorbar(aes(ymin = rmse.test.mean - rmse.test.ci, ymax = rmse.test.mean + rmse.test.ci, y = rmse.test.mean, col = "95% C.I."),
                  width = 0.2,
                  color = "black", show.legend = T,
                  position = position_dodge(0.9)) +
    scale_color_manual(name = "95% C.I.", values=c("95% C.I." = "black"), guide = "legend") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim=c(6, 1.02*(max(df.group$rmse.test.mean) + max(df.group$rmse.test.ci)))) +
    ggtitle(title) + labs(y = ylab) +
    theme(plot.title = element_text(size = 22, face = "bold"),
          panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1),
          axis.text.x = ggplot2::element_text(face = "bold", color = "black", angle = 30, size = 15, hjust = 1),
          axis.text.y = ggplot2::element_text(face = "bold", color = "black", size = 15),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(face = "bold", size = 20, margin=margin(0,15,0,0)),
          panel.grid.minor = element_line(color = "gray", size = 0.1, linetype = 2),
          panel.grid.major = element_line(color = "gray", size = 0.2, linetype = 2),
          panel.background = element_rect(fill = "#F5F5F5"),
          aspect.ratio = 9/20,
          legend.title = element_text(size=20),
          legend.text = element_text(size=15)
    )

  if(!is.na(output.file)){
    ggsave(output.file, plot, width = 16, height = 9, dpi = 300, units = "in")
  }
  plot
}


#' Plots the RMSE
#'
#' For the two given groups (see \link[CovCoExpNets]{groupMetrics}),
#' it plots the RMSE for the conditions and its 95% confidence
#' interval.
#'
#' @param df.group.A data.frame obtained from the function \link[CovCoExpNets]{groupMetrics}
#' @param df.group.B data.frame obtained from the function \link[CovCoExpNets]{groupMetrics}
#' @param output_file path the output file to save the plot. If not provided, the plot will
#'   not be stored
#' @param title string, title of the plot
#' @param ylab string, y-axis title
#' @param legend.title string, title of the legend
#' @param first.arg string, group name of df.group.A
#' @param second.arg string, group name of df.group.B
#'
#' @return ggplot object containing the RMSE comparison of the two groups
#' @export
plotComparisonRMSE <- function(df.group.A, df.group.B, output.file = NA, title= "", ylab = "", legend.title = "Algorithm", first.arg = "GLMNET", second.arg = "CovCoExpNets"){
  require(ggplot2)
  df.rmse.combined = rbind(cbind(df.group.A, Group = first.arg), cbind(df.group.B, Group = second.arg))

  plot = ggplot(df.rmse.combined, aes(x = reorder(Condition, rmse.test.mean), y = rmse.test.mean, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", show.legend = T) +
    geom_errorbar(aes(ymin = rmse.test.mean - rmse.test.ci, ymax = rmse.test.mean + rmse.test.ci, y = rmse.test.mean, col = "95% C.I."),
                  width = 0.2,
                  color = "black", show.legend = T,
                  position = position_dodge(0.9)) +
    scale_color_manual(name = "95% C.I.", values=c("95% C.I." = "black"), guide = "legend") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim=c(6, 1.02*(max(df.rmse.combined$rmse.test.mean) + max(df.rmse.combined$rmse.test.ci)))) +
    ggtitle(title) + labs(y = ylab) +
    theme(plot.title = element_text(size = 22, face = "bold"),
          panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1),
          axis.text.x = ggplot2::element_text(face = "bold", color = "black", angle = 30, size = 15, hjust = 1),
          axis.text.y = ggplot2::element_text(face = "bold", color = "black", size = 15),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(face = "bold", size = 20, margin=margin(0,15,0,0)),
          panel.grid.minor = element_line(color = "gray", size = 0.1, linetype = 2),
          panel.grid.major = element_line(color = "gray", size = 0.2, linetype = 2),
          panel.background = element_rect(fill = "#F5F5F5"),
          aspect.ratio = 9/20,
          legend.title = element_text(size=20),
          legend.text = element_text(size=15)
    ) + guides(fill=guide_legend(title=legend.title))

  if(!is.na(output.file)){
    ggsave(output.file, plot, width = 16, height = 9, dpi = 300, units = "in")
  }
  plot
}

#' Plots the adjusted R2
#'
#' For the two given groups (see \link[CovCoExpNets]{groupMetrics}),
#' it plots the adjusted R2 for the conditions and its 95% confidence
#' interval.
#'
#' @param df.group.A data.frame obtained from the function \link[CovCoExpNets]{groupMetrics}
#' @param df.group.B data.frame obtained from the function \link[CovCoExpNets]{groupMetrics}
#' @param output_file path the output file to save the plot. If not provided, the plot will
#'   not be stored
#' @param title string, title of the plot
#' @param ylab string, y-axis title
#' @param legend.title string, title of the legend
#' @param first.arg string, group name of df.group.A
#' @param second.arg string, group name of df.group.B
#'
#' @return ggplot object containing the adjusted R2 comparison of the two groups
#' @export
plotComparisonR2_adj <- function(df.group.A, df.group.B, output.file = NA, title= "", ylab = "", legend.title = "Algorithm", first_arg = "GLMNET", second_arg = "CovCoExpNets"){
  require(ggplot2)
  df.rmse.combined = rbind(cbind(df.group.A, Group = first_arg), cbind(df.group.B, Group = second_arg))

  #df.rmse.combined =  rbind(cbind(df.group.glmnet, Group = "GLMNET"), cbind(df.group.covco, Group = "CovCoExpNets"))
  plot.order = df.rmse.combined %>%
    filter(Group == "GLMNET") %>%
    arrange(r2_adj.train.mean) %>%
    mutate(labels = factor(Condition))
  df.rmse.combined = df.rmse.combined %>% mutate(labels = factor(Condition, levels = plot.order$labels, ordered = T))

  plot = ggplot(df.rmse.combined, aes(x = labels, y = r2_adj.train.mean, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", show.legend = T) +
    geom_errorbar(aes(ymin = r2_adj.train.mean - r2_adj.train.ci, ymax = r2_adj.train.mean + r2_adj.train.ci, y = r2_adj.train.mean, col = "95% C.I."),
                  width = 0.2,
                  color = "black", show.legend = T,
                  position = position_dodge(0.9)) +
    scale_color_manual(name = "95% C.I.", values=c("95% C.I." = "black"), guide = "legend") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim=c(0.99*min(df.rmse.combined$r2_adj.train.mean)-max(df.rmse.combined$r2_adj.train.ci), 1)) +
    ggtitle(title) + labs(y = ylab) +
    theme(plot.title = element_text(size = 22, face = "bold"),
          panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1),
          axis.text.x = ggplot2::element_text(face = "bold", color = "black", angle = 30, size = 15, hjust = 1),
          axis.text.y = ggplot2::element_text(face = "bold", color = "black", size = 15),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(face = "bold", size = 20, margin=margin(0,15,0,0)),
          panel.grid.minor = element_line(color = "gray", size = 0.1, linetype = 2),
          panel.grid.major = element_line(color = "gray", size = 0.2, linetype = 2),
          panel.background = element_rect(fill = "#F5F5F5"),
          aspect.ratio = 9/20,
          legend.title = element_text(size=20),
          legend.text = element_text(size=15),
          plot.margin = margin(0, 0, 0, 2.2, "cm")
    ) + guides(fill=guide_legend(title=legend.title))

  if(!is.na(output.file)){
    ggsave(output.file, plot, width = 16, height = 9, dpi = 300, units = "in")
  }
  plot
}

#' Plots the returned predictors
#'
#' For the two given groups (see \link[CovCoExpNets]{groupMetrics}),
#' it plots the returned predictors for the conditions and its 95% confidence
#' interval.
#'
#' @param df.group.A data.frame obtained from the function \link[CovCoExpNets]{groupMetrics}
#' @param df.group.B data.frame obtained from the function \link[CovCoExpNets]{groupMetrics}
#' @param output_file path the output file to save the plot. If not provided, the plot will
#'   not be stored
#' @param title string, title of the plot
#' @param ylab string, y-axis title
#' @param legend.title string, title of the legend
#' @param first.arg string, group name of df.group.A
#' @param second.arg string, group name of df.group.B
#'
#' @return ggplot object containing the returned predictors comparison of the two groups
#' @export
plotComparisonReturnedPredictors <- function(df.group.A, df.group.B, output.file = NA, title= "", ylab = "", legend.title = "Algorithm", first_arg = "GLMNET", second_arg = "CovCoExpNets"){
  require(ggplot2)
  df.rmse.combined = rbind(cbind(df.group.A, Group = first_arg), cbind(df.group.B, Group = second_arg))

  #df.rmse.combined =  rbind(cbind(df.group.glmnet, Group = "GLMNET"), cbind(df.group.covco, Group = "CovCoExpNets"))
  plot.order = df.rmse.combined %>%
    filter(Group == "GLMNET") %>%
    arrange(Returned.predictors.mean) %>%
    mutate(labels = factor(Condition))
  df.rmse.combined = df.rmse.combined %>% mutate(labels = factor(Condition, levels = plot.order$labels, ordered = T))

  plot = ggplot(df.rmse.combined, aes(x = labels, y = Returned.predictors.mean, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", show.legend = T) +
    geom_errorbar(aes(ymin = Returned.predictors.mean - Returned.predictors.ci, ymax = Returned.predictors.mean + Returned.predictors.ci, y = Returned.predictors.mean, col = "95% C.I."),
                  width = 0.2,
                  color = "black", show.legend = T,
                  position = position_dodge(0.9)) +
    scale_color_manual(name = "95% C.I.", values=c("95% C.I." = "black"), guide = "legend") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim=c(0,#0.97*min(df.rmse.combined$Returned.predictors.mean) - max(df.rmse.combined$Returned.predictors.ci),
                           1.01*max(df.rmse.combined$Returned.predictors.mean) + max(df.rmse.combined$Returned.predictors.ci)
    )) +
    ggtitle(title) + labs(y = ylab) +
    theme(plot.title = element_text(size = 22, face = "bold"),
          panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1),
          axis.text.x = ggplot2::element_text(face = "bold", color = "black", angle = 30, size = 15, hjust = 1),
          axis.text.y = ggplot2::element_text(face = "bold", color = "black", size = 15),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(face = "bold", size = 20, margin=margin(0,15,0,0)),
          panel.grid.minor = element_line(color = "gray", size = 0.1, linetype = 2),
          panel.grid.major = element_line(color = "gray", size = 0.2, linetype = 2),
          panel.background = element_rect(fill = "#F5F5F5"),
          aspect.ratio = 9/20,
          legend.title = element_text(size=20),
          legend.text = element_text(size=15)
    ) + guides(fill=guide_legend(title=legend.title))

  if(!is.na(output.file)){
    ggsave(output.file, plot, width = 16, height = 9, dpi = 300, units = "in")
  }
  plot
}


#' ToDo
#'
#' ToDo
#'
#' @param genes.freq data.frame obtained from the function \link[CovCoExpNets]{groupMetrics}
#' @param heatplots data.frame obtained from the function \link[CovCoExpNets]{groupMetrics}
#' @param diag (optinal) which parameter to register in the diagional. Only
#'   other possible value is "genes". By default, the jaccard index.
#' @param show.text string, title of the plot
#' @param output.path string, y-axis title
#' @param draw.plot string, title of the legend
#'
#' @return ToDo
#' @export
plotAllSimilarity = function(genes.freq, heatplots, diag = "default", show.text = T, output.path = NA, draw.plot = T){
  require(ggplot2)

  if(missing(heatplots)) heatplots = measureStabilityJaccard(genes.freq, diag = diag)
  return.list = is(genes.freq, "list")
  if(!return.list){
    genes.freq = list(genes.freq)
    heatplots = list(heatplots)
  }

  r = foreach(i = 1:length(genes.freq), .combine = "c") %do%{
    genes.freq.i = genes.freq[[i]]
    max.iter = max(genes.freq.i$iter)
    heatmap = heatplots[[i]]

    df.heatmap = expand.grid(X = as.character(1:max.iter), Y = as.character(1:max.iter))
    df.heatmap$text = c(heatmap)
    df.heatmap$similarity = c(heatmap)

    heatplot = ggplot(df.heatmap, aes(x = X, y = Y, fill = similarity)) +
      geom_tile(color = "black", lwd = 0.25, linetype = 1) +
      xlab("Repetition") + ylab("Repetition") +
      scale_fill_gradientn(colors = hcl.colors(50, "RdYlGn")[1:47], name = paste0("Jaccard index\n"), limits = c(0, 1)) +
      theme(aspect.ratio = 1) +
      theme(plot.title = element_text(size = 22, face = "bold"),
            panel.border = element_rect(colour = "black", fill = NA, size = 1),
            axis.text.x = element_text(face = "bold", colour = "black", size = 15),
            axis.text.y = element_text(face = "bold", colour = "black", size = 15),
            axis.title.x = element_text(face = "bold", size = 20),
            axis.title.y = element_text(face = "bold", size = 20, margin=margin(0,15,0,0)),
            panel.grid.minor = element_line(color = "gray", size = 0.1, linetype = 2),
            panel.grid.major = element_line(color = "gray", size = 0.2, linetype = 2),
            panel.background = element_rect(fill = "#F5F5F5"),
            legend.key.size = unit(0.5, 'cm'), #change legend key size
            legend.key.height = unit(2, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=15, face = "bold"), #change legend title font size
            legend.text = element_text(size=12))
    if(return.list) heatplot = heatplot + ggtitle(paste0("Tissue: ", names(genes.freq)[i]))
    if(show.text) heatplot = heatplot + geom_text(aes(label = round(text, 2)), size = 6)

    if(draw.plot) heatplot %>% print()

    if(!is.na(output.path)){
      output.file = if(return.list) paste0(output.path, gsub(" ", "_", names(genes.freq)[i]), ".png") else output.path
      ggsave(output.file, heatplot, width = 10, height = 9, dpi = 300, units = "in")
    }

    heatplot
  }
}

plotConditionCoincidences <- function(genes.relevant, title = "Titulo", diag = "default", output.file = NA){
  heatmap = matrix(0, nrow=length(genes.relevant), ncol=length(genes.relevant))
  r = foreach(i = 1:length(genes.relevant)) %:%
    foreach(j = i:length(genes.relevant)) %do% {
      common.genes = length(intersect(genes.relevant[[i]], genes.relevant[[j]]))

      if(i == j){
        if(is.na(diag)){
          heatmap[i, j] = NA
        }else if(diag == "default"){
          heatmap[i, j] = common.genes
        }else{
          heatmap[i, j] = diag
        }
      }else{
        heatmap[i, j] = common.genes
        heatmap[j, i] = common.genes
      }
    };

  #heatmap[upper.tri(heat_map)] = NA
  plot.coin = expand.grid(X = names(genes.relevant), Y = names(genes.relevant))
  plot.coin$text = c(heatmap)
  diag(heatmap) = NA
  plot.coin$color = c(heatmap)
  plot.coin = plot.coin %>% filter(!is.na(text))

  library(pheatmap)

  plot.coincidences = ggplot(plot.coin, aes(x = X, y = Y, fill = color)) +
    geom_tile(color = "black", lwd = 0.3, linetype = 1) +
    geom_text(aes(label = round(text, 1)), size = 5) +
    theme_bw() +
    ggtitle(title) +
    scale_fill_gradientn(colors = hcl.colors(30, "oslo")[1:20], name = "Common genes", guide = "colourbar", limits = c(0, max(plot.coin$color, na.rm = T))) + #, breaks=seq(0.80, 1, 0.05), limits = c(0.8, 1), na.value = "#D3D3D3") +
    coord_fixed() +
    scale_x_discrete(limits = rev) +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(2, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=12),
          panel.border = element_blank(),
          plot.title = element_text(size = 20, face = "bold"),
          axis.text.x = ggplot2::element_text(angle = 30, size = 11, hjust = 1, vjust = 1),
          axis.text.y = ggplot2::element_text(angle = 30, size = 11, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank()) +
    guides(fill = guide_colourbar(frame.linewidth = 0.5,frame.colour = "black")); plot.coincidences
  if(!is.na(output.file)){
    ggsave(output.file, width = 10, height = 9, dpi = 300, units = "in")
  }
  plot.coincidences
}
