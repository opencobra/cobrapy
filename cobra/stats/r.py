#cobra.r.py
#NOTE: This code is experimental.  Do not use.
from os import listdir
from rpy2 import robjects
import rpy2.robjects.numpy2ri
r = robjects.r
#Load the necessary R libraries and create the helper functions.
import cobra
print "WARNING cobra.stats.r is experimental and it's use is strongly "+\
      "discouraged.  Use rpy2 instead."
r_path = cobra.__path__[0] + '/stats/r_scripts/'
if 'tools.R' in listdir(r_path):
    r('''source("%stools.R")'''%r_path)

    
r( '''
   library(gplots)
   library(CarbonEL)
   library(RColorBrewer)
   ps.options(pointsize=10)
   # Generate red-to-green colorscale
   redgreen <- function(n) colorpanel(n, "red", "black", "green")
   greenred <- function(n) colorpanel(n, "green", "black", "red")
   bluered  <- function(n) colorpanel(n, "blue","white","red")
   redblue  <- function(n) colorpanel(n, "red","white","blue")
   whitered  <- function(n) colorpanel(n, "white","red")
   whiteblack  <- function(n) colorpanel(n, "white","black")
   whiteblue  <- function(n) colorpanel(n, "white","blue")
   blackred  <- function(n) colorpanel(n, "black", "red")
   greenblack  <- function(n) colorpanel(n, "green", "black")
   update_names <- function( the_vector, the_names) {
      names( the_vector) = the_names
      return( the_vector)
   }
   calculate_growth <- function(the_data, x_variable = "Time", start_row = 1, finish_row = 5) {
      return( apply( the_data[ start_row:finish_row, ], 2, function( y) lm(y ~ the_data[ start_row:finish_row, x_variable ])$coeff[ 2 ]))
   
   }
   read_table <- function(row_names=1,...) {
      return(read.table(row.names=row_names,...))
    }
   col_names <- function(the_matrix, the_names) {
      colnames(the_matrix) = the_names
      return(the_matrix)
   }
   delete_columns <- function(the_matrix, the_columns){
      return(the_matrix[,-the_columns])
      }
   row_names <- function(the_matrix, the_names) {
      rownames(the_matrix) = the_names
      return(the_matrix)
   }
   ncaHeatmap <- function(the_data, the_names, the_scores, file_name = FALSE,
                           scale_data = FALSE) {
      require(gplots)
      if (file_name == FALSE) {
         quartz()
      } else {
         pdf(file_name)
      }
      if (scale_data) the_data = scaleData(the_data)
      heatmap.2(as.matrix(the_data[, the_names]),   dendrogram = "none",Rowv=FALSE,Colv=FALSE,scale="none",density.info="none", trace = "none",col =bluered(31),lwid=c(2,5),lhei=c(10,3),lmat = rbind(c(2,1),c(4,3)), labCol=paste(the_scores, the_names, sep =  "-"))
      if (file_name != FALSE) {
          dev.off()
      }
  }

  build_dendrogram <- function(the_data, distance_method = "pearson",
                                linkage_method = "complete") {
      #Builds a dendrogram based on the distances between the rows in a matrix
      if (distance_method == "pearson") {
           the_dist = pearsonDistance(t(the_data))
      } else {
           the_dist = dist(the_data, method = distance_method)
      }
      return(as.dendrogram(hclust(the_dist, method = linkage_method))) 

  }
  write_table <-function(the_data, the_file, row_names = TRUE, col_names = TRUE) {
    write.table(the_data, file = the_file, row.names = row_names, sep="\t",
                 quote = FALSE, col.names = col_names)
  }
  heatmap2 <- function(the_data, file_name = FALSE,
                           dendrogram_type = "none", distance_method = "euclidian",
                           linkage_method = "complete",
                           scale_data = FALSE, order_on_distance_from = "",
                           color_panel = NULL, number_of_colors = 31,pointsize=6,...) {
      require(CarbonEL)
      require(gplots)
      the_rowv = FALSE
      the_colv = FALSE
      if (order_on_distance_from != "") {
          if(distance_method == "euclidian" || any(the_data == NA)) {
             the_order = names(sort(as.matrix(dist(t(the_data)))[ order_on_distance_from,]))
          }
          the_data = the_data[,the_order]
      } else {
         #Calculate the distances and build the dendrogram, if specified.
         if (dendrogram_type == "column"  ||
            dendrogram_type == "both") {
            the_colv = build_dendrogram(t(the_data), distance_method,
                                        linkage_method)
         }
         if (dendrogram_type == "row"  ||
            dendrogram_type == "both") {
            the_rowv = build_dendrogram(the_data, distance_method,
                                      linkage_method)
         }
      }
      if (file_name == FALSE) {
         quartz()
      } else {
         pdf(file_name,pointsize=pointsize)
      }
      if (is.null(color_panel) ) {
          color_panel = bluered(number_of_colors)
      } 
      if (scale_data) the_data = scaleData(the_data)
      heatmap.2(as.matrix(the_data), dendrogram = dendrogram_type,
          Rowv=the_rowv,Colv=the_colv,scale="none",density.info="none",
          trace = "none",col = color_panel,lwid=c(1,1),lhei=c(3,15),# mar=c(5,15),
          #lmat = rbind(c(4,3),c(2,1)),
	  na.color = "grey",...)
      if (file_name != FALSE) {
          dev.off()
      }
  }
  histogram_of_categorical <- function(the_categorical_list, show_plot = TRUE, ...) {
      the_categories = levels(as.factor(the_categorical_list))
      frequency_vector = c()
      for (the_category in the_categories) {
          frequency_vector = c(frequency_vector, length(which(the_categorical_list == the_category)))
      }
      names(frequency_vector) = the_categories
      frequency_vector = sort(frequency_vector, decreasing = TRUE)
      if (show_plot) {
          quartz()
          barplot(frequency_vector, ...)
      }
      return(frequency_vector)

  }
  save_data<-function(the_data, the_file){ save(the_data, file=the_file) }
  adjust_p <- function(the_data, the_method) { return(p.adjust(the_data, the_method)) }
'''
)
