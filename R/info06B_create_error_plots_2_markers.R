# create 2x2 plot of 2 markers by 2 accuracy metrics (% incorrect, % correct) 
  # and overlay species from a location on them (if sp_data is not NULL)

library(patchwork)
library(dplyr)
library(ggplot2)

create_error_plots_2_markers <- function(eplots, 
                                         rubric = "rdp90",
                                         marker_pair = c("MiMammalU_12S","Vences_16S"),
                                         sp_data = ct_data){
  
  # set up the data
  eplot1 <- eplots[[rubric]][[marker_pair[1]]]
  eplot2 <- eplots[[rubric]][[marker_pair[2]]]
  if(!is.null(sp_data)){
    sp_data1 <- sp_data %>% select(contains(paste0(marker_pair[1], "_", rubric))) %>%
      dplyr::select(nnd = contains("nnd"), n_seqs = contains("n_seqs"))
    sp_data2 <- sp_data %>% select(contains(paste0(marker_pair[2], "_", rubric))) %>%
      dplyr::select(nnd = contains("nnd"), n_seqs = contains("n_seqs"))
  }
  
  # create top left plot
  top_left <- eplot1$both$i + ggtitle(paste0(marker_pair[1], " marker"))
  
  if(!is.null(sp_data)){top_left <- top_left + 
    geom_point(aes(x = jitter(nnd), 
                   y = jitter(n_seqs+0.5)),
               data = sp_data1, shape = 21, color = "gray", fill = "transparent") 
  }
  
  # create top right plot
  top_right <- eplot1$both$c + ggtitle("")
  
  if(!is.null(sp_data)){top_right <- top_right + 
    geom_point(aes(x = jitter(nnd), 
                   y = jitter(n_seqs+0.5)),
               data = sp_data1, shape = 21, color = "gray", fill = "transparent") 
  }
  
  # create bottom left plot
  bottom_left <- eplot2$both$i + ggtitle(paste0(marker_pair[2], " marker"))
  
  if(!is.null(sp_data)){bottom_left <- bottom_left + 
    geom_point(aes(x = jitter(nnd), 
                   y = jitter(n_seqs+0.5)),
               data = sp_data2, shape = 21, color = "gray", fill = "transparent") 
  }
  
  # create bottom right plot
  bottom_right <- eplot2$both$c + ggtitle("")
  
  if(!is.null(sp_data)){bottom_right <- bottom_right + 
    geom_point(aes(x = jitter(nnd), 
                   y = jitter(n_seqs+0.5)),
               data = sp_data2, shape = 21, color = "gray", fill = "transparent") 
  }
  
  
  top_markers_eplot <- (top_left | top_right) / (bottom_left | bottom_right)
  
  return(top_markers_eplot)
  
}