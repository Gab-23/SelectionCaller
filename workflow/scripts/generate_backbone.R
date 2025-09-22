generate_backbone <- function(windows){
  
  `%>%` <- magrittr::`%>%`

  windows <- windows %>% 
              dplyr::group_by(chr) %>% 
              dplyr::mutate(binID = paste0(gsub(pattern = "chr", replacement = "", x = chr), 
                                    "_", 
                                    1:dplyr::n()
                                          )
                           )
  return(windows)

}