generate_backbone <- function(windows){                                                       # windows is a df with chr start end
  
  `%>%` <- magrittr::`%>%`

  windows <- windows %>% 
              dplyr::group_by(chr) %>%                                                        # for each chromosome
              dplyr::mutate(binID = paste0(gsub(pattern = "chr", replacement = "", x = chr),  # generate binID
                                    "_", 
                                    1:dplyr::n()
                                          )
                           )
  return(windows)

}