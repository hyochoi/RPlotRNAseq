split.char <- function(barcode){
  tmp.split <- strsplit(barcode,"-")[[1]] ;
  return(paste0(tmp.split[1],"-",tmp.split[2],"-",tmp.split[3],"-",tmp.split[4]))
}
