underscore <- function(species){
  
  species_with_underscores <- gsub(" ", "_", species)
  return(species_with_underscores)
  
}

ununderscore <- function(species){
  
  species_without_underscores <- gsub("_", " ", species)
  return(species_without_underscores)
  
}