process.group <- function(group, group.weight) {
  
  # Validation and correction of group and group.weight
  gf <- factor(group)
  label <- levels(gf)
  if (is.numeric(group) | is.integer(group)) {
    label <- paste0("Group ", label)
  }
  group <- as.integer(gf)
  
  if (any(levels(gf)=='Group 0')) {
    stop('The group name "0" is not allowed: use Z to fit a model with unpenalized variables', call.=FALSE)
  } 
  if (any(levels(gf)=='Unpenalized')) {
    stop('The group name "Unpenalized" is not allowed: the name is reserved for the variables in Z', call.=FALSE)
  } 
  
  if (missing(group.weight)) {
    group.weight <- rep(1, length(label))
    names(group.weight) <- label
  } else {
    if (any(group.weight < 0)) stop('A negative group.weight is not allowed', call.=FALSE)
    if (any(group.weight == 0)) stop('A group.weight of 0 is not allowed: use Z to fit a model with unpenalized variables', call.=FALSE)
    if (length(group.weight) != length(label)) stop("Group.weight should contain one value per group", call.=FALSE)
    if (storage.mode(group.weight) != "double") storage.mode(group.weight) <- "double"
  }
  return(structure(group, levels=label, group.weight=group.weight))
}
