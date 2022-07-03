process.penalty <- function(penalty,pvar,pgr,vargamma,grgamma,vartau,grtau,alpha) {
  
  # Validation and correction of penalty related input
  if(!is.null(penalty)){
    switch (penalty,
            "sgl" = pvar <- pgr <- "lasso",
            "sgs" = pvar <- pgr <- "scad",
            "sgm" = pvar <- pgr <- "mcp",
            "sge" = pvar <- pgr <- "exp")
  }
  
  if ((vargamma <= 1 | grgamma <= 1) & (pvar=="mcp" | pgr=="mcp")){
    warning("Gamma must be > 1 for MCP and was set to its default value.")
    vargamma <- grgamma <- 3
  } 
  if ((vargamma <= 2 | grgamma <= 2) & (pvar=="scad" | pgr=="scad")){
    warning("Gamma must be > 2 for SCAD and was set to its default value.")
    vargamma <- 4
  } 
  if (alpha > 1 | alpha < 0){
    warning("alpha must be in [0, 1] and was set to its default value.")
    alpha    <- 1/3
  } 
  
  pvar_int   <- 1
  pgr_int    <- 1
  
  if(pvar == "mcp")   pvar_int <- 2
  if(pvar == "exp")   pvar_int <- 3
  if(pvar == "scad")  pvar_int <- 4
  
  if(pgr == "mcp")   pgr_int <- 2
  if(pgr == "exp")   pgr_int <- 3
  if(pgr == "scad")  pgr_int <- 4
  
  return(list(pvar_int=pvar_int,pgr_int=pgr_int))
}
