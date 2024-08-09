setGeneric("DropHaploLocus", function(trio, haplo.obs, geno.vec ) standardGeneric("DropHaploLocus"))

setMethod("DropHaploLocus", signature( trio = "Trio", haplo.obs = "matrix", geno.vec = "numeric"), function( trio, haplo.obs, geno.vec ){
  g1 <- geno.vec[trio@id]
  g2 <- geno.vec[trio@spouse]
  h1 <- haplo.obs[trio@id,]
  h2 <- haplo.obs[trio@spouse,]
  for (i in 1:length(trio@offspring)){
    if( is.character(trio@offspring[[i]]) ){
      if( any(is.na(haplo.obs[ trio@offspring[[i]],]))){
        hoff <- drop.haplo.locus.fn(g1,g2, geno.vec[trio@offspring[[i]]],h1,h2)
        haplo.obs[ trio@offspring[[i]],] <- hoff
      }
    }else{
      # this is the case where the offspring is a trio object
      if( any(is.na(haplo.obs[ trio@offspring[[i]]@id,]))){
        hoff <- drop.haplo.locus.fn(g1,g2, geno.vec[trio@offspring[[i]]@id],h1,h2)
        haplo.obs[ trio@offspring[[i]]@id,] <- hoff
        haplo.obs <- DropHaploLocus(trio@offspring[[i]], haplo.obs, geno.vec)
      }
    }
  }
  return( haplo.obs )
})

setGeneric("DropHaplo", function(trio, haplo.obs ) standardGeneric("DropHaplo"))

setMethod("DropHaplo", signature( trio = "Trio", haplo.obs = "matrix"), function( trio, haplo.obs){
  h1 <- haplo.obs[trio@id,]
  h2 <- haplo.obs[trio@spouse,]
  for (i in 1:length(trio@offspring)){
    if( is.character(trio@offspring[[i]]) ){
      if( any(is.na(haplo.obs[ trio@offspring[[i]],]))){
        hoff <- drop.haplo.fn(h1,h2)
        haplo.obs[ trio@offspring[[i]],] <- hoff
      }
    }else{
      # this is the case where the offspring is a trio object
      if( any(is.na(haplo.obs[ trio@offspring[[i]]@id,]))){
        hoff <- drop.haplo.fn(h1,h2)
        haplo.obs[ trio@offspring[[i]]@id,] <- hoff
        haplo.obs <- DropHaplo(trio@offspring[[i]], haplo.obs)
      }
    }
  }
  return( haplo.obs )
})
