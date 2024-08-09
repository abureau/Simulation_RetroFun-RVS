res <- list()
n_sim <- 35
for (i in 1:n_sim){
  res <- append(res, readRDS(paste0("sim_object_arcturus_", i, ".rds")))
}
saveRDS(res, "sim_object.rds")
