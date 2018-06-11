# Generate example data for metaglomr.
set.seed(1000)
features <- matrix(
  rnbinom(n=1024, size=1, mu=100),
  nrow = 16
)

samples <- data.frame(
  sample_id = paste0("S_", c(1:16)),
  study_group = c(rep("case", 8), rep("control", 8)),
  sample_type = rep(c("A", "B", "C", "D"), 4)
)

taxa <- read.csv("data-raw/taxa.csv", row.names = 1, header = T)

rownames(features) <- samples$sample_id
colnames(features) <- rownames(taxa)

devtools::use_data(features, samples, taxa, overwrite = TRUE)
