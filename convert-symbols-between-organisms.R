library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
drosphilia = useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
converted = getLDS(attributes = c("external_gene_name"),
                   mart = drosphilia,
                   attributesL = c("hgnc_symbol"), martL = human,
                   uniqueRows=T)
colnames(converted) = c("drosophila", "symbol")

hg = readr::read_csv("~/Downloads/Single-cell RNA-seq cell cycle markers - hsapiens.csv")

dros = hg %>%
  left_join(converted) %>%
  dplyr::select(phase, drosophila) %>%
  unique() %>%
  complete.cases()
colnames(dros) = c("phase", "symbol")
readr::write_csv(dros, "~/Downloads/dros_cycle.csv")
# Get a list of the worksheets (ws)
ws <- gs_ws_ls(gs)
print(ws)
