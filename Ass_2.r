data_dir <- file.path("R/CGS_data/data","SRP073813")
data_file <- file.path(data_dir,"SRP073813.tsv")
metadata_file <- file.path(data_dir,"metadata_SRP073813.tsv")
# Read Data ---------------------------------------------------------------
library(org.Hs.eg.db)
library(magrittr)
metadata<- readr::read_tsv(metadata_file)
print(metadata)
expression_df <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("Gene")
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

# mapping -----------------------------------------------------------------

mapped_list <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = expression_df$Gene,
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "ENTREZID", # The type of gene identifiers you would like to map to
  multiVals = "list"
)
head(mapped_list)

mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "Entrez") %>%
  tidyr::unnest(cols = Entrez)

multi_mapped <- mapped_df %>%
  # Let's count the number of times each Ensembl ID appears in `Ensembl` column
  dplyr::count(Ensembl, name = "entrez_id_count") %>%
  # Arrange by the genes with the highest number of Entrez IDs mapped
  dplyr::arrange(desc(entrez_id_count))
head(multi_mapped)
collapsed_mapped_df <- mapped_df %>%
  # Group by Ensembl IDs
  dplyr::group_by(Ensembl) %>%
  # Collapse the Entrez IDs `mapped_df` into one column named `all_entrez_ids`
  dplyr::summarize(all_entrez_ids = paste(Entrez, collapse = ";"))
collapsed_mapped_df %>%
  # Filter `collapsed_mapped_df` to include only the rows where
  # `all_entrez_ids` values include the ";" character --
  # these are the rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_entrez_ids, ";")) %>%
  # We only need a preview here
  head()
final_mapped_df <- data.frame(
  "HUGO" = mapIds(
    org.Hs.eg.db, # Replace with annotation package for your organism
    keys = expression_df$Gene,
    keytype = "ENSEMBL", # Replace with the gene identifiers used in your data
    column = "SYMBOL", # The type of gene identifiers you would like to map to
    multiVals = "first" # Keep only the first mapped value for each Ensembl ID
  )
) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  # Add the multiple mappings data from `collapsed_mapped_df` using Ensembl IDs
  dplyr::inner_join(collapsed_mapped_df, by = "Ensembl") %>%
  # Now let's add on the rest of the expression data
  dplyr::inner_join(expression_df, by = c("Ensembl" = "Gene"))
final_mapped_df %>%
  # Filter `final_mapped_df` to rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_entrez_ids, ";")) %>%
  head()
readr::write_tsv(final_mapped_df, file.path(
  data_dir,
  "SRP073813_Entrez_IDs.tsv" # Replace with a relevant output file name
))




