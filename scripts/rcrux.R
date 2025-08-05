setwd("/Users/bfj994/Documents/barcodeMiner/")

library(rCRUX)
library(taxonomizr)

accession_taxa_sql_path <- "accessionTaxa.sql"
prepareDatabase(accession_taxa_sql_path)
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "ncbi-blast-2.10.1+/bin/", sep = ":"))


forward_primer_seq = "GGGCAATCCTGAGCCAA"

reverse_primer_seq =  "CCATTGAGTCTCTGCACCTATC"

output_directory_path <- "output" # path to desired output directory

metabarcode_name <- "trnL" # desired name of metabarcode locus

accession_taxa_sql_path <- "accessionTaxa.sql" # path to taxonomizr sql database

blast_db_path <- "ncbi-blast-2.10.1+/DB/viri_full_chlo"  # path to blast formatted database

seeds_output_path <- "output/get_seeds_local/trnL_filtered_get_seeds_local_output_with_taxonomy.csv"

get_seeds_local(forward_primer_seq,
                reverse_primer_seq,
                metabarcode_name,
                output_directory_path,
                accession_taxa_sql_path,
                blast_db_path, evalue = 300)

blast_seeds(seeds_output_path,
            blast_db_path,
            accession_taxa_sql_path,
            output_directory_path,
            metabarcode_name)    
