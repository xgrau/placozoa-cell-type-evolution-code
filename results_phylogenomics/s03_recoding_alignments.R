# libraries
library("Biostrings")

args = commandArgs(trailingOnly=TRUE)
ali_fn = args[1]

# read alignment
ali   = Biostrings::readAAMultipleAlignment(ali_fn)

# recoding dictionaries
recoding_dictionary = list(
  "SR4" = list(
    "A" = "AGNPST",
    "C" = "CHWY",
    "G" = "DEKQR",
    "T" = "FILMV"
  ),
  "Dayhoff6" = list(
    "0" = "AGPST",
    "1" = "DENQ", 
    "2" = "HKR", 
    "3" = "ILMV",
    "4" = "FWY",
    "5" = "C"
  ),
  "SR6" = list(
    "0" = "APST",
    "1" = "DENG",
    "2" = "QKR",
    "3" = "MIVL",
    "4" = "WC",
    "5" = "FYH"
  ),
)

# references:
# Dayhoff6: M.O. Dayhoff, R.M. Schwartz, B.C. Orcutt, A model of evolutionary change in proteins
# SR6: On reduced amino acid alphabets for phylogenetic inference, Mol. Biol. Evol., 24 (2007), pp. 2139-2150

for (dic in names(recoding_dictionary)) {
  
  # where does the output go?
  out_fn = gsub(".fasta", sprintf(".rec%s.fasta", dic), ali_fn)
  
  # apply recoding
  ali_r = ali
  for (i in 1:length(recoding_dictionary[[dic]])) {
    n = names(recoding_dictionary[[dic]])[i]
    message(sprintf("%s recoding, change to %s", dic, n))
    for (a in strsplit(recoding_dictionary[[dic]][[i]], split = "")[[1]] ) {
      ali_r = Biostrings::chartr(old = a, new = n, x = ali_r)
    }
  }
  
  # encode ambiguous characters as gaps
  ali_r = Biostrings::chartr(old = "X", new = "-", x = ali_r)
  ali_r = Biostrings::chartr(old = "N", new = "-", x = ali_r)
  
  # save as string set
  ali_r = Biostrings::AAStringSet(ali_r)
  
  # save as fasta
  Biostrings::writeXStringSet(ali_r, filepath = out_fn)
  
}
