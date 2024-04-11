args <- commandArgs(trailingOnly = TRUE)

# Leggi il percorso del file di testo dai parametri della riga di comando
file_pacchetti <- args[1]

# Leggi i nomi dei pacchetti dal file di testo
pacchetti <- readLines(file_pacchetti)

# Installa i pacchetti
for (pacchetto in pacchetti) {
  if (!require(pacchetto, quietly = TRUE)) {
    install.packages(pacchetto, repos = "https://cran.r-project.org", dependencies = TRUE)
  }
}

# pak::pak('r-lib/ragg', dependencies = TRUE)
pak::pak('tidyverse', dependencies = TRUE)
pak::pak('spatstat', dependencies = TRUE)
download.file(url = "https://cran.r-project.org/src/contrib/Archive/spatstat.core/spatstat.core_2.4-4.tar.gz", destfile =  "spatstat.core_2.4-4.tar.gz")
install.packages(pkgs='spatstat.core_2.4-4.tar.gz', type='source', repos=NULL)
BiocManager::install("EBImage")
remotes::install_github('r-lib/textshaping')
remotes::install_github('akoyabio/phenoptr')
remotes::install_github('akoyabio/phenoptrReports')
