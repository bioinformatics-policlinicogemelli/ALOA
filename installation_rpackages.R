args <- commandArgs(trailingOnly = TRUE)

# Leggi il percorso del file di testo dai parametri della riga di comando
file_pacchetti <- args[1]

# Leggi i nomi dei pacchetti dal file di testo
pacchetti <- readLines(file_pacchetti)

# Installa i pacchetti
for (pacchetto in pacchetti) {
  install.packages(pacchetto, repos = "https://cran.r-project.org", dependencies = TRUE)
}