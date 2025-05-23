#Copyright 2024 bioinformatics-policlinicogemelli

#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

args <- commandArgs(trailingOnly = TRUE)

file_pacchetti <- args[1]

# Read packages name from text
pacchetti <- readLines(file_pacchetti)

# Install packages
for (pacchetto in pacchetti) {
  if (!require(pacchetto, quietly = TRUE)) {
    install.packages(pacchetto, repos = "https://cran.r-project.org", dependencies = TRUE)
  }
}

pak::pak('tidyverse', dependencies = TRUE)
pak::pak('spatstat', dependencies = TRUE)
download.file(url = "https://cran.r-project.org/src/contrib/Archive/spatstat.core/spatstat.core_2.4-4.tar.gz", destfile =  "spatstat.core_2.4-4.tar.gz")
install.packages(pkgs='spatstat.core_2.4-4.tar.gz', type='source', repos=NULL)
BiocManager::install("EBImage")
remotes::install_github('r-lib/textshaping')
remotes::install_version("rjson", "0.2.20")
remotes::install_github('akoyabio/rtree')
remotes::install_github('akoyabio/phenoptr')
remotes::install_github('akoyabio/phenoptrReports')
