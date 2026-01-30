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

# logger.R - Python-friendly version
# Copyright 2024 bioinformatics-policlinicogemelli

# This logger prints messages to stdout.
# Python (rpy2) captures stdout and writes it to the log file.

log4r_info <- function(message) {
  cat(sprintf("[%s] INFO: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
}

log4r_warn <- function(message) {
  cat(sprintf("[%s] WARNING: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
}

log4r_error <- function(message) {
  cat(sprintf("[%s] ERROR: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
}


