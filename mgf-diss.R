# This script parses and groups MGF (mascot generic format) MS/MS spectra from multiple files
# into a single file. The reordering and grouping of particular entries (between BEGIN IONS
# and END IONS statements) is based on their parent ion m/z and retention time values (the
# grouping is done by applying a selected window of m/z and RT values, respectively).
# A minimalistic user interface is provided through svDialogs package.

# load the libraries
library(readr)
library(stringr)
library(data.table)
library(svDialogs)

# rename files
files <- choose.files()
newfiles <- gsub("\\.mgf", "\\.txt", files)
file.rename(files, newfiles)
# ...and read them in
texty <- gsub("\r\n", "\n", lapply(newfiles, read_file)) #  the write.table function below introduces an extra \r at each new line...
write.table(texty, "merged.txt", quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)

test <- read_file("merged.txt")

# rename the files back to .mgf
file.rename(newfiles, files)

# create a list of individual AutoMS mgf spectra, split on "END IONS"
q <- unlist(str_split(test, "END IONS\r\n")) # str_split returns a 1 element list - use unlist...

# "END IONS" was stripped, add it to the end of every element
s <- lapply(q, paste, "END IONS\r\n", sep="")

# extract values from strings like this...
# precursor mass
p <- function(i) {
  x <- str_extract(s[[i]],"PEPMASS=\\d{1,4}\\.\\d{1,}") #  use regex for PEPMASS
  y <- as.numeric(unlist(str_split(x, pattern = "="))[2]) # take the second part, i.e. the numeric value
  return(y)
}

# rt in seconds
t <- function(i) {
  x <- str_extract(s[[i]],"RTINSECONDS=\\d{1,4}(\\.\\d{1,})?") # in some cases RTINSECS might be without the decimal part
  y <- as.numeric(unlist(str_split(x, pattern = "="))[2]) # take the second part, i.e. the numeric value
  return(y)
}

# sum of intensities
sumi <- function(i) {
  sum_i <- str_extract_all(s[[i]], "\\t\\d{1,}\\.\\d{1,}") # the intensity value comes always after a Tab
  sumec <-  sum(as.numeric(unlist(lapply(sum_i, str_split, boundary("word")))))
  return(sumec)
}
# create a list of lists, sth like this...(change the number of iterations based on the length of the list)
el <- list()
for (i in 1:length(s)) { 
  el[[i]] <- list(x = p(i), y = t(i), w = sumi(i), z = s[[i]])
}

# convert the el list to a data table through a dataframe
dt <- setDT(data.frame(matrix(unlist(el), nrow = length(el), byrow = TRUE), stringsAsFactors = FALSE))

# change character columns to numbers
dt$X1 <- as.numeric(dt$X1)
dt$X2 <- as.numeric(dt$X2)
dt$X3 <- as.numeric(dt$X3)

# rename columns
oldColNames <- c("X1", "X2", "X3", "X4")
newColNames <- c("parent", "time", "sumint", "mgf")
setnames(dt, old = oldColNames, new = newColNames)

# delete rows where sumint == 0
dt <- dt[-which(sumint == 0)]
# order dt by retention time
setorderv(dt, cols = "time", order = 1L)

# clean up...
rm(texty, test, q, s, el)
file.remove("merged.txt")

# create a new data.table from groups by same parent and time
# the grouping limits (criteria) are "ppm" for parent ions in ppm and
# "timeDelta" in seconds for retention times - defaults are 60 ppm and 60 s
ppm <- as.numeric(dlg_input(message = "Enter mass precission (ppm)", default = "60")$res)
if (identical(ppm, character(0))) {
  dlg_message("The script will terminate", "ok")
  q(save = "no")
} 

timeDelta <- as.numeric(dlg_input(message = "Enter RT tolerance (s)", default = "60")$res)
if (identical(ppm, character(0))) {
  dlg_message("The script will terminate", "ok")
  q(save = "no")
} 

while (nrow(dt) > 1) {
  
  temp_dt <- dt[1]
  del_vec <- c(1)
  
  j <- 2
  
  while ((dt$time[j] - dt$time[1]) < timeDelta ) {
    
    # ppm mass defect formula
    delta <- 1000000 * abs(dt$parent[j] - dt$parent[1]) / dt$parent[1]
    
    if (delta < ppm) {
      temp_dt <- rbind(temp_dt, dt[j])
      del_vec <- append(del_vec, j) # construct a vector of row numbers to delete them from dt later
    }
    if (j < nrow(dt)) {
      j <- j + 1
    } else
        break
  }
  
  dt <- dt[-del_vec] # remove the just snipped out rows from dt
  # order the temp_dt by sumint highest value first
  setorderv(temp_dt, cols = "sumint", order = -1L)
  
  # write output to text file
  for (k in 1:nrow(temp_dt)) {
    cat(str_remove_all(temp_dt$mgf[k], "\r"), file = "output.txt", sep = "", append = TRUE) # remove the carriage return "\r" because it introduces an extra empty line
  }
}

# clean up...
rm(dt, temp_dt)

# save the resulting merged .txt as .mgf in the same location as the original .mgfs
brokenC <- unlist(str_split(files[1], "\\\\"))
brokenC1 <- brokenC[1:(length(brokenC)-1)]

pasted <- ""

for (i in 1:length(brokenC1)) {
  pasted <- paste0(pasted, brokenC1[i], "\\")
}

file.rename("output.txt", paste0(pasted, "merged_A.mgf"))

