library("miniCRAN")
library ("dplyr")

# use Revolution Analytics CRAN mirror
#revolution <- c(CRAN = "http://cran.microsoft.com")
#revolution <- c(CRAN = " https://www.icesi.edu.co/CRAN/")
mirrorCRAN <- c(CRAN = "http://cloud.r-project.org/")



# Specify list of packages to download
pkgs = c("packman", "rrBLUP", "parallel","config","dplyr","stringi","qqman",
		  "VennDiagram","RColorBrewer","circlize","gplots", "rmarkdown",
		   "kableExtra" ,"doParallel", "ldsep", "yaml", "BiocManager", "multtest")

pkgList <- pkgDep (pkgs, 
				   repos = mirrorCRAN, 
				   type = "source", 
				   suggests = FALSE)
                   #Rversion = "3.6")

ml_pkg <- pkgDep(c("foreach"), 
                 repos = mirrorCRAN, type = "source",
                 suggests = TRUE)

pkgList <- pkgList %>% union(ml_pkg)
localPath <- normalizePath("/opt/miniCRAN", winslash = "/")

# Create temporary folder for miniCRAN
dir.create(localPath, showWarnings=F)

# Make repo for source and win.binary
makeRepo(pkgList, path = localPath, repos = mirrorCRAN, type = c("source"))

# List all files in miniCRAN
list.files(localPath, recursive = TRUE, full.names = FALSE)

# Check for available packages
pkgAvail(repos = localPath, type = "souce")[, c(1:3, 5)]


.libPaths ("/opt/rlibs")
