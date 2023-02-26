#!/usr/bin/Rscript
#
# Install R libraries
MULTIGWAS_HOME = Sys.getenv ("MULTIGWAS_HOME")

libPath    = paste0  (MULTIGWAS_HOME, "/opt/Rlibs")
libRepo    = paste0  (MULTIGWAS_HOME, "/install/repo")
if (!dir.exists (libPath))
	dir.create (libPath)

.libPaths (libPath)
message ("\n\nInstalling R libraries for MultiGWAS into: ", libPath)

multigwas_packages = c("rrBLUP", "parallel","config","dplyr", "stringi","qqman", 
					   "VennDiagram", "RColorBrewer","circlize", "gplots", "rmarkdown", 
					   "kableExtra" ,"doParallel", "ldsep", "yaml", "BiocManager", 
					   "rlang", "ggplot2", "gtable", "scam", "tidyr")
packagesInstalled    = rownames(installed.packages())
packagesNotInstalled = setdiff(multigwas_packages, rownames(installed.packages()))  

message ("\nInstalled packages: ", paste (packagesInstalled, collapse=", "))
message ("\nPackages to install: ", paste (packagesNotInstalled, collapse=", "));Sys.sleep (3)
install.packages (packagesNotInstalled, repos=paste0("file://", libRepo))

#if (!"multtest" %in% packagesInstalled)
#	BiocManager::install("multtest", site_repository=libRepo)

if (!"GWASpoly" %in% packagesInstalled)
	install.packages(paste0(libRepo, "/src/contrib",'/GWASpoly23.tgz'), lib=libPath, repos=NULL, type="source", dependencies=T) 
#------------------------------------------------------------------------------

