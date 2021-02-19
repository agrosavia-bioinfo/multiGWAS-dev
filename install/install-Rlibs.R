#!/usr/bin/Rscript
#
# Create multiGWAS profile according to current directory
# 

MULTIGWAS_HOME = strsplit (getwd (), "/install")[[1]][1]
message ("MULTIGWAS_HOME=",MULTIGWAS_HOME)

message ("Creating multiGWAS profile...")
sink ("multiGWAS_profile.sh", append=F)
writeLines ("\n#------------------- multiGWAS.R profile ---------------------")
writeLines (paste0 ("export MULTIGWAS_HOME=", MULTIGWAS_HOME))
writeLines ("MULTIGWAS_TOOLS=$MULTIGWAS_HOME/opt/tools")
writeLines ("MULTIGWAS_MAIN=$MULTIGWAS_HOME/main")
writeLines ("export PATH=$PATH:$MULTIGWAS_TOOLS:$MULTIGWAS_MAIN")
sink()

# Write into .bashrc
profileFile = paste0 (path.expand ("~"), "/.bashrc")
sink (profileFile, append=T)
writeLines ("\n#------------------- multiGWAS.R tool profile ---------------------")
writeLines (paste0 (". ", MULTIGWAS_HOME, "/multiGWAS_profile.sh"))
sink ()

message ("\nMultiGWAS is ready to use, right after installed!\n")

# Install R libraries
libpath = paste0 (MULTIGWAS_HOME, "/opt/Rlibs")
if (!dir.exists (libpath))
	dir.create (libpath)

message (libpath)
message ("\n\nInstalling R libraries...\n\n")

.libPaths (libpath)

if (!require("pacman")) 
	install.packages('pacman', lib=libpath, repos='http://cran.us.r-project.org')

pacman::p_load("rrBLUP", "parallel","config","dplyr","stringi","qqman",
			   "VennDiagram","RColorBrewer","circlize","gplots", "rmarkdown",
			   "kableExtra" ,"doParallel", "ldsep", "yaml",
			   "BiocManager",
			   "multtest")

if(!'GWASpoly'%in% installed.packages()[,"Package"])
	install.packages(paste0(MULTIGWAS_HOME,'/opt/tools/GWASpoly_1.3.tar.gz'),
					 lib=libpath, repos=NULL, type="source") 
# Copy multiGWAS profile to multiGWAS home
x=file.copy ("multiGWAS_profile.sh", "..", overwrite=T)

message ("\n------------------------------------------\n")
message ("Close the terminal to finish the installation process")
message ("Then, open a new terminal and write: ")
message ("multiGWAS")
message ("\n------------------------------------------\n")




