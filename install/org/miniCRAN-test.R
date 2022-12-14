
.libPaths ("/opt/miniCRAN/libs")


multiGWAS_packages = c("packman", "rrBLUP", "parallel","config","dplyr","stringi","qqman",
			   "VennDiagram","RColorBrewer","circlize","gplots", "rmarkdown",
			   "kableExtra" ,"doParallel", "ldsep", "yaml",
			   "BiocManager",
			   "multtest")


install.packages (multiGWAS_packages, repos="file:///opt/miniCRAN")
BiocManager::install("multtest", site_repository="file:///opt/miniCRAN")
