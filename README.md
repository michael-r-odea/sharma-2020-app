# sharma-2020-app

This is a Shiny application for exploring the P0 and P5 mouse DRG scRNA-sequencing data from the Sharma et al., 2020 paper.

To use this app, you must have R (version ≥4.0.0) installed. If you do not, install at https://cran.r-project.org/bin/macosx/

To use this app on a Mac, follow these steps:
	
	1. Open the Terminal app (which can be found in the Applications folder). 
	2. Enter the following lines into the command line, one by one, hitting 'Enter' after each line:
	
```
R
install.packages('shiny')
install.packages('dplyr')
install.packages('ggplot2')
install.packages('Seurat')
install.packages('gridExtra')
install.packages('DT')
shiny::runUrl("https://github.com/michael-r-odea/sharma-2020-app/raw/main/app.zip")
```

The app should then appear in your web browser window shortly. 

