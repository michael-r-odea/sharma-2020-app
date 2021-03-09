# Sharma et al. 2020 P0 & P5 DRG scRNA-seq Data Visualization App

This is a Shiny application for exploring the P0 and P5 mouse DRG scRNA-sequencing data from [Sharma et al. 2020](https://www.nature.com/articles/s41586-019-1900-1).

## Usage
To use this app, you must have R (version ≥4.0.0) installed. If you do not, install from [CRAN](https://cran.r-project.org/bin/macosx/). If you're unsure if you have the latest version of R, run 'R' in the Terminal and check the version number. 

To use this app on a Mac/Linux, follow these steps: 

1. Open the Terminal app (which can be found in your Applications folder) and enter the following line: 
```
R
```

Next, install the required dependencies: (note: this step only needs to be performed once, the first time you use this app) 
```
install.packages(c('shiny', 'dplyr', 'ggplot2', 'Seurat', 'gridExtra', 'DT'), repos = "https://cloud.r-project.org")
```

Finally, launch the app using the following line:
```
shiny::runUrl("https://github.com/michael-r-odea/sharma-2020-app/raw/main/app.zip")
```

The app should then appear in your web browser window shortly. 


To use this app in future instances, simply run the following lines from the Terminal:
```
R
shiny::runUrl("https://github.com/michael-r-odea/sharma-2020-app/raw/main/app.zip")
```
