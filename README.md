# CDGnet

Package repository for the `CDGnet` application. `CDGnet` is described in the following manuscript: https://ascopubs.org/doi/full/10.1200/CCI.19.00097



## Installation

Install from github using the `remotes` package:

```r
remotes::install_github('SiminaB/cdgnet-package', dependencies=TRUE)
```

## Preparation

`CDGnet` uses pathway information that needs to be download and processed before using the application. This only needs to be done once per installation.

```r
CDGnet::download_and_process_KEGG()
```

## Usage

`CDGnet` is a shiny app. To run the app we provide a helper function

```r
runCDGnet()
```

Once the app is running a molecular profile can be uploaded from a CSV file as directed in the app.
