#!/usr/bin/env Rscript

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  install.packages(
    "rmarkdown",
    repos = "https://cloud.r-project.org"
  )
}

input <- "ANR_TF_gallery.Rmd"

if (!file.exists(input)) {
  stop(
    "Missing RMarkdown report: ",
    input
  )
}

rmarkdown::render(
  input = input,
  output_file = "ANR_TF_gallery.html",
  clean = TRUE,
  envir = new.env(
    parent = globalenv()
  )
)

message(
  "Gallery rendered: ",
  normalizePath(
    "ANR_TF_gallery.html"
  )
)
