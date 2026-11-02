# config.R
# Setting environment variables for input/output filesystem roots
#   CARE_DATA_ROOT    extractionScripts pipeline outputs
#   CARE_PLOT_ROOT    figure panels
#   CARE_GENOME_ROOT  cellranger references and FASTAs for Figure2's polyA etc
#   CARE_PUBLIC_ROOT  downloaded public datasets (Collora et al, Frouard et al)
#

data.root <- Sys.getenv("CARE_DATA_ROOT",
                        "/projects/b1042/GoyalLab/egrody/extractedData/")
plot.root <- Sys.getenv("CARE_PLOT_ROOT",
                        "/projects/b1042/GoyalLab/egrody/plotData/")
genome.root <- Sys.getenv("CARE_GENOME_ROOT",
                          "/projects/b1042/GoyalLab/egrody/genomes/")
public.root <- Sys.getenv("CARE_PUBLIC_ROOT",
                          "/projects/b1042/GoyalLab/egrody/publicData/")

for (root in c("data.root", "plot.root", "genome.root", "public.root")) {
  value <- get(root)
  if (!grepl("/$", value)) {
    assign(root, paste0(value, "/"))
  }
}
rm(root, value)

# Where this file lives, so functions.R can reach the pipeline library beside it
# without depending on the working directory. source() records the path it was
# given in its own frame, so scanning the call stack finds it whether config.R
# was sourced as "config.R" from plotScripts or by absolute path from a job
# script; the working directory is only a fallback for an interactive paste.
care.scripts.dir <- local({
  for (i in seq_len(sys.nframe())) {
    o <- get0("ofile", envir = sys.frame(i), inherits = FALSE)
    if (!is.null(o) && is.character(o) && length(o) == 1 && nzchar(o))
      return(dirname(normalizePath(o)))
  }
  normalizePath(getwd())
})

