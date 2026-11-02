# packages.R
# Installing and loading all required packages with version tracking


# Packages are split into a core set attached at every startup and optional
# groups that are installed and version-tracked but left unattached until a
# script asks for them with load_packages().
#
# Attaching every package costs most of a session's memory in the Bioconductor
# genomics stack, which no Paper figure uses, so that stack is left optional.
#
# Optional packages are still installed at startup, so a fresh system ends up
# with an identical environment either way -- only the attach step is deferred.

# Progress chatter is off by default, so sourcing this file from a notebook
# prints nothing on a run where everything is already installed -- and nothing
# but failures on a run that installs. Set options(care.packages.verbose = TRUE)
# or CARE_PACKAGES_VERBOSE=1 to see install progress and the installer's own
# output. Failures are reported either way.
packages_verbose <- function() {
  isTRUE(getOption("care.packages.verbose",
                   nzchar(Sys.getenv("CARE_PACKAGES_VERBOSE"))))
}
say <- function(...) if (packages_verbose()) cat(...)

# Runs an installer with its output and messages captured, and replays them
# only if it fails -- where they are the only record of why. Under
# packages_verbose() it runs with output attached instead.
run_quiet <- function(code, what) {
  if (packages_verbose()) return(invisible(isTRUE(tryCatch({ code(); TRUE },
                                                   error = function(e) FALSE))))
  log_file <- tempfile()
  con <- file(log_file, open = "wt")
  sink(con)
  sink(con, type = "message")
  ok <- tryCatch({ code(); TRUE },
                 error = function(e) { message("error: ", conditionMessage(e)); FALSE })
  sink(type = "message")
  sink()
  close(con)
  if (!ok) {
    cat("ERROR: ", what, " failed:\n", sep = "")
    captured <- tryCatch(readLines(log_file, warn = FALSE), error = function(e) character())
    if (length(captured)) cat(paste0("  ", captured, collapse = "\n"), "\n", sep = "")
  }
  unlink(log_file)
  invisible(ok)
}

core_packages <- c(
  "data.table",
  "dplyr",
  "ggforce",
  "ggnewscale",
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "ggsignif",
  "ggtext",
  "ggVennDiagram",
  "grid",
  "gridExtra",
  "MASS",
  "matrixStats",
  "parallel",
  "patchwork",
  "pheatmap",
  "purrr",
  "reshape2",
  "rlang",
  "scales",
  "Seurat",
  "stringr",
  "tibble",
  "tidyr"
)

# Groups are a convenience for load_packages(); it also accepts bare package
# names, so the grouping does not have to be exhaustive or exclusive.
optional_packages <- list(
  genomics     = c("Biostrings", "Rsamtools", "rtracklayer", "GenomicAlignments",
                   "rhdf5", "GenomicRanges", "IRanges"),
  sce          = c("SingleCellExperiment", "zellkonverter", "glmGamPoi"),
  # miloR: KNN-neighbourhood differential abundance, used for the cluster-free
  # splice+ enrichment test. edgeR/statmod/BiocNeighbors come with it.
  neighborhood = c("miloR", "edgeR", "statmod", "BiocNeighbors"),
  hashing      = c("demuxmix"),
  # vegan/cluster: PERMANOVA (adonis2) and silhouette width for the host-state
  # dispersion test. vegan compiles, so it is deferred rather than attached.
  dispersion   = c("vegan", "cluster"),
  # presto: fast Wilcoxon rank-sum. Seurat's FindMarkers/FindAllMarkers pick it
  # up automatically when it is installed, so it needs no call of its own.
  # leidenbase: the backend Seurat 5.5 runs FindClusters(algorithm = 4) through.
  # It is the default leiden_method, so without it every Leiden call errors.
  seurat_extra = c("scCustomize", "harmony", "clustree", "presto", "leidenbase"),
  enrichment   = c("fgsea", "msigdbr"),
  lineage      = c("GEMLI", "igraph"),
  utils        = c("openxlsx", "transport", "philentropy", "stringdist",
                   "doParallel", "foreach", "RANN", "R.utils", "svglite",
                   "reticulate", "UpSetR")
)

# Bioconductor packages are not on CRAN and need BiocManager to resolve.
bioc_packages <- c("fgsea", "glmGamPoi", "zellkonverter", "SingleCellExperiment",
                   "Biostrings", "Rsamtools", "rtracklayer", "GenomicAlignments",
                   "rhdf5", "GenomicRanges", "IRanges", "demuxmix",
                   "miloR", "edgeR", "statmod", "BiocNeighbors")

# GEMLI and presto install from GitHub via install_gemli() / install_presto()
# below, not from CRAN.
optional_cran_packages <- setdiff(unlist(optional_packages, use.names = FALSE),
                                  c(bioc_packages, "GEMLI", "presto"))

# packages_versions.csv sits next to this script, so its path is resolved from
# the script's own location rather than from the caller's working directory.
# source() records the sourced file in the calling frame's `ofile`; outside a
# source() call (Rscript -e, console paste) that is absent and the working
# directory is the only thing left to fall back on.
script_dir <- local({
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  getwd()
})

versions_file <- file.path(script_dir, "packages_versions.csv")

# Clean up stale 00LOCK-* directories left behind by previous failed
# installs. These would otherwise cause "failed to lock directory" errors
# the next time install.packages tries the same package.
local({
  for (lib in .libPaths()) {
    locks <- list.files(lib, pattern = "^00LOCK", full.names = TRUE)
    if (length(locks) > 0) {
      say("Removing", length(locks), "stale install lock(s) from", lib, ":",
          paste(basename(locks), collapse = ", "), "\n")
      unlink(locks, recursive = TRUE, force = TRUE)
    }
  }
})

# Disk-based "is package installed?" check. requireNamespace() returns FALSE
# whenever a package's own dependencies fail to load — so if dplyr is missing,
# Seurat (which depends on dplyr) appears "missing" too, even though it's
# fully installed. Using installed.packages() looks at DESCRIPTION files on
# disk and avoids that cascade.
is_pkg_installed <- function(pkgs) {
  pkgs %in% rownames(installed.packages())
}

# Base-priority packages (grid, parallel, stats, ...) ship inside R itself and
# are not distributed on CRAN, so their version is R's version and no install
# can ever change it. They are excluded from version pinning and from version
# tracking; pinning one produces a warning on every startup and an install that
# cannot succeed.
is_base_pkg <- function(pkgs) {
  ip <- installed.packages()
  priority <- ip[match(pkgs, rownames(ip)), "Priority"]
  !is.na(priority) & priority == "base"
}

# Subprocess installer: launches the SAME R binary the parent is running and
# inherits the parent's .libPaths() so installs land where the parent will
# find them. Uses install.packages with an explicit CRAN mirror; the
# subprocess starts with no packages loaded, so dependency upgrades that
# would error with "Updating loaded packages" in the parent will succeed.
# CRAN-only by default; Bioc-dependent packages should be passed with
# use_bioc = TRUE so BiocManager resolves the Bioc deps. Output is streamed
# to stdout so install errors are visible.
install_in_subprocess <- function(pkgs, label = "packages",
                                  use_bioc = FALSE) {
  if (length(pkgs) == 0) return(invisible(0))
  say("Installing ", length(pkgs), " missing ", label, " in subprocess: ",
      paste(pkgs, collapse = ", "), "\n", sep = "")
  rscript <- file.path(R.home("bin"), "Rscript")
  pkg_vec <- paste0('"', pkgs, '"', collapse = ", ")
  lib_vec <- paste0('"', .libPaths(), '"', collapse = ", ")
  body <- if (use_bioc) {
    sprintf(
      'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE); BiocManager::install(c(%s), update = FALSE, ask = FALSE, dependencies = TRUE, quiet = TRUE)',
      pkg_vec
    )
  } else {
    sprintf(
      'install.packages(c(%s), dependencies = TRUE, repos = "https://cloud.r-project.org", quiet = TRUE)',
      pkg_vec
    )
  }
  install_expr <- sprintf('.libPaths(c(%s)); %s', lib_vec, body)

  # Streaming the installer's output would put every compilation line in the
  # caller's notebook, so it is captured and replayed only when the install
  # fails -- where it is the only record of why.
  if (packages_verbose()) {
    say("Rscript -e ", shQuote(install_expr), "\n", sep = "")
    exit_code <- system2(rscript, args = c("-e", shQuote(install_expr)),
                         stdout = "", stderr = "")
  } else {
    output <- suppressWarnings(
      system2(rscript, args = c("-e", shQuote(install_expr)),
              stdout = TRUE, stderr = TRUE))
    exit_code <- attr(output, "status")
    exit_code <- if (is.null(exit_code)) 0L else exit_code
    if (exit_code != 0) {
      cat("ERROR: installing ", paste(pkgs, collapse = ", "),
          " failed (exit ", exit_code, "):\n", sep = "")
      cat(paste0("  ", output, collapse = "\n"), "\n", sep = "")
    }
  }
  invisible(exit_code)
}

# attach = FALSE installs and version-tracks without calling library(), which is
# how optional groups are handled: a fresh system still gets every package, but
# the session only pays the memory cost of what it actually attaches.
load_packages_with_versions <- function(packages, versions_file,
                                        update_versions_file = TRUE,
                                        attach = TRUE) {

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", quiet = TRUE)
  }

  # Load existing versions file
  if (file.exists(versions_file)) {
    versions_df <- read.csv(versions_file, stringsAsFactors = FALSE)
  } else {
    versions_df <- data.frame(package = character(), version = character())
  }

  # ---- Phase 0: Install any missing packages (latest) BEFORE loading -------
  # All installs happen up front so that no library() call later in this
  # function can fail because a transitive dep isn't yet on disk (e.g.
  # library(clustree) -> ggraph -> dplyr).
  #
  # Run install.packages in an Rscript subprocess. Calling it directly from
  # the parent R session fails with "Updating loaded packages" whenever a
  # transitive dependency (Matrix, rlang, etc.) is already loaded — the
  # subprocess starts with a clean namespace, so upgrades are unblocked.

  missing_pkgs <- packages[!is_pkg_installed(packages)]
  if (length(missing_pkgs) > 0) {
    say("Installing", length(missing_pkgs), "missing CRAN package(s):",
        paste(missing_pkgs, collapse = ", "), "\n")

    # Sanity-check CRAN reachability up front. If we can't even fetch the
    # index, install.packages will silently produce zero installs; better
    # to fail loudly here so the user knows it's an internet/module issue
    # not a packages.R bug.
    cran_ok <- tryCatch({
      idx <- available.packages(repos = "https://cloud.r-project.org")
      say("CRAN index reached: ", nrow(idx), " packages available.\n", sep = "")
      TRUE
    }, error = function(e) {
      cat("ERROR: cannot reach CRAN: ", conditionMessage(e), "\n",
          "       Run `module load openssl/1.0.2k-fips-gcc-12.3.0` from your\n",
          "       shell BEFORE launching R, then re-source this file.\n",
          sep = "")
      FALSE
    }, warning = function(w) {
      cat("WARNING reaching CRAN: ", conditionMessage(w), "\n",
          "       If this mentions libssl.so.10, run\n",
          "       `module load openssl/1.0.2k-fips-gcc-12.3.0` BEFORE R,\n",
          "       then re-source this file.\n", sep = "")
      FALSE
    })

    if (cran_ok) {
      # dependencies = NA (default) installs Depends/Imports/LinkingTo only —
      # avoids "Updating loaded packages" errors that fire when TRUE pulls
      # in Suggests touching already-loaded base packages. quiet = FALSE
      # so install errors are visible.
      run_quiet(function() {
        install.packages(missing_pkgs, dependencies = NA,
                         repos = "https://cloud.r-project.org", quiet = TRUE)
      }, paste("install of", paste(missing_pkgs, collapse = ", ")))
    }
  }

  # Re-check after install — anything still missing on disk should warn
  # loudly so the user can act, but we continue so as much as possible loads.
  still_missing <- missing_pkgs[!is_pkg_installed(missing_pkgs)]
  if (length(still_missing) > 0) {
    cat("WARNING: still missing after install attempt: ",
        paste(still_missing, collapse = ", "), "\n",
        "If you see 'Updating loaded packages' errors above, ",
        "restart R and re-source. If you see 'internet routines cannot be ",
        "loaded' or 'cannot reach CRAN', load the network module before R ",
        "(e.g. `module load curl openssl`). If you see \"cannot find -l<name>\" ",
        "the package compiles against a system library the linker cannot see; ",
        "see the ggtext note at the end of this file.\n", sep = "")
  }

  # ---- Phase 1: Install version-mismatched packages in a subprocess --------
  # A subprocess has a clean namespace environment, so unloadNamespace() will
  # never fail due to reverse-dependency imports from the main session.

  pkgs_to_install <- character()
  vers_to_install <- character()

  for (pkg in packages) {
    if (is_base_pkg(pkg)) next
    target_version <- versions_df$version[versions_df$package == pkg]
    if (length(target_version) == 0 || is.na(target_version)) next

    installed_version <- tryCatch(
      as.character(packageVersion(pkg)),
      error = function(e) NA_character_
    )

    if (is.na(installed_version) || installed_version != target_version) {
      # Only ever install upwards. An installed version NEWER than the pin is
      # reported, not replaced: something upgraded it on purpose (a dependency
      # resolution, a manual install), and silently downgrading it every
      # session is a loop, not reproducibility. Re-pin with
      #   write.csv(data.frame(package = <pkg>, version = <installed>), ...)
      # or delete the row to accept the newer version.
      if (!is.na(installed_version) &&
          utils::compareVersion(installed_version, target_version) > 0) {
        say("NOTE: ", pkg, " ", installed_version, " is newer than the pin (",
            target_version, "); keeping the installed version.\n", sep = "")
        next
      }
      # A package whose namespace is already loaded cannot be swapped under a
      # running session: the install would succeed on disk and the following
      # library() would fail with "version X cannot be unloaded" because other
      # loaded packages import it. Restarting R is the only way to apply it.
      if (pkg %in% loadedNamespaces()) {
        say("NOTE: ", pkg, " is already loaded at ", installed_version,
            "; restart R to install the pinned ", target_version, ".\n", sep = "")
        next
      }
      pkgs_to_install <- c(pkgs_to_install, pkg)
      vers_to_install <- c(vers_to_install, target_version)
    }
  }

  if (length(pkgs_to_install) > 0) {
    say("Installing", length(pkgs_to_install),
        "package(s) with version mismatches...\n")

    # Each package is installed in its OWN Rscript process. This prevents
    # namespace contamination: installing CoGAPS in subprocess A may load
    # fgsea/data.table as transitive dependencies during dependency resolution.
    # If data.table were then installed in the same subprocess, the loaded
    # fgsea namespace would block unloadNamespace("data.table"). Separate
    # processes guarantee a clean namespace environment per install.
    #
    # dependencies = FALSE because these packages are already installed (just
    # at the wrong version). Their dependency trees are already satisfied on
    # disk. Using TRUE would trigger dependency resolution that can load
    # namespaces and recreate the original conflict.

    for (i in seq_along(pkgs_to_install)) {
      pkg <- pkgs_to_install[i]
      ver <- vers_to_install[i]
      say("  [", i, "/", length(pkgs_to_install), "] ", pkg, " ", ver, " ...\n", sep = "")

      install_expr <- sprintf(
        'options(repos = c(CRAN = "https://cloud.r-project.org")); if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", quiet = TRUE); remotes::install_version("%s", version = "%s", dependencies = FALSE, quiet = TRUE, upgrade = "never")',
        pkg, ver
      )

      if (packages_verbose()) {
        exit_code <- system2("Rscript", args = c("-e", shQuote(install_expr)),
                             stdout = "", stderr = "")
      } else {
        output <- suppressWarnings(
          system2("Rscript", args = c("-e", shQuote(install_expr)),
                  stdout = TRUE, stderr = TRUE))
        exit_code <- attr(output, "status")
        exit_code <- if (is.null(exit_code)) 0L else exit_code
      }
      if (exit_code != 0) {
        cat("    WARNING: install of", pkg, ver, "exited with code", exit_code, "\n")
      }
    }
  }

  # ---- Phase 2: Load all packages and track new ones -----------------------
  # By this point every package in `packages` should be installed (Phase 0
  # handled missing, Phase 1 handled version mismatches). library() here
  # only loads; no installs happen inside the loop.

  new_packages <- character()
  new_versions <- character()

  for (pkg in packages) {
    installed_version <- tryCatch(
      as.character(packageVersion(pkg)),
      error = function(e) NA_character_
    )

    if (is.na(installed_version)) {
      cat("WARNING:", pkg, "still missing after install phase; skipping load\n")
      next
    }

    if (attach) {
      load_ok <- tryCatch({
        # quietly/warn.conflicts do not propagate to packages attached through
        # Depends, so the Bioconductor stack still announces itself and its
        # masks without these wrappers.
        suppressPackageStartupMessages(suppressMessages(
          library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
        ))
        TRUE
      }, error = function(e) {
        cat("WARNING: failed to load ", pkg, ": ", conditionMessage(e), "\n",
            sep = "")
        FALSE
      })
      if (!load_ok) next
    }

    if (is_base_pkg(pkg)) next

    target_version <- versions_df$version[versions_df$package == pkg]
    if (length(target_version) == 0 || is.na(target_version)) {
      # Newly installed (no pin) — record in versions file
      new_packages <- c(new_packages, pkg)
      new_versions <- c(new_versions, installed_version)
    } else if (target_version != installed_version) {
      cat("WARNING:", pkg, "loaded as", installed_version,
          "but pin requests", target_version, "\n")
    }
  }

  # ---- Update versions file with newly tracked packages --------------------

  if (update_versions_file && length(new_packages) > 0) {
    new_entries <- data.frame(package = new_packages, version = new_versions)
    versions_df <- rbind(versions_df, new_entries)
    versions_df <- versions_df[order(versions_df$package), ]
    write.csv(versions_df, versions_file, row.names = FALSE)
  }

  invisible(versions_df)
}

# Core CRAN packages: installed, version-tracked and attached.
load_packages_with_versions(core_packages, versions_file,
                            update_versions_file = TRUE, attach = TRUE)

# Optional CRAN packages: installed and version-tracked, but not attached.
load_packages_with_versions(optional_cran_packages, versions_file,
                            update_versions_file = TRUE, attach = FALSE)

# Bioconductor packages are installed after the CRAN ones so that any CRAN
# dependencies they need (e.g. dplyr for fgsea helpers) are already on disk.
# None are attached here; every one of them belongs to an optional group.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", quiet = TRUE)
}
missing_bioc <- bioc_packages[!is_pkg_installed(bioc_packages)]
install_in_subprocess(missing_bioc, label = "Bioconductor package(s)",
                      use_bioc = TRUE)

# ---- GitHub-only packages ------------------------------------------------
# GEMLI (github.com/UPSUTER/GEMLI) is not on CRAN or Bioconductor and lives in
# the GEMLI_package_v0 subdirectory of its repo.
#
# Its DESCRIPTION lists HiClimR in Imports, which pulls in ncdf4 and therefore
# the netCDF C library. That toolchain is not available on Quest, so a plain
# install_github fails. HiClimR supplies exactly one symbol to GEMLI —
# fastCor(), called only in the fast = TRUE branch of calculate_correlations().
# quantify_clusters_iterative() hardcodes fast = FALSE, so no GEMLI code path
# used for lineage prediction ever reaches it. The install therefore clones the
# repo, drops HiClimR from Imports, and installs from the local copy. Calling
# calculate_correlations(fast = TRUE) directly will error under this install.

install_gemli <- function() {
  if (is_pkg_installed("GEMLI")) return(invisible(TRUE))
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  }

  clone_dir <- file.path(tempdir(), "GEMLI_repo")
  unlink(clone_dir, recursive = TRUE, force = TRUE)

  clone_ok <- tryCatch(
    system2("git", args = c("clone", "--depth", "1",
                            "https://github.com/UPSUTER/GEMLI.git", clone_dir),
            stdout = FALSE, stderr = FALSE) == 0,
    error = function(e) FALSE
  )
  if (!clone_ok) {
    cat("WARNING: could not clone GEMLI from GitHub; skipping install.\n")
    return(invisible(FALSE))
  }

  pkg_dir  <- file.path(clone_dir, "GEMLI_package_v0")
  desc     <- file.path(pkg_dir, "DESCRIPTION")
  desc_txt <- readLines(desc)
  desc_txt <- desc_txt[!grepl("HiClimR", desc_txt)]
  writeLines(desc_txt, desc)

  run_quiet(function() {
    remotes::install_local(pkg_dir, upgrade = "never", quiet = TRUE)
  }, "GEMLI install")
  invisible(is_pkg_installed("GEMLI"))
}

# Installed but not attached; it belongs to the "lineage" optional group.
install_gemli()

# presto (github.com/immunogenomics/presto) is not on CRAN or Bioconductor.
# Unlike GEMLI it installs straight from the repo root with no DESCRIPTION
# surgery, so install_github is enough. remotes provides install_github and is
# already a dependency here, so this does not pull in devtools as well.
install_presto <- function() {
  if (is_pkg_installed("presto")) return(invisible(TRUE))
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  }

  run_quiet(function() {
    remotes::install_github("immunogenomics/presto", upgrade = "never",
                            quiet = TRUE)
  }, "presto install")
  invisible(is_pkg_installed("presto"))
}

# Installed but not attached; it belongs to the "seurat_extra" optional group.
install_presto()

# ---- On-demand attaching --------------------------------------------------
# Attaches an optional group (or bare package names) in the script that needs
# it, e.g. load_packages("lineage") before a GEMLI chunk. Anything already
# attached is skipped, so repeated calls are free.

load_packages <- function(...) {
  requested <- unlist(list(...), use.names = FALSE)
  pkgs <- unlist(lapply(requested, function(x) {
    if (x %in% names(optional_packages)) optional_packages[[x]] else x
  }), use.names = FALSE)
  pkgs <- unique(pkgs)

  missing_pkgs <- pkgs[!is_pkg_installed(pkgs)]
  if (length(missing_pkgs) > 0) {
    if ("GEMLI" %in% missing_pkgs) {
      install_gemli()
      missing_pkgs <- setdiff(missing_pkgs, "GEMLI")
    }
    if ("presto" %in% missing_pkgs) {
      install_presto()
      missing_pkgs <- setdiff(missing_pkgs, "presto")
    }
    install_in_subprocess(intersect(missing_pkgs, bioc_packages),
                          label = "Bioconductor package(s)", use_bioc = TRUE)
    install_in_subprocess(setdiff(missing_pkgs, bioc_packages),
                          label = "CRAN package(s)")
  }

  for (pkg in pkgs) {
    if (paste0("package:", pkg) %in% search()) next
    tryCatch(
      suppressPackageStartupMessages(suppressMessages(
        library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
      )),
      error = function(e) {
        cat("WARNING: failed to load ", pkg, ": ", conditionMessage(e), "\n",
            sep = "")
      }
    )
  }
  invisible(pkgs)
}

# Resolve namespace conflicts — ensure dplyr wins over MASS and plyr
if ("dplyr" %in% loadedNamespaces()) {
  conflict_fns <- c("select", "filter", "count", "rename", "mutate", "arrange", "summarise", "summarize")
  for (fn in conflict_fns) {
    if (exists(fn, envir = asNamespace("dplyr")))
      assign(fn, get(fn, envir = asNamespace("dplyr")), envir = .GlobalEnv)
  }
}

# ---- Leiden clustering backend (Seurat FindClusters algorithm = 4) ----
# Seurat 5.5 runs Leiden through the R package leidenbase, which is its default
# `leiden_method`. It is declared in optional_packages$seurat_extra above;
# install it with install.packages("leidenbase") if a run reports
#   Package 'leidenbase' is required for leiden_method = 'leidenbase'
#
# The package is installed per R version, so a library tree for one R minor
# version having it says nothing about another. Every pipeline script loads
# R/4.4.0; run the notebooks under the same R, or the clustering backend and
# the package versions behind every cluster label differ from the ones that
# produced the cached objects.
#
# Older Seurat reached leidenalg in Python through reticulate instead. That
# path is unused here, and installing leidenalg does not satisfy Seurat 5.5.

# ---- Compiled dependencies that need a system library at build time ----
# ggtext pulls in gridtext, which pulls in jpeg, which links against libjpeg.
# /usr/lib64 carries only the runtime libjpeg.so.62, with no libjpeg.so for the
# linker to resolve -ljpeg against, so the build fails with
#   ld: cannot find -ljpeg
# even though CRAN is reachable and every other package installs. Loading the
# libjpeg-turbo module is not enough on its own: it leaves LIBRARY_PATH pointing
# at the previous spack tree. Give the linker and compiler the path outright,
# matching the compiler R itself was built with (R 4.5.1 is gcc-12.4.0):
#
#   JPEG=/gpfs/software/2025/spack/opt/spack/linux-rhel8-x86_64/gcc-12.4.0/\
#   libjpeg-turbo-3.0.3-s2cnlywiuswtzefilmwyx7jwwrml67cq
#   export LIBRARY_PATH="$JPEG/lib64:$LIBRARY_PATH"
#   export CPATH="$JPEG/include:$CPATH"
#   Rscript -e 'install.packages(c("jpeg","gridtext","ggtext"))'
#
# Build time only. Once installed, ggtext loads with none of these set.
