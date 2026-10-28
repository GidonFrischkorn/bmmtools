# Run the setup call of a grid, or one slice of it.
#
#   Rscript grid-slice.R grid-args.R info      cores per slice, cells, dir
#   Rscript grid-slice.R grid-args.R setup     the setup call
#   Rscript grid-slice.R grid-args.R <k> <n>   slice k of n
#
# launch.sh calls it; it can also be run by hand, or from a scheduler's
# job script with k from the array index. launch.sh adds one argument
# after setup and after <k> <n>, a mark by which it finds the process
# with ps; it is not read here.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("usage: Rscript grid-slice.R <grid-args.R> info|setup|<k> <n>")
}
source(args[[1L]])
if (!exists("grid_args") || !is.list(grid_args)) {
  stop(args[[1L]], " must define a list called `grid_args`.")
}
reps <- if (is.null(grid_args$reps)) 1L else grid_args$reps

if (identical(args[[2L]], "info")) {
  # the cores one slice uses: its chains, times the threads of each chain
  chains <- if (is.null(grid_args$chains)) 4L else grid_args$chains
  threads <- grid_args$threads
  threads <- if (is.null(threads)) {
    1L
  } else if (is.list(threads)) {
    threads$threads
  } else {
    threads
  }
  # whole numbers in full: launch.sh compares them, and 3e+05 is no number.
  # The dir as R reads it, so that a `~` in it is the same folder to the
  # shell. On a line of its own that launch.sh finds by its marker, even
  # after a grid-args.R that printed without a newline.
  cat(
    "\nlaunch-info:",
    format(chains * threads, scientific = FALSE),
    format(nrow(grid_args$grid) * reps, scientific = FALSE),
    normalizePath(path.expand(grid_args$dir), mustWork = FALSE), "\n"
  )
  quit(save = "no", status = 0L)
}

cat("bmmtools", as.character(utils::packageVersion("bmmtools")), "\n")
print(as.data.frame(bmmtools::run_info()))

if (identical(args[[2L]], "setup")) {
  cells <- do.call(
    bmmtools::recovery_grid, c(grid_args, list(cells = integer(0)))
  )
  cat("setup done:", nrow(cells), "cells\n")
  quit(save = "no", status = 0L)
}

k <- as.integer(args[[2L]])
n <- as.integer(args[[3L]])
cells <- bmmtools::grid_tasks(grid_args$grid, reps, n)[[k]]
cat("slice", k, "of", n, ": cells", min(cells), "to", max(cells), "\n")
out <- do.call(bmmtools::recovery_grid, c(grid_args, list(cells = cells)))
status <- attr(out, "cells")$status
cat(
  "slice", k, "done:", length(status), "cells,",
  sum(status == "ok"), "ok,", sum(status == "error"), "failed\n"
)
# a slice with failed cells ends with EXIT 1, so --status shows it
if (any(status == "error")) {
  quit(save = "no", status = 1L)
}
