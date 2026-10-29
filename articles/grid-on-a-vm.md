# Running a grid on a cloud VM

A recovery grid fits one model per cell, and the cells do not depend on
each other. On a laptop they run one after another, each on as many
cores as it has chains. A cloud virtual machine with 32 or 64 cores can
run several cells at once. This article shows how to split one
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
call over several R processes on one such machine, with the launcher
that ships with bmmtools. Nothing here needs a batch scheduler. The same
R side works under one, and a sketch for SLURM is at the end.

The split changes where cells run, not what they are. Every cell keeps
the seed and the file names that one sequential run would give it. The
processes together therefore write the cells that a single call writes,
and
[`collect_grid()`](https://www.gfrischkorn.org/bmmtools/reference/collect_grid.md)
reads them as one grid.

## How a split run works

A split run has three steps, all in one grid directory:

1.  **The setup call** is your
    [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
    call with `cells = integer(0)`. It runs once. It writes `grid.rds`,
    runs the preflight fit, fits the priors of `prior_sd` if you asked
    for them, and simulates the data of every cell. It fits no cell.
2.  **The slices** are the same call, each with
    `cells = grid_tasks(grid, reps, n)[[k]]`, one per process.
    [`grid_tasks()`](https://www.gfrischkorn.org/bmmtools/reference/grid_tasks.md)
    cuts the cells into `n` runs of consecutive cells. Each slice first
    checks that its call matches the one the setup call recorded, then
    fits its own cells. The check covers the design, the seed, the
    prior, the model and formula, and the fitter’s arguments. A slice of
    another grid stops before it fits anything.
3.  **Collecting.** `collect_grid(dir)` reads every cell’s results from
    its files. `attr(x, "run_info")` has one row per call that reached
    its cells: which cells it ran, for how many it found results from an
    earlier run, how many failed, when it started and ended, and the
    machine and package versions. A call that stopped with an error
    before its first fit leaves no row.

``` r

args <- list(
  model = bmm::mixture2p(resp_error = "y"),
  grid = expand.grid(n_subjects = c(20, 50), n_trials = c(30, 100)),
  pars = c(kappa = log(8), thetat = qlogis(0.75)),
  sds = c(kappa = 0.3, thetat = 0.5),
  dir = "fits/mixture2p", reps = 10, seed = 1,
  chains = 4, cores = 4, iter = 2000, backend = "cmdstanr"
)

# once
do.call(recovery_grid, c(args, list(cells = integer(0))))

# in process k of 4
tasks <- grid_tasks(args$grid, args$reps, 4)
do.call(recovery_grid, c(args, list(cells = tasks[[k]])))

# afterwards, here or on another machine
out <- collect_grid("fits/mixture2p")
summary(out)
attr(out, "run_info")
```

Run each slice in its own R process. A parallel backend inside one
session can change the random number generator its workers use, and then
the fits differ from a sequential run. The launcher below starts one
`Rscript` per slice.

## The launcher

bmmtools ships three files for this, which you copy and edit:

``` r

file.copy(
  system.file("templates", "vm", package = "bmmtools"), ".",
  recursive = TRUE
)
```

- `vm/grid-args.R` is your grid, written once as a list called
  `grid_args`: everything you would pass to
  [`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md)
  except `cells`. Every process reads this file, so the setup call and
  the slices cannot drift apart.
- `vm/grid-slice.R` runs the setup call (`setup`) or one slice
  (`<k> <n>`). You can run it by hand.
- `vm/launch.sh` runs the setup call and then starts the slices. It also
  reports on them and restarts one.

From the directory that `dir` in `grid-args.R` is relative to:

``` bash
bash vm/launch.sh vm/grid-args.R
```

The launcher reads the machine’s cores and the cores one slice uses: its
`chains`, times the `threads` of each chain. It starts as many slices as
fit, so 32 cores and 4 chains give 8 slices. `-k` sets the number
yourself. The launcher refuses a number that would put more chains than
cores on the machine unless you add `--force`.

The setup call runs in the foreground and logs to
`<dir>/logs/setup.log`. If it fails, nothing else starts, and stopping
the launcher stops the setup call too. It also stops when the SSH
session ends, and simulating a large grid can take a while, so start the
launcher inside `tmux` or `screen`. The slices run in the background
under `nohup`, so they go on after the SSH session ends. Each one writes
`<dir>/logs/slice-<k>.log`, which ends with a line `EXIT <status>` when
the process finishes: `EXIT 0` when every cell has a fit, `EXIT 1` when
some failed. To see how far they are:

``` bash
bash vm/launch.sh --status vm/grid-args.R
```

`--status` prints the state of the setup call, the number of cells with
a result and, for each slice, its exit line, `running` with its process
id, or `stopped without an exit line` for a process that is gone without
one (the machine restarted, or something killed it hard). A recorded
process id that another process has taken since, as after a restart or
by a slice of another grid, does not count as running. This needs a
working `ps`, as does stopping every process a slice’s R started
(below). Without one the launcher says so: it then counts a slice as
running while a process group with its recorded id is live, and the
setup call while a process with its recorded id is. A `kill` of a slice
still stops its R and the processes R started in the slice’s process
group. Only processes in a session of their own survive it, such as the
Stan chains that cmdstanr starts through processx.

The launcher runs one launch at a time. It refuses to start while
another launch sets up or starts the same grid, and while a slice of its
last launch still runs, so two launches do not fit the same cells at
once. It cannot see slices that you start by hand.

To stop a slice, send its process id from `<dir>/logs/pids` a plain
`kill`. The slice’s R process and the processes R started, such as the
Stan chains, stop with it, and the log gets its exit line. A `kill -9`
stops only the process you name and leaves no exit line. Its R fits on,
so `--status` still shows the slice as running, with a note that the
process is gone and the `kill` command, with the ids of its R and of
what R started, that stops the rest. Without a working `ps` the note
gives `kill -- -<id>` instead, which stops the processes of the slice’s
process group, its R among them, but not those in a session of their
own. A `kill -9` on the launcher during the setup call leaves the setup
running: a new launch is refused until it ends, and no slices start from
it. `--status` then shows the `kill` command that stops the setup and
what it started. Without a working `ps` it names only the setup’s own
process, and what the setup started in a session of its own outlives
that `kill`.

If the slices die when you log out, the machine ends a user’s processes
at logout (`KillUserProcesses=yes` in systemd’s `logind.conf`), which
`nohup` does not prevent; check that setting first.

## When a slice stops

A cell’s results are written the moment it finishes, each file under a
temporary name first and then renamed. A slice that is killed leaves the
cells it finished and nothing half-written. Start it again:

``` bash
bash vm/launch.sh --rerun 3 vm/grid-args.R
```

This runs slice 3 of the last launch in the background, as the launcher
first did, and does not repeat the setup call. The slice skips the cells
that have results and fits the rest. A refitted cell draws the same
chains it would have drawn the first time: a split run always has a
`seed`.

## Bringing the results home

The fits are large. The results are small.
[`collect_grid()`](https://www.gfrischkorn.org/bmmtools/reference/collect_grid.md)
needs `grid.rds`, `prior-sd.rds` (if the grid has one), the `run-info/`
folder, and each cell’s `-sim.rds` and `-est.rds` files. It does not
need the fits:

``` bash
rsync -av --include='*/' --include='grid.rds' --include='prior-sd.rds' \
  --include='run-info/*' --include='*-sim.rds' --include='*-est.rds' \
  --exclude='*' vm:project/fits/mixture2p/ fits/mixture2p/
```

## What to set up on the machine

The launcher installs nothing. The machine needs R, CmdStan with a C++
toolchain, and the versions of bmm, brms, cmdstanr and bmmtools you mean
to report. Each process compiles the model once at its start.
[`run_info()`](https://www.gfrischkorn.org/bmmtools/reference/run_info.md)
on the machine shows what is installed, and each slice’s record keeps
it.

Leave `cmdstanr_write_stan_file_dir` unset, or give each process its own
directory. Processes that start together and share one directory would
compile the same executable at the same time.

## Several machines

Run the setup call on one machine and copy the grid directory to the
others. Install the same versions of bmmtools, bmm and brms on every
machine: `grid.rds` records these three, and a slice on a machine with
other versions of them stops. Keep R and CmdStan the same too. The check
does not cover them, but each slice’s run record shows what it ran on.
Then give each machine its share of the slices, for example slices 1–8
of 16 on the first and 9–16 on the second. The launcher always runs all
of its slices, so here you start them by hand, one `Rscript` per slice,
with the same total on every machine:

``` bash
nohup Rscript vm/grid-slice.R vm/grid-args.R 9 16 > slice-9.log 2>&1 &
```

The slices write disjoint files as long as every machine splits the grid
into the same 16. Copying the directories back into one gives the whole
grid.

## Under a batch scheduler

On a cluster with a scheduler, each array task runs one slice. A SLURM
sketch:

``` bash
#!/bin/bash
#SBATCH --array=1-16
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
Rscript vm/grid-slice.R vm/grid-args.R "$SLURM_ARRAY_TASK_ID" 16
echo "EXIT $?"
```

Run the setup call before you submit the array, on the login node or as
a job of its own that the array waits for:

``` bash
Rscript vm/grid-slice.R vm/grid-args.R setup
```

The run records fill in `scheduler`, `job_id`, `array_task_id` and
`cores_granted` from the scheduler’s environment (SLURM, PBS, LSF or
SGE).
