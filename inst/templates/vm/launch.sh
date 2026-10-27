#!/usr/bin/env bash
# Run a recovery grid on one machine as K processes, each a slice.
#
#   bash launch.sh grid-args.R              the setup call, then K slices
#   bash launch.sh -k 6 grid-args.R         K = 6 slices
#   bash launch.sh --force -k 20 grid-args.R
#   bash launch.sh --status grid-args.R     how far the slices are
#   bash launch.sh --rerun 3 grid-args.R    slice 3 of the last launch again
#
# K defaults to the cores of the machine divided by the cores one slice
# uses (chains x threads), and is refused above that without --force.
# The setup call runs in the foreground; the slices run in the
# background under nohup, so they go on when the SSH session ends. Each
# writes <dir>/logs/slice-<k>.log, which ends with "EXIT <status>".
# A launch is refused while a slice of the last one still runs; to stop
# a slice, kill the pid that <dir>/logs/pids lists for it.
# Run it from the directory that grid-args.R's `dir` is relative to.

set -u
here="$(cd "$(dirname "$0")" && pwd)"
slice="$here/grid-slice.R"

if [ "${1:-}" = "--run-slice" ]; then
  # one background slice: launch.sh calls itself with this under nohup.
  # R runs as a child, so a TERM sent here reaches it and the log still
  # gets its exit line.
  Rscript "$slice" "$2" "$3" "$4" > "$5" 2>&1 &
  child=$!
  trap 'kill -TERM "$child" 2>/dev/null' TERM INT
  wait "$child"
  code=$?
  if kill -0 "$child" 2>/dev/null; then
    wait "$child"
    code=$?
  fi
  echo "EXIT $code" >> "$5"
  exit 0
fi

whole() {
  case "${1:-}" in
    '' | *[!0-9]*) return 1 ;;
  esac
  [ "$1" -ge 1 ] 2>/dev/null
}

k=""
rerun=""
force=0
status=0
args_file=""
while [ $# -gt 0 ]; do
  case "$1" in
    -k | --rerun)
      if [ $# -lt 2 ] || ! whole "$2"; then
        echo "launch.sh: $1 needs a whole number of 1 or more" >&2
        exit 2
      fi
      if [ "$1" = "-k" ]; then k="$2"; else rerun="$2"; fi
      shift 2
      ;;
    --force) force=1; shift ;;
    --status) status=1; shift ;;
    -h | --help) sed -n '2,18p' "$0"; exit 0 ;;
    *) args_file="$1"; shift ;;
  esac
done
if [ -z "$args_file" ] || [ ! -f "$args_file" ]; then
  echo "launch.sh: give the grid-args.R file (see --help)" >&2
  exit 2
fi

# the last line only: grid-args.R may print
info=$(Rscript "$slice" "$args_file" info | tail -n 1)
read -r per_slice n_cells dir <<< "$info"
if ! whole "${per_slice:-}" || ! whole "${n_cells:-}" || [ -z "${dir:-}" ]; then
  echo "launch.sh: could not read $args_file" >&2
  exit 2
fi
logs="$dir/logs"

alive() {
  [ -n "$1" ] && kill -0 "$1" 2>/dev/null
}
# the pid last recorded for slice $1
pid_of() {
  [ -f "$logs/pids" ] || return 0
  awk -v k="$1" '$1 == k { p = $2 } END { print p }' "$logs/pids"
}
start_slice() {
  # emptied here, so that the old exit line is gone before the slice starts
  : > "$logs/slice-$1.log"
  nohup bash "$here/launch.sh" --run-slice "$args_file" "$1" "$2" \
    "$logs/slice-$1.log" > /dev/null 2>&1 &
  echo "$1 $!" >> "$logs/pids"
}
n_last=$(cat "$logs/slices" 2>/dev/null || true)

if [ "$status" = 1 ]; then
  if [ ! -d "$dir" ]; then
    echo "no grid in $dir yet"
    exit 0
  fi
  done_cells=$(find "$dir" -maxdepth 1 -name 'cell-*-est.rds' | wc -l |
    tr -d ' ')
  echo "cells with a result: $done_cells of $n_cells"
  i=1
  while [ "$i" -le "${n_last:-0}" ]; do
    end=$(grep -E '^EXIT [0-9]+$' "$logs/slice-$i.log" 2>/dev/null | tail -n 1)
    pid=$(pid_of "$i")
    if [ -n "$end" ]; then
      state=$end
    elif alive "$pid"; then
      state="running (pid $pid)"
    else
      state="stopped without an exit line"
    fi
    echo "slice-$i.log: $state"
    i=$((i + 1))
  done
  exit 0
fi

if [ -n "$rerun" ]; then
  if [ -z "$n_last" ]; then
    echo "launch.sh: no earlier launch in $logs to rerun a slice of" >&2
    exit 2
  fi
  if [ "$rerun" -gt "$n_last" ]; then
    echo "launch.sh: the last launch had $n_last slices" >&2
    exit 2
  fi
  pid=$(pid_of "$rerun")
  if alive "$pid"; then
    echo "launch.sh: slice $rerun still runs (pid $pid)" >&2
    exit 2
  fi
  start_slice "$rerun" "$n_last"
  echo "started slice $rerun of $n_last again; log in $logs/slice-$rerun.log"
  exit 0
fi

busy=""
i=1
while [ "$i" -le "${n_last:-0}" ]; do
  pid=$(pid_of "$i")
  if alive "$pid"; then busy="$busy slice $i (pid $pid)"; fi
  i=$((i + 1))
done
if [ -n "$busy" ]; then
  echo "launch.sh: these slices of the last launch still run:$busy." \
    "Wait for them, or stop them first." >&2
  exit 2
fi

# OpenMP's settings would make nproc count fewer cores than the machine has
cores=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc 2>/dev/null ||
  getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null)
if ! whole "$cores"; then
  echo "launch.sh: could not read the number of cores; give -k" \
    "and --force" >&2
  [ -n "$k" ] && [ "$force" = 1 ] || exit 2
  cores=$((k * per_slice))
fi
k_max=$((cores / per_slice))
if [ -z "$k" ]; then
  k=$k_max
fi
echo "cores $cores, cores per slice $per_slice (chains x threads):" \
  "at most $k_max slices at once"
if [ "$k" -gt "$k_max" ] && [ "$force" != 1 ]; then
  echo "launch.sh: $k slices x $per_slice cores is more than $cores cores;" \
    "use -k $k_max, or --force" >&2
  exit 2
fi
if [ "$k" -lt 1 ]; then
  echo "launch.sh: one slice needs $per_slice cores and there are $cores;" \
    "use -k 1 --force" >&2
  exit 2
fi
if [ "$k" -gt "$n_cells" ]; then
  k=$n_cells
fi

mkdir -p "$logs"
echo "setup: $logs/setup.log"
Rscript "$slice" "$args_file" setup > "$logs/setup.log" 2>&1
code=$?
echo "EXIT $code" >> "$logs/setup.log"
if [ "$code" != 0 ]; then
  echo "launch.sh: the setup call failed; see $logs/setup.log" >&2
  exit 1
fi

: > "$logs/pids"
echo "$k" > "$logs/slices"
for i in $(seq 1 "$k"); do
  start_slice "$i" "$k"
done
echo "started $k slices over $n_cells cells; logs in $logs"
printf 'check with: bash %q --status %q\n' "$here/launch.sh" "$args_file"
