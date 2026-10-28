#!/usr/bin/env bash
# Run a recovery grid on one machine as K processes, each a slice.
#
#   bash launch.sh grid-args.R              the setup call, then K slices
#   bash launch.sh -k 6 grid-args.R         K = 6 slices
#   bash launch.sh --force -k 20 grid-args.R
#   bash launch.sh --status grid-args.R     how far the setup and slices are
#   bash launch.sh --rerun 3 grid-args.R    slice 3 of the last launch again
#
# K defaults to the cores of the machine divided by the cores one slice
# uses (chains x threads), and is refused above that without --force.
# The setup call runs in the foreground, so start the launcher inside
# tmux or screen when the setup is long. The slices run in the
# background under nohup, so they go on when the SSH session ends. Each
# writes <dir>/logs/slice-<k>.log, which ends with "EXIT <status>".
# One launch at a time: a launch is refused while another one sets up
# the same grid or while a slice of the last one still runs. To stop a
# slice, kill the pid that <dir>/logs/pids lists for it (a plain kill,
# not kill -9). It needs a working ps to tell its processes apart.
# Run it from the directory that grid-args.R's `dir` is relative to.

set -u
# an exported CDPATH can send cd to another directory of the same name,
# and makes cd print the directory, so $(cd ... && pwd) gives two lines
unset CDPATH
here="$(cd "$(dirname "$0")" && pwd)"
slice="$here/grid-slice.R"

# ps checks that a recorded pid is still ours and finds what R started;
# without a working ps, a slice counts as running while a process group
# with its recorded id is live, and the setup call while a process with
# its recorded id is. -ww: the whole command line, however wide the
# terminal.
if ps -ww -p $$ -o pid= > /dev/null 2>&1; then have_ps=1; else have_ps=0; fi

# appends "EXIT <status $1>" to log $2 on a line of its own, also after
# output that ended without a newline
log_exit() {
  if [ -s "$2" ] && [ -n "$(tail -c 1 "$2")" ]; then echo >> "$2"; fi
  echo "EXIT $1" >> "$2"
}

# $1 and every process below it. processx starts Stan's chains in a
# session of their own, which a signal to the process group misses.
tree() {
  if [ "$have_ps" = 0 ]; then echo "$1"; return; fi
  ps -ww -A -o pid= -o ppid= | awk -v top="$1" '
    { parent[$1] = $2 }
    END {
      keep[top] = 1
      grew = 1
      while (grew) {
        grew = 0
        for (p in parent) {
          if (!(p in keep) && (parent[p] in keep)) { keep[p] = 1; grew = 1 }
        }
      }
      for (p in keep) print p
    }'
}

if [ "${1:-}" = "--run-slice" ]; then
  # one background slice: launch.sh starts this under nohup, as the
  # leader of a process group of its own. R and what it starts run in
  # that group or below R, so a TERM sent here reaches all of them, and
  # the log still gets its exit line. The mark $6 goes to R as well,
  # unused there: it tells both processes apart as this grid's slice $3
  # (see mark below).
  Rscript "$slice" "$2" "$3" "$4" "$6" > "$5" 2>&1 &
  child=$!
  trap 'trap "" TERM INT
    kill -TERM $(tree "$child") 2>/dev/null
    kill -TERM -- "-$$" 2>/dev/null' TERM INT
  wait "$child"
  code=$?
  if kill -0 "$child" 2>/dev/null; then
    wait "$child"
    code=$?
  fi
  log_exit "$code" "$5"
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
    -h | --help) sed -n '2,/^$/p' "$0"; exit 0 ;;
    -*)
      echo "launch.sh: unknown option $1 (see --help)" >&2
      exit 2
      ;;
    *)
      if [ -n "$args_file" ]; then
        echo "launch.sh: give one grid-args.R file, not $args_file and $1" >&2
        exit 2
      fi
      args_file="$1"
      shift
      ;;
  esac
done
if [ -z "$args_file" ] || [ ! -f "$args_file" ]; then
  echo "launch.sh: give the grid-args.R file (see --help)" >&2
  exit 2
fi
if [ "$have_ps" = 0 ]; then
  echo "launch.sh: ps does not work here, so a recorded pid that another" \
    "process has taken counts as running, and a kill of a slice misses" \
    "what its R started in a session of its own, such as the Stan chains" \
    "that cmdstanr starts through processx" >&2
fi

# grid-args.R may print, so the answer is the line that starts with
# grid-slice.R's marker
out=$(Rscript "$slice" "$args_file" info 2>&1)
info=$(printf '%s\n' "$out" | sed -n 's/^launch-info: //p' | tail -n 1)
if [ -z "$info" ]; then
  echo "launch.sh: 'grid-slice.R $args_file info' gave no launch-info line;" \
    "it printed:" >&2
  printf '%s\n' "$out" | tail -n 20 >&2
  exit 2
fi
read -r per_slice n_cells dir <<< "$info"
if ! whole "${per_slice:-}" || ! whole "${n_cells:-}" || [ -z "${dir:-}" ]; then
  echo "launch.sh: could not read the launch-info line: $info" >&2
  exit 2
fi
logs="$dir/logs"
lock="$logs/launch.lock"

# the logs folder with links resolved, the same on every call
real_logs() { (cd "$logs" 2> /dev/null && pwd -P); }
# the mark of this grid's process $1 (setup, or slice-<k>): an argument
# that the setup call, and a slice's wrapper and R, get last. It is made
# from the resolved logs folder, but of ASCII alone, which ps prints as
# it is in any locale; a path with a non-ASCII letter in it, ps may print
# escaped. It is compared as a whole argument, so slice-1 is not
# slice-10.
mark() {
  local real
  real=$(real_logs)
  [ -n "$real" ] || return 1
  printf 'bmmtools-grid=%s-%s\n' \
    "$(printf %s "$real" | cksum | awk '{ print $1 "-" $2 }')" "$1"
}
# whether pid $1 has this grid's mark for $2 as an argument; without a
# working ps, whether pid $1 is live
marked() {
  local m
  [ -n "$1" ] || return 1
  if [ "$have_ps" = 0 ]; then
    kill -0 "$1" 2>/dev/null
    return
  fi
  m=$(mark "$2") || return 1
  # the mark goes through the environment, so that no command line here
  # has it
  ps -ww -p "$1" -o args= 2>/dev/null | m="$m" awk '
    { for (i = 1; i <= NF; i++) if ($i == ENVIRON["m"]) found = 1 }
    END { exit !found }'
}
# whether pid $1 runs a command line that contains the text $2
runs() {
  [ -n "$1" ] || return 1
  if [ "$have_ps" = 0 ]; then
    kill -0 "$1" 2>/dev/null
    return
  fi
  ps -ww -p "$1" -o args= 2>/dev/null | grep -qF -- "$2"
}
# whether pid $1 is a launch.sh that sets up or starts a grid, and not
# the wrapper of a slice (bash launch.sh --run-slice ...), of this grid
# or another one
launch_runs() {
  runs "$1" launch.sh || return 1
  [ "$have_ps" = 0 ] && return 0
  ! runs "$1" --run-slice
}
# how slice $2, whose process group is $1, runs: "slice" while its
# wrapper does, "orphan <pid of its R>" while only its R does (the
# wrapper was killed with -9), or nothing. Both have this grid's mark
# for slice $2 as an argument, so a pid that another process has taken
# since (after a reboot, or another grid's slice) does not count.
# Without a working ps: "slice" while the group and pid $1 are live,
# "group" while only the group is (the wrapper was killed with -9).
slice_state() {
  local m
  [ -n "$1" ] || return 0
  if [ "$have_ps" = 0 ]; then
    if kill -0 -- "-$1" 2>/dev/null; then
      if kill -0 "$1" 2>/dev/null; then echo slice; else echo group; fi
    fi
    return 0
  fi
  m=$(mark "slice-$2") || return 0
  ps -ww -A -o pid= -o pgid= -o args= | m="$m" awk -v g="$1" '
    $2 == g {
      for (i = 3; i <= NF; i++) if ($i == ENVIRON["m"]) {
        if ($1 == g) w = 1; else r = $1
      }
    }
    END { if (w) print "slice"; else if (r) print "orphan", r }'
}
# the pid last recorded for slice $1, which is also its process group
pid_of() {
  [ -f "$logs/pids" ] || return 0
  awk -v k="$1" '$1 == k { p = $2 } END { print p }' "$logs/pids"
}
running_text() {
  case "$2" in
    orphan*)
      # R's processes are listed while R runs: once it is gone, the
      # chains it started have no parent to be found by
      r=${2#orphan }
      echo "running (pid $1 is gone but its R, pid $r, is not;" \
        "kill $(tree "$r" | tr '\n' ' ')stops R and what it started)"
      ;;
    group)
      # without ps, R's pid cannot be found, but it is in the group
      echo "running (pid $1 is gone but processes of its group are not;" \
        "kill -- -$1 stops them)"
      ;;
    *) echo "running (pid $1)" ;;
  esac
}
start_slice() {
  # emptied here, so that the old exit line is gone before the slice starts
  : > "$logs/slice-$1.log"
  # set -m gives the slice a process group of its own, whose id is the
  # pid recorded here
  set -m
  nohup bash "$here/launch.sh" --run-slice "$args_file" "$1" "$2" \
    "$(real_logs)/slice-$1.log" "$(mark "slice-$1")" \
    < /dev/null > /dev/null 2>&1 &
  echo "$1 $!" >> "$logs/pids"
  set +m
}
# whether the lock is stale: neither its launch.sh nor the setup call
# that one started still runs (a launch.sh killed with -9 leaves its
# setup running). The setup call has this grid's setup mark as an
# argument, so a setup-pid that another process has taken since (another
# grid's setup or slice, or a tail -f of setup.log) does not count. A
# lock with no pid yet is not stale: another launch made it a moment ago.
lock_stale() {
  holder=$(cat "$lock/pid" 2>/dev/null || true)
  [ -n "$holder" ] && ! launch_runs "$holder" &&
    ! marked "$(cat "$lock/setup-pid" 2>/dev/null || true)" setup
}
# stops the setup call and what it started, when the launcher is stopped
stop_setup() {
  if [ -n "${setup_pid:-}" ]; then
    # shellcheck disable=SC2046 # one pid per word
    kill -TERM $(tree "$setup_pid") 2>/dev/null
  fi
}
# held from before the setup until the slices are started, so that two
# launches never set up or start the same grid at once
take_lock() {
  mkdir -p "$logs"
  got=0
  if mkdir "$lock" 2>/dev/null; then
    got=1
  elif lock_stale && mkdir "$lock.takeover" 2>/dev/null; then
    # a stale lock is taken over under a second lock, and checked again
    # there, so that of two launches that find it stale only one takes it
    if lock_stale; then
      rm -rf "$lock"
      if mkdir "$lock" 2>/dev/null; then got=1; fi
    fi
    rmdir "$lock.takeover"
  fi
  if [ "$got" = 0 ]; then
    holder=$(cat "$lock/pid" 2>/dev/null || true)
    sp=$(cat "$lock/setup-pid" 2>/dev/null || true)
    who="another launch.sh${holder:+ (pid $holder)} is"
    if ! launch_runs "$holder" && marked "$sp" setup; then
      who="the setup call (pid $sp) of a launch.sh that is gone is"
    fi
    also=""
    if [ -d "$lock.takeover" ]; then also=" and $lock.takeover"; fi
    echo "launch.sh: $who setting up or starting this grid. If you know" \
      "it is gone, remove $lock$also" >&2
    exit 2
  fi
  echo $$ > "$lock/pid"
  trap 'rm -rf "$lock"' EXIT
  trap 'stop_setup; exit 129' HUP
  trap 'stop_setup; exit 130' INT
  trap 'stop_setup; exit 143' TERM
}
n_last=$(cat "$logs/slices" 2>/dev/null || true)

if [ "$status" = 1 ]; then
  if [ ! -d "$dir" ]; then
    echo "no grid in $dir yet"
    exit 0
  fi
  end=$(grep -E '^EXIT [0-9]+$' "$logs/setup.log" 2>/dev/null | tail -n 1)
  holder=$(cat "$lock/pid" 2>/dev/null || true)
  sp=$(cat "$lock/setup-pid" 2>/dev/null || true)
  if [ ! -f "$logs/setup.log" ]; then
    state="none yet"
  elif [ -n "$end" ]; then
    state=$end
  elif launch_runs "$holder"; then
    state="running (launch.sh pid $holder)"
  elif marked "$sp" setup; then
    state="running (pid $sp), but its launch.sh is gone: no slices will"
    if [ "$have_ps" = 1 ]; then
      state="$state start; kill $(tree "$sp" | tr '\n' ' ')stops it and what it"
      state="$state started"
    else
      # without ps, what the setup started in a session of its own cannot
      # be found, and outlives this kill
      state="$state start; kill $sp stops it"
    fi
  else
    state="stopped without an exit line"
  fi
  echo "setup.log: $state"
  if [ -n "$n_last" ] && [ "$logs/setup.log" -nt "$logs/slices" ]; then
    echo "(the slices below are of the launch before this setup)"
  fi
  done_cells=$(find "$dir" -maxdepth 1 -name 'cell-*-est.rds' | wc -l |
    tr -d ' ')
  echo "cells with a result: $done_cells of $n_cells"
  i=1
  while [ "$i" -le "${n_last:-0}" ]; do
    end=$(grep -E '^EXIT [0-9]+$' "$logs/slice-$i.log" 2>/dev/null | tail -n 1)
    pid=$(pid_of "$i")
    run=$(slice_state "$pid" "$i")
    if [ -n "$end" ]; then
      state=$end
    elif [ -n "$run" ]; then
      state=$(running_text "$pid" "$run")
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
  take_lock
  pid=$(pid_of "$rerun")
  run=$(slice_state "$pid" "$rerun")
  if [ -n "$run" ]; then
    echo "launch.sh: slice $rerun still runs (pid $pid). Stop it first" \
      "(see --status); if you know it is gone, remove its line from" \
      "$logs/pids" >&2
    exit 2
  fi
  start_slice "$rerun" "$n_last"
  echo "started slice $rerun of $n_last again; log in $logs/slice-$rerun.log"
  exit 0
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

take_lock
# read again under the lock: a launch may have ended since
n_last=$(cat "$logs/slices" 2>/dev/null || true)
busy=""
i=1
while [ "$i" -le "${n_last:-0}" ]; do
  pid=$(pid_of "$i")
  if [ -n "$(slice_state "$pid" "$i")" ]; then
    busy="$busy slice $i (pid $pid)"
  fi
  i=$((i + 1))
done
if [ -n "$busy" ]; then
  echo "launch.sh: these slices of the last launch still run:$busy." \
    "Wait for them, or stop them first (see --status). If you know" \
    "they are gone, remove $logs/pids" >&2
  exit 2
fi

echo "setup: $logs/setup.log"
# waited for in the background, so that its pid can go in the lock: a
# launch.sh killed with -9 leaves the lock held while the setup runs.
# The mark goes to R as well, unused there: it tells the process apart
# as this grid's setup call.
Rscript "$slice" "$args_file" setup "$(mark setup)" \
  > "$logs/setup.log" 2>&1 &
setup_pid=$!
echo "$setup_pid" > "$lock/setup-pid"
wait "$setup_pid"
code=$?
setup_pid=""
log_exit "$code" "$logs/setup.log"
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
