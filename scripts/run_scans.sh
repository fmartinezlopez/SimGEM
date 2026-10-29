#!/usr/bin/env bash
# run_scans.sh -- run SimGEM scan files one after another (a simple queue).
#
#   scripts/run_scans.sh [options] scans/a.json [scans/b.json ...]
#
#   -j JOBS    workers per point (default: $SIMGEM_JOBS, else all cores)
#   -e EVENTS  events per point (overrides the scan files)
#   -n         dry run: list points, field maps and work for every scan
#   -d         detach: keep running after you log out (nohup); log in runs/scans/
#   -x         stop the queue at the first scan that fails
#   -h         this help
#
# Scans are resumable, so re-running a queue only does what is missing.
# Environment: SIMGEM_PYTHON = python that has gmsh (default: python3),
#              SIMGEM_JOBS, SIMGEM_RUNS, SIMGEM_FIELDS (see README).
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${SIMGEM_PYTHON:-python3}"
JOBS="" EVENTS="" DRY=0 DETACH=0 STOP=0

usage() { sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
while getopts "j:e:ndxh" opt; do
  case $opt in
    j) JOBS=$OPTARG ;;
    e) EVENTS=$OPTARG ;;
    n) DRY=1 ;;
    d) DETACH=1 ;;
    x) STOP=1 ;;
    h) usage 0 ;;
    *) usage 2 ;;
  esac
done
shift $((OPTIND - 1))
[ $# -ge 1 ] || { echo "no scan files given" >&2; usage 2; }

for f in "$@"; do
  [ -f "$f" ] || { echo "no such scan file: $f" >&2; exit 2; }
done
if ! "$PY" -c "import gmsh" 2>/dev/null; then
  echo "'$PY' cannot import gmsh; set SIMGEM_PYTHON to the right python" >&2
  echo "(e.g. SIMGEM_PYTHON=~/miniconda3/envs/garfield/bin/python)" >&2
  exit 2
fi

if [ $DETACH -eq 1 ]; then
  logdir="${SIMGEM_RUNS:-runs}/scans"
  mkdir -p "$logdir"
  log="$logdir/queue-$(date +%Y%m%d-%H%M%S).log"
  args=()
  [ -n "$JOBS" ] && args+=(-j "$JOBS")
  [ -n "$EVENTS" ] && args+=(-e "$EVENTS")
  [ $STOP -eq 1 ] && args+=(-x)
  [ $DRY -eq 1 ] && args+=(-n)
  nohup "$0" "${args[@]+"${args[@]}"}" "$@" > "$log" 2>&1 < /dev/null &
  echo "queue running in the background: PID $!"
  echo "log: $log      (follow with: tail -f $log)"
  exit 0
fi

opts=()
[ -n "$JOBS" ] && opts+=(--jobs "$JOBS")
[ -n "$EVENTS" ] && opts+=(--events "$EVENTS")
[ $DRY -eq 1 ] && opts+=(--dry-run)

declare -a names results times
t_start=$(date +%s)
for f in "$@"; do
  echo
  echo "################ $(date '+%F %T')  $f"
  t0=$(date +%s)
  "$PY" "$ROOT/gem.py" scan "$f" "${opts[@]+"${opts[@]}"}"
  rc=$?
  case $rc in
    0) res="ok" ;;
    3) res="done, some points FAILED" ;;
    *) res="FAILED (exit $rc)" ;;
  esac
  names+=("$f"); results+=("$res"); times+=($(( $(date +%s) - t0 )))
  if [ $rc -ne 0 ] && [ $STOP -eq 1 ]; then
    echo "stopping the queue (-x)"
    break
  fi
done

echo
echo "################ queue summary ($(( ($(date +%s) - t_start) / 60 )) min)"
fail=0
for i in "${!names[@]}"; do
  printf "  %-40s %-28s %5d min\n" "${names[$i]}" "${results[$i]}" $(( ${times[$i]} / 60 ))
  [ "${results[$i]}" = "ok" ] || fail=1
done
exit $fail
