#!/usr/bin/env bash
#
# Start, pause and resume the HR-WSI download (scripts/1b-download_data-hrwsi.py)
# as a detached background job, so it keeps running after Claude, the terminal,
# or the login session goes away.
#
#   scripts/hrwsi-download.sh start [max_workers]   begin, or resume where it left off
#   scripts/hrwsi-download.sh stop                  pause: stop cleanly, keep what is on disk
#   scripts/hrwsi-download.sh status                is it running, and how far along
#   scripts/hrwsi-download.sh log                   follow the live log (Ctrl-C to detach)
#
# Pausing is safe at any moment. The downloader writes each file to a .part and
# renames it only once complete, so every file on disk is whole; `start` rescans
# the archive (a few minutes) and skips everything already downloaded.
#
# `start 2` (or 1) halves the parallel transfers if you want the download to
# continue while leaving bandwidth for other things.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

out_dir="derived/eo/hrwsi"
log_file="logs/hrwsi-download.log"
pid_file="$out_dir/.download.pid"

running_pid() {
  [[ -f "$pid_file" ]] || return 1
  local pid
  pid="$(cat "$pid_file" 2>/dev/null || true)"
  [[ -n "$pid" ]] || return 1
  kill -0 "$pid" 2>/dev/null || return 1
  echo "$pid"
}

on_disk() {
  local files bytes
  files="$(find "$out_dir" -type f -name '*.tif' 2>/dev/null | wc -l | tr -d ' ')"
  bytes="$(find "$out_dir" -type f ! -name '*.part' -exec stat -f%z {} + 2>/dev/null \
           | awk '{s+=$1} END {printf "%.2f", s/1e9}')"
  echo "${files:-0} raster files, ${bytes:-0.00} GB on disk"
}

case "${1:-status}" in
  start)
    if pid="$(running_pid)"; then
      echo "Already running (pid $pid). Use 'stop' first, or 'status'."
      exit 0
    fi
    workers="${2:-}"
    mkdir -p "$out_dir" logs
    # Drop partial transfers from any previous stop; they are re-fetched whole.
    find "$out_dir" -name '*.part' -delete 2>/dev/null || true
    args=(uv run python scripts/1b-download_data-hrwsi.py)
    [[ -n "$workers" ]] && args+=(--max-workers "$workers")
    echo "--- started $(date '+%Y-%m-%d %H:%M:%S') ${workers:+(max_workers=$workers)} ---" >> "$log_file"
    nohup "${args[@]}" >> "$log_file" 2>&1 &
    echo $! > "$pid_file"
    sleep 2
    if pid="$(running_pid)"; then
      echo "Started (pid $pid). It keeps running if you close Claude or this terminal."
      echo "  progress : scripts/hrwsi-download.sh status"
      echo "  pause    : scripts/hrwsi-download.sh stop"
    else
      echo "Failed to start; last log lines:"; tail -20 "$log_file"; exit 1
    fi
    ;;

  stop)
    if ! pid="$(running_pid)"; then
      echo "Not running. $(on_disk)"
      rm -f "$pid_file"
      exit 0
    fi
    kill -TERM "$pid" 2>/dev/null || true
    for _ in $(seq 1 20); do
      kill -0 "$pid" 2>/dev/null || break
      sleep 0.5
    done
    kill -0 "$pid" 2>/dev/null && kill -KILL "$pid" 2>/dev/null || true
    rm -f "$pid_file"
    echo "--- paused $(date '+%Y-%m-%d %H:%M:%S') ---" >> "$log_file"
    echo "Paused. $(on_disk)"
    echo "Resume any time with: scripts/hrwsi-download.sh start"
    ;;

  status)
    if pid="$(running_pid)"; then
      echo "RUNNING (pid $pid)"
    elif [[ -f "$log_file" ]] \
         && [[ "$(grep -n -e '^Done:' -e '^--- started' "$log_file" | tail -1)" == *Done:* ]]; then
      # The last thing the log recorded was a completed run, not an interruption.
      echo "COMPLETE"
      grep '^Done:' "$log_file" | tail -1
    else
      echo "PAUSED / not running"
    fi
    echo "$(on_disk)"
    if [[ -f "$log_file" ]]; then
      echo "--- last log lines ---"
      tail -5 "$log_file"
    fi
    ;;

  log)
    tail -f "$log_file"
    ;;

  *)
    sed -n '3,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 1
    ;;
esac
