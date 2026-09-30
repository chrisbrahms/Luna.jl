#!/bin/sh
# Overnight queue of the fixed-versus-adaptive study.
#
#     benchmark/modal_rules/queue.sh <plan file> <R deadline> <R kill> <deadline> <kill>
#
# Times are HH:MM (the next occurrence). The plan file has one run per line,
# `lane est_seconds case/rule/prtol`. Lanes `R1` ... `R4` hold the accuracy runs and run
# concurrently (one serial process each); their run times are indicative only. A run in an
# R lane is not started if its estimate would take it past <R deadline>, and whatever is
# still running at <R kill> is killed. Then lane `T` -- the timed runs, written to
# fields_timed/ and results/timed.csv -- runs alone, with
# <deadline> and <kill> in the same roles. The kills are by PID, from this queue's own list.
# run.jl skips runs whose result file exists, so the queue can be restarted as it is.
set -u
PLAN=$1; RDEADLINE=$2; RKILL=$3; DEADLINE=$4; KILLAT=$5
cd "$(dirname "$0")/../.."
LOGDIR=benchmark/modal_rules/logs; mkdir -p $LOGDIR
RPIDS=$LOGDIR/pids_R; TPIDS=$LOGDIR/pids_T
: > $RPIDS; : > $TPIDS

epoch() { date -j -f "%Y-%m-%d %H:%M" "$(date +%Y-%m-%d) $1" +%s; }
tomorrow_if_past() { t=$(epoch "$1"); [ "$t" -lt "$(date +%s)" ] && t=$((t + 86400)); echo $t; }
RDL=$(tomorrow_if_past "$RDEADLINE"); RKT=$(tomorrow_if_past "$RKILL")
DL=$(tomorrow_if_past "$DEADLINE"); KT=$(tomorrow_if_past "$KILLAT")

watchdog() { # kill time, pid file, name
    ( now=$(date +%s); [ $1 -gt $now ] && sleep $(($1 - now))
      echo "$(date +%H:%M) watchdog $3: kill time reached"
      while read -r p; do kill "$p" 2>/dev/null && echo "$(date +%H:%M) watchdog $3: killed $p"; done < $2 ) \
      >> $LOGDIR/watchdog.log 2>&1 &
    echo $!
}
WR=$(watchdog $RKT $RPIDS R)
WT=$(watchdog $KT $TPIDS T)

lane() { # name, deadline, pid file
    grep "^$1 " "$PLAN" | while read -r _ est spec; do
        now=$(date +%s)
        if [ $((now + est)) -gt $2 ]; then
            echo "$(date +%H:%M) $1 skip (deadline) $spec"; continue
        fi
        echo "$(date +%H:%M) $1 start $spec (est ${est}s)"
        # the timed runs repeat accuracy runs; they go to fields_timed/ and results/timed.csv
        if [ "$1" = T ]; then F=timed; else F=""; fi
        if [ -n "$F" ]; then
            FIELDS=$F OUT=$F RUNS="$spec" julia --project=benchmark -t 1 \
                benchmark/modal_rules/run.jl >> $LOGDIR/$1.log 2>&1 &
        else
            RUNS="$spec" julia --project=benchmark -t 1 benchmark/modal_rules/run.jl \
                >> $LOGDIR/$1.log 2>&1 &
        fi
        jp=$!; echo $jp >> $3; wait $jp; rc=$?
        echo "$(date +%H:%M) $1 end $spec (exit $rc)"
    done
}

LANES=""
for L in R1 R2 R3 R4; do
    lane $L $RDL $RPIDS &
    LANES="$LANES $!"; echo $! >> $RPIDS
done
for p in $LANES; do wait $p; done
echo "$(date +%H:%M) R lanes done"
lane T $DL $TPIDS
echo "QUEUE_DONE $(date)"
kill $WR $WT 2>/dev/null
