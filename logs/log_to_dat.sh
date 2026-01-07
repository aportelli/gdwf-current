#!/usr/bin/env bash

set -euo pipefail

if (( $# != 1 )); then
    echo "usage: $(basename "$0") <log file>" 1>&2
    exit 1
fi
log=$1

eps=$(sed -n 's/.*epsilon= *\([0-9.]\+\).*/\1/p' "${log}" | head -n1)
awk -v eps="${eps}" '
/1st order numerical derivative/ { mode=1 }
/2nd order numerical derivative/ { mode=2 }
/Diff --/ {
    if (mode==1) {
        a=$(NF-2); b=$NF
    } else if (mode==2) {
        c=$(NF-2); d=$NF
        print eps, a, b, c, d 
    }
}
' "${log}"
