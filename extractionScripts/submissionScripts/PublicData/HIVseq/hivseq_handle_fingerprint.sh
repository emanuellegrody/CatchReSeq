#!/bin/bash
# hivseq_handle_fingerprint.sh  (Figure3, HIVseq library chemistry)
# Counts how many reads carry the HIV capture sequences immediately after the RT
# primer handle, per library. Only Primer libraries do; NOTHING libraries do not.
# Figure3.Rmd cites this as the evidence behind its Primer vs NOTHING comparison
# in the ART-suppressed draws.
#
# Run: bash hivseq_handle_fingerprint.sh
set -uo pipefail

# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../../pipelineScripts"

# Filesystem root for the downloaded HIVseq FASTQs. Override to run against a
# different copy:
#   PUBLIC_ROOT=/my/publicData bash hivseq_handle_fingerprint.sh
PUBLIC_ROOT="${PUBLIC_ROOT:-/projects/b1042/GoyalLab/egrody/publicData}"

cd "${PUBLIC_ROOT}/HIVseq/fastq"
PY="${PIPELINE_DIR}/publicData/hivseq_handle_fingerprint.py"
N=4000000
while read label run; do
  f=$(ls ${run}_4.fastq.gz ${run}_2.fastq.gz 2>/dev/null | head -1)
  python3 "$PY" "$label" "$f" "$N" &
done <<LIST
PID1052w0r1_NOTHING SRR28878447
PID1052w0r3_Primer SRR28878441
PID1052w0r5_Primer SRR28878435
PID1052w70r2_NOTHING SRR28878444
PID1052w70r4_Primer SRR28878438
PID1052w70r6_Primer SRR28878432
PID8027w0r1_NOTHING SRR28878388
PID8027w0r2_NOTHING SRR28878384
PID8027w0r5_Primer SRR28878396
PID8027w0r6_Primer SRR28878392
PID8027w45r3_NOTHING SRR28878382
PID8027w45r4_NOTHING SRR28878400
PID8027w45r7_Primer SRR28878379
PID8027w45r8_Primer SRR28878431
PID8026w0r5_Primer SRR28878409
PID8026w24r1_Primer SRR28878424
LIST
wait
