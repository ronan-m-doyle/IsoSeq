#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Sequencing Master Pipeline v0.9
# ----------------------------------------------------------------------
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# -------------------------- ARG PARSE --------------------------
THREADS=1
SAMPLESHEET=""
RUN_DIR=""

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --samplesheet) SAMPLESHEET="$2"; shift; shift ;;
        --threads) THREADS="$2"; shift; shift ;;
        --rundir) RUN_DIR="$2"; shift; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$SAMPLESHEET" || -z "$THREADS" || -z "$RUN_DIR" ]]; then
    echo "❌ ERROR: Must specify --samplesheet, --threads, --rundir"
    exit 1
fi

RUN_DIR=$(realpath "$RUN_DIR")
[[ -d "$RUN_DIR" ]] || { echo "Run directory not found"; exit 1; }

SAMPLESHEET=$(realpath "$SAMPLESHEET")
[[ -f "$SAMPLESHEET" ]] || { echo "Samplesheet not found"; exit 1; }

RUN_NAME=$(basename "$RUN_DIR")

# ---------------------- PIPELINE LOG DIR -----------------------
mkdir -p "/data/IsoSeq_results/pipeline_logs"
RUN_TIMESTAMP=$(date '+%Y%m%d_%H%M%S')

touch "/data/IsoSeq_results/pipeline_logs/pipeline_summary_${RUN_NAME}_${RUN_TIMESTAMP}.txt"
MASTER_SUMMARY="/data/IsoSeq_results/pipeline_logs/pipeline_summary_${RUN_NAME}_${RUN_TIMESTAMP}.txt"

echo -e "Sample\tStatus\tFailed_Step\tTotal_Time(s)" > "$MASTER_SUMMARY"

echo "Monitoring sequencing run..."

STAGE1_FLAG="/data/IsoSeq_results/pipeline_logs/.8h_complete_${RUN_NAME}"
STAGE2_FLAG="/data/IsoSeq_results/pipeline_logs/.24h_complete_${RUN_NAME}"
STAGE3_FLAG="/data/IsoSeq_results/pipeline_logs/.48h_complete_${RUN_NAME}"
STAGE4_FLAG="/data/IsoSeq_results/pipeline_logs/.72h_complete_${RUN_NAME}"

while true; do

    # ---------------------------------
    # Re-resolve summary file each iteration.
    # At end-of-run the sequencer deletes the .tmp file and writes a
    # permanent sequencing_summary*.txt into $RUN_DIR. Re-resolving here
    # ensures we always point at whichever file currently exists.
    # Priority: finalised file in RUN_DIR > mid-run .tmp in shared temp tree.
    # ---------------------------------
    CURRENT_SUMMARY=$(find "$RUN_DIR" -name 'sequencing_summary*.txt' 2>/dev/null | head -n 1)
    if [[ -z "$CURRENT_SUMMARY" ]]; then
        CURRENT_SUMMARY=$(find /data/reads/tmp -path "*/${RUN_NAME}/*" -name 'sequencing_summary*.txt.tmp' 2>/dev/null | head -n 1)
    fi
 
    if [[ ! -f "$CURRENT_SUMMARY" ]]; then
        echo "Waiting for sequencing data..."
        sleep 600
        continue
    fi
 
    SUMMARY_FILE="$CURRENT_SUMMARY"
    echo "Summary file: $SUMMARY_FILE"
 
    # ---------------------------------
    # Safety check: ensure at least 2 data rows exist
    # ---------------------------------
    if [[ $(wc -l < "$SUMMARY_FILE") -lt 3 ]]; then
        echo "Waiting for sequencing data..."
        sleep 600
        continue
    fi

    # ---------------------------------
    # Calculate active sequencing time
    # Column position of "start_time" varies between MinKNOW/Guppy versions
    # (e.g. R9 summaries lack the filename_bam column that R10 has, shifting
    # every later column left by one), so look the column up by header name
    # each time rather than assuming a fixed index.
    # ---------------------------------
    last_start=$(awk -F'\t' '
        NR==1 {
            for (i=1; i<=NF; i++) { if ($i=="start_time") col=i }
            if (!col) { print "NO_START_TIME_COL"; exit 1 }
            next
        }
        { if ($col>max) max=$col }
        END { print int(max) }
    ' "$SUMMARY_FILE")

    if [[ "$last_start" == "NO_START_TIME_COL" || -z "$last_start" ]]; then
        echo "❌ Could not find 'start_time' column in $SUMMARY_FILE - skipping this iteration"
        sleep 600
        continue
    fi

    elapsed_hours=$((last_start / 3600))

    echo "[$(date)] Active sequencing time: ${elapsed_hours} hours"

    # ---------------------------------
    # 8h Stage (≥8 hour)
    # ---------------------------------
    if (( last_start >= 28800 )) && [[ ! -f "$STAGE1_FLAG" ]]; then
        echo "▶ Running 8h analysis"
        if bash "${SCRIPT_DIR}/isoseq_pipeline_8h.sh" "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"; then
            touch "$STAGE1_FLAG"
        else
            echo "❌ 8h analysis failed"
        fi
    fi

    #---------------------------------
    #24h Stage (≥24 hours)
    #---------------------------------
    if (( last_start >= 86400 )) && [[ ! -f "$STAGE2_FLAG" ]]; then
        echo "▶ Running 24h analysis"
        if bash "${SCRIPT_DIR}/isoseq_pipeline_24h.sh" "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"; then
            touch "$STAGE2_FLAG"
        else
            echo "❌ 24h analysis failed"
        fi
        echo "▶ Running phylogenetic analysis (24h)"
        if ! bash "${SCRIPT_DIR}/isoseq_pipeline_tree.sh" "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"; then
            echo "❌ Phylogenetic analysis (24h) failed"
        fi
    fi

    #---------------------------------
    #48h Stage (≥48 hours)
    #---------------------------------
    if (( last_start >= 172800 )) && [[ ! -f "$STAGE3_FLAG" ]]; then
        echo "▶ Running 48h analysis"
        if bash "${SCRIPT_DIR}/isoseq_pipeline_48h.sh" "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"; then
            touch "$STAGE3_FLAG"
        else
            echo "❌ 48h analysis failed"
        fi
        echo "▶ Running phylogenetic analysis (48h)"
        if ! bash "${SCRIPT_DIR}/isoseq_pipeline_tree.sh" "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"; then
            echo "❌ Phylogenetic analysis (48h) failed"
        fi
    fi

    # ---------------------------------
    # Detect End of Run (Final summary file)
    # ---------------------------------
    FINAL_SUMMARY_FILE=$(find "$RUN_DIR" -name "final_summary_*.txt" | head -n 1)

    if [[ -n "$FINAL_SUMMARY_FILE" ]]; then

        echo "[$(date)] End-of-run detected"
        echo "Final summary file: $FINAL_SUMMARY_FILE"

   	# ---------------------------------
  	# 72h Stage (End of run)
	# ---------------------------------
        if [[ ! -f "$STAGE4_FLAG" ]]; then
        	echo "▶ Running 72h analysis"
        	if bash "${SCRIPT_DIR}/isoseq_pipeline_72h.sh" "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"; then
        		touch "$STAGE4_FLAG"
        	else
        		echo "❌ 72h analysis failed"
        	fi
    	fi

        # ---------------------------------
    	# Tree building Stage - Will always re-run if pipeline is started
    	# ---------------------------------
        echo "▶ Running phylogenetic analysis (End of run)"
        if ! bash "${SCRIPT_DIR}/isoseq_pipeline_tree.sh" "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"; then
            echo "❌ Phylogenetic analysis (End of run) failed"
        fi

        echo "All stages complete. Exiting monitor."
        break
    fi

    # ---------------------------------
    # Sleep 1 hour
    # ---------------------------------
    echo "Sleeping for 1 hour..."
    sleep 3600

done

echo "Pipeline complete. Run log saved in $MASTER_SUMMARY"