#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Sequencing Master Pipeline v0.8
# ----------------------------------------------------------------------
IFS=$'\n\t'

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
    # ---------------------------------
    last_start=$(awk -F'\t' 'NR>1 {if ($11>max) max=$11} END {print int(max)}' "$SUMMARY_FILE")

    elapsed_hours=$((last_start / 3600))

    echo "[$(date)] Active sequencing time: ${elapsed_hours} hours"

    # ---------------------------------
    # 8h Stage (≥8 hour)
    # ---------------------------------
    if (( last_start >= 28800 )) && [[ ! -f "$STAGE1_FLAG" ]]; then
        echo "▶ Running 8h analysis"
        bash isoseq_pipeline_8h.sh "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"
        touch "$STAGE1_FLAG"
    fi

    #---------------------------------
    #24h Stage (≥24 hours)
    #---------------------------------
    if (( last_start >= 86400 )) && [[ ! -f "$STAGE2_FLAG" ]]; then
        echo "▶ Running 24h analysis"
        bash isoseq_pipeline_24h.sh "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"
        echo "▶ Running phylogenetic analysis (24h)"
        bash isoseq_pipeline_tree.sh "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"
        touch "$STAGE2_FLAG"
    fi

    #---------------------------------
    #48h Stage (≥48 hours)
    #---------------------------------
    if (( last_start >= 172800 )) && [[ ! -f "$STAGE3_FLAG" ]]; then
        echo "▶ Running 48h analysis"
        bash isoseq_pipeline_48h.sh "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"
        echo "▶ Running phylogenetic analysis (48h)"
        bash isoseq_pipeline_tree.sh "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"
        touch "$STAGE3_FLAG"
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
        	bash isoseq_pipeline_72h.sh "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"
        	touch "$STAGE4_FLAG"
    	fi

        # ---------------------------------
    	# Tree building Stage - Will always re-run if pipeline is started
    	# ---------------------------------
        echo "▶ Running phylogenetic analysis (End of run)"
        bash isoseq_pipeline_tree.sh "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"

        echo "All stages complete. Exiting monitor."
        break
    fi

    # ---------------------------------
    # Sleep 30 minutes
    # ---------------------------------
    echo "Sleeping for 30 minutes..."
    sleep 1800

done

echo "Pipeline complete. Run log saved in $MASTER_SUMMARY"
