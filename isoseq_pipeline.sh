#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Sequencing Master Pipeline v0.7
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

SUMMARY_FILE=$(find "$RUN_DIR" -name "sequencing_summary*.txt" | head -n 1)

if [[ ! -f "$SUMMARY_FILE" ]]; then
    echo "❌ sequencing_summary.txt not found"
    exit 1
fi

# ---------------------- PIPELINE LOG DIR -----------------------
RUN_NAME=$(basename "$RUN_DIR")
mkdir -p "/data/IsoSeq_results/pipeline_logs"
RUN_TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
<<<<<<< HEAD
<<<<<<< HEAD
touch "pipeline_logs/pipeline_summary_${RUN_NAME}_${RUN_TIMESTAMP}.tsv"
MASTER_SUMMARY="pipeline_logs/pipeline_summary_${RUN_NAME}_${RUN_TIMESTAMP}.tsv"
=======
touch "pipeline_logs/pipeline_summary_${RUN_NAME}_${RUN_TIMESTAMP}.txt"
MASTER_SUMMARY="pipeline_logs/pipeline_summary_${RUN_NAME}_${RUN_TIMESTAMP}.txt"
>>>>>>> 5697c77 (	modified:   isoseq_pipeline.sh)
=======
touch "/data/IsoSeq_results/pipeline_logs/pipeline_summary_${RUN_NAME}_${RUN_TIMESTAMP}.txt"
MASTER_SUMMARY="/data/IsoSeq_results/pipeline_logs/pipeline_summary_${RUN_NAME}_${RUN_TIMESTAMP}.txt"
>>>>>>> 9312397 (Updates for version 0.7 of the pipeline)
echo -e "Sample\tStatus\tFailed_Step\tTotal_Time(s)" > "$MASTER_SUMMARY"

echo "Monitoring sequencing run..."
echo "Summary file: $SUMMARY_FILE"

<<<<<<< HEAD
STAGE1_FLAG="pipeline_logs/.8h_complete_${RUN_NAME}"
STAGE2_FLAG="pipeline_logs/.24h_complete_${RUN_NAME}"
STAGE3_FLAG="pipeline_logs/.48h_complete_${RUN_NAME}"
STAGE4_FLAG="pipeline_logs/.72h_complete_${RUN_NAME}"
<<<<<<< HEAD
TREE_FLAG="pipeline_logs/.tree_complete_${RUN_NAME}"
=======
>>>>>>> 5697c77 (	modified:   isoseq_pipeline.sh)
=======
STAGE1_FLAG="/data/IsoSeq_results/pipeline_logs/.8h_complete_${RUN_NAME}"
STAGE2_FLAG="/data/IsoSeq_results/pipeline_logs/.24h_complete_${RUN_NAME}"
STAGE3_FLAG="/data/IsoSeq_results/pipeline_logs/.48h_complete_${RUN_NAME}"
STAGE4_FLAG="/data/IsoSeq_results/pipeline_logs/.72h_complete_${RUN_NAME}"
>>>>>>> 9312397 (Updates for version 0.7 of the pipeline)

while true; do

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
        touch "$STAGE2_FLAG"
    fi

    #---------------------------------
    #48h Stage (≥48 hours)
    #---------------------------------
    if (( last_start >= 172800 )) && [[ ! -f "$STAGE3_FLAG" ]]; then
        echo "▶ Running 48h analysis"
        bash isoseq_pipeline_48h.sh "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"
        touch "$STAGE3_FLAG"
    fi

    # ---------------------------------
    # Detect End of Run (Final summary file)
    # ---------------------------------
<<<<<<< HEAD
    if (( last_start >= 258000 )) && [[ ! -f "$STAGE4_FLAG" ]]; then
        echo "▶ Running 72h analysis"
        bash isoseq_pipeline_72h.sh "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"
        touch "$STAGE4_FLAG"
=======
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
    
    	# Tree building Stage - Will always re-run if pipeline is started
    	# ---------------------------------
    #    echo "▶ Running phylogenetic analysis"
    #    bash isoseq_pipeline_tree.sh "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"

>>>>>>> 5697c77 (	modified:   isoseq_pipeline.sh)
        echo "All stages complete. Exiting monitor."
        break
    fi
    
    # ---------------------------------
    # Tree building Stage
    # ---------------------------------
    if (( last_start >= 258000 )) && [[ ! -f "$TREE_FLAG" ]]; then
        echo "▶ Running phylogenetic analysis"
        bash isoseq_pipeline_tree.sh "$SAMPLESHEET" "$THREADS" "$RUN_DIR" "$MASTER_SUMMARY" "$RUN_NAME"
        touch "$TREE_FLAG"
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
