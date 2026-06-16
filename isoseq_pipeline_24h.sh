#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Sequencing Pipeline v0.8 - 24h stage - (Assembly, QC, MLST, AMR)
# ----------------------------------------------------------------------
IFS=$'\n\t'

source "$(dirname "$0")/isoseq_pipeline_config.sh"

SAMPLESHEET="$1"
THREADS="$2"
RUN_DIR="$3"
MASTER_SUMMARY="$4"
RUN_NAME="$5"

# ---------------------- READ SAMPLESHEET -----------------------
SAMPLES=(); BARCODES=(); ORGANISMS=(); COLLECTIONDATE=()
while IFS=, read -r barcode sample organism collection; do
    [[ "$barcode" == "barcode" || -z "$barcode" ]] && continue
    if [[ -z "${species_dict_size[$organism]:-}" ]]; then
        echo "❌ Unknown organism: $organism"; exit 1
    fi
    SAMPLES+=("$sample")
    BARCODES+=("$barcode")
    ORGANISMS+=("$organism")
    COLLECTIONS+=("$collection")
done < "$SAMPLESHEET"

# ============================================================== 
# PER-SAMPLE PIPELINE 
# ============================================================== 

for idx in "${!SAMPLES[@]}"; do 
{
    sample="${SAMPLES[$idx]}"
    barcode="${BARCODES[$idx]}"
    organism="${ORGANISMS[$idx]}"
    collection="${COLLECTIONS[$idx]}"

    sample_dir="/data/IsoSeq_results/${sample}/24h"
    mkdir -p "${sample_dir}/logs"
    LOGFILE="${sample_dir}/logs/pipeline.log"

    mkdir -p "${sample_dir}/samplesheet"
    cp -v "$SAMPLESHEET" "${sample_dir}/samplesheet/samplesheet.csv"
    touch "${sample_dir}/samplesheet/${RUN_NAME}"

    log() { echo -e "[$(date '+%F %T')] $*"; }
    exec 3>&1 4>&2
    exec > >(tee -a "$LOGFILE") 2>&1

    if [[ -f "/data/IsoSeq_results/${sample}/8h/kraken2_top5_taxa.csv" ]]; then
        top_hit=$(awk -F, 'NR>1 && $1 {print $2; exit}' "/data/IsoSeq_results/${sample}/8h/kraken2_top5_taxa.csv")
        if [[ "$top_hit" != "$organism" ]]; then
            log "8h analysis failed for $sample - sample skipped"; continue
            fi
    fi

    log "========== Processing sample: $sample ($organism) =========="
    start_total=$(date +%s)
    status="OK"; failed_step="-"

    # ---------------- STEP 1: Concatenate ----------------
    step="concat"; step_start=$(date +%s)
    log "▶ Step $step"
    src_pattern="${RUN_DIR}/**/fastq_pass/barcode${barcode}/*.fastq.gz"
    shopt -s globstar nullglob
    files=( $src_pattern )
    if (( ${#files[@]} == 0 )); then
        log "⚠️ No FASTQ found for barcode${barcode}"
        status="FAILED"; failed_step="$step"
        echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"
        continue
    fi
    log "Concatenating ${#files[@]} fastq files → ${sample_dir}/${sample}.fastq.gz"
    if ! cat "${files[@]}" > "${sample_dir}/${sample}_raw.fastq.gz"; then
       log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    ontime --to 24h -o "${sample_dir}/${sample}.fastq.gz" "${sample_dir}/${sample}_raw.fastq.gz"  # Use ontime to filter for only reads for this step
    rm -vf "${sample_dir}/${sample}_raw.fastq.gz"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 2: Porechop & NanoStat ----------------
    step="porechop_nanostat"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p ${sample_dir}/trimmed ${sample_dir}/qc
    fq="${sample_dir}/${sample}.fastq.gz"; out_trim="${sample_dir}/trimmed/${sample}.fastq.gz"
    if ! porechop -i "$fq" -o "$out_trim" -t "$THREADS" --no_split > /dev/tty; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    NanoStat --fastq "$out_trim" -n "${sample_dir}/qc/${sample}_nanostat.txt" -t "$THREADS"
    rm -vf "$fq"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 3: Assembly ----------------
    step="assembly"; step_start=$(date +%s)
    log "▶ Step $step"
    outdir="${sample_dir}/assemblies"; mkdir -p "$outdir"
    size="${species_dict_size[$organism]}"
    if ! bash /data/IsoSeq/scripts/autoautocycler.sh -o "$outdir/" -t "$THREADS" -c "2" -s "$size" -a "metamdbg myloasm" "${sample_dir}/trimmed/${sample}.fastq.gz"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 4: Coverage depth  ----------------
    step="mapping_depth"; step_start=$(date +%s)
    log "▶ Step $step"
    asm="${sample_dir}/assemblies/${sample}.fasta"; depth="${sample_dir}/qc/genome_coverage_depth.txt"
    if ! minimap2 -x map-ont --secondary no -a -t "$THREADS" "$asm" "${sample_dir}/trimmed/${sample}.fastq.gz" | samtools sort -O BAM | samtools coverage - | cut -f 7 | head -n 2 > "$depth"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 5: QUAST ----------------
    step="quast"; step_start=$(date +%s)
    log "▶ Step $step"
    if ! quast -o "${sample_dir}/qc/quast" -t "$THREADS" -m 100 -l "$sample" "$asm"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 6: Bandage ----------------
    step="bandage"; step_start=$(date +%s)
    log "▶ Step $step"
    gfa="${sample_dir}/assemblies/${sample}/autocycler_out/consensus_assembly.gfa"
    if ! Bandage image "$gfa" "${sample_dir}/qc/assembly_image.svg"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 7: BUSCO ----------------
    step="busco"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p "${sample_dir}/qc/busco_results"
    lineage="${species_dict_busco[$organism]}"
    if ! busco -i "$asm" -f -m genome -l "$lineage" -c "$THREADS" --out_path "${sample_dir}/qc/busco_results" -o "$sample"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 8: abritamr (AMRFINDER) ----------------
    step="amrfinder"; step_start=$(date +%s)
    log "▶ Step $step"
    amr_species="${species_dict_amrfinder[$organism]}"
    if ! conda run -n "$ENV_ABRITAMR" abritamr run -c "${sample_dir}/assemblies/${sample}.fasta" -px "$sample_dir" -j "$THREADS" -sp "$amr_species"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 9: MLST contigs ----------------
    step="mlst_contigs"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p ${sample_dir}/mlst_contigs
    scheme="${species_dict_mlst_contigs[$organism]}"
    if ! conda run -n "$ENV_MLST" mlst --scheme "$scheme" --full --blastdb "$MLST_BLAST_DB" --datadir "$MLST_DB" "${sample_dir}/assemblies/${sample}.fasta" > "${sample_dir}/mlst_result.tsv"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- Stage complete ----------------
    end_total=$(date +%s)
    total_time=$((end_total-start_total))
    status="24h COMPLETE"
    echo -e "${sample}\t${status}\t${failed_step}\t${total_time}" >> "${MASTER_SUMMARY}"
    log "✅ Sample $sample 24h analysis completed in ${total_time}s"
    exec 1>&3 2>&4
    exec 3>&- 4>&-
}; done
