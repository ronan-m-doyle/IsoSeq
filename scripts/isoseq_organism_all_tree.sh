#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Sequencing Pipeline v0.9 - Organism-wide all-sample phylogenetic analysis
#
# Takes ONLY an organism and a thread count.
#
# Samples for an organism live in per-ST sub-folders. This script:
#   1. Collates every ST's sample_collection_dates.csv into a single
#      organism-level CSV (adding an "st" column) at:
#          ${species}_trees/sample_collection_dates.csv
#   2. Picks ONE reference for the whole organism.
#   3. Maps/calls/builds consensus for samples not already processed,
#      ONE SAMPLE AT A TIME, with each sample using the full --threads
#      value (minimap2 -t, clair3 --threads).
#   4. Builds one combined tree + one combined SNP matrix across all
#      STs, once, at the end - also using the full --threads value.
#      Nothing is written back to the original per-sample directories;
#      all outputs stay in ${species}_trees/.
# ----------------------------------------------------------------------

IFS=$'\n\t'

source "$(dirname "$0")/../isoseq_pipeline_config.sh"

log() { echo -e "[$(date '+%F %T')] $*"; }
export -f log

usage() {
    cat <<EOF
Usage: $0 -o|--organism "<organism>" -t|--threads <n>

<organism> must exactly match one of the following (case-sensitive):
$(printf '  - %s\n' "${!species_dict_size[@]}" | sort)
EOF
    exit 1
}

# ---------------------------- ARG PARSING ----------------------------
ORGANISM=""
THREADS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--organism) ORGANISM="$2"; shift 2 ;;
        -t|--threads)  THREADS="$2"; shift 2 ;;
        -h|--help)     usage ;;
        *) echo "❌ Unknown option: $1"; usage ;;
    esac
done

[[ -z "$ORGANISM" || -z "$THREADS" ]] && usage

if [[ -z "${species_dict_size[$ORGANISM]:-}" ]]; then
    echo "❌ Unknown organism: $ORGANISM"
    usage
fi

species="${species_dict_read_mlst[$ORGANISM]}"
species_tree_dir="${TREE_BASE}/${species}_trees"

if [[ ! -d "$species_tree_dir" ]]; then
    echo "❌ Tree folder does not exist: $species_tree_dir"
    exit 1
fi

mkdir -p "${species_tree_dir}/logs"
LOGFILE="${species_tree_dir}/logs/pipeline.log"
exec > >(tee -a "$LOGFILE") 2>&1

log "========== Organism-wide tree: $ORGANISM ($species) =========="
log "Tree directory: $species_tree_dir"

# ---------------------------- STEP 0: Collate per-ST CSVs ----------------------------
step="collate_csvs"; step_start=$(date +%s)
combined_csv="${species_tree_dir}/sample_collection_dates.csv"
tmp_csv=$(mktemp)
echo "sample,collection_date,st" > "$tmp_csv"

shopt -s nullglob
ST_DIRS=("${species_tree_dir}"/ST*/)
if [[ "${#ST_DIRS[@]}" -eq 0 ]]; then
    log "❌ No ST sub-folders found under $species_tree_dir"
    rm -f "$tmp_csv"
    exit 1
fi

for st_dir in "${ST_DIRS[@]}"; do
    st_name=$(basename "$st_dir")   # e.g. ST5
    st_num="${st_name#ST}"
    st_csv="${st_dir}sample_collection_dates.csv"
    if [[ ! -f "$st_csv" ]]; then
        log "⚠️ No sample_collection_dates.csv in $st_dir - skipping"
        continue
    fi
    tail -n +2 "$st_csv" | while IFS=, read -r sample_name collection_date; do
        [[ -z "$sample_name" ]] && continue
        echo "${sample_name},${collection_date},${st_num}" >> "$tmp_csv"
    done
done
mv "$tmp_csv" "$combined_csv"
step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s - combined CSV: $combined_csv"

# ---------------------------- LOAD SAMPLES ----------------------------
SAMPLES=(); declare -A COLLECTION_OF; declare -A ST_OF
while IFS=, read -r sample_name collection_date st_num; do
    [[ "$sample_name" == "sample" || -z "$sample_name" ]] && continue
    fasta="${species_tree_dir}/ST${st_num}/${sample_name}.fasta"
    fastq="${species_tree_dir}/ST${st_num}/${sample_name}.fastq.gz"
    if [[ ! -f "$fasta" || ! -f "$fastq" ]]; then
        log "⚠️ Missing fasta/fastq for $sample_name in ST${st_num} - skipping"
        continue
    fi
    SAMPLES+=("$sample_name")
    COLLECTION_OF["$sample_name"]="$collection_date"
    ST_OF["$sample_name"]="$st_num"
done < "$combined_csv"

if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
    log "❌ No usable samples found for $ORGANISM"
    exit 1
fi

log "Found ${#SAMPLES[@]} sample(s) for $ORGANISM across ${#ST_DIRS[@]} ST folder(s): ${SAMPLES[*]}"

# ---------------------------- STEP 1: Pick reference ----------------------------
step="tree_reference_pick"; step_start=$(date +%s)
REF_CIRCULAR=false
REF_SAMPLE=""

while IFS=, read -r sample_name collection_date st_num; do
    [[ "$sample_name" == "sample" || -z "$sample_name" ]] && continue
    fasta="${species_tree_dir}/ST${st_num}/${sample_name}.fasta"
    [[ ! -f "$fasta" ]] && continue

    if head -n 1 "$fasta" | grep -q "circular=true"; then
        REF_SAMPLE="$sample_name"
        REF_ST="$st_num"
        REF_CIRCULAR=true
        break
    fi

    if [[ -z "$REF_SAMPLE" ]]; then
        REF_SAMPLE="$sample_name"
        REF_ST="$st_num"
    fi
done < <(tail -n +2 "$combined_csv" | sort -t, -k2,2n)

if [[ -z "$REF_SAMPLE" ]]; then
    log "❌ No samples available to use as reference"
    exit 1
fi

grep -A1 "^>1 " "${species_tree_dir}/ST${REF_ST}/${REF_SAMPLE}.fasta" > "${species_tree_dir}/ref_full_genome.fasta"
REF_FASTA="${species_tree_dir}/ref_full_genome.fasta"

if $REF_CIRCULAR; then
    log "$REF_SAMPLE (ST${REF_ST}) chosen as circular reference sample for mapping as root of tree."
else
    log "⚠️ No circular complete genomes found. $REF_SAMPLE (ST${REF_ST}, oldest sample) chosen as non-circular reference for mapping as root of tree."
fi

if ! samtools faidx "$REF_FASTA"; then
    log "❌ Step $step failed"
    exit 1
fi
step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

# ---------------------------- STEP 1b: Detect reference change since last run ----------------------------
ref_marker_file="${species_tree_dir}/.reference_used"
if [[ -f "$ref_marker_file" ]]; then
    previous_ref=$(cat "$ref_marker_file")
    if [[ "$previous_ref" != "$REF_SAMPLE" ]]; then
        log "⚠️ Reference changed since last run (was: $previous_ref, now: $REF_SAMPLE)."
        log "⚠️ All existing per-sample consensus fastas are against the OLD reference - purging them to force a full re-run."
        for sample in "${SAMPLES[@]}"; do
            rm -vf "${species_tree_dir}/${sample}_consensus.fasta"
        done
        # Downstream combined outputs are now stale too, regardless of per-sample state
        rm -vf "${species_tree_dir}/alignment.fasta" \
               "${species_tree_dir}/${species}"*.tre \
               "${species_tree_dir}/${species}"*.tree \
               "${species_tree_dir}/${species}.recombination_predictions.gff" \
               "${species_tree_dir}/${species}.filtered_polymorphic_sites.fasta" \
               "${species_tree_dir}/${species}.filtered_polymorphic_sites_core.fasta" \
               "${species_tree_dir}/${species}.snp_distances.tsv" \
               "${species_tree_dir}/${species}.node_labelled.final_tree.png"
    else
        log "Reference unchanged since last run ($REF_SAMPLE) - existing consensus fastas will be reused."
    fi
else
    log "No previous reference record found (first run for this organism) - proceeding fresh."
fi

# ---------------------------- STEP 2: Map / call / consensus per sample (serial, full threads each) ----------------------------
step="sample_processing"; step_start=$(date +%s)
for sample in "${SAMPLES[@]}"; do
    st="${ST_OF[$sample]}"
    sample_fastq="${species_tree_dir}/ST${st}/${sample}.fastq.gz"
    consensus_fasta="${species_tree_dir}/${sample}_consensus.fasta"
    bam="${species_tree_dir}/${sample}.bam"

    if [[ -f "$consensus_fasta" ]]; then
        log "↷ $sample already has a consensus fasta - skipping"
        continue
    fi

    log "-------- Processing sample: $sample (ST${st}) --------"

    # ---- Mapping ----
    log "Mapping $sample to $REF_SAMPLE"
    if ! minimap2 -x lr:hq -a -t "$THREADS" "$REF_FASTA" "$sample_fastq" \
        | samtools view -h -F 4 \
        | samtools sort -O BAM -o "$bam"; then
        log "❌ Mapping failed for $sample - skipping sample"
        continue
    fi
    samtools index "$bam"

    # ---- Variant calling ----
    log "Calling variants for $sample against $REF_SAMPLE"
    if ! conda run -n "$ENV_CLAIR3" run_clair3.sh \
            --bam_fn="$bam" \
            --ref_fn="$REF_FASTA" \
            --threads="$THREADS" \
            --platform="ont" \
            --model_path="$CLAIR3_MODEL" \
            --output="${species_tree_dir}/clair3/${sample}" \
            --sample_name="$sample" \
            --min_coverage="10" \
            --include_all_ctgs \
            --haploid_precise \
            --no_phasing_for_fa \
            --enable_long_indel; then
        log "❌ Variant calling failed for $sample - skipping sample"
        continue
    fi
    mv -v "${species_tree_dir}/clair3/${sample}/merge_output.vcf.gz" "${species_tree_dir}/${sample}_full.vcf.gz"
    bcftools norm -a -m - "${species_tree_dir}/${sample}_full.vcf.gz" |
    bcftools norm -a -d "none" |
    bcftools view -v "snps" -f "PASS" -O "z" |
    bcftools filter -i "QUAL>=50" -O "z" > "${species_tree_dir}/${sample}_filtered.vcf.gz"
    bcftools index -f "${species_tree_dir}/${sample}_filtered.vcf.gz"

    # ---- Consensus ----
    log "Generating consensus sequence for $sample"
    bedtools genomecov -ibam "$bam" -bga | awk '$4 < 10' > "${species_tree_dir}/${sample}.bed"
    if ! bcftools consensus -H A -f "$REF_FASTA" -m "${species_tree_dir}/${sample}.bed" "${species_tree_dir}/${sample}_filtered.vcf.gz" > "$consensus_fasta"; then
        log "❌ Consensus generation failed for $sample"
        continue
    fi
    sed -i "s/>.*$/>${sample}/" "$consensus_fasta"

    # Only the consensus fasta is needed downstream (alignment.fasta) - clean up
    # everything else for this sample now to save space.
    log "Removing intermediate files for $sample"
    rm -vf "$bam" "${bam}.bai" \
           "${species_tree_dir}/${sample}_full.vcf.gz" \
           "${species_tree_dir}/${sample}_filtered.vcf.gz" \
           "${species_tree_dir}/${sample}_filtered.vcf.gz.csi" \
           "${species_tree_dir}/${sample}.bed"
    rm -vrf "${species_tree_dir}/clair3/${sample}"

    log "✅ $sample done"
done
step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

# Mapping/variant-calling/consensus for all samples completed against $REF_SAMPLE -
# record this now, since recombination/tree building below always rerun regardless
# and shouldn't gate whether this reference is considered "used".
echo "$REF_SAMPLE" > "${species_tree_dir}/.reference_used"

cat "${species_tree_dir}"/*_consensus.fasta > "${species_tree_dir}/alignment.fasta"

# ---------------------------- STEP 3: Recombination filtering + tree ----------------------------
step="recombination"; step_start=$(date +%s)
rm -vrf ${PWD}/tmp*
if ! conda run -n $ENV_GUBBINS run_gubbins.py --threads "$THREADS" \
        -p "${species_tree_dir}/${species}" \
        --tree-builder "iqtree" \
        -o "$REF_SAMPLE" \
        --best-model \
        -v \
        -f "90" \
        -m "2" \
        --p-value "0.5" \
        --trimming-ratio "1.5" \
        "${species_tree_dir}/alignment.fasta"; then
    log "❌ Step $step failed"
    exit 1
fi
coresnpfilter -c "1.0" "${species_tree_dir}/${species}.filtered_polymorphic_sites.fasta" > "${species_tree_dir}/${species}.filtered_polymorphic_sites_core.fasta"
snp-dists -j 1 "${species_tree_dir}/${species}.filtered_polymorphic_sites_core.fasta" > "${species_tree_dir}/${species}.snp_distances.tsv"
step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

# ---------------------------- STEP 4: Tree image ----------------------------
step="tree"; step_start=$(date +%s)
if ! conda run -n $ENV_GUBBINS plot_gubbins.R -t "${species_tree_dir}/${species}.node_labelled.final_tree.tre" \
    -r "${species_tree_dir}/${species}.recombination_predictions.gff" \
    -o "${species_tree_dir}/${species}.node_labelled.final_tree.png" \
    --taxon-label-size "2" \
    --tree-width "10" \
    --show-taxa \
    --tree-axis-expansion "100"; then
    log "❌ Step $step failed"
    exit 1
fi
step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

log "Final tree image: ${species_tree_dir}/${species}.node_labelled.final_tree.png"
log "Final SNP matrix: ${species_tree_dir}/${species}.snp_distances.tsv"
log "✅ Organism-wide phylogenetic analysis complete for $ORGANISM ($species)"