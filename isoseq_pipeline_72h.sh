#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Sequencing Pipeline - 72h stage - (Assembly, QC, MLST, AMR, Trees)
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

    sample_dir="results/${sample}/72h"
    mkdir -p "${sample_dir}/logs"
    LOGFILE="${sample_dir}/logs/pipeline.log"

    mkdir -p "${sample_dir}/samplesheet"
    cp -v "$SAMPLESHEET" "${sample_dir}/samplesheet/samplesheet.csv"
    touch "${sample_dir}/samplesheet/${RUN_NAME}"

    log() { echo -e "[$(date '+%F %T')] $*"; }
    exec 3>&1 4>&2
    exec > >(tee -a "$LOGFILE") 2>&1

    if [[ -f "${sample}/8h/kraken2_top5_taxa.csv" ]]; then
        top_hit=$(awk -F, 'NR>1 && $1 {print $2; exit}' "${sample}/8h/kraken2_top5_taxa.csv")
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
    ontime --to 72h -o "${sample_dir}/${sample}.fastq.gz" "${sample_dir}/${sample}_raw.fastq.gz"  # Use ontime to filter for only reads for this step
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
    if ! bash scripts/autoautocycler.sh -o "$outdir/" -t "$THREADS" -s "$size" -a "flye raven" "${sample_dir}/trimmed/${sample}.fastq.gz"; then
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
    gfa="${sample_dir}/assemblies/${sample}.gfa"
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
    if ! conda run -n "$ENV_MLST" mlst --scheme "$scheme" "${sample_dir}/assemblies/${sample}.fasta" > "${sample_dir}/mlst_result.tsv"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 10: Copy genomes for ST specifc tree ----------------
    step="copy_genomes"; step_start=$(date +%s)
    log "▶ Step $step"
    species="${species_dict_read_mlst[$organism]}"
    species_tree_dir="${TREE_BASE}/${species}_trees"
    # Adds a step in MLST result exists
    if [[ -f "${sample_dir}/mlst_result.tsv" ]]; then
        log "✅ MLST result exists."
        else
        log "⚠️ MLST result does not exist, ST set to \"none\"."
        sequence_type=`echo "_NONE"`
    fi
    mlst_tree_dir="${species_tree_dir}/ST${sequence_type}"
    mkdir -p "$species_tree_dir"
    mkdir -p "$mlst_tree_dir"
    cp -v "${sample_dir}/assemblies/${sample}.fasta" "${mlst_tree_dir}/${sample}.fasta" || log "⚠️ copy fasta failed"
    cp -v "${sample_dir}/trimmed/${sample}.fastq.gz" "${mlst_tree_dir}/${sample}.fastq.gz" || log "⚠️ copy fastq failed"

    collection_file="${mlst_tree_dir}/sample_collection_dates.csv"

    # Create header if file does not exist
    if [[ ! -f "$collection_file" ]]; then
        echo "sample,collection_date" > "$collection_file"
    fi

    # Append current sample and collection date
    echo "${sample},${collection}" >> "$collection_file"

    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"
    
    # ---------------- STEP 11: Picking Reference for ST specifc tree ----------------
    step="tree_reference_pick"; step_start=$(date +%s)
    REF_SAMPLE=""

    # Sort by collection date (oldest first)
    while IFS=, read -r sample_name collection_date; do
        [[ "$sample_name" == "sample" ]] && continue

        fasta="${mlst_tree_dir}/${sample_name}.fasta"

        if [[ -f "$fasta" ]] && \
            grep ">1" "$fasta" | grep -q "circular=true"; then
            REF_SAMPLE="$sample_name"
            break
        fi
    done < <(tail -n +2 "$collection_file" | sort -t, -k2,2n)

    if [[ -z "$REF_SAMPLE" ]]; then
        log "❌ No circular complete genomes available for reference"
        status="FAILED"
        failed_step="$step"
        echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"
        continue
    fi

    log "$REF_SAMPLE chosen as circular reference sample for mapping as root of tree."
    
    grep -a1 ">1" "${mlst_tree_dir}/${REF_SAMPLE}.fasta" > "${mlst_tree_dir}/ref_full_genome.fasta"
    REF_FASTA="${mlst_tree_dir}/ref_full_genome.fasta"
    samtools faidx $REF_FASTA
    
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 12: Mapping against tree reference ----------------
    step="reference_mapping"; step_start=$(date +%s)
    bam="${sample}.bam"
    log "Mapping $sample to $REF_SAMPLE and outputting $bam";
    if ! minimap2 -x map-ont -a -t "$THREADS" "$REF_FASTA" "$mlst_tree_dir/${sample}.fastq.gz" \
        | samtools view -h -F 4 \
        | samtools sort -O BAM -o "${mlst_tree_dir}/${bam}"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    samtools index "${mlst_tree_dir}/${bam}"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 13: Variant calling with Clair3 ----------------
    step="variant_calling"; step_start=$(date +%s)
    log "Calling variants for $sample against $REF_SAMPLE and outputting VCF";
    if ! conda run -n $ENV_CLAIR3 run_clair3.sh \
            --bam_fn="${mlst_tree_dir}/${bam}" \
            --ref_fn="$REF_FASTA" \
            --threads="$THREADS" \
            --platform="ont" \
            --model_path="$CLAIR3_MODEL" \
            --output="${mlst_tree_dir}/clair3/${sample}" \
            --sample_name="$sample" \
            --min_coverage="10" \
            --include_all_ctgs \
            --haploid_precise \
            --no_phasing_for_fa \
            --enable_long_indel; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    mv -v "${mlst_tree_dir}/clair3/${sample}/merge_output.vcf.gz" "${mlst_tree_dir}/${sample}_full.vcf.gz"
    # Filter variants to split MNVs to SNVs and keep only varints that pass filter
    bcftools norm -a -m - "${mlst_tree_dir}/${sample}_full.vcf.gz" |
    bcftools norm -a -d "none" |
    bcftools view -v "snps" -f "PASS" -O "z" > "${mlst_tree_dir}/${sample}_filtered.vcf.gz"
    bcftools index -f "${mlst_tree_dir}/${sample}_filtered.vcf.gz"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 14: Generate consensus and multi-alignment FASTA ----------------
    step="consensus"; step_start=$(date +%s)
    log "Generating consensus sequence for $sample";
    bedtools genomecov -ibam "${mlst_tree_dir}/${bam}" -bga | awk '$4 < 10' > "${mlst_tree_dir}/${sample}.bed" # Generate coverage less than 10 bed file
    if ! bcftools consensus -f "$REF_FASTA" -m "${mlst_tree_dir}/${sample}.bed" "${mlst_tree_dir}/${sample}_filtered.vcf.gz" > "${mlst_tree_dir}/${sample}_consensus.fasta"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi # Generate consensus, masking low cov regions
    sed -i "s/>.*$/>${sample}/" "${mlst_tree_dir}/${sample}_consensus.fasta" # Change name of sequences in consensus fasta file to sample name
    cat ${mlst_tree_dir}/*_consensus.fasta > ${mlst_tree_dir}/alignment.fasta
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 15: Filter recombination and generate tree and SNP distance matrix ----------------
    step="recombination"; step_start=$(date +%s)
    # Remove intermediate folders that are made by a failed gubbins run
    rm -vrf ${PWD}/tmp*
    if ! conda run -n $ENV_GUBBINS run_gubbins.py --threads "$THREADS" \
            -p "${mlst_tree_dir}/ST${sequence_type}" \
            --first-tree-builder "iqtree-fast"  \
            --tree-builder "iqtree" \
            -o "$REF_SAMPLE" \
            --first-model "JC" \
            --model "GTRGAMMA" \
            -v \
            -f "90" \
            ${mlst_tree_dir}/alignment.fasta; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    snp-dists -j 20 "${mlst_tree_dir}/gubbins.filtered_polymorphic_sites.fasta" > snp_distances.tsv
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 16: Generate tree image ----------------
    step="tree"; step_start=$(date +%s)
    conda run -n $ENV_GUBBINS plot_gubbins.R -t "${mlst_tree_dir}/ST${sequence_type}.node_labelled.final_tree.tre" \
        -r "${mlst_tree_dir}/ST${sequence_type}.recombination_predictions.gff" \
        -o "${sample_dir}/ST${sequence_type}.node_labelled.final_tree.png" \
        --taxon-label-size "2" \
        --tree-width "10" \
        --show-taxa \
        --tree-axis-expansion "100"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- Stage complete ----------------
    end_total=$(date +%s)
    total_time=$((end_total-start_total))
    status="72h COMPLETE"
    echo -e "${sample}\t${status}\t${failed_step}\t${total_time}" >> "${MASTER_SUMMARY}"
    log "✅ Sample $sample 72h analysis completed in ${total_time}s"
    exec 1>&3 2>&4
    exec 3>&- 4>&-
}; done
