#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Sequencing Pipeline v0.8 - Phylogenetic analysis stage - (Trees)
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

    # Set sample_dir to the latest available timepoint
    sample_dir=""
    for timepoint in 72h 48h 24h; do
        candidate="/data/IsoSeq_results/${sample}/${timepoint}"
        if [[ -d "$candidate" ]]; then
            sample_dir="$candidate"
            break
        fi
    done

    if [[ -z "$sample_dir" ]]; then
        log "❌ No timepoint directory found for $sample (checked 72h, 48h, 24h) - skipping"
        continue
    fi

    mkdir -p "${sample_dir}/logs"
    LOGFILE="${sample_dir}/logs/pipeline.log"

    log() { echo -e "[$(date '+%F %T')] $*"; }
    exec 3>&1 4>&2
    exec > >(tee -a "$LOGFILE") 2>&1

    if [[ -f "${sample_dir}/mlst_result.tsv" ]]; then
        log "Assembly complete for $sample - proceed"
        else 
        log "Assembly failed for $sample - sample skipped"; continue
    fi

    log "========== Processing sample: $sample ($organism) =========="
    start_total=$(date +%s)
    status="OK"; failed_step="-"

    # ---------------- STEP 1: Copy genomes for ST specifc tree ----------------
    step="copy_genomes"; step_start=$(date +%s)
    log "▶ Step $step"
    species="${species_dict_read_mlst[$organism]}"
    species_tree_dir="${TREE_BASE}/${species}_trees"

    # Add a step that if contig MLST result is ND, then use krocus mlst result instead
    if [[ -f "${sample_dir}/mlst_result.tsv" ]]; then
        contig_result=$(awk 'NR==2{print $3; exit}' "${sample_dir}/mlst_result.tsv")
        if ! [[ "$contig_result" =~ ^[0-9]+$ ]]; then
            log "❌ Contig MLST result is ST \"ND\" so using krocus mlst result instead"
            sequence_type=$(awk '{print $1; exit}' "/data/IsoSeq_results/${sample}/8h/mlst_result.tsv")
            log "✅ MLST result has succesfully changed from ST${contig_result} to ST${sequence_type}"
            else
                log "✅ Contig result is ST${contig_result} and present so do not need to use krocus mlst result"
                sequence_type="$contig_result"
        fi
        else
        log "⚠️ Contig MLST result does not exist, skip sample."
        status="FAILED"
        failed_step="$step"
        echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"
        continue
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
    
    # ---------------- STEP 2: Picking Reference for ST specifc tree ----------------
    step="tree_reference_pick"; step_start=$(date +%s)
    REF_CIRCULAR=false

    # Sort by collection date (oldest first), prefer circular=true
    while IFS=, read -r sample_name collection_date; do
        [[ "$sample_name" == "sample" ]] && continue

        fasta="${mlst_tree_dir}/${sample_name}.fasta"
        [[ ! -f "$fasta" ]] && continue

        if head -n 1 "$fasta" | grep -q "circular=true"; then
            REF_SAMPLE="$sample_name"
            REF_CIRCULAR=true
            break
        fi

        # First non-circular fallback (oldest due to sort order)
        [[ -z "$REF_SAMPLE" ]] && REF_SAMPLE="$sample_name"

    done < <(tail -n +2 "$collection_file" | sort -t, -k2,2n)

    # Fail only if no samples found at all
    if [[ -z "$REF_SAMPLE" ]]; then
        log "❌ No samples available to use as reference"
        status="FAILED"
        failed_step="$step"
        echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"
        continue
    fi

    grep -A1 "^>1 " "${mlst_tree_dir}/${REF_SAMPLE}.fasta" > "${mlst_tree_dir}/ref_full_genome.fasta"
    REF_FASTA="${mlst_tree_dir}/ref_full_genome.fasta"

    if $REF_CIRCULAR; then
        log "$REF_SAMPLE chosen as circular reference sample for mapping as root of tree."
    else
        log "⚠️ No circular complete genomes found. $REF_SAMPLE (oldest sample) chosen as non-circular reference for mapping as root of tree."
    fi

    if ! samtools faidx "$REF_FASTA"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi

    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 3: Mapping against tree reference ----------------
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

    # ---------------- STEP 4: Variant calling with Clair3 ----------------
    step="variant_calling"; step_start=$(date +%s)
    log "Calling variants for $sample against $REF_SAMPLE and outputting VCF";
    if ! conda run -n $ENV_CLAIR3 run_clair3.sh \
            --bam_fn="${mlst_tree_dir}/${bam}" \
            --ref_fn="$REF_FASTA" \
            --threads=1 \
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
    bcftools view -v "snps" -f "PASS" -O "z" |
    bcftools filter -i "QUAL>=50" -O "z" > "${mlst_tree_dir}/${sample}_filtered.vcf.gz"
    bcftools index -f "${mlst_tree_dir}/${sample}_filtered.vcf.gz"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 5: Generate consensus and multi-alignment FASTA ----------------
    step="consensus"; step_start=$(date +%s)
    log "Generating consensus sequence for $sample";
    bedtools genomecov -ibam "${mlst_tree_dir}/${bam}" -bga | awk '$4 < 10' > "${mlst_tree_dir}/${sample}.bed" # Generate coverage less than 10 bed file
    if ! bcftools consensus -H A -f "$REF_FASTA" -m "${mlst_tree_dir}/${sample}.bed" "${mlst_tree_dir}/${sample}_filtered.vcf.gz" > "${mlst_tree_dir}/${sample}_consensus.fasta"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi # Generate consensus, masking low cov regions
    sed -i "s/>.*$/>${sample}/" "${mlst_tree_dir}/${sample}_consensus.fasta" # Change name of sequences in consensus fasta file to sample name
    cat ${mlst_tree_dir}/*_consensus.fasta > ${mlst_tree_dir}/alignment.fasta
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 6: Filter recombination and generate tree and SNP distance matrix ----------------
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
            --model-fitter "raxmlng" \
            --recon-model "GTR" \
            -v \
            -f "90" \
            -m "2" \
            --p-value "0.5" \
            --trimming-ratio "1.5" \
            ${mlst_tree_dir}/alignment.fasta; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    coresnpfilter -c "1.0" "${mlst_tree_dir}/ST${sequence_type}.filtered_polymorphic_sites.fasta" > "${mlst_tree_dir}/ST${sequence_type}.filtered_polymorphic_sites_core.fasta"
    snp-dists -j 1 "${mlst_tree_dir}/ST${sequence_type}.filtered_polymorphic_sites_core.fasta" > ${mlst_tree_dir}/ST${sequence_type}.snp_distances.tsv
    grep "$sample" ${mlst_tree_dir}/ST${sequence_type}.snp_distances.tsv > ${sample_dir}/snp_distances.tsv
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 7: Generate tree image ----------------
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
    status="Phylogenetic analysis COMPLETE"
    echo -e "${sample}\t${status}\t${failed_step}\t${total_time}" >> "${MASTER_SUMMARY}"
    log "✅ Sample $sample Phylogenetic analysis completed in ${total_time}s"
    exec 1>&3 2>&4
    exec 3>&- 4>&-
}; done
