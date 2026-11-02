#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=subsample
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=2G
#SBATCH -t 1:00:00
#SBATCH --output=logs/subsample_%j.txt
#SBATCH --error=logs/subsample_%j.err

module purge
module load seqkit

# Check if the correct number of arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: sbatch $0 <sample_name> <input_path> <transcriptomepath> <outputpath>"
    exit 1
fi

# Assign command line arguments to variables
sample_name="$1"
input_path="${2%/}"
transcriptomepath="${3%/}"
outputpath="${4%/}"

# Define subsampling depths (in number of reads)
subsample_depths=(150000000 200000000 300000000 400000000)
subsample_labels=("150M" "200M" "300M" "400M")

# Number of threads for seqkit
THREADS=4

# Find R1 files for this sample
r1_files=($(find "$input_path" -maxdepth 1 -name "${sample_name}*_R1_*.fastq.gz" | sort))

if [ ${#r1_files[@]} -eq 0 ]; then
    echo "Error: No R1 FASTQ files found matching sample name '$sample_name' in '$input_path'."
    exit 1
fi

echo "Found ${#r1_files[@]} R1 file(s) for sample '$sample_name':"
for file in "${r1_files[@]}"; do
    echo "  - $(basename "$file")"
done

# Create output directory
subsample_dir="${input_path}/subsample"
mkdir -p "$subsample_dir"

# Get the maximum depth needed
max_depth=${subsample_depths[-1]}
max_label=${subsample_labels[-1]}

# Process each R1/R2 pair
for r1_file in "${r1_files[@]}"; do
    # Find corresponding R2 file
    r2_file="${r1_file/_R1_/_R2_}"
    
    if [ ! -f "$r2_file" ]; then
        echo "ERROR: R2 file not found for $(basename "$r1_file")"
        echo "Expected: $(basename "$r2_file")"
        exit 1
    fi
    
    r1_basename=$(basename "$r1_file" .fastq.gz)
    r2_basename=$(basename "$r2_file" .fastq.gz)
    
    echo "========================================="
    echo "Processing paired files:"
    echo "  R1: $(basename "$r1_file")"
    echo "  R2: $(basename "$r2_file")"
    echo "========================================="
    
    # STEP 1: Count reads from R1 (R2 should have same count)
    echo "Step 1: Counting reads..."
    start_time=$(date +%s)
    total_reads=$(seqkit stats -j $THREADS "$r1_file" 2>/dev/null | grep -v "^file" | awk '{print $4}' | tr -d ',')
    end_time=$(date +%s)
    
    # Verify we got a valid number
    if [ -z "$total_reads" ] || ! [[ "$total_reads" =~ ^[0-9]+$ ]]; then
        echo "  ERROR: Failed to count reads. Got: '$total_reads'"
        echo "  Trying alternative method..."
        total_reads=$(zcat "$r1_file" 2>/dev/null | wc -l | awk '{print int($1/4)}')
        if [ -z "$total_reads" ] || ! [[ "$total_reads" =~ ^[0-9]+$ ]]; then
            echo "  ERROR: Alternative counting also failed. Exiting."
            exit 1
        fi
    fi
    
    echo "  Total reads: $total_reads (counted in $((end_time - start_time)) seconds)"
    
    # Calculate fraction for maximum depth
    if (( $(echo "$max_depth > $total_reads" | bc -l) )); then
        echo "  Warning: Max depth ($max_depth) exceeds total reads. Using all reads."
        max_fraction=1.0
    else
        max_fraction=$(echo "scale=10; $max_depth / $total_reads" | bc -l)
    fi
    
    # Validate fraction
    if [ -z "$max_fraction" ]; then
        echo "  ERROR: Failed to calculate fraction"
        exit 1
    fi
    
    echo "  Calculated fraction: $max_fraction"
    
    # STEP 2: Subsample R1 and R2 together to maximum depth
    r1_output_max="${subsample_dir}/${r1_basename/${sample_name}/${sample_name}_${max_label}}.fastq.gz"
    r2_output_max="${subsample_dir}/${r2_basename/${sample_name}/${sample_name}_${max_label}}.fastq.gz"
    
    echo "Step 2: Subsampling paired files to ${max_label} reads..."
    start_time=$(date +%s)
    
    # CRITICAL: Use seqkit with 2-pass mode to handle paired files
    # Process R1
    seqkit sample -p $max_fraction --rand-seed 100 --two-pass -j $THREADS "$r1_file" -o "$r1_output_max"
    # Process R2 with SAME seed to maintain pairing
    seqkit sample -p $max_fraction --rand-seed 100 --two-pass -j $THREADS "$r2_file" -o "$r2_output_max"
    
    end_time=$(date +%s)
    echo "  Completed in $((end_time - start_time)) seconds"
    
    # Verify both files were created
    if [ ! -f "$r1_output_max" ] || [ ! -f "$r2_output_max" ]; then
        echo "  ERROR: Failed to create subsampled files"
        exit 1
    fi
    
    # STEP 3: Downsample to lower depths
    echo "Step 3: Creating lower depth subsamples..."
    for ((i=${#subsample_depths[@]}-2; i>=0; i--)); do
        depth="${subsample_depths[$i]}"
        label="${subsample_labels[$i]}"
        
        echo "  Downsampling to ${label}..."
        start_time=$(date +%s)
        
        downsample_fraction=$(echo "scale=10; $depth / $max_depth" | bc -l)
        
        r1_output="${subsample_dir}/${r1_basename/${sample_name}/${sample_name}_${label}}.fastq.gz"
        r2_output="${subsample_dir}/${r2_basename/${sample_name}/${sample_name}_${label}}.fastq.gz"
        
        # Use same seed for both R1 and R2 to maintain pairing
        seed=$((100 + i))
        seqkit sample -p $downsample_fraction --rand-seed $seed --two-pass -j $THREADS "$r1_output_max" -o "$r1_output"
        seqkit sample -p $downsample_fraction --rand-seed $seed --two-pass -j $THREADS "$r2_output_max" -o "$r2_output"
        
        end_time=$(date +%s)
        echo "    Completed in $((end_time - start_time)) seconds"
        
        # Verify both files were created
        if [ ! -f "$r1_output" ] || [ ! -f "$r2_output" ]; then
            echo "    ERROR: Failed to create downsampled files"
            exit 1
        fi
    done
    
    echo "  All subsampling completed for this pair"
done

# STEP 4: Submit cellranger jobs
echo ""
echo "========================================="
echo "Submitting cellranger jobs..."
echo "========================================="
for label in "${subsample_labels[@]}"; do
    subsampled_name="${sample_name}_${label}"
    echo "Submitting job for ${subsampled_name}..."
    sbatch mkcounts.sh "$subsampled_name" "$subsample_dir" "$transcriptomepath" "$outputpath"
done

echo ""
echo "========================================="
echo "Complete!"
echo "========================================="
echo "Subsampled FASTQ files are in: $subsample_dir"
echo "End time: $(date)"
echo "Total job time: $SECONDS seconds"