#!/bin/bash

# Check if the input file argument is provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Input file containing sample information
input_file="$1"
threads=30
genome="$2"
pgk="$3"

workdir="$2"

cd $workdir

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file not found."
    exit 1
fi

# Iterate over each line in the input file

# Initialize line number
line_num=0


while IFS= read -r line; do
    ((line_num++))

    if [ $line_num -eq 1 ]; then
        continue
    fi

    # Extract sample path and parental path
    family_id=$(echo "$line" | awk '{print $1}'| tr -d '\r')
    sample_id=$(echo "$line" | awk '{print $2}'| tr -d '\r')
    sample_path=$(echo "$line" | awk '{print $7}' | tr -d '\r')
    paternal_id=$(echo "$line" | awk '{print $3}'| tr -d '\r')
    paternal_path=$(echo "$line" | awk '{print $8}'| tr -d '\r')
    maternal_id=$(echo "$line" | awk '{print $4}'| tr -d '\r')
    maternal_path=$(echo "$line" | awk '{print $9}'| tr -d '\r')



    # Check if the sample directory exists

    cd $workdir
    # Perform actions on the target file
    echo "Processing target sample: ${sample_id}"
    
    # Check if the directory named after the sample ID exists in the current working directory
    cd "family${family_id}"

    targetdir="${workdir}/family${family_id}/1.preprocess/"
    ptrim1="family${family_id}_${paternal_id}_r1.trim.fastq.gz"
    ptrim2="family${family_id}_${paternal_id}_r2.trim.fastq.gz"

    mtrim1="family${family_id}_${maternal_id}_r1.trim.fastq.gz"
    mtrim2="family${family_id}_${maternal_id}_r2.trim.fastq.gz"

    ctrim1="family${family_id}_${sample_id}_r1.trim.fastq.gz"
    ctrim2="family${family_id}_${sample_id}_r2.trim.fastq.gz"

    #-- RUN ETCHING
    cd etching_call

    if [ ! -d parental_call ]; then
        mkdir parental_call
    fi

    cd parental_call

    if [ ! -d "paternal_${paternal_id}" ]; then
        mkdir "paternal_${paternal_id}"
    fi

    cd "paternal_${paternal_id}"

    outputfile="paternal_${paternal_id}"
    checkoutputfile="${outputfile}.scored.filtered.typed.vcf"
    outputdir="${workdir}/family${family_id}/etching_call/parental_call/paternal_${paternal_id}/"
    if [ ! -f $checkoutputfile ]; then
        petching_cmd="/home/hyunwoo/programs/etching/ETCHING_v1.4.2_release_20230421/bin/etching -1 ${workdir}/family${family_id}/1.preprocess/${ptrim1} -2 ${workdir}/family${family_id}/1.preprocess/${ptrim2} -g ${genome} -f ${pgk} -t 30 --keep-kmc -o ${outputfile} --output-dir ${outputdir} --work-dir ${outputdir} >& ${outputdir}/log.etching_run &"
        echo $petching_cmd
        eval $petching_cmd
    fi

    cd ../

    if [ ! -d  "maternal_${maternal_id}" ]; then
        mkdir "maternal_${maternal_id}"
    fi

    cd "maternal_${maternal_id}"

    outputfile="maternal_${maternal_id}"
    checkoutputfile="${outputfile}.scored.filtered.typed.vcf"
    outputdir="${workdir}/family${family_id}/etching_call/parental_call/maternal_${maternal_id}/"
    if [ ! -f $checkoutputfile ]; then
        metching_cmd="/home/hyunwoo/programs/etching/ETCHING_v1.4.2_release_20230421/bin/etching -1 ${workdir}/family${family_id}/1.preprocess/${mtrim1} -2 ${workdir}/family${family_id}/1.preprocess/${mtrim2} -g ${genome} -f ${pgk} -t 30 --keep-kmc -o ${outputfile} --output-dir ${outputdir} --work-dir ${outputdir} >& ${outputdir}/log.etching_run &"
        echo $metching_cmd
        eval $metching_cmd
    fi

    wait

done < "$input_file"
