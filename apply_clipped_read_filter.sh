#!/bin/bash


# Check if the input file argument is provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <inputdir> <outputfile>"
    exit 1
fi

inputfile="$1"
workdir="$2"
bamdir="$3"
cohort="$4"
src="/home/hyunwoo/src/python/svProject/dnsv_workflow/etching_trio/clipped_read_filter.py"

cd $workdir

# Check if the input directory exists
if [ ! -d "$workdir" ]; then
    echo "Error: Input directory not found."
    exit 1
fi

line_num=0
while IFS= read -r line; do
    ((line_num++))

    if [ $line_num -eq 1 ]; then
        continue
    fi

    family_id=$(echo "$line" | awk '{print $1}' | tr -d '\r')
    child_id=$(echo "$line" | awk '{print $2}' | tr -d '\r')
    paternal_id=$(echo "$line" | awk '{print $3}' | tr -d '\r')
    maternal_id=$(echo "$line" | awk '{print $4}' | tr -d '\r')

    if [ $cohort == "CEPH" ]; then
        child_regex="SRR[0-9]+_${child_id}_[^Child].+\.BQSR.Apply.bam$"
        paternal_regex="SRR[0-9]+_${paternal_id}_[^Child].+\.BQSR.Apply.bam$"
        maternal_regex="SRR[0-9]+_${maternal_id}_[^Child].+\.BQSR.Apply.bam$"
    elif [ $cohort == "ABC" ] || [ $cohort == 'CTC' ]; then 
        child_regex=".+_${child_id}.final.bam$"
        paternal_regex=".+_${paternal_id}.final.bam$"
        maternal_regex=".+_${maternal_id}.final.bam$"
        #child_bai_regex=".+_${child_id}.final.bam.bai$"
        #paternal_bai_regex=".+_${paternal_id}.final.bam.bai$"
        #maternal_bai_regex=".+_${maternal_id}.final.bam.bai$"
    elif [ $cohort == "RD" ]; then
        child_regex=".+_${child_id}*.bam$"
        paternal_regex=".+_${paternal_id}*.bam$"
        maternal_regex=".+_${maternal_id}*.bam$"
    fi

    #echo $child_regex

    child_bam=$(find "$bamdir" -type l | grep -E "$child_regex")
    paternal_bam=$(find "$bamdir" -type l | grep -E "$paternal_regex")
    maternal_bam=$(find "$bamdir" -type l | grep -E "$maternal_regex")


#    child_bai=$(find "$bamdir" -type l | grep -E "$child_bai_regex")
#    paternal_bai=$(find "$bamdir" -type l | grep -E "$paternal_bai_regex")
#    maternal_bai=$(find "$bamdir" -type l | grep -E "$maternal_bai_regex")
    child_bai="${child_bam}.bai"
    paternal_bai="${paternal_bam}.bai"
    maternal_bai="${paternal_bam}.bai"


    echo $family_id
    echo $child_bam 
    echo $paternal_bam 
    echo $maternal_bam

    if [ ! -f $child_bam ] || [ ! -f $paternal_bam ] || [ ! -f $maternal_bam ]; then
        echo "Error : input bam file not found"
        if [ ! -f $child_bam ]; then
            echo "${child_bam} NOT FOUND"
        fi
        if [ ! -f $paternal_bam ]; then
            echo "${paternal_bam} NOT FOUND"
        fi
        if [ ! -f $maternal_bam ]; then
            echo "${maternal_bam} NOT FOUND"
        fi
        #exit 1
        continue
    fi

    cd $bamdir
#    echo $child_bai
    if [ ! -f $child_bai ] || [ ! -f $paternal_bai ] || [ ! -f $maternal_bai ]; then
        echo "Error : bam file index not found"
        if [ ! -f $child_bai ]; then
            echo $child_bai not found... trying indexing
#            samtools index -@ 30 $child_bam
        fi
        if [ ! -f $paternal_bai ]; then
            echo $paternal_bai not found... trying indexing
#            samtools index -@ 30 $paternal_bam
        fi
        if [ ! -f $maternal_bai ]; then
            echo $maternal_bai not found... trying indexing
#            samtools index -@ 30 $maternal_bam
        fi
    fi

    cd $workdir

    family="family${family_id}"
    targetdir="${workdir}/${family}/etching_call_v2/subject${child_id}/"
#    targetdir="${workdir}/${family}/etching_call_v2/${child_id}/"
#    targetdir="${workdir}/${family}/etching_call/${child_id}/"

    targetfile="${targetdir}/${family}_subject${child_id}.scored.filtered.typed.contigFilt_50bp_PASS_nonLCR_nonIntraBND.addVAF.dnsv.nonBND.parentVCFfilter.vcf"
#    targetfile="${targetdir}/${family}_${child_id}.scored.filtered.typed.contigFilt_50bp_PASS_nonLCR_nonIntraBND.addVAF.dnsv.nonBND.parentVCFfilter.vcf"
#    targetfile="${targetdir}/${family}_${child_id}.scored.filtered.typed.contigFilt_50bp_PASS_nonLCR_nonBND.addVAF.dnsv.vcf"
    outputfile="${targetdir}/$(basename ${targetfile%.vcf}).addClipRead.vcf"
#    if [ ! -f $outputfile ]; then
#        cmd="python ${src} -i ${child_bam} -m ${maternal_bam} -p ${paternal_bam} -v ${targetfile} -o ${outputfile}"
    cmd="python ${src} -i ${child_bam} -m ${maternal_bam} -p ${paternal_bam} -v ${targetfile}"
    echo $cmd
    eval $cmd
#    fi
done < "$inputfile"


#for family in "$workdir"/*/; do
#    # Check if the path is a directory
#    if [ -d "$family" ]; then
#        echo "Processing subdirectory: $family"
#        targetdir="${family}/etching_call_v2/"
##        targetdir="${family}/etching_call/"
#        for subject in "$targetdir"/*/; do
#            if [ -d "$subject" ]; then
#                familynum=$(basename ${family})
#                subjectnum=$(basename ${subject})
#                #targetfile="${subject}/${familynum}_${subjectnum}.scored.filtered.typed.contigFilt_50bp_PASS_nonLCR_nonIntraBND.addVAF.dnsv.vcf"
#                if [[ "$subjectnum" == subject* ]]; then
#                    echo $subjectnum
##                    targetfile="${subject}/${familynum}_${subjectnum}.scored.filtered.typed.contigFilt_50bp_PASS_nonLCR_nonIntraBND.addVAF.dnsv.nonBND.parentVCFfilter.vcf"
#                    targetfile="${subject}/${familynum}_${subjectnum}.scored.filtered.typed.contigFilt_50bp_PASS_nonLCR_nonBND.addVAF.dnsv.nonBND.parentVCFfilter.vcf"
#                    outputfile="$(basename ${targetfile%.vcf}).clip_read_filter.vcf"
#
#                    child_bam=""
#                    maternal_bam=""
#                    paternal_bam=""
#
#                    echo $(basename $targetfile)
#                    echo $outputfile
#                    cmd="python ${src} -i ${} -m ${} -p ${} -v ${targetfile} -o ${outputfile}"
#                    #echo $cmd
#                    #eval $cmd
#                fi
#            fi
#        done
#    fi
#done
