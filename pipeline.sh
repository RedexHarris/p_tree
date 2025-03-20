#!/bin/bash

## To use this script for tree construction, the following conditions must be met:
## 1. Create a miniconda environment named 'biotools' and install the required tools(parallel, fastp, spades, quast, snippy, gubbins)
##    in the environment. The best way is to use my dockerfile to build a container and then run the process within the container.
## 2. Ensure that the 'input' folder contains .fastq.gz files.
## 3. Ensure that the 'reference' folder contains a reference genome .fna file, renamed as 'reference.fasta'.
## 4. The naming and paths mentioned above can be customized as needed.


# Work options
your_input_file_type_is="fastq"

###################################### Please do not modify any of the following code ! #######################################
# Setup the working directory
INPUT_DIR="/data/input"
OUTPUT_DIR="/data/output"
REFERENCE_DIR="/data/reference"

# Activate conda environment
export PATH="~/miniconda/bin:$PATH"
source activate biotools38

if [[ "$your_input_file_type_is" == "fastq" ]]
then
    # Determine whether input data is paired-end sequencing and create a datalist for each .fastq.gz file
    cd ${INPUT_DIR}
    if [ -s *_R2.fastq.gz ]
    then
        # Run fastp on 6 sequencing at the same time for quality control
        echo "####################### Running fastp for quality control #######################"
        mkdir -p ${OUTPUT_DIR}/fastp
        for r1seq in *_R1.fastq.gz
        do
            keyword=$(echo ${r1seq%%_*})
            r2seq=${keyword}"_R2.fastq.gz"
            echo "fastp -i $r1seq -I $r2seq -o ${OUTPUT_DIR}/fastp/clean_$r1seq -O ${OUTPUT_DIR}/fastp/clean_$r2seq -j ${OUTPUT_DIR}/fastp/$keyword.json -h ${OUTPUT_DIR}/fastp/$keyword.html" >> fastp.cmd.list
        done
        parallel -j 6 < fastp.cmd.list
        cd ${OUTPUT_DIR}/fastp

        # Run spades on each clean .fastq.gz file for genome assembly
        echo "###################### Running SPAdes for genome assembly #######################"
        mkdir -p ${OUTPUT_DIR}/spades
        for r1seq in clean_*_R1.fastq.gz
        do
            keyword=$(echo ${r1seq} | awk -F'_' '{print $2}')
            r2seq="clean_"${keyword}"_R2.fastq.gz"
            echo "spades.py -1 $r1seq -2 $r2seq -o ${OUTPUT_DIR}/spades/$keyword -t 4 --careful" >> spades.cmd.list
        done
        parallel -j 1 < spades.cmd.list

        # Run fastp and spades for single-end sequencing
    else
        echo "####################### Running fastp for quality control #######################"
        mkdir -p ${OUTPUT_DIR}/fastp
        for seq in *.fastq.gz
        do
            echo "fastp -i $seq -o ${OUTPUT_DIR}/fastp/clean_$seq -j ${OUTPUT_DIR}/fastp/$seq.json -h ${OUTPUT_DIR}/fastp/$seq.html" >> fastp.cmd.list
        done
        parallel -j 6 < fastp.cmd.list
        cd ${OUTPUT_DIR}/fastp

        echo "###################### Running SPAdes for genome assembly #######################"
        mkdir -p ${OUTPUT_DIR}/spades
        for seq in clean_*.fastq.gz
        do
            keyword=$(echo ${r1seq} | awk -F'_.' '{print $2}')
            echo "spades.py -s $seq -o ${OUTPUT_DIR}/spades/$keyword -t 16" >> spades.cmd.list
        done
        parallel -j 3 < spades.cmd.list
    fi
    wait

    # Collect scaffolds.fasta in each folder
    mkdir -p ${OUTPUT_DIR}/fasta
    cd ${OUTPUT_DIR}/spades
    for keyword in $(ls)
    do
        cp ${keyword}/scaffolds.fasta ${OUTPUT_DIR}/fasta/${keyword}.fasta
    done
    FASTA_DIR="/data/output/fasta"

elif [[ "$your_input_file_type_is" == "fasta" ]]
then
    FASTA_DIR="/data/input"
fi

cd ${FASTA_DIR}
ls *.fasta > fasta.list
conda activate biotools37

# Run quast for assembly quality evaluate
echo "################## Running Quast for assembly quality evaluate ##################"
mkdir -p ${OUTPUT_DIR}/quast
ls *.fasta > fasta.list
mapfile -t fastas < fasta.list
quast.py "${fastas[@]}" -r ${REFERENCE_DIR}/reference.fasta -o ${OUTPUT_DIR}/quast -t 12

: << EOF
####################### This code uses refseq_masher to match the best reference genome and download it #######################

echo "Running refseq_masher & Extracting the top match..."
mkdir -p ${OUTPUT_DIR}/refseq_masher
refseq_masher matches ${OUTPUT_DIR}/spades/contigs.fasta \
      -o ${OUTPUT_DIR}/refseq_masher/matches.tsv;
REF_GENOME=$(awk 'NR==2 {print $1}' ${OUTPUT_DIR}/refseq_masher/matches.tsv)

echo "Downloading reference genome..."
datasets download genome accession ${REF_GENOME} --filename ${OUTPUT_DIR}/refseq_masher/${REF_GENOME}.zip
unzip -o ${OUTPUT_DIR}/refseq_masher/${REF_GENOME}.zip \
      -d ${OUTPUT_DIR}/refseq_masher/${REF_GENOME}
mv ${OUTPUT_DIR}/refseq_masher/${REF_GENOME}/ncbi_dataset/data/${REF_GENOME}/*.fna \
   ${OUTPUT_DIR}/reference.fasta

## The above code currently has some issues: it fails to capture the best-matching reference genome squence id from the .tsv ##
## result file and connot perform batch matching to identify the best match. It is recommended to avoid using it. #############


# Run snippy for variant calling analysis
echo "######################## Running snippy for SNP calling #########################"
for contigs in $(< fasta.list)
do
    keyword=$(echo ${contigs%%.*})
    echo "snippy --cpus 4 --contigs $contigs --ref ${REFERENCE_DIR}/reference.fasta --outdir ${OUTPUT_DIR}/snippy/$keyword" >> snippy.cmd.list
done
parallel -j 3 < snippy.cmd.list
mapfile -t lines < fasta.list
declare -a args=()
for line in "${lines[@]}"
do
    args+=("${line%.fasta}")
done
cd ${OUTPUT_DIR}/snippy
snippy-core --prefix core --ref ${REFERENCE_DIR}/reference.fasta "${args[@]}"
snippy-clean_full_aln core.full.aln > clean.full.aln
EOF

# Run snippy for SNP calling
echo "######################## Running snippy for SNP calling #########################"
mkdir -p ${OUTPUT_DIR}/snippy
for file in $(< fasta.list)
do
    filename=$(echo ${file%%.*})
    echo -e "$filename\t$FASTA_DIR/$file" >> snippy_run.tab
done
snippy-multi snippy_run.tab --cpus 4 --ref ${REFERENCE_DIR}/reference.fasta > snippy_run.sh
mv snippy_run.sh ${OUTPUT_DIR}/snippy
cd ${OUTPUT_DIR}/snippy
sh ./snippy_run.sh
snippy-clean_full_aln core.full.aln > clean.full.aln

# Run gubbins for phylogenetic tree construction
echo "############## Running gubbins for phylogenetic tree construction ###############"
cd ${OUTPUT_DIR}
mkdir gubbins
run_gubbins.py --prefix gubbins ${OUTPUT_DIR}/snippy/clean.full.aln
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln
FastTree -gtr -nt clean.core.aln > clean.core.tre

# Organize intermediate files
cd ${OUTPUT_DIR}
if [-d "/fasta/"]
then
    cd ${OUTPUT_DIR}/fasta && rm *.list && cd ..
fi
cd fastp && rm *.list && cd ..
mv gubbins.* ${OUTPUT_DIR}/gubbins/ && mv clean.* ${OUTPUT_DIR}/gubbins/
cd ${INPUT_DIR} && rm *.list && rm *.tab && mv fastp.html ${OUTPUT_DIR}/fastp/

echo "############################## Pipeline completed ###############################"