#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Process RNAseq data for ONNV-infected mosquitoes
# ----------------------------------------------------------------------------------------

# --- Download reads

mkdir data
cd data/
wget -c -i ../URL.txt
cd ..

# --- Clone repository

git clone https://github.com/cahende/RNAseq.git

mkdir -p RNAseq/data/rawData/

# --- Rename weirdly named files

cd data

mv TX1D2Mplus_S4_L001_R1_001.fastq.gz TX1_D2_Mplus_S4_L001_R1_001.fastq.gz
mv TX1D2Mplus_S4_L001_R2_001.fastq.gz TX1_D2_Mplus_S4_L001_R2_001.fastq.gz
mv TX1D2Mplus_S4_L002_R1_001.fastq.gz TX1_D2_Mplus_S4_L002_R1_001.fastq.gz
mv TX1D2Mplus_S4_L002_R2_001.fastq.gz TX1_D2_Mplus_S4_L002_R2_001.fastq.gz

mv TX1_D7Fplus_S14_L001_R1_001.fastq.gz TX1_D7_Fplus_S14_L001_R1_001.fastq.gz
mv TX1_D7Fplus_S14_L001_R2_001.fastq.gz TX1_D7_Fplus_S14_L001_R2_001.fastq.gz
mv TX1_D7Fplus_S14_L002_R1_001.fastq.gz TX1_D7_Fplus_S14_L002_R1_001.fastq.gz
mv TX1_D7Fplus_S14_L002_R2_001.fastq.gz TX1_D7_Fplus_S14_L002_R2_001.fastq.gz

cd ..

# --- Combine reads from same sample

# Format should be:
# RNAseq/data/rawData/{{sample}}_{read}.fastq.gz

cd data

RAW_PATH=../RNAseq/data/rawData/

for REP in TX1 TX2 TX3; do

    for D in 2 7; do

        for SEX in M F; do

            for INF in "plus" ""; do

                IND=${REP}_D${D}_${SEX}${INF}_

                NEWINF=$INF
                if [[ "$NEWINF" == "" ]]; then
                    NEWINF=neg
                fi

                NEW_IND=${REP}_D${D}_${SEX}_${NEWINF}

                echo "Copying $NEW_IND..."

                for READ in R1 R2; do
                    cat $IND*${READ}* > $RAW_PATH/${NEW_IND}_${READ}.fastq.gz
                done
            done
        done
    done
done

ls $RAW_PATH/*R1* | sed -e "s/.*\///" -e "s/_R.*//" > \
    ../RNAseq/data/individual_list.txt

cd ../RNAseq

# --- General variables

IND_FILE=data/individual_list.txt
QUEUE=open

# --- Trim adapter sequences and low quality bases from reads

module load trimmomatic

for IND in `cat $IND_FILE`; do

    CMD="cd `pwd`; module load trimmomatic; \
        java -jar $TRIMMOMATIC PE -phred33 \
            -threads 8 \
            data/rawData/${IND}_R1.fastq.gz \
            data/rawData/${IND}_R2.fastq.gz \
            data/${IND}_R1_paired.fastq.gz \
            data/${IND}_R1_unpaired.fastq.gz \
            data/${IND}_R2_paired.fastq.gz \
            data/${IND}_R2_unpaired.fastq.gz \
            ILLUMINACLIP:`dirname $TRIMMOMATIC`/adapters/TruSeq3-PE.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

    echo $CMD | qsub -V -A $QUEUE -M cxb585@psu.edu -m abe \
        -l nodes=1:ppn=8,walltime=12:00:00,mem=8gb -N t_$IND

done

# --- Download reference genome and GTF (An. stephensi) and make index for STAR aligner

VECTOR_BASE=https://www.vectorbase.org/download

mkdir -p genomes
# wget -O genomes/AsteI2.fa.gz \
#     $VECTOR_BASE/anopheles-stephensi-indianscaffoldsastei2fagz
# gunzip -c genomes/AsteI2.fa.gz > genomes/AsteI2.fa
#
# wget -O genomes/AsteI2.gtf.gz \
#     $VECTOR_BASE/anopheles-stephensi-indianbasefeaturesastei23gtfgz
# gunzip -c genomes/AsteI2.gtf.gz > genomes/AsteI2.gtf

wget -O genomes/AgamP4.fa.gz \
    $VECTOR_BASE/anopheles-gambiae-pestchromosomesagamp4fagz
gunzip -c genomes/AgamP4.fa.gz > genomes/AgamP4.fa

wget -O genomes/AgamP4.gtf.gz \
    $VECTOR_BASE/anopheles-gambiae-pestbasefeaturesagamp412gtfgz
gunzip -c genomes/AgamP4.gtf.gz > genomes/AgamP4.gtf

module load star

CMD="cd `pwd`; module load star; \
    STAR \
        --runMode genomeGenerate \
        --genomeDir genomes \
        --genomeFastaFiles genomes/AgamP4.fa \
        --sjdbGTFfile genomes/AgamP4.gtf \
        --sjdbOverhang 74 \
        --genomeSAindexNbases 13"

QUEUE=open
echo $CMD | qsub -V -A $QUEUE -M cxb585@psu.edu -m abe \
    -l nodes=1:ppn=8,walltime=12:00:00,mem=8gb -N idx

# --- Align reads to reference genome

module load star

for IND in `cat $IND_FILE`; do

    mkdir -p results/$IND

    CMD="cd `pwd`; module load star; \
        for R in R1 R2; do \
            gunzip -c data/${IND}_\${R}_paired.fastq.gz > \
                data/${IND}_\${R}_paired.fastq; \
        done; \
        STAR \
            --genomeDir genomes \
            --runThreadN 8 \
            --readFilesIn \
                data/${IND}_R1_paired.fastq \
                data/${IND}_R2_paired.fastq \
            --sjdbGTFfile genomes/AgamP4.gtf \
            --sjdbOverhang 74 \
            --outFileNamePrefix results/$IND/$IND"

    QUEUE=open
    echo $CMD | qsub -V -A $QUEUE -M cxb585@psu.edu -m abe \
        -l nodes=1:ppn=8,walltime=12:00:00,mem=8gb -N a_$IND

done

# --- Sort and index BAM files

module load samtools

for IND in `cat $IND_FILE`; do

    CMD="cd `pwd`; module load samtools; \
        samtools view -Sb results/$IND/${IND}Aligned.out.sam | \
            samtools sort -o results/$IND/${IND}.sort.bam; \
        samtools index results/$IND/${IND}.sort.bam \
            results/$IND/${IND}.sort.bam.bai"

    echo $CMD | qsub -V -A $QUEUE -M cxb585@psu.edu -m abe \
        -l nodes=1:ppn=1,walltime=4:00:00,mem=8gb -N b$IND

done

# --- Fix mate pairs

module load picard

for IND in `cat $IND_FILE`; do

    CMD="cd `pwd`; module load picard; \
        java -jar $PICARD/picard.jar FixMateInformation \
            INPUT=results/$IND/${IND}.sort.bam \
            OUTPUT=results/$IND/${IND}.sort.fix.bam \
            SORT_ORDER=coordinate \
            VALIDATION_STRINGENCY=LENIENT \
            CREATE_INDEX=true"

    echo $CMD | qsub -V -A $QUEUE -M cxb585@psu.edu -m abe \
        -l nodes=1:ppn=1,walltime=4:00:00,mem=8gb -N m$IND

done

# --- Filter for mapped and paired reads

module load bamtools

for IND in `cat $IND_FILE`; do

    CMD="cd `pwd`; module load bamtools; \
        bamtools filter -isMapped true \
            -in results/$IND/${IND}.sort.fix.bam \
            -out results/$IND/${IND}.sort.fix.map.bam"

    echo $CMD | qsub -V -A $QUEUE -M cxb585@psu.edu -m abe \
        -l nodes=1:ppn=1,walltime=12:00:00,mem=8gb -N f$IND

done

# --- Remove duplicate reads

module load picard

for IND in `cat $IND_FILE`; do

    CMD="cd `pwd`; module load picard; \
        java -jar $PICARD/picard.jar MarkDuplicates \
            INPUT=results/$IND/${IND}.sort.fix.map.bam \
            OUTPUT=results/$IND/${IND}.sort.fix.map.dup.bam \
            VALIDATION_STRINGENCY=SILENT \
            REMOVE_DUPLICATES=true \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 \
            METRICS_FILE=results/$IND/${IND}.duplog.txt"

    echo $CMD | qsub -V -A $QUEUE -M cxb585@psu.edu -m abe \
        -l nodes=1:ppn=1,walltime=12:00:00,mem=8gb -N d$IND

done

# --- Filter for mapping quality

module load bamtools

for IND in `cat $IND_FILE`; do

    CMD="cd `pwd`; module load bamtools; \
        bamtools filter \
            -mapQuality '>=20' \
            -length '>=50' \
            -in results/$IND/${IND}.sort.fix.map.dup.bam \
            -out results/$IND/${IND}.final.bam"

    echo $CMD | qsub -V -A $QUEUE -M cxb585@psu.edu -m abe \
        -l nodes=1:ppn=1,walltime=12:00:00,mem=8gb -N q$IND

done

mkdir -p final_bams
cp results/*/*.final.bam final_bams

# ----------------------------------------------------------------------------------------
# --- Backup to AWS
# ----------------------------------------------------------------------------------------

aws s3api create-bucket --bucket onnv-anopheles-response-bam-files --region us-east-1
aws s3 sync final_bams s3://onnv-anopheles-response-bam-files/final_bams
