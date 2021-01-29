#### Samples link https://www.ncbi.nlm.nih.gov/sra?term=SRP136364

#control	run assession
#SRX3841241	SRR6890763
#SRX3841256	SRR6890748
#treated	
#SRX3841235	SRR6890769
#SRX3841255	SRR6890749


#Step 1: check the Quality
#control1_left.fastq control1_right.fastq control2_left.fastq control2_right.fastq treated1_left.fastq treated1_right.fastq  treated2_left.fastq treated2_right.fastq


# step 2: trimming and cleaning using cutadapt

########  syntex for paired end : cutadapt -b <Adapter seq> -B <Adapter seq> -q 30,30 -m 20 -o <outputfile_left.fastqc> -p <outputfile_right.fastqc> input_left.fastq input_right.fastq ############################

#cutadapt -b AGATCGGAAGAG -B AGATCGGAAGAG -q 30,30 -m 20 -o trim_c1_left.fastq -p trim_c1_right.fastq control1_left.fastq control1_right.fastq
#cutadapt -b AGATCGGAAGAG -B AGATCGGAAGAG -q 30,30 -m 20 -o trim_c2_left.fastq -p trim_c2_right.fastq control2_left.fastq control2_right.fastq
#cutadapt -b AGATCGGAAGAG -B AGATCGGAAGAG -q 30,30 -m 20 -o trim_t1_left.fastq -p trim_t1_right.fastq treated1_left.fastq treated1_right.fastq
#cutadapt -b AGATCGGAAGAG -B AGATCGGAAGAG -q 30,30 -m 20 -o trim_t2_left.fastq -p trim_t2_right.fastq treated2_left.fastq treated2_right.fastq

#step3: Again check the Quality

#fastqc trim_c1_left.fastq trim_c1_right.fastq trim_c2_left.fastq trim_c2_right.fastq trim_t1_left.fastq trim_t1_right.fastq trim_t2_left.fastq trim_t2_right.fastq

####  step4: to combine the all left reads and all right reads as allreads.left.fastq and allreads.right.fastq for assembly

#cat trim_c1_left.fastq trim_c2_left.fastq trim_t1_left.fastq trim_t2_left.fastq > new_allreads.left.fastq

#cat trim_c1_right.fastq trim_c2_right.fastq trim_t1_right.fastq trim_t2_right.fastq >> new_allreads.right.fastq


#### step5: convert fastq to fasta 

#sed '/^@/!d;s//>/;N' allreads.left.fastq > allreads.left.fa
#sed '/^@/!d;s//>/;N' allreads.right.fastq > allreads.right.fa

#########    step6: Trinity denovo assembly ##################

#/home/Documents/NGS_tools/trinityrnaseq-Trinity-v2.6.6/Trinity --seqType fa \
#--left /home/Documents/local/Samarth/RNASeq_denovo_sampledata/allreads.left.fa \
#--right /home/Documents/local/Samarth/RNASeq_denovo_sampledata/allreads.right.fa \
#--SS_lib_type RF \
#--max_memory 1G \
#--CPU 1 \
#--output /home/Documents/local/Samarth/RNASeq_denovo_sampledata/TRINITY_outdir

##NOTE: assembly output dir must have "trinity" word in its name


#to check the Stast of Assembly:

#/home/Documents/NGS_tools/trinityrnaseq-Trinity-v2.6.6/util/TrinityStats.pl /home/Documents/local/Samarth/RNASeq_denovo_sampledata/TRINITY_outdir/Trinity.fasta


#### step7 : Abundance Estimate using RSEM method 


#control1

#/home/Documents/NGS_tools/trinityrnaseq-Trinity-v2.6.6/util/align_and_estimate_abundance.pl --seqType fq \
#--transcripts /home/Documents/local/Samarth/RNASeq_denovo_sampledata/TRINITY_outdir/Trinity.fasta \
#--left /home/Documents/local/Samarth/RNASeq_denovo_sampledata/trim_c1_left.fastq \
#--right /home/Documents/local/Samarth/RNASeq_denovo_sampledata/trim_c1_right.fastq \
#--SS_lib_type RF \
#--est_method RSEM \
#--aln_method bowtie \
#--trinity_mode \
#--prep_reference \
#--output_dir  /home/Documents/local/Samarth/RNASeq_denovo_sampledata/RSEM_control1

#control2

#/home/Documents/NGS_tools/trinityrnaseq-Trinity-v2.6.6/util/align_and_estimate_abundance.pl --seqType fq \
#--transcripts  /home/Documents/local/Samarth/RNASeq_denovo_sampledata/TRINITY_outdir/Trinity.fasta \
#--left /home/Documents/local/Samarth/RNASeq_denovo_sampledata/trim_c2_left.fastq \
#--right /home/Documents/local/Samarth/RNASeq_denovo_sampledata/trim_c2_right.fastq \
#--SS_lib_type RF \
#--est_method RSEM \
#--aln_method bowtie \
#--trinity_mode \
#--prep_reference \
#--output_dir /home/Documents/local/Samarth/RNASeq_denovo_sampledata/RSEM_control2

#treated1

#/home/Documents/NGS_tools/trinityrnaseq-Trinity-v2.6.6/util/align_and_estimate_abundance.pl --seqType fq \
#--transcripts /home/Documents/local/Samarth/RNASeq_denovo_sampledata/TRINITY_outdir/Trinity.fasta \
#--left /home/Documents/local/Samarth/RNASeq_denovo_sampledata/trim_t1_left.fastq \
#--right /home/Documents/local/Samarth/RNASeq_denovo_sampledata/trim_t1_right.fastq \
#--SS_lib_type RF \
#--est_method RSEM \
#--aln_method bowtie \
#--trinity_mode \
#--prep_reference \
#--output_dir /home/Documents/local/Samarth/RNASeq_denovo_sampledata/RSEM_treated1

#treated2

#/home/Documents/NGS_tools/trinityrnaseq-Trinity-v2.6.6/util/align_and_estimate_abundance.pl --seqType fq \
#--transcripts /home/Documents/local/Samarth/RNASeq_denovo_sampledata/TRINITY_outdir/Trinity.fasta \
#--left /home/Documents/local/Samarth/RNASeq_denovo_sampledata/trim_t2_left.fastq \
#--right /home/Documents/local/Samarth/RNASeq_denovo_sampledata/trim_t2_right.fastq \
#--SS_lib_type RF \
#--est_method RSEM \
#--aln_method bowtie \
#--trinity_mode \
#--prep_reference \
#--output_dir /home/Documents/local/Samarth/RNASeq_denovo_sampledata/RSEM_treated2



### Step8 : Generating the count matrices for DE Analysis

#/home/Documents/NGS_tools/trinityrnaseq-Trinity-v2.6.6/util/abundance_estimates_to_matrix.pl \
#--est_method RSEM \
#--out_prefix abundance_counts \
#--name_sample_by_basedir \
#--gene_trans_map none \
#/home/Documents/local/Samarth/RNASeq_denovo_sampledata/RSEM_control1/RSEM.isoforms.results /home/Documents/local/Samarth/RNASeq_denovo_sampledata/RSEM_control2/RSEM.isoforms.results /home/Documents/local/Samarth/RNASeq_denovo_sampledata/RSEM_treated1/RSEM.isoforms.results /home/Documents/local/Samarth/RNASeq_denovo_sampledata/RSEM_treated2/RSEM.isoforms.results

### step 9: Differential Expression Profiling Using EdgeR

## Biological Analysis Plan 

#c1 c2->C  t1 t2 ->T

#T vs C


#/home/Documents/NGS_tools/trinityrnaseq-Trinity-v2.6.6/Analysis/DifferentialExpression/run_DE_analysis.pl \
#--matrix /home/Documents/local/Samarth/RNASeq_denovo_sampledata/abundance_counts.isoform.counts.matrix \
#--method edgeR \
#--samples_file sample_file_TvsC \
#--dispersion 0.1 \
#--output /home/Documents/local/Samarth/RNASeq_denovo_sampledata/EDGER_Tvs_C



######   For ANNOTATION   ###############



### Step10: to predict all long orfs in trinity.fasta(transcritopme assembly)

#/home/Documents/NGS_tools/TransDecoder-TransDecoder-v5.3.0/TransDecoder.LongOrfs -t /home/Documents/local/Samarth/RNASeq_denovo_sampledata/TRINITY_outdir/Trinity.fasta

#note: we will use longest_orf.pep as query for BLAST
####  to do annotaion we have to done BLAST using longest_orf.pep as Query and have to take https://www.uniprot.org/downloads take reviewed swiss-prot fasta seq as database #######
#a) to create Database

#makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprotdb

#b) to search Query in database
#blastp -query longest_orfs.pep -db uniprotdb -outfmt '6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen' -out blast_result.txt

#%quer coverage=(exact match/Qlen)*100
#exact match=length-mismatch-gaps



