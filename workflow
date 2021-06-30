## Reads pre-processing and alignment ##
The strand orientation of the original RNA was preserved during the process of library construction and reads 2 yields sequences sense to the original RNA. Thus, only reads 2 was used for m 6 A/m 6 Am signal identiﬁcation in our study. Raw sequencing reads were ﬁrst subjected to Trim_galore.

trim_galore -q 20 --phred33 –fastqc --stringency 3 -o ./reads2_trimed  ${file}

All reads that mapped to human rRNA by bowtie2 were removed. Processed reads were mapped to genome.

bowtie2 --end-to-end -p 5 -x rRNA_index -U ${file}_trimed -S ${file}.bwt2-rRNA.sam --un-gz ${file} _de-rRNA.fq.gz

echo “hisat2 -p 5 -k 1 -x hg19_index -U ${file}_de-rRNA.fq.gz | samtools view -@ 10 -h -q 1 -Sbo ${file}.bam -" > ${file}.map_genome.sh
sh ${file}.map_genome.sh


## Identiﬁcation of 5’-UTR peaks in m6A-IP sample ##
For genome-base peak caller MACS2 , the effective genome size was set to 2.7 * 109 for human, under the option of -nomodel and P-value cutoff 0.01.

macs2 callpeak -t ${m6A_IP.bam} -c ${m6A_input.bam} -f BAM -g hs --nomodel -q 0.01 -n  ${macs2_m6A_peak}

macs2 callpeak -t ${m7G_IP.bam} -c ${m7G_input.bam} -f BAM -g hs --nomodel -q 0.01 -n  ${macs2_m7G_peak}

bedtools intersect -wa -a ${macs2_m6A_peak} -b ${macs2_m7G_peak} > ${macs2_m6A_peak_UTR5}

Compare the peak intensity of ${macs2_m6A_peak_UTR5} and ${macs2_m7G_peak}, m6A peak intensity higher than m7G peak was used for m6Am signal identification.
Calculation of the demethylase-sensitivity peaks.
${macs2_m6A_peak_UTR5} were used for analysis. Secondly, differential methylation peak between FTO (+) and FTO (-) sample were identified using exomePeak (P<0.01).
gtf=c("NCBI_RefSeq_hg19.GTF")
demethylase_sensitivity_peaks <- exomepeak(GENE_ANNO_GTF=gtf, IP_BAM=c(${nonFTO_m6AIP_bam}),INPUT_BAM=c(${nonFTO_Input_bam}),TREATED_IP_BAM=c(${FTO_m6AIP_bam}),TREATED_INPUT_BAM=c(${FTO_input_bam}), EXPERIMENT_NAME = c("${emethylase_sensitivity_peaks}") )

An m6Am peak was identified when: (1) peak was high-confidence in m6A-IP sample; (2) peak was demethylase-sensitivity between the FTO (+) and FTO (-) sample.


## Identification of m6Am site ##
To identify the potential m6Am site, we first need to calculate the start reads and start rate for each locus on the genome. 
samtools sort -o ${FTO_bam_sort} ${FTO_bam}
samtools sort -o ${nonFTO_bam_sort} ${nonFTO_bam}
samtools mpileup -BQ0 -d 10000000 -o ${FTO_bam_sort.pileup} -f hg19_genome.fa ${FTO_bam_sort}
samtools mpileup -BQ0 -d 10000000 -o ${nonFTO_bam_sort.pileup} -f hg19_genome.fa ${nonFTO_bam_sort}
perl get_started_reads.pl  ${FTO_bam_sort.pileup}  ${FTO_bam_sort.pileup_ started}
perl get_started_reads.pl ${nonFTO_bam_sort.pileup}  ${nonFTO_bam_sort.pileup_ started}

We defined a “start rate difference” score (SRD score) for each nucleotide within an m6Am peak.
(1)	m1=[FTO (-) start reads] / [FTO (-) depth]
(2)	m2=[FTO (-) depth - FTO (+) depth] / [FTO (-) depth]
(3)	SRD score = m1 + m2

The SRD score were used to identify the m6Am site. The defined a position i was defined to be m6Am when the following criteria were met: (1) position i must be located within m6Am peak and carry adenosine residue; (2) the start reads of position i are not less than 20 in FTO (-) sample; (3) the start rate of position i in FTO (-) sample is greater than that in m7G-IP sample; (4) the start reads coverage in m7G-IP sample is greater than that in input; (5) the position i possesses the top 3 of SRD scores (SRD score>1).

