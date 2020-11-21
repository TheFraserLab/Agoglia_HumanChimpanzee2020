# the module command actually needs to edit the shell environment, and is not a real executable
module = """module () {
        eval $($LMOD_CMD bash "$@")
        [ $? = 0 ] && eval $(${LMOD_SETTARG_CMD:-:} -s sh)
        }"""

# Separate data from code!
configfile: 'config.json'

species = [ "chimp"]
samples = config['samples']
read_length = 150

rule all:
    input: 
#expand('analysis/{species}/{sample}/{sample}_ase_by_reads_merged.txt', species=species, sample=samples)
        expand('analysis/{species}/{sample}/{sample}_ase_by_reads_merged.chr.txt', species=species, sample=samples)

def seqprep_inputs(wildcards):
    return config["samples"][wildcards.sample]


rule seqprep:
    input: 
        seqprep_inputs
    output: 
        R1="analysis/{sample}/r1_trimmed.fq.gz", 
        R2="analysis/{sample}/r2_trimmed.fq.gz"
    shell: """
    /home/users/ragoglia/bin/SeqPrep/SeqPrep \
            -f {input[0]} -r {input[1]} \
            -1 {output.R1} -2 {output.R2} \
			-A GATCGGAAGAGCACACGTCT -B GATCGGAAGAGCGTCGTGTA
    """
 
overhang = read_length - 1

rule star_genome:
    input: 
        fasta="/scratch/users/ragoglia/reference/{species}/{species}.fasta",
        gff="/scratch/users/ragoglia/reference/{species}/{species}.gff"
    threads: 12
    output:
        "/scratch/users/ragoglia/reference/{species}/STAR/Genome"
    shell: """ {module}; module load STAR/2.5.1b;
    rm -rf /scratch/users/ragoglia/reference/{wildcards.species}/STAR/tmp
    STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir /scratch/users/ragoglia/reference/{wildcards.species}/STAR \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gff} \
            --sjdbOverhang {overhang} \
            --outTmpDir /scratch/users/ragoglia/reference/{wildcards.species}/STAR/tmp
    """

rule star_align:
    input: 
        R1="analysis/{sample}/r1_trimmed.fq.gz", R2="analysis/{sample}/r2_trimmed.fq.gz",
        genome="/scratch/users/ragoglia/reference/{species}/STAR/Genome",
    output: 
        bam="analysis/{species}/{sample}/STAR1/Aligned.out.bam",
        sj="analysis/{species}/{sample}/STAR1/SJ.out.tab"
    threads: 10
    params: mem="64G"
    shell: """
    {module}
    module load STAR/2.5.1b
    rm -rf analysis/{wildcards.species}/{wildcards.sample}/STAR1/STARtmp
    STAR \
            --genomeDir /scratch/users/ragoglia/reference/{wildcards.species}/STAR \
            --outFileNamePrefix analysis/{wildcards.species}/{wildcards.sample}/STAR1/ \
            --outSAMattributes MD NH \
            --outSAMtype BAM Unsorted \
            --outTmpDir analysis/{wildcards.species}/{wildcards.sample}/STAR1/STARtmp \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --readFilesIn {input.R1} {input.R2} \
			--outFilterMultimapNmax 1
    """
 
rule star_align2:
    input: 
        R1="analysis/{sample}/r1_trimmed.fq.gz", R2="analysis/{sample}/r2_trimmed.fq.gz",
        genome="/scratch/users/ragoglia/reference/{species}/STAR/Genome",
        sj="analysis/{species}/{sample}/STAR1/SJ.out.tab"
    output: 
        "analysis/{species}/{sample}/STAR2/Aligned.out.bam"
    threads: 10
    params: mem="64G"
    shell: """ {module}; module load STAR/2.5.1b;
    rm -rf analysis/{wildcards.species}/{wildcards.sample}/STAR2/STARtmp
    STAR \
            --genomeDir /scratch/users/ragoglia/reference/{wildcards.species}/STAR \
            --outFileNamePrefix analysis/{wildcards.species}/{wildcards.sample}/STAR2/ \
            --outSAMattributes MD NH \
            --outSAMtype BAM Unsorted \
            --outTmpDir analysis/{wildcards.species}/{wildcards.sample}/STAR2/STARtmp \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --readFilesIn {input.R1} {input.R2} \
			--outFilterMultimapNmax 1 \
			--sjdbFileChrStartEnd {input.sj}
    """
#
#rule rmdup:
#    input: "analysis/{species}/{sample}/STAR2/Aligned.sortedByCoord.out.bam"
#    output: "analysis/{species}/{sample}/rmdup.bam"
#    log: "logs/{sample}_rmdups.log"
#    shell: """
#    java -jvm-args -jar /scratch/users/ragoglia/bin/picard.jar MarkDuplicates I={input} O={output} M={log} \
#		REMOVE_DUPLICATES=true \
#		DUPLICATE_SCORING_STRATEGY=RANDOM
#        """
#
#FOR TESTING ONLY
rule sort:
    input: 
        "analysis/{genome}/{sample}/STAR2/Aligned.out.bam"
    output:
        "analysis/{genome}/{sample}/STAR2/Aligned.out.sort.bam"
    shell: """ {module}; module load samtools/1.3;
    samtools sort \
        {input} -o {output}
    """


#FOR TESTING ONLY 
rule rmdup:
     input: "analysis/{species}/{sample}/STAR2/Aligned.out.sort.bam"
     output: "analysis/{species}/{sample}/rmdup.bam"
     shell: """
     python ~/bin/Genomics/rmdup_for_ase.py \
         -i {input} -o {output} -p -s /home/users/ragoglia/bin/samtools-0.1.19/samtools
     """
rule index:
    input:
        "analysis/{genome}/{sample}/rmdup.bam"
    output:
        "analysis/{genome}/{sample}/rmdup.bam.bai"
        
    shell: " {module}; module load samtools/1.3; samtools index {input}"

rule wasp_find_snps:
    input: 
        bam="analysis/{genome}/{sample}/rmdup.bam",
        snpdir="/home/users/ragoglia/SNPs/ase/{genome}/wasp_split_species"
    output: 
        "analysis/{genome}/{sample}/rmdup.remap.fq1.gz",
        "analysis/{genome}/{sample}/rmdup.remap.fq2.gz",
        "analysis/{genome}/{sample}/rmdup.keep.bam",
        "analysis/{genome}/{sample}/rmdup.to.remap.bam",

    shell: """{module}; module load fraserconda/5.0;
    python find_intersecting_snps.py \
        -p -P -s\
        {input.bam} {input.snpdir}
    """

rule wasp_remap:
    input:
        R1="analysis/{genome}/{sample}/rmdup.remap.fq1.gz",
        R2="analysis/{genome}/{sample}/rmdup.remap.fq2.gz",
        genome="/scratch/users/ragoglia/reference/{genome}/STAR/Genome",
        sj="analysis/{genome}/{sample}/STAR1/SJ.out.tab"
    output:
        "analysis/{genome}/{sample}/rmdup.remap.bam",
    threads: 10
    shell: """{module}; module load STAR/2.5.1b;
    rm -rf analysis/{wildcards.genome}/{wildcards.sample}/STARtmp
    STAR \
            --genomeDir /scratch/users/ragoglia/reference/{wildcards.genome}/STAR \
            --outFileNamePrefix analysis/{wildcards.genome}/{wildcards.sample}/remap \
            --outSAMattributes MD NH \
            --outSAMtype BAM Unsorted \
            --outTmpDir analysis/{wildcards.genome}/{wildcards.sample}/STARtmp \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --outFilterMultimapNmax 1 \
            --sjdbFileChrStartEnd {input.sj} \
            --readFilesIn {input.R1} {input.R2}
    mv analysis/{wildcards.genome}/{wildcards.sample}/remapAligned.out.bam {output}
            """
rule sort_wasp:
    input: 
        "analysis/{genome}/{sample}/rmdup.remap.bam"
    output:
        "analysis/{genome}/{sample}/rmdup.remap.sort.bam"
    shell: """ {module}; module load samtools/1.3; 
    samtools sort \
        -n \
        {input} -o {output}
    """


rule wasp_filter:
    input: 
        toremap="analysis/{genome}/{sample}/rmdup.to.remap.bam",
        remapped="analysis/{genome}/{sample}/rmdup.remap.sort.bam",
    output:
        "analysis/{genome}/{sample}/rmdup.remap.kept.bam",

    shell:""" {module}; module load fraserconda/5.0;
    python /home/users/ragoglia/bin/Hornet/mapping/filter_remapped_reads.py \
        -p \
        {input.toremap} {input.remapped} \
        {output}
    """

            
rule sort_ase:
    input: 
        "analysis/{genome}/{sample}/rmdup.remap.kept.bam"
    output:
        "analysis/{genome}/{sample}/rmdup.remap.kept.sort.bam"
    shell: """ {module}; module load samtools/1.3; 
    samtools sort \
        {input} -o {output}
    """

rule sort_keep:
    input: 
        "analysis/{genome}/{sample}/rmdup.keep.bam"
    output:
        "analysis/{genome}/{sample}/rmdup.keep.sort.bam"
    shell: """ {module}; module load samtools/1.3; 
    samtools sort \
        {input} -o {output}
    """
 


rule wasp_merge:
    input:
        "analysis/{genome}/{sample}/rmdup.remap.kept.sort.bam",
        "analysis/{genome}/{sample}/rmdup.keep.sort.bam",
    output:
        "analysis/{genome}/{sample}/rmdup.remap.kept.merged.bam"
    shell:
        "{module}; module load samtools/1.3; samtools merge {output} {input}"


rule index_ase:
    input:
        "analysis/{genome}/{sample}/rmdup.remap.kept.merged.bam"
    output:
        "analysis/{genome}/{sample}/rmdup.remap.kept.merged.bam.bai"
        
    shell: " {module}; module load samtools/1.3; samtools index {input}"

rule count_ase:
    input:
        snps="/home/users/ragoglia/SNPs/ase/{genome}/ASE_SNPs.FILTER.SPLIT_SPECIES.bed",      
        gff="/scratch/users/ragoglia/reference/{genome}/{genome}.chr.gff",
        reads="analysis/{genome}/{sample}/rmdup.remap.kept.merged.bam",
        readsidx="analysis/{genome}/{sample}/rmdup.remap.kept.merged.bam.bai",
    conda:"/home/users/ragoglia/envs/test.yml" 
    output:
        "analysis/{genome}/{sample}/{sample}_ase_by_reads_merged.chr.txt"
#shell: """{module}; module load conda;
    shell: """python ~/bin/ASEr/bin/GetGeneASEbyReads.py \
        {input.snps} {input.gff} {input.reads}\
        -o {output}
    """

#ENSEMBL = "89"
#GRCh = "GRCh38"
#ChimpV = "2.1.4"
#
#rule get_human_fasta:
#    input: 
#    output: "reference/human.fasta"
#    version: ENSEMBL
#    shell: """{module}; module load wget;
#    rm -f {output}.gz
#    wget -O {output}.gz ftp://ftp.ensembl.org/pub/release-{ENSEMBL}/fasta/homo_sapiens/dna/Homo_sapiens.{GRCh}.dna.toplevel.fa.gz
#    gunzip {output}.gz
#    """
#
#rule get_chimp_fasta:
#    input: 
#    output: "reference/chimp.fasta"
#    version: ENSEMBL
#    shell: """{module}; module load wget;
#    rm -f {output}.gz
#    wget -O {output}.gz ftp://ftp.ensembl.org/pub/release-{ENSEMBL}/fasta/pan_troglodytes/dna/Pan_troglodytes.CHIMP{ChimpV}.dna_rm.toplevel.fa.gz
#    gunzip {output}.gz
#    """
#
#rule get_human_gff:
#    input: 
#    output: "reference/human.gff"
#    version: ENSEMBL
#    shell: """{module}; module load wget;
#    rm -f {output}.gz
#    wget -O {output}.gz ftp://ftp.ensembl.org/pub/release-{ENSEMBL}/gff3/homo_sapiens/Homo_sapiens.{GRCh}.{ENSEMBL}.chr.gff3.gz
#    gunzip {output}.gz
#    """
#
#rule get_chimp_gff:
#    input: 
#    output: "reference/chimp.gff"
#    version: ENSEMBL
#    shell: """{module}; module load wget;
#    rm -f {output}.gz
#    wget -O {output}.gz ftp://ftp.ensembl.org/pub/release-{ENSEMBL}/gff3/pan_troglodytes/Pan_troglodytes.CHIMP{ChimpV}.{ENSEMBL}.chr.gff3.gz
#    gunzip {output}.gz
#    """
#
#rule get_snps:
#    input:
#        "reference/{reference}.fasta",
#    params:
#        alternate = lambda wildcards: ["chimp", "human"][wildcards.reference=="chimp"]
#    output:
#        "analysis/{reference}/snps.bed"
#    shell: "echo {wildcards.reference} {params.alternate}> {output}"
