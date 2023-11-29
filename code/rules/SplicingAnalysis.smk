rule ExtractJuncs:
    input:
        bam = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai",
    output:
        "SplicingAnalysis/juncfiles/{sample}.junc",
    params:
        strand = "0"
    conda:
        "../envs/regtools.yml"
    log:
        "logs/ExtractJuncs/{sample}.log"
    shell:
        """
        (regtools junctions extract -m 20 -s {params.strand} {input.bam} > {output}) &> {log}
        """

rule Download_hg38_basic_gtf:
    output:
        "Gencode/hg38_Primary_basic_v44.gtf"
    shell:
        """
        wget -O- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.basic.annotation.gtf.gz | zcat > {output.chr_gtf}
        """


rule annotate_juncfiles:
    input:
        fa = FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.fa"),
        fai = FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai"),
        gtf = FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"),
        juncs = "SplicingAnalysis/juncfiles/{sample}.junc",
    output:
        counts = "SplicingAnalysis/juncfiles/{sample}.junccounts.tsv.gz"
    log:
        "logs/annotate_juncfiles/{sample}.log"
    conda:
        "../envs/regtools.yml"
    shell:
        """
        (regtools junctions annotate {input.juncs} {input.fa} {input.gtf} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{$4=$1"_"$2"_"$3"_"$6; print $4, $5}}' | gzip - > {output.counts} ) &> {log}
        """


rule ConcatJuncFilesAndKeepUniq:
    input:
        ExpandAllSamplesInFormatStringFromGenomeNameWildcard("SplicingAnalysis/juncfiles/{sample}.junc")
    output:
        "SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.junc"
    log:
        "logs/ConcatJuncFilesAndKeepUniq/{GenomeName}.log"
    conda:
        "../envs/r_2.yml"
    shell:
        """
        Rscript scripts/Collapse_Juncsfiles.R {output} {input} &> {log}
        """

