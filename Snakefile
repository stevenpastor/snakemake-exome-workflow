configfile: "config.yaml"

rule all:
    input:
        "report.html"


rule bwa_map:
    input:
        fa = expand("{refdir}/{ref}.fa.gz", refdir=config["reference_dir"], 
            ref=config["references"]),
        r1 = expand("{sampdir}/{sample}_R1.fq.gz", sampdir=config["samples_dir"], 
            sample=config["samples"]),
        r2 = expand("{sampdir}/{sample}_R2.fq.gz", sampdir=config["samples_dir"], 
            sample=config["samples"])
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: config["threads"]
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input.fa} {input.r1} {input.r2} | "
        "samtools view -@ {threads} -Sb - > {output}) 2> {log}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.bam")
    threads: config["threads"]
    shell:
        "samtools sort -@ {threads} -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"


rule bcftools_call:
    input:
        fa = expand("{refdir}/{ref}.fa.gz", refdir=config["reference_dir"], 
            ref=config["references"]),
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "variant_calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"


rule report:
    input:
        "variant_calls/all.vcf"
    output:
        "report.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as vcf:
            n_calls = sum(1 for l in vcf if not l.startswith("#"))

        report("""
        An example variant calling workflow
        ===================================

        Reads were mapped to the Human (hg38)
        reference genome and variants were called jointly with
        SAMtools/BCFtools.

        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], T1=input[0])
