import yaml

def load_configfile(configfile, verbose, info='Config'):
   with open(configfile, "r") as f:
       config = yaml.load(f)

# extract the base and group names from the sampleinfo
def get_sampleinfo_info(sampleinfo, column):
    with open(sampleinfo) as f:
        lines = f.read().split('\n')[:-1]
        data = ""
        for i, line in enumerate(lines):
            if i > 0: # header
                data += line.split()[column] + " "
    return data

def multiqc_input_check(return_value):
    infiles = []
    indir = ""
    infiles.append( expand("01_fastq/adaptor_trimmed/fastqc/{reads}_fastqc.html", reads=reads) )
    indir += " 01_fastq/adaptor_trimmed/fastqc "
    infiles.append( expand("01_fastq/adaptor_trimmed/{reads}.fastq.gz", reads=reads) )
    indir += "01_fastq/adaptor_trimmed "
    infiles.append( expand("03_mapping"+spike_prefix+"/{base}.bam", base=base) )
    indir += "03_mapping"+spike_prefix
    if return_value == "infiles":
        return(infiles)
    else:
        return(indir)

groups = get_sampleinfo_info(config["sampleinfo"], 2).split()
base = get_sampleinfo_info(config["sampleinfo"], 1).split()
reads = config["reads"]
spike_prefix = config["spikeIn_prefix"] if config["map_spike"] == True else ""
fold_change_prefix = "tss_foldch_"+config["fold_change"]+"_" #changed after running

def spike_mapping():
    output = []
    if config["map_spike"]:
        output.append([expand("04_removedup{spike}/CSobject.Rdata", spike = spike_prefix),
                      expand("04_removedup{spike}/{base}.filtered.bam", base = base, spike = spike_prefix),
                      expand("multiQC{spike}/multiqc_report.html", spike = spike_prefix)
                      ])
    return(output)


rule all:
    input:
        expand("04_removedup/CSobject.Rdata"),
        expand("04_removedup/{base}.filtered.bam", base = base),
        expand("05_bigwigs/with_duplicates/{base}.bw", base = base),
        expand("05_bigwigs/without_duplicates/{base}.bw", base = base),
        "06_tss_calling/CSobject.Rdata",
        expand("07_tss_annotation/{fold_change}{group}.annotated.tsv", fold_change=fold_change_prefix, group=groups),
        expand("08_plots/"),
        expand("multiQC/multiqc_report.html"),
        spike_mapping()

## FASTQ linking
rule FASTQ:
    input:
        config["indir"]+"/{read}.fastq.gz"
    output:
        "01_fastq/raw/{read}.fastq.gz"
    shell:
        "( [ -f {output} ] || ln -s -r {input} {output} ) && touch -h {output}"

## TRIMMING
rule TrimGalore:
    input:
        r1 = "01_fastq/raw/"+reads[0]+".fastq.gz",
        r2 = "01_fastq/raw/"+reads[1]+".fastq.gz"
    output:
        r1 = "01_fastq/adaptor_trimmed/"+reads[0]+".fastq.gz",
        r2 = "01_fastq/adaptor_trimmed/"+reads[1]+".fastq.gz"
    params:
        tmp1 = "01_fastq/adaptor_trimmed/"+reads[0]+"_val_1.fq.gz",
        tmp2 = "01_fastq/adaptor_trimmed/"+reads[1]+"_val_2.fq.gz",
        opts = str(config["trim_options"] or ''),
        output_dir = "01_fastq/adaptor_trimmed/"
    log:
        "01_fastq/adaptor_trimmed/Logs/TrimGalore.log"
    shell:
        config["trim_galore_path"]+"trim_galore "
        "--path_to_cutadapt "+config["cutadapt_path"]+"cutadapt "
        "--output_dir {params.output_dir} "
        "--paired "
        "--stringency 3 "
        "{params.opts} "
        "{input.r1} {input.r2} "
        "&> {log} "
        "&& (mv {params.tmp1} {output.r1}; mv {params.tmp2} {output.r2})"

## FASTQC
rule FastQC_on_trimmed:
    input:
        "01_fastq/adaptor_trimmed/{read}.fastq.gz"
    output:
        "01_fastq/adaptor_trimmed/fastqc/{read}_fastqc.html"
    params:
        output_dir = "01_fastq/adaptor_trimmed/fastqc"
    log:
        "01_fastq/adaptor_trimmed/fastqc/Logs/FastQC_trimmed.{read}.log"
    threads: 2
    shell:
        "module load FastQC && fastqc -o {params.output_dir} {input} &> {log}"

## DEMULTIPLEXING
rule demultiplex_fastq:
    input:
        r1 = rules.TrimGalore.output.r1,
        r2 = rules.TrimGalore.output.r2,
        sampleinfo = config["sampleinfo"]
    output:
        "02_splitting/CSobject.Rdata",
        expand("02_splitting/{base}_R1.fastq.gz", base = base),
        expand("02_splitting/{base}_R2.fastq.gz", base = base)
    params:
        expname = config["expname"],
        outdir = "02_splitting",
        rscript = config["demult"]
    log:
        "02_splitting/Logs/demultiplex.log"
    threads: 20
    shell:
        "export R_LIBS_USER="+config["R_libs_path"] +
        " && "+os.path.join(config["R_path"],"Rscript") +
        " {params.rscript} {params.expname}"
        " {input.r1} {input.r2} {input.sampleinfo} {params.outdir} {threads} &> {log}"

## MAPPING
if config["mapping_prg"] == "subread":
    rule subread:
        input:
            r1 = "02_splitting/{base}_R1.fastq.gz",
            r2 = "02_splitting/{base}_R2.fastq.gz"
        output:
            "03_mapping"+spike_prefix+"/{base}.bam"
        params:
            gtf = lambda wildcards: "" if config["map_spike"] else "-a "+ config["genes_gtf"],
            index = lambda wildcards: config["spike_index_subread"] if config["map_spike"] else config["subread_index"]
        log:
            "03_mapping{}/Logs/subread.log".format(spike_prefix)
        threads: 15
        shell:
            "module load samtools && module load subread &&"
            " subjunc -D 500 -d 20 -T {threads}"
            " {params.gtf}"
            " -i {params.index}"
            " -r {input.r1}"
            " -R {input.r2} |"
            " samtools sort"
            " -O BAM -@ {threads}"
            " -o {output}"
            " &> {log}"

elif config["mapping_prg"] == "STAR":
    if config["paired"]:
        rule STAR:
            input:
                r1 = "02_splitting/{base}_R1.fastq.gz",
                r2 = "02_splitting/{base}_R2.fastq.gz"
            output:
                 "03_mapping"+spike_prefix+"/{base}.bam"
            params:
                star_options = str(config["star_options"] or ''),
                gtf = lambda wildcards: "" if config["map_spike"] else "--sjdbGTFfile "+ config["genes_gtf"],
                index = lambda wildcards: config["spike_index_star"] if config["map_spike"] else config["star_index"],
                sample_dir = "03_mapping"+spike_prefix+"/{base}",
                prefix = "03_mapping"+spike_prefix+"/{base}/{base}."
            benchmark:
                "03_mapping"+spike_prefix+"/.benchmark/STAR.{base}.benchmark"
            threads: 12
            shell:
                "( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} ) &&"
                " module load STAR &&"
                " STAR"
                " --runThreadN {threads}"
                " {params.star_options}"
                " --limitBAMsortRAM 14806948926"
                " --sjdbOverhang 100"
                " --readFilesCommand zcat"
                " --outSAMunmapped Within"
                " --outSAMtype BAM SortedByCoordinate"
                " {params.gtf}"
                " --genomeDir {params.index}"
                " --readFilesIn {input.r1} {input.r2}"
                " --outFileNamePrefix {params.prefix}"
                " && mv {params.prefix}Aligned.sortedByCoord.out.bam {output}"

## INDEXING
rule bam_index:
    input:
        "03_mapping"+spike_prefix+"/{base}.bam"
    output:
        "03_mapping"+spike_prefix+"/{base}.bam.bai"
    log:
        "03_mapping{}/Logs/bam.index.log".format(spike_prefix)
    shell:
        "module load samtools && "
        "samtools index {input} &> {log}"

## DUPLICATE REMOVAL
rule remove_dup:
    input:
        expand("03_mapping"+spike_prefix+"/{base}.bam.bai", base = base)
    output:
        expand("04_removedup"+spike_prefix+"/CSobject.Rdata"),
        expand("04_removedup"+spike_prefix+"/{base}.filtered.bam", base = base)
    params:
        output_dir = "04_removedup"+spike_prefix,
        input_dir = "03_mapping"+spike_prefix,
        last_cs = "02_splitting/CSobject.Rdata",
        rscript = config["removedup"]
    log:
        "04_removedup"+spike_prefix+"/Log/removedup.log"
    threads: 10
    shell:
        "export R_LIBS_USER="+config["R_libs_path"]+" && "+os.path.join(config["R_path"],"Rscript")+ " {params.rscript} "
        "{params.last_cs} {params.input_dir} "
        "{params.output_dir} {threads} &> {log}"

## INDEXING DUPremoved BAMs
rule bam_index_rmdup:
    input:
        "04_removedup"+spike_prefix+"/{base}.filtered.bam"
    output:
        "04_removedup"+spike_prefix+"/{base}.filtered.bam.bai"
    log:
        "04_removedup{}/Logs/bam.index.log".format(spike_prefix)
    shell:
        "module load samtools && "
        "samtools index {input} &> {log}"

## BAMCOVERAGE
# with dup
rule bam_coverage:
    input:
        bam = "03_mapping/{base}.bam",
        bai = "03_mapping/{base}.bam.bai"
    output:
        "05_bigwigs/with_duplicates/{base}.bw"
    log:
        "05_bigwigs/with_duplicates/Logs/bam.coverage.log"
    benchmark:
        "05_bigwigs/with_duplicates/Logs/bam.coverage.benchmark"
    threads: 20
    shell:
        "module load deeptools/3.0.2 && "
        "bamCoverage "
        "--Offset 1 5 "
        "--normalizeUsing CPM -p 20 -bs 1 "
        "-b {input.bam} "
        "-o {output}  &> {log}"

# without dup
rule bam_coverage_withoutdup:
    input:
        bam = "04_removedup/{base}.filtered.bam",
        bai = "04_removedup/{base}.filtered.bam.bai"
    output:
        "05_bigwigs/without_duplicates/{base}.bw"
    log:
        "05_bigwigs/without_duplicates/Logs/bam.coverage.log"
    benchmark:
        "05_bigwigs/without_duplicates/Logs/bam.coverage.benchmark"
    threads: 20
    shell:
        "module load deeptools/3.0.2 && "
        "bamCoverage "
        "--Offset 1 5 "
        "--normalizeUsing CPM -p 20 -bs 1 "
        "-b {input.bam} "
        "-o {output}  &> {log}"

rule tss_calling:
    input:
        "04_removedup/CSobject.Rdata"
    output:
        outfile = "06_tss_calling/CSobject.Rdata",
        bedfiles = expand("06_tss_calling/{foldchange}{group}.bed", foldchange=fold_change_prefix, group=groups)
    params:
        last_cs = "04_removedup/CSobject.Rdata",
        sampleinfo = config["sampleinfo"],
        output_dir = "06_tss_calling/",
        fold_change = config["fold_change"],
        prefix = "tss_foldch_{}".format(config["fold_change"]),
        rscript = config["tssdetection"]
    log:
        "06_tss_calling/tss.log"
    threads: 10
    shell:
        "export R_LIBS_USER="+config["R_libs_path"]+" && "+os.path.join(config["R_path"],"Rscript")+ " {params.rscript} "
        "{params.last_cs} {params.sampleinfo} "
        "{params.output_dir} {threads} {params.fold_change} {params.prefix} &> {log}"

rule tss_annotation:
    input:
        "06_tss_calling/CSobject.Rdata",
        expand("06_tss_calling/{foldchange}{group}.bed", foldchange=fold_change_prefix, group=groups)
    output:
        "07_tss_annotation/{tss}.annotated.tsv"
    params:
        last_cs = "06_tss_calling/CSobject.Rdata",
        input_bed = "06_tss_calling/{tss}.bed",
        plotFile = "07_tss_annotation/{tss}.pdf",
        GTFfile = config["GTFfile"],
        enhancer_file = config["enhancer_file"],
        repeat_file = config["repeat_file"],
        dhs_file = "NA",
        rscript = config["tss_annotation"]
    log:
        "07_tss_annotation/tss.log"
    threads:
        10
    shell:
        "export R_LIBS_USER="+config["R_libs_path"]+" && "+os.path.join(config["R_path"],"Rscript")+ " {params.rscript} "
        "{params.input_bed} {output} {params.plotFile} "
        "{params.GTFfile} {params.enhancer_file} {params.repeat_file} {params.dhs_file} &> {log}"

rule plot_stats:
    input:
        "06_tss_calling/CSobject.Rdata"
    output:
        "08_plots/"
    params:
        prefix = "08_plots/tss_foldch_{}".format(config["fold_change"]),
        rscript = config["plotscripts"],
        last_cs = "06_tss_calling/CSobject.Rdata"
    log:
        "08_plots/plots.log"
    shell:
        "export R_LIBS_USER="+config["R_libs_path"]+" && "+os.path.join(config["R_path"],"Rscript")+ " {params.rscript} "
        "{params.last_cs} {params.prefix} &> {log}"

rule multiQC:
    input:
        multiqc_input_check(return_value = "infiles")
    output:
        "multiQC"+spike_prefix+"/multiqc_report.html"
    params:
        indirs = multiqc_input_check(return_value = "indir")
    log:
        "multiQC"+spike_prefix+"/multiQC.log"
    shell:
        "module load MultiQC && "
        "multiqc -o {output} -f {params.indirs} &> {log}"
