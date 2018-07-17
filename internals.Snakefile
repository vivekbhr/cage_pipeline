import yaml


# ------- define functions --------------

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


#------------ define some globals -------------------

# check whether annotation exists
if config["knownAnnotations"]:
    annotation_folder = config["knownAnnotations"]
else:
    annotation_folder = "gtf_annotations"

# define groups
groups = get_sampleinfo_info(config["sampleinfo"], 2).split()
groups.append('merged')
base = get_sampleinfo_info(config["sampleinfo"], 1).split()
reads = config["reads"]

spike_prefix = config["spikeIn_prefix"] if config["map_spike"] is True else ""
fold_change_prefix = "tss_foldch_"+config["fold_change"] #changed after running
