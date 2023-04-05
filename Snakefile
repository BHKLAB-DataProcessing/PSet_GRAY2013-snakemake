from os import path
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"],
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)

prefix = config["prefix"]
filename = config["filename"]
rna_tool = config["rna_tool"]
rna_ref = config["rna_ref"]
is_filtered = config["filtered"]
filtered = 'filtered' if is_filtered is not None and is_filtered == 'filtered' else ''

basePath = "https://orcestradata.blob.core.windows.net/gray"

rna_tool_dir = rna_tool.replace('-', '_')
rnaseq_dir = path.join(prefix, "processed",
                       rna_tool_dir, rna_tool_dir + '_' + rna_ref)
rna_ref_file = rna_ref.replace('_', '.') + '.annotation.RData'

rule get_pset:
    input:
        S3.remote(prefix + "download/drugs_with_ids.csv"),
        S3.remote(prefix + "download/cell_annotation_all.csv"),
        S3.remote(prefix + "download/" + rna_tool_dir + '.tar.gz'),
        S3.remote(prefix + 'download/' + rna_ref_file),
        S3.remote(prefix + 'download/JRGraySRRMapping.csv'),
        S3.remote(prefix + "processed/drug_norm_post.RData"),
        S3.remote(prefix + "processed/profiles.RData"),
        S3.remote(prefix + "data/gb-2013-14-10-r110-s1.xlsx")
    output:
        S3.remote(prefix + filename)
    shell:
        """
        Rscript scripts/GRAY.R {prefix} {filename} {rna_tool} {rna_ref} {filtered}
        """

rule recalculate_and_assemble_slice:
    input:
        S3.remote(prefix + "processed/raw_sense_slices.zip")
    output:
        S3.remote(prefix + "processed/profiles.RData")
    shell:
        """
        Rscript scripts/recalculateAndAssembleSlice.R {prefix}
        """

rule get_sens_data:
    input:
        S3.remote(prefix + "download/drugs_with_ids.csv"),
        S3.remote(prefix + "download/cell_annotation_all.csv"),
        S3.remote(prefix + "data/gb-2013-14-10-r110-s9.txt"),
    output:
        S3.remote(prefix + "processed/drug_norm_post.RData"),
        S3.remote(prefix + "processed/raw_sense_slices.zip")
    shell:
        """
        Rscript scripts/downloadSensData.R {prefix}
        """

rule download_annotation:
    output:
        S3.remote(prefix + "download/drugs_with_ids.csv"),
        S3.remote(prefix + "download/cell_annotation_all.csv"),
        S3.remote(prefix + 'download/' + rna_ref_file),
        S3.remote(prefix + 'download/JRGraySRRMapping.csv')
    shell:
        """
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/drugs_with_ids.csv' \
            -O {prefix}download/drugs_with_ids.csv
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/cell_annotation_all.csv' \
            -O {prefix}download/cell_annotation_all.csv
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/{rna_ref_file}' \
            -O {prefix}download/{rna_ref_file}
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/JRGraySRRMapping.csv' \
            -O {prefix}download/JRGraySRRMapping.csv
        """

rule download_data:
    output:
        S3.remote(prefix + "download/" + rna_tool_dir + '.tar.gz'),
        S3.remote(prefix + "data/gb-2013-14-10-r110-s9.txt"),
        S3.remote(prefix + "data/gb-2013-14-10-r110-s1.xlsx")
    shell:
        """
        wget '{basePath}/RNA-seq/{rna_tool_dir}.tar.gz' -O {prefix}download/{rna_tool_dir}.tar.gz
        wget 'https://github.com/BHKLAB-DataProcessing/getGRAY2013/raw/master/data/gb-2013-14-10-r110-s9.txt' \
            -O {prefix}data/gb-2013-14-10-r110-s9.txt
        wget 'https://github.com/BHKLAB-DataProcessing/getGRAY2013/raw/master/data/gb-2013-14-10-r110-s1.xlsx' \
            -O {prefix}data/gb-2013-14-10-r110-s1.xlsx
        """
