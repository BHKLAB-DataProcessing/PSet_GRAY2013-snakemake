from os import path
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
# S3 = S3RemoteProvider(
#     access_key_id=config["key"],
#     secret_access_key=config["secret"],
#     host=config["host"],
#     stay_on_remote=False
# )

prefix = config["prefix"]
rna_tool = 'Kallisto-0.46.1'
rna_ref = 'Gencode_v33'
basePath = "https://orcestradata.blob.core.windows.net/gray"

rna_tool_dir = rna_tool.replace('-', '_')
rnaseq_dir = path.join(prefix, "processed",
                       rna_tool_dir, rna_tool_dir + '_' + rna_ref)
rna_ref_file = rna_ref.replace('_', '.') + '.annotation.RData'

rule get_pset:
    input:
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv",
        prefix + "download/" + rna_tool_dir + '.tar.gz',
        prefix + 'download/' + rna_ref_file,
        prefix + 'download/JRGraySRRMapping.csv',
        prefix + "processed/drug_norm_post.RData",
        prefix + "processed/profiles.RData"
    output:
        prefix + "GRAY2013.rds"
    shell:
        """
        Rscript {prefix}scripts/GRAY.R {prefix} {rna_tool} {rna_ref}
        """

rule recalculate_and_assemble_slice:
    input:
        prefix + "processed/raw_sense_slices.zip"
    output:
        prefix + "processed/profiles.RData"
    shell:
        """
        Rscript {prefix}scripts/recalculateAndAssembleSlice.R {prefix}
        """

rule get_sens_data:
    input:
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv",
        prefix + "data/gb-2013-14-10-r110-s9.txt",
    output:
        prefix + "processed/drug_norm_post.RData",
        prefix + "processed/raw_sense_slices.zip"
    shell:
        """
        Rscript {prefix}scripts/downloadSensData.R {prefix}
        """

rule download_annotation:
    output:
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv",
        prefix + 'download/' + rna_ref_file,
        prefix + 'download/JRGraySRRMapping.csv'
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
        prefix + "download/" + rna_tool_dir + '.tar.gz'
    shell:
        """
        wget '{basePath}/RNA-seq/{rna_tool_dir}.tar.gz' -O {prefix}download/{rna_tool_dir}.tar.gz
        """
