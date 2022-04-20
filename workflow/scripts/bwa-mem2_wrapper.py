# Any Python script in the scripts folder will be able to import from this module.
__author__ = "Sebastian Röner"
__copyright__ = (
    "Copyright 2022, Sebastian Röner"
)
__email__ = "sebastian.roener@bih-charite.de"
__license__ = "MIT"


from os import path
import sys
from wsgiref.handlers import BaseCGIHandler

from snakemake.shell import shell

inputs = snakemake.input
params = snakemake.params

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Check inputs/arguments.
if not isinstance(snakemake.input.reads, str) and len(snakemake.input.reads) not in {
    1,
    2,
}:
    raise ValueError("input must have 1 (single-end) or 2 (paired-end) elements")


non_merged_cmd = ""
singleton_cmd = ""

if snakemake.input.get("noadapter_R1") and snakemake.input.get("noadapter_R2"):
    non_merged_cmd = "bwa-mem2 mem -t {snakemake.threads} -R \"{snakemake.params.RG}\" {snakemake.input.ref} {snakemake.input.noadapter_R1} {snakemake.input.noadapter_R2}| grep -v \"^@\" || true ; "

if snakemake.input.get("single_reads"):
    singleton_cmd = "bwa-mem2 mem -t {snakemake.threads} -R \"{snakemake.params.RG}\" {snakemake.input.ref} {snakemake.input.single_reads} | grep -v \"^@\" || true ; "

base_cmd = f"((bwa-mem2 mem -t {{snakemake.threads}} -R \"{{snakemake.params.RG}}\" {{snakemake.input.ref}} {{snakemake.input.reads}}; \
{non_merged_cmd}\
{singleton_cmd}\
) | samtools view -b -o {{snakemake.output.mapped_reads}} - ) 2>{{snakemake.log}}"

print(f"Executing command:\n {base_cmd}", file=sys.stderr)

shell(base_cmd)

