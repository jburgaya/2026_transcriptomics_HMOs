from pathlib import Path
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("reads_path", type=Path, help="Path containing fastq.gz files")
args = parser.parse_args()

pattern = re.compile(r"(.+)_L1_(R[12])_001_.*\.fastq\.gz$")

for f in args.reads_path.glob("*.fastq.gz"):
    m = pattern.match(f.name)
    if m:
        new_name = f"{m.group(1)}_{m.group(2)}.fastq.gz"
        f.rename(f.with_name(new_name))
