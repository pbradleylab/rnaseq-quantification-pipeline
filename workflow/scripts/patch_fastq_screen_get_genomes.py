"""Patch FastQ Screen's bundled genome downloader to use HTTPS URLs."""

from pathlib import Path
import sys


script = Path(sys.argv[1])
text = script.read_text()

replacements = {
    "www.bioinformatics.babraham.ac.uk/projects/fastq_screen/genome_locations.txt": (
        "https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/genome_locations.txt"
    ),
    'wget --no-check-certificate -r --no-parent -R \'index.html*\' $download_folder': (
        'wget --no-check-certificate -r --no-parent -R \'index.html*\' https://$download_folder'
    ),
}

for old, new in replacements.items():
    if old not in text:
        raise SystemExit(f"Expected FastQ Screen downloader text not found: {old}")
    text = text.replace(old, new, 1)

script.write_text(text)
