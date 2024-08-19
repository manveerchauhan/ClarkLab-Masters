# Adapted from: https://www.reddit.com/r/bioinformatics/comments/1c57c28/whats_the_fastest_pythonic_way_to_read_fastq/?rdt=49974&onetap_auto=true&one_tap=true
# To use: python extractReadsForBarcodeList.py barcodes.txt output_pref sample1.fastq.gz sample2.fastq.gz

import sys
import gzip
import multiprocessing

def readfq(fastq):
    state = -1
    lines = ["", "", "", ""]
    for line in fastq:
        if state < 1 and line.startswith("@"):
            lines[0] = line
            state = 1
        elif state < 4:
            lines[state] = line
            state += 1
            if 4 == state:
                state = 0
                yield lines

def load_barcodes(barcode_file):
    with open(barcode_file, "r") as f:
        barcodes = set(line.strip() for line in f)
    return barcodes

def barcode_filter(lines, valid_barcodes):
    # Extract barcode from the header line
    header = lines[0].strip()
    barcode = header.split()[0][1:17]  # Extract the first 16 nucleotides after '@'

    # Check if the barcode is in the valid barcode set
    if barcode in valid_barcodes:
        return True
    return False

def processFastq(fn, valid_barcodes, output_fn):
    with gzip.open(fn, "rt") if fn.endswith('.gz') else open(fn, "rt") as fastq:
        with open(output_fn, "w") as output_fastq:
            for lines in readfq(fastq):
                if barcode_filter(lines, valid_barcodes):
                    output_fastq.write("".join(lines))

def main(barcode_file, output_prefix, *fastq_files):
    valid_barcodes = load_barcodes(barcode_file)
    
    with multiprocessing.Pool() as pool:
        pool.starmap(processFastq, [(fn, valid_barcodes, f"{output_prefix}_{fn.split('/')[-1]}") for fn in fastq_files])

if __name__ == "__main__":
    # Usage: python script.py barcode.txt output_prefix file1.fastq.gz file2.fastq.gz ...
    main(sys.argv[1], sys.argv[2], *sys.argv[3:])
