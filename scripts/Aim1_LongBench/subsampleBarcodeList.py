# adapted from https://www.geeksforgeeks.org/python-random-sample-function/
# Usage: python subsampleBarcodeList.py barcode.txt subsampleNum outputFileName

import sys
from random import sample

def load_barcodes(barcode_file):
    with open(barcode_file, "r") as f:
        barcodes = list(line.strip() for line in f)
    return barcodes

def subsampleToDesiredLength(barcodes, n):
    subsampledList = sample(barcodes, n)
    print("Successful downsampled barcode list to length: ", n)
    return subsampledList

def exportSubsampledBarcodes(BC_List_Subset, fileName):
    outputFile = fileName
    # Open the file in write mode
    with open(outputFile, "w") as file:
        # Iterate over each element in the list
        for element in BC_List_Subset:
            # Write the element to the file, followed by a newline character
            file.write(f"{element}\n")
            
def main(originalBarcodeList, desiredBarcodeNum, outputPrefix):
    OriginalBarcodes = load_barcodes(originalBarcodeList)
    subsampledBarcodes = subsampleToDesiredLength(OriginalBarcodes, int(desiredBarcodeNum))
    
    outputFileName = f"{outputPrefix}_{desiredBarcodeNum}BCs.txt"
    exportSubsampledBarcodes(subsampledBarcodes, outputFileName)
    print("Successfully exported subsampled list to: ", outputFileName)
    

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])