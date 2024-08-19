I'm using the annotations that Yupei provided me here @ /data/gpfs/projects/punim2251/LongBench_data/ont_sc/cell_line_annotation

I've split samples up into batches based on the number of cells identified by BLAZE

Batch 1:
H69 - 449               4,419,109 total reads
H146 - 486              4,865,297 total reads
SHP77 - 473             6,440,158 total reads
HCC827 - 514

Batch 2:
H211 - 352
H526 - 265

Batch 3:
H1975 - 722
H2228 - 795


Batch 4 - Subsampled to 473 barcodes:
H2228             11,989,938 reads            15k reads/cell            
H1975             12,460,762 reads            16k reads/cell            
SHP77             6,440,158 reads             12k reads/cell            
HCC827            14,018,210 reads            28k reads/cell            ** Start with this


Rarefaction breakpoints:
1) 100%
2) 75%
3) 50%
4) 25%
5) 10%

Pipeline for Nextflow:
subsampleBarcodeList.py -> extractReadsForBarcodeList.py -> rarefyMatchedFastQs.py