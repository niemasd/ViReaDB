# vireadb
vireadb: Viral Read Database

# Notes

## Convert BAM to CRAM
```
samtools view --output-fmt-option lossy_names=1 --output-fmt-option archive=1 --output-fmt-option use_lzma=1 --output-fmt-option version=3.0 --output-fmt-option level=9 -@ 8 -F 4 -T NC_045512.2.fas -C -o example_untrimmed_sorted.cram example_untrimmed_sorted.bam
```
