# ViReaDB
ViReaDB is a user-friendly database for storing reference-compressed viral sequence data and computing consensus genome sequences.

## Installation
ViReaDB can be installed using `pip`:


```bash
sudo pip install vireadb
```

If you are using a machine on which you lack administrative powers, ViReaDB can be installed locally using `pip`:

```bash
pip install --user vireadb
```

### Dependencies
ViReaDB is a Python package that depends on the [Pysam](https://github.com/pysam-developers/pysam) and [NumPy](https://numpy.org/) packages.

If your read mappings are not already in the CRAM format (i.e., they are in the SAM or BAM format), you will also need to install [Samtools](https://github.com/samtools/samtools) for conversion to CRAM. If your reads are in the FASTQ format (meaning they need to be mapped to the reference genome), you will also need to install [Minimap2](https://github.com/lh3/minimap2). Both tools must be in your `PATH`. The [Dockerfile](Dockerfile) in this repository may be helpful in setting up those tools.

## Usage
Typical usage should be as follows:

1. Import the `vireadb` package
2. Use `vireadb.create_db` to create a new database, or use `vireadb.load_db` to load an existing database
3. Add new datasets to the database, and/or compute/query consensus sequences of items in the database

```python
import vireadb
db = vireadb.load_db("example.db")
seq = db.get_consensus("my_sample")
```

Full documentation can be found at [https://niemasd.github.io/ViReaDB/](https://niemasd.github.io/ViReaDB/), and more examples can be found in the [ViReaDB Wiki](https://github.com/niemasd/ViReaDB/wiki).

## Citing ViReaDB
If you use ViReaDB in your work, please cite this GitHub repository (a manuscript is in preparation):

> Moshiri N (2022). "ViReaDB: A user-friendly database for compactly storing viral sequence data and rapidly computing consensus genome sequences." *bioRxiv*. [doi:10.1101/2022.10.21.513318](https://doi.org/10.1101/2022.10.21.513318)

# Tips and Tricks
The following are some tips or tricks when using ViReaDB.

## Convert to CRAM While Mapping
In many viral sequencing workflows, the first step is to trim reads (e.g. using [fastp](https://github.com/OpenGene/fastp)), and the second step is to map reads (e.g. using [Minimap2](https://github.com/lh3/minimap2)). Most often, folks pipe the output of their read mapper to [samtools](https://github.com/samtools/samtools) to convert to the [BAM format](https://en.wikipedia.org/wiki/Binary_Alignment_Map). However, the storage savings achieved by ViReaDB when storing raw reads is due to converting to the [CRAM format](https://en.wikipedia.org/wiki/CRAM_(file_format)), which ViReaDB performs using samtools internally. However, when ViReaDB is given a CRAM file, it simply stores the CRAM data as-is, without any additional conversion (which is much faster).

Thus, rather than converting your read mapping output to BAM on-the-fly and then having ViReaDB separately convert the BAM to CRAM, it is faster to simply convert the read mapping output to CRAM on-the-fly and then add the resulting CRAM to ViReaDB directly:

```bash
minimap2 -x sr -a <REF_GENOME_MMI> -t <THREADS> <R1_FASTQ> <R2_FASTQ> | samtools view --output-fmt-option version=3.0 --output-fmt-option use_lzma=1 --output-fmt-option archive=1 --output-fmt-option level=9 -C -T <REF_GENOME_FASTA> -@ <THREADS> --output-fmt-option lossy_names=1 -F 4 -o <OUTPUT_CRAM>
```
