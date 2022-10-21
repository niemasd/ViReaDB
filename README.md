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
ViReaDB is a Python package that depends on the [pysam](https://github.com/pysam-developers/pysam) package.

If your read mappings are not already in the CRAM format (i.e., they are in the SAM or BAM format), you will also need to install [samtools](https://github.com/samtools/samtools) for conversion to CRAM. If your reads are in the FASTQ format (meaning they need to be mapped to the reference genome), you will also need to install [minimap2](https://github.com/lh3/minimap2). Both tools must be in your `PATH`. The [Dockerfile](Dockerfile) in this repository may be helpful in setting up those tools.

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

> https://github.com/niemasd/ViReaDB
