# vireadb
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
In addition to installing ViReaDB itself (which is a Python package), you will also need to install its dependencies, [samtools](https://github.com/samtools/samtools) and [minimap2](https://github.com/lh3/minimap2), both of which must be in your `PATH`. The [Dockerfile](Dockerfile) in this repository may be helpful in setting up those tools.
