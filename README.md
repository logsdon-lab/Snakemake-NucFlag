# Snakemake-NucFlag
A workflow to map PacBio HiFi reads back to an assembly and check for misassemblies at specific regions via [`NucFlag`](https://github.com/logsdon-lab/NucFlag).


### Getting Started
```bash
git clone git@github.com:logsdon-lab/Snakemake-NucFlag.git
```

#### Configuration
Files can be passed multiple ways:

##### Assemblies
By path.
```yaml
samples:
-   name: "1"
    asm_fa: "1.fa"
```

By `fofn`.
```yaml
samples:
-   name: "1"
    asm_fofn: "1.fofn"
```

By directory and file extension.
```yaml
samples:
-   name: "1"
    asm_dir: "1/"
    asm_ext: "fa.gz"
```

##### Reads
By `fofn`.
```yaml
samples:
-   name: "1"
    read_fofn: "1.fofn"
```

By directory and file extension.
```yaml
samples:
-   name: "1"
    read_dir: "1/"
    read_ext: "bam"
```

### Usage
```bash
snakemake -np -c 1 --configfile config/config.yaml
```

### Module
To incorporate this into a workflow.

```python
SAMPLE_NAMES = ["sample_1"]
NUCFLAG_CFG = {
    "samples": [
        {
            "name": sm,
            "asm_fa": f"{sm}.fa",
            "read_dir": f"reads/{sm}/",
            "read_ext": "bam",
            "region_bed": f"regions/{sm}_region.bed"
        }
        for sm in SAMPLE_NAMES
    ],
    # Other nucflag parameters
    **config["nucflag"]
}


module NucFlag:
    snakefile:
        github(
            "logsdon-lab/Snakemake-NucFlag",
            path="workflow/Snakefile",
            branch="main"
        )
    config: NUCFLAG_CFG

use rule * from NucFlag as *

rule all:
    input:
        expand(rules.nucflag.input, sm=SAMPLE_NAMES),
```

### Output

|Path|Description|
|-|-|
|`./{sample}/{contig}.png`|Per-base coverage graph plot with misassemblies highlighted.|
|`./{sample}_cen_misassemblies.bed`|Bed file with misassemblies and their coordinates per contig.|
|`./{sample}_cen_status.bed`|Bed file with each centromeric contig, coordinates, and status. Either `good` or `misassembled`.|
