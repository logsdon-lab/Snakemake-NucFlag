# Snakemake-NucFlag
[![CI](https://github.com/logsdon-lab/Snakemake-NucFlag/actions/workflows/main.yml/badge.svg)](https://github.com/logsdon-lab/Snakemake-NucFlag/actions/workflows/main.yml)

A workflow to map PacBio HiFi reads back to an assembly and check for misassemblies at specific regions via [`NucFlag`](https://github.com/logsdon-lab/NucFlag).


### Getting Started
```bash
git clone git@github.com:logsdon-lab/Snakemake-NucFlag.git --recursive
cd Snakemake-NucFlag
```

### Configuration
Configuration is handled via `nucflag.toml`.

To read more, refer to the [`NucFlag` documentation](https://github.com/logsdon-lab/NucFlag/wiki/2.-Configuration).

#### Input
Files can be passed multiple ways in the `samples` section of `config.yaml`:

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

#### Configuration

##### General
General configuration can be filled in `config.yaml`:
```yaml
# Output directory
output_dir: "results/nucflag"
# Output pileup signals to {output_dir}/{sm}_pileup.
output_pileup: false
# Output ideogram
output_ideogram: false
# Output breakdown
output_breakdown: false
# Log directory
logs_dir: "logs/nucflag"
# Benchmarks directory
benchmarks_dir: "benchmarks/nucflag"
# Job resources. Memory in GB.
threads_aln: 8
mem_aln: 30G
processes_nucflag: 12
mem_nucflag: 50G
# samtools view filter flag.
samtools_view_flag: 2308
# Minimum length to display for ideogram.
ideogram_filter_length: 10000000
# Height of individual track.
ideogram_track_height: 2.0
# Minimum length to display for breakdown plot.
breakdown_filter_length: 10000000
# Show breakdown by percent or length.
breakdown_type: percent
```

##### By Sample
The following optional are optional per sample:
```yaml
samples: [
    {
        # Regions to check. If omitted, checks entire assembly.
        region_bed: "",
        # nucflag configuration
        config: "",
        # Regions to ignore.
        ignore_bed: "",
        # Regions to overlap.
        overlay_beds: [],
        # Cytobands
        # Only applicable with output_ideogram
        # Expects BED file (chrom, chromStart, chromEnd, name, bandType) used in pyideogram.
        cytobands: ""
    }
]
```

### Output
|Path|Description|
|-|-|
|`./{output_dir}/{sample}/{contig}.png`|Per-base coverage, mismatch, insertion, and deletion pileup plot with calls highlighted.|
|`./{output_dir}/{sample}_misassemblies.bed`|BED9 file with potential misassemblies with their coordinates per contig.|
|`./{output_dir}/{sample}_status.bed`|BED file with each chromosome, coordinates, status (Either `correct` or `misassembled`), and percent of each call within region.|
|`./{output_dir}/{sample}_pileup/{contig}.bw`|(Optional) BigWigs of above signals. Use `bigtools` to merge.|
|`./{output_dir}/{sample}_ideogram.(pdf/png)`|(Optional) Ideogram.|
|`./{output_dir}/{sample}_breakdown.(pdf/png)`|(Optional) Breakdown.|


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
            "region_bed": f"regions/{sm}_region.bed",
            "overlay_beds": []
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
