# Global multi-omics reveals that phage–host arms race drives resistome dynamics in livestock manure composting

This repository contains the bioinformatics pipeline and custom scripts for the research article: **"Global multi-omics reveals that phage–host arms race drives resistome dynamics in livestock manure composting"**.

We integrated **420 metagenomes （9.1T）** and **42 viromes （840Gb）** to construct the **Global Viral Dataset (GVD_LMC)** and **Global Antiviral Defense System Catalogue (GADSC_LMC)**, revealing the mechanisms governing resistome persistence.

---

## 📂 Repository Structure

The workflow is modularized into 6 steps, corresponding to the Materials & Methods section of the manuscript.

```text
.
├── scripts/                  # Custom Python/R helper scripts
├── 01_read_profiling.sh      # Step 1: Reads QC, ARGs & Taxonomy Profiling
├── 02_assembly_binning.sh    # Step 2: Assembly, MAG Recovery & Annotation
├── 03_viral_discovery.sh     # Step 3: Viral Identification (GVD_LMC)
├── 04_viral_taxonomy.sh      # Step 4: Viral Taxonomy & Lifestyle
├── 05_host_ads_network.sh    # Step 5: Host Assignment & Immune Network
├── 06_visualization.R        # Step 6: Figure Generation & Statistics
└── README.md                 # Documentation

```

---

## 🛠 Prerequisites

This pipeline is designed for Linux environments. We strongly recommend using **Conda** to manage dependencies.

### Core Tools Required

Ensure the following tools are installed and accessible in your `$PATH`:

* **Reads & MAGs:** `metaWRAP`, `Diamond`, `MetaPhlAn`, `Prodigal`, `GTDB-Tk`, `dRep`, `CoverM`.
* **Viral Analysis:** `VirSorter2`, `VirFinder`, `geNomad`, `CheckV`, `vConTACT2`.
* **Host & Defense:** `BLAST+`, `PADLOC`, `CRT`.
* **Visualization:** R (v4.0+) with `ggplot2`, `pheatmap`, `vegan`, `UpSetR`.

---

## 🚀 Pipeline Usage

### Step 1: Read-based Profiling

Quality control of raw reads and profiling of ARGs (CARD) and microbial communities (MetaPhlAn).


### Step 2: Genome-resolved Analysis (MAGs)

Assembly, hybrid binning (metaBAT2/MaxBin2/CONCOCT), dereplication, and taxonomic classification of MAGs.


### Step 3: Viral Discovery (GVD_LMC)

This module implements a **multi-tool ensemble approach** (VirSorter2 + VirFinder + geNomad) followed by strict filtering and CheckV quality control to construct the GVD_LMC.


### Step 4: Viral Taxonomy & Lifestyle

Assigns taxonomy to vOTUs and predicts viral lifestyles (Lytic vs. Temperate).


### Step 5: Host Assignment & Immune Network

Constructs the phage-host interaction network (CRISPR/tRNA/Homology) and profiles Antiviral Defense Systems (ADS) using PADLOC.


### Step 6: Visualization & Statistics

Reproduces the main figures (Global Map, Upset Plot, PCoA, Heatmaps) and statistical tests (PERMANOVA).


---

## 📊 Data Availability

* **Newly Generated Raw Sequencing Data:**  Available at  NGDC (National Genomics Data Center, https://ngdc.cncb.ac.cn/) under CRA034768.
* **Processed Data (GVD_LMC & GADSC_LMC):** Available at figshare **(10.6084/m9.figshare.31272853)**.


---

## 📧 Contact

For questions regarding the code, please open an issue or contact:
**[Junya Zhang]** (jyzhang@rcees.ac.cn)
