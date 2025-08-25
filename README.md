

#  README

## 1. Data sources

This task is based on on publicly available sequencing data from the study "“Transcriptional landscape of the cell cycle in a model thermoacidophilic archaeon reveals similarities to eukaryotes”which focused on understanding the gene expression of Saccharolobus islandicus relative to its different cell-cycle stages. Cultures of S. islandicus were made to begin the cell cycle at once by a 6-hour treatment with acetic acid. After the removal of the acetic acid, the cells resumed their cell cycle. Total RNA was extracted from 3 samples at three specific time points after synchronization: sample 1 (2 hours and 30 minutes (in M-G1 phase)), sample 2 (4 hours (in S phase)), sample 3 (6 hours( in G2 phase)).The samples were sequenced using an Illumina NextSeq 2000 sequencer.
The subsampled and cleaned FASTQs are stored in `data/` and are used as the inputs for the workflow.

---

## 2. How to download

The samples were obtained by first fetching the sample’s SRA record (SRX) on GEO using geofetch, then downloaded the fastq of the selected samples using their corresponding SRR accession.
### Code for downloading

```bash
geofetch -i GSE296035 --just-metadata
SRRS=("SRX28623476" "	SRX28623471" "SRX28623484")

for SRR in "${SRRS[@]}"; do
    echo "Downloading $SRR ..."
    prefetch "$SRR"
    fastq-dump --gzip --split-files "$SRR"
done
```


---

## 3. Pre-processing / subsampling

INCLUDE THE METHOD YOU USED TO SUBSAMPLE, MINATURIZE, OR TRIM DOWN

1. **STEP 1** ...

Example:

```bash
CODE TO SUBSAMPLE
```


---

## 4. How the workflow works
The workflow files is stored in workflow/ and it is divided into different steps:
DESCRIBE THE WORKFLOW HERE - NOTE THE BELOW ARE JUST EXAMPLES, REPLACE WITH YOUR OWN - YOURS CAN TAKE A VERY DIFFERENT FORMAT
The workflow files is stored in `workflow/`.

---

### Step 1 – Quality Control (example)

**Purpose:** Remove low-quality reads and adapter sequences
**Tools:** `fastp`, `cutadapt`, `trimmomatic`
**Inputs:** Subsampled FASTQ files (from `data/fastq_subsampled/`)
**Outputs:** Cleaned FASTQs, QC reports (`.html`, `.json`, or `.txt`)
**Command:**

```bash
fastp --in1 sample.fastq.gz --out1 cleaned.fastq.gz ...
```

---

### Step 2 ...

**Purpose:** ...
**Tools:** ...
**Inputs:** ...
**Outputs:** ...
**Command:**


---

### Step X – Analysis (e.g., DESeq2, variant calling, etc.)

**Purpose:** ...
**Tools:** ...
**Inputs:** ...
**Outputs:** ...
**Command:**

