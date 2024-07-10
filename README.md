# crosslink-by-sequence

## 1. Setups on Unix systems

### 1.1 Install from Github

```bash
pip install https://github.com/Gabaldonlab/crosslink-by-sequence/archive/main.zip
```

### 1.2 Install from source

```bash
git clone https://github.com/gabaldonlab/crosslink-by-sequence
cd crosslink-by-sequence
make install
```

### 1.3 Install in a fresh virtual environment

```bash
git clone https://github.com/gabaldonlab/crosslink-by-sequence
cd crosslink-by-sequence
source scripts/install-with-fresh-env.sh
```

### 1.4 Uninstall

```bash
make uninstall
```

---

## 2. Examples

### 2.1 Command xample

```bash
crosslink-by-sequence \
		--target_reference_species_fasta_gzip_file ./test_data/input_data/reference_proteomes/8.7165.faa.gz \
		--target_fasta_gzip_files ./test_data/input_data/target_proteomes/0.7165.fasta.gz \
		--output_directory ./test_data/output_data \
		--tmp_directory ./test_data/output_data/tmp \
		--max_threads 4 \
		--minimum_coverage 0.5 \
		--minimum_identity 0.5 \
        --verbose
```

### 2.2 Output example

```bash
dbid	extid	version	protid	score
0	tr|A0NAB3|A0NAB3_ANOGA	0	41948896	1
0	tr|A0NAZ1|A0NAZ1_ANOGA	0	41948784	1
0	tr|A0NBF5|A0NBF5_ANOGA	0	41951081	1
0	tr|A0NC40|A0NC40_ANOGA	0	41952308	1
0	tr|A0NEB5|A0NEB5_ANOGA	0	41949551	1
0	tr|A0NEL3|A0NEL3_ANOGA	0	41950072	1
0	tr|A0NH68|A0NH68_ANOGA	0	41956125	1
0	tr|A7URP6|A7URP6_ANOGA	0	41950891	1
0	tr|A7URS0|A7URS0_ANOGA	0	41950804	1
etc...
etc...
etc...
```

---

## 3. Singularity image

## 3.1. Build Singularity image for HPC cluster

**NOTE**: Inside the HPC you won't have sudo permissions to build the image. So you will have to build it locally on your working machine, then copy it to the remote one.

```bash
sudo make build-singularity-image
```

## 3.2. Run Singularity image example

```bash
singularity run --cleanenv crosslink_by_sequence_singularity.sif crosslink-by-sequence \
		--target_reference_species_fasta_gzip_file ./test_data/input_data/reference_proteomes/8.7165.faa.gz \
		--target_fasta_gzip_files ./test_data/input_data/target_proteomes/0.7165.fasta.gz \
		--output_directory ./test_data/output_data \
		--tmp_directory ./test_data/output_data/tmp \
		--max_threads 4 \
		--minimum_coverage 0.5 \
		--minimum_identity 0.5 \
		--verbose
```

---

## 4. Acknowledgments

### 4.1. Diamond v2.1.9.163 (C) Max Planck Society for the Advancement of Science, Benjamin Buchfink, University of Tuebingen

This tool uses the Diamond software for sequence alignment.
Please cite the following when using Diamond in your work:

-   (Documentation, support and updates available at)[http://www.diamondsearch.org]
-   Please cite: (Nature Methods (2021))[http://dx.doi.org/10.1038/s41592-021-01101-x]
