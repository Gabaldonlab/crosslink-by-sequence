# crosslink-by-sequence

## 1. Setups on Unix systems

### 1.1 Install from source

```bash
git clone https://github.com/gabaldonlab/crosslink-by-sequence
cd crosslink-by-sequence
make install
```

### 1.2 Install in a fresh virtual environment
```bash
source scripts/install-with-fresh-env.sh
```

### 1.2 Uninstall

```bash
make uninstall
```

## 2. Command example

```bash
crosslink-by-sequence \
    --target_reference_species_fasta_gzip_file ./crosslink-by-sequence/test_data/input_data/reference_proteomes/8.9612.faa.gz \
    --target_fasta_gzip_files ./crosslink-by-sequence/test_data/input_data/target_proteomes/0.9615.fasta.gz \
    --output_directory ./crosslink-by-sequence/test_data/output_data \
    --tmp_directory ./crosslink-by-sequence/output_data/tmp \
    --max_threads 4 \
    --minimum_coverage 0.5 \
    --minimum_identity 0.5
```

## 3. Singularity image

## 3.1. Build Singularity image for HPC cluster

**NOTE**: Inside the HPC you won't have sudo permissions to build the image. So you will have to build it locally on your working machine, then copy it to the remote one.

```bash
sudo make build-singularity-image
```

## 3.2. Run Singularity image

```bash
singularity run --cleanenv crosslink_by_sequence_singularity.sif crosslink-by-sequence < ...rest_of_the_arguments... >
```
