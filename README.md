# crosslink-by-sequence

## 1. Install
```bash
pip install -e ./crosslink_by_sequence
```

## 2. Command example
```bash
python3 /gpfs/projects/bsc40/current/avlasova/projects/metaphors/scripts_pipeline/xlink_analyze_blat.py \
    --output_directory /gpfs/projects/bsc40/dfuentes/projects/qfo_2022/2022_vs_metaphorsIDs_0.5_0.5/ \
    --tmp_directory /gpfs/projects/bsc40/current/dfuentes/projects/qfo_2022/fasta.xlink/tmp_5psl/ \
    --max_threads 4 \
    --refSpecie /gpfs/projects/bsc40/project/pipelines/metaphors/metaphors-db-2019/data/fasta.reference/8.9612.faa.gz \
    --minimum_coverage 0.5 \
    --minimum_identity 0.5 \
    --target_fasta_gzip_files /gpfs/projects/bsc40/current/dfuentes/projects/qfo_2022/2022_proteomes/0.9615.fasta.gz
```
