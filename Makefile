install:
	pip install -e .

uninstall:
	pip uninstall crosslink-by-sequence && rm -rf ./crosslink_by_sequence.egg-info ./crosslink_by_sequence/bin

build-singularity-image:
	singularity build crosslink_by_sequence_singularity.sif crosslink_by_sequence_singularity.def

run-dummy:
	crosslink-by-sequence \
		--target_reference_species_fasta_gzip_file ./test_data/input_data/reference_proteomes/8.7165.faa.gz \
		--target_fasta_gzip_files ./test_data/input_data/target_proteomes/0.7165.fasta.gz \
		--output_directory ./test_data/output_data \
		--tmp_directory ./test_data/output_data/tmp \
		--max_threads 4 \
		--minimum_coverage 0.5 \
		--minimum_identity 0.5

run-dummy-as-module:
	python3 -m crosslink_by_sequence \
		--target_reference_species_fasta_gzip_file ./test_data/input_data/reference_proteomes/8.7165.faa.gz \
		--target_fasta_gzip_files ./test_data/input_data/target_proteomes/0.7165.fasta.gz \
		--output_directory ./test_data/output_data \
		--tmp_directory ./test_data/output_data/tmp \
		--max_threads 4 \
		--minimum_coverage 0.5 \
		--minimum_identity 0.5

