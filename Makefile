install:
	pip install -e .

uninstall:
	pip uninstall crosslink-by-sequence && rm -rf ./crosslink_by_sequence.egg-info ./crosslink_by_sequence/bin

build-singularity-image:
	singularity build crosslink_by_sequence_singularity.sif crosslink_by_sequence_singularity.def