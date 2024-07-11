#!/usr/bin/env python3.10
"""
  phylomizer - automated phylogenetic reconstruction pipeline - it resembles the
  steps followed by a phylogenetist to build a gene family tree with error-
  control of every step

  Copyright (C) 2024 - Toni Gabaldon, Dániel Májer, Diego Fuentes, Anna Vlasova

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from __future__ import annotations

import gzip
import hashlib
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import Callable
from typing import Generator
from typing import IO

from crosslink_by_sequence.diamond_operations import (
    crosslink_with_diamond,
    make_diamond_db,
)
from crosslink_by_sequence.utils import (
    get_md5_fasta,
    run_command_with_return,
    time_it,
)
import pandas as pd
from Bio import SeqIO  # type: ignore

from crosslink_by_sequence.argument_parser import CrosslinkBySequenceArgs


VERBOSE: bool = False


def crosslink_with_md5(
    target_hashes: dict[str, list[str]], reference_hashes: dict[str, list[str]]
) -> pd.DataFrame:
    """
    Crosslinks the reference and target sequences by their MD5 hashes.

    @return pd.DataFrame (target to reference)
        columns=["extid", "score", "protid"]
    """
    target_to_reference: pd.DataFrame = pd.DataFrame(
        columns=["extid", "score", "protid"]
    )
    shared_md5_keys: list[str] = [
        md5 for md5 in target_hashes if md5 in reference_hashes
    ]
    for md5 in shared_md5_keys:
        # metaID is a list, can be a few elements
        metaID = reference_hashes[md5]
        for extID in target_hashes[md5]:
            # add ext2meta entry
            for metIDEl in metaID:
                target_to_reference = pd.concat(
                    [
                        target_to_reference,
                        pd.DataFrame(
                            [[extID, "1", metIDEl]],
                            columns=target_to_reference.columns,
                        ),
                    ]
                )
    return target_to_reference


def write_log_file(
    output_directory: str,
    file_prefix: str,
    number_all_result: int,
    number_matched: int,
    orphans_number: int,
    orphans_percentage: float,
) -> str:
    """
    Writes to log file as a .tsv line.

    @return str
        Path of the created log file.
    """
    log_file_path: str = os.path.join(output_directory, f"{file_prefix}.log")
    with open(log_file_path, "wt", encoding="UTF-8") as opened_log_file:
        opened_log_file.write(
            f"{file_prefix}\t{number_all_result}\t{number_matched}\t{orphans_number}\t{orphans_percentage:.2f}\n"
        )
    return log_file_path


def process_taxid(
    output_directory: str,
    tmp_dir: str,
    reference_file: str,
    target_file: str,
    reference_hashes: dict[str, list[str]],
    minimum_coverage: float,
    minimum_identity: float,
    verbose: bool,
) -> str:
    db_id: int = 0
    version: int = 0
    file_name: str = os.path.basename(target_file)
    file_prefix: str = os.path.splitext(file_name)[0]
    if VERBOSE:
        print(f"[INFO] Starting to process: {file_prefix}")

    target_hashes: dict[str, list[str]] = get_md5_fasta(target_file)
    target_to_reference_tmp: pd.DataFrame = crosslink_with_md5(
        target_hashes, reference_hashes
    )

    number_all: list[str] = run_command_with_return(
        f"zcat {target_file} | grep -c '>'"
    )
    number_all_result: int = int(number_all[0])
    number_matched: int = target_to_reference_tmp["extid"].nunique()

    leftover: int = number_all_result - number_matched
    if VERBOSE:
        print(
            f"Number of leftovers from crosslinking by MD5 hashes: {leftover}"
        )

    has_not_matched_by_md5: bool = leftover > 0
    if has_not_matched_by_md5:
        target_to_reference_tmp = crosslink_with_diamond(
            target_to_reference_tmp,
            tmp_dir,
            reference_file,
            target_file,
            minimum_coverage,
            minimum_identity,
            1,
            verbose,
        )

    target_to_reference_tmp["dbid"] = db_id
    target_to_reference_tmp["version"] = version
    target_to_reference_tmp = target_to_reference_tmp[
        ["dbid", "extid", "version", "protid", "score"]
    ]
    target_to_reference_tmp.drop_duplicates(inplace=True)

    result_file_path: str = os.path.join(
        output_directory, f"{file_prefix}.target_to_reference.tbl.gz"
    )
    target_to_reference_tmp.to_csv(
        result_file_path,
        compression="gzip",
        index=False,
        sep="\t",
        header=True,
    )

    # Calculate number of matched and non-matched proteins.
    number_matched = target_to_reference_tmp["extid"].nunique()
    orphans_number: int = number_all_result - number_matched
    orphans_percentage: float = orphans_number * 100 / number_all_result

    write_log_file(
        output_directory,
        file_prefix,
        number_all_result,
        number_matched,
        orphans_number,
        orphans_percentage,
    )

    return file_prefix


def print_parallelized_processes_logs(
    target_fasta_files: list[str], parallelized_processes_results: list[str]
) -> None:
    for idx, file_name_prefix in enumerate(parallelized_processes_results, 1):
        if not file_name_prefix:
            continue
        if VERBOSE:
            print(f" {idx} / {len(target_fasta_files)}  {file_name_prefix}")
            print(f"[INFO] done {file_name_prefix} ")


def compute(
    target_fasta_files: list[str],
    reference_file: str,
    output_directory: str,
    tmp_directory: str,
    max_threads: int,
    minimum_coverage: float,
    minimum_identity: float,
    run_as_sync: bool,
) -> None:
    reference_hashes: dict[str, list[str]] = get_md5_fasta(reference_file)

    make_diamond_db(tmp_directory, reference_file)

    process_taxid_args: list[
        tuple[str, str, str, str, dict[str, list[str]], float, float, bool]
    ] = [
        (
            output_directory,
            tmp_directory,
            reference_file,
            fasta_file,
            reference_hashes,
            minimum_coverage,
            minimum_identity,
            VERBOSE,
        )
        for fasta_file in target_fasta_files
    ]

    if run_as_sync:
        parallelized_processes_results: list[str] = []
        for args in process_taxid_args:
            result: str = process_taxid(
                args[0],
                args[1],
                args[2],
                args[3],
                args[4],
                args[5],
                args[6],
                args[7],
            )
            parallelized_processes_results.append(result)
    else:
        with ProcessPoolExecutor(max_workers=max_threads) as executor:
            parallelized_processes_results: list[str] = list(
                executor.map(process_taxid, process_taxid_args)
            )

    print_parallelized_processes_logs(
        target_fasta_files, parallelized_processes_results
    )


@time_it
def main() -> int:
    args: CrosslinkBySequenceArgs = CrosslinkBySequenceArgs.get_arguments()
    VERBOSE = args.verbose
    if VERBOSE:
        print(f"Options: {str(args)}")

    # Create out and tmp directory if they dont exists.
    Path(args.output_directory).mkdir(exist_ok=True, parents=True)
    Path(args.tmp_directory).mkdir(exist_ok=True, parents=True)

    # Crosslink proteomes.
    compute(
        args.target_fasta_gzip_files,
        args.target_reference_species_fasta_gzip_file,
        args.output_directory,
        args.tmp_directory,
        args.max_threads,
        args.minimum_coverage,
        args.minimum_identity,
        args.run_as_sync,
    )

    if not args.keep_tmp_directory:
        shutil.rmtree(args.tmp_directory)

    return 0
