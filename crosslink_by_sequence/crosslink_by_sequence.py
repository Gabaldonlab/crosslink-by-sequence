#!/usr/bin/env python3.10
"""
  crosslink-by-sequence - A tool to crosslink proteins between two proteomes
  from different sources using their sequences, aiming to avoid the inaccuracies
  associated with source-specific IDs.

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

from dataclasses import dataclass
import os
import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

from crosslink_by_sequence.argument_parser import CrosslinkBySequenceArgs
from crosslink_by_sequence.diamond_operations import crosslink_with_diamond
from crosslink_by_sequence.diamond_operations import make_diamond_db
from crosslink_by_sequence.utils import get_md5_fasta
from crosslink_by_sequence.utils import run_command_with_return
from crosslink_by_sequence.utils import time_it


VERBOSE: bool = False


@dataclass
class ProcessTaxidArgs:
    output_directory: str
    tmp_dir: str
    reference_file: str
    target_file: str
    reference_hashes: dict[str, list[str]]
    db_file_path: str
    minimum_coverage: float
    minimum_identity: float
    verbose: bool


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
        reference_id = reference_hashes[md5]
        for target_id in target_hashes[md5]:
            for metIDEl in reference_id:
                target_to_reference = pd.concat(
                    [
                        target_to_reference,
                        pd.DataFrame(
                            [[target_id, "1", metIDEl]],
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


def process_taxid(process_taxid_args: ProcessTaxidArgs) -> str:
    db_id: int = 0
    version: int = 0
    file_name: str = os.path.basename(process_taxid_args.target_file)
    file_prefix: str = os.path.splitext(file_name)[0]
    if VERBOSE:
        print(f"[INFO] Starting to process: {file_prefix}")

    target_hashes: dict[str, list[str]] = get_md5_fasta(process_taxid_args.target_file)
    target_to_reference_tmp: pd.DataFrame = crosslink_with_md5(
        target_hashes, process_taxid_args.reference_hashes
    )

    number_all: list[str] = run_command_with_return(
        f"zcat {process_taxid_args.target_file} | grep -c '>'"
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
            process_taxid_args.tmp_dir,
            process_taxid_args.reference_file,
            process_taxid_args.target_file,
            process_taxid_args.db_file_path,
            process_taxid_args.minimum_coverage,
            process_taxid_args.minimum_identity,
            1,
            process_taxid_args.verbose,
        )

    target_to_reference_tmp["dbid"] = db_id
    target_to_reference_tmp["version"] = version
    target_to_reference_tmp = target_to_reference_tmp[
        ["dbid", "extid", "version", "protid", "score"]
    ]
    target_to_reference_tmp.drop_duplicates(inplace=True)

    result_file_path: str = os.path.join(
        process_taxid_args.output_directory, f"{file_prefix}.target_to_reference.tbl.gz"
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
        process_taxid_args.output_directory,
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
    max_threads: int,
    run_as_sync: bool,
    process_taxid_args: list[ProcessTaxidArgs],
) -> list[str]:
    processes_results: list[str] = []
    if run_as_sync:
        for args in process_taxid_args:
            result: str = process_taxid(args)
            processes_results.append(result)
    else:
        with ProcessPoolExecutor(max_workers=max_threads) as executor:
            processes_results = list(
                executor.map(process_taxid, process_taxid_args)
            )
    return processes_results


@time_it
def main() -> int:
    args: CrosslinkBySequenceArgs = CrosslinkBySequenceArgs.get_arguments()
    VERBOSE = args.verbose
    if VERBOSE:
        print(f"Options: {str(args)}")

    Path(args.output_directory).mkdir(exist_ok=True, parents=True)
    Path(args.tmp_directory).mkdir(exist_ok=True, parents=True)

    reference_hashes: dict[str, list[str]] = get_md5_fasta(
        args.target_reference_species_fasta_gzip_file
    )
    db_file_path: str = make_diamond_db(
        args.tmp_directory, args.target_reference_species_fasta_gzip_file, VERBOSE
    )

    process_taxid_args: list[ProcessTaxidArgs] = [
        ProcessTaxidArgs(
            output_directory=args.output_directory,
            tmp_dir=args.tmp_directory,
            reference_file=args.target_reference_species_fasta_gzip_file,
            target_file=fasta_file,
            reference_hashes=reference_hashes,
            db_file_path=db_file_path,
            minimum_coverage=args.minimum_coverage,
            minimum_identity=args.minimum_identity,
            verbose=VERBOSE,
        )
        for fasta_file in args.target_fasta_gzip_files
    ]

    processes_results = compute(
        args.max_threads, args.run_as_sync, process_taxid_args
    )

    print_parallelized_processes_logs(
        args.target_fasta_gzip_files, processes_results
    )

    if not args.keep_tmp_directory:
        shutil.rmtree(args.tmp_directory)

    return 0
