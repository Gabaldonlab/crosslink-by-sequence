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

from dataclasses import dataclass
import gzip
import hashlib
import os
import subprocess as sp
import sys
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path
from typing import IO, Any, Generator
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from Bio import SeqIO  # type: ignore

from crosslink_by_sequence.argument_parser import CrosslinkBySequenceArgs


VERBOSE: bool = False


@dataclass
class DiamondOutputFileLine:
    query_id: str  # Query ID
    subject_id: str  # Subject ID
    identity: float  # Percentage of identical matches
    alignment_length: int  # Alignment length
    mismatches_number: str  # Number of mismatches
    gap_openings_number: str  # Number of gap openings
    query_alignment_start: str  # Start of alignment in query
    query_alignment_end: str  # End of alignment in query
    subject_alignment_start: str  # Start of alignment in subject
    subject_alignment_end: str  # End of alignment in subject
    expected_value: str  # Expected value
    bit_score: str  # Bit score

    @classmethod
    def yield_output_file_lines(
        cls,
        output_file_path: str,
    ) -> Generator[DiamondOutputFileLine, None, None]:
        with open(
            output_file_path, "r", encoding="UTF-8"
        ) as opened_diamond_file:
            for line in opened_diamond_file:
                stripped_line: str = line.strip()
                if not stripped_line or stripped_line.startswith("#"):
                    continue
                    # skip lines with errors
                split_line: list[str] = stripped_line.split("\t")
                if len(split_line) != 12:
                    if VERBOSE:
                        print(
                            f"  Error: blat2xlinks:"
                            " Expected 12 elements in line,"
                            f" got {len(split_line)} {str(split_line)}!",
                            file=sys.stderr,
                        )
                    continue

                mapped_line = DiamondOutputFileLine(
                    query_id=split_line[0],
                    subject_id=split_line[1],
                    identity=float(split_line[2]) / 100,
                    alignment_length=int(split_line[3]),
                    mismatches_number=split_line[4],
                    gap_openings_number=split_line[5],
                    query_alignment_start=split_line[6],
                    query_alignment_end=split_line[7],
                    subject_alignment_start=split_line[8],
                    subject_alignment_end=split_line[9],
                    expected_value=split_line[10],
                    bit_score=split_line[11],
                )
                yield mapped_line


def get_md5_fasta(fasta_file_path: str) -> dict[str, list[str]]:
    """
    @return dict[str, list[str]]
        { FASTA_md5: [ protid1, protid2, ... ] }
    """
    hash_to_protein: dict[str, list[str]] = {}
    with gzip.open(fasta_file_path, mode="rt") as opened_gzip_file:
        for fasta_record in SeqIO.parse(opened_gzip_file, "fasta"):
            external_id: str = str(fasta_record.id)
            seq: str = str(fasta_record.seq)
            md5: str = hashlib.md5(seq.encode("utf-8")).hexdigest()
            hash_to_protein.setdefault(md5, [])
            hash_to_protein[md5].append(external_id)
    return hash_to_protein


def make_diamond_db(tmp_dir: str, reference_file: str) -> None:
    db_file_name: str = f"{os.path.basename(reference_file)}db"
    output_file_path: str = os.path.join(tmp_dir, db_file_name)
    diamond_bin_path: str = str(Path(__file__).parent / "bin" / "diamond")
    cmd: str = (
        f"{diamond_bin_path} makedb -d {output_file_path} --in {reference_file} "
    )
    if not os.path.isfile(output_file_path + ".dmnd"):
        run_shell_command(cmd, False)
        print(cmd)
    else:
        print("[INFO] Diamond DB file exists for ", reference_file)


def get_sequence_lengths(fasta_file_path: str) -> dict[str, int]:
    """
    @return dict[str, int]
        { external_id: <sequence_length: int> }
    """
    with gzip.open(fasta_file_path, mode="rt") as handle:
        sequence_lengths: dict[str, int] = {
            str(fasta_record.id): len(str(fasta_record.seq))
            for fasta_record in SeqIO.parse(handle, "fasta")
        }
    return sequence_lengths


def run_diamond(
    tmp_dir: str, missing_file_name: str, threads: int, reference_file: str
) -> str:
    """
    Runs Diamond as a process with the given parameters.

    @return str - The output file's path.
    """
    db_file_name: str = f"{os.path.basename(reference_file)}db"
    database_file_path: str = os.path.join(tmp_dir, db_file_name)
    output_file_name: str = f"{os.path.basename(missing_file_name)}diamond"
    output_file_path: str = os.path.join(tmp_dir, output_file_name)

    # -tileSize=8 caused `Internal error genoFind.c 2225` error
    diamond_bin_path: str = str(Path(__file__).parent / "bin" / "diamond")
    cmd: str = (
        f"{diamond_bin_path} blastp "
        f" --db {database_file_path}"
        f" --query {missing_file_name}"
        f" --out {output_file_path}"
        " --outfmt 6"
        " --header"
        " --strand both"
        f" --threads {threads}"
    )

    # TODO: check existance of the file
    run_shell_command(cmd=cmd, skip_error=False)
    return output_file_path


def crosslink_with_diamond(
    ext2meta: pd.DataFrame,
    tmp_dir: str,
    reference_file: str,
    missing_file_name: str,
    minimum_coverage: float,
    minimum_identity: float,
    threads: int,
) -> pd.DataFrame:
    output_file_path: str = run_diamond(
        tmp_dir, missing_file_name, threads, reference_file
    )
    missing_sequence_lengths: dict[str, int] = get_sequence_lengths(
        missing_file_name
    )
    reference_sequence_lengths: dict[str, int] = get_sequence_lengths(
        reference_file
    )

    # parse diamond output
    unique_query_ids: set[str] = set()
    unique_subject_ids: set[str] = set()
    query_to_target: dict[str, list[tuple[str, float]]] = {}
    output_file_lines: Generator[DiamondOutputFileLine, None, None] = (
        DiamondOutputFileLine.yield_output_file_lines(output_file_path)
    )
    for line in output_file_lines:
        length: int = reference_sequence_lengths[line.subject_id]
        if (
            missing_sequence_lengths[line.query_id]
            > reference_sequence_lengths[line.subject_id]
        ):
            length = missing_sequence_lengths[line.query_id]

        unique_query_ids.add(line.query_id)
        unique_subject_ids.add(line.subject_id)

        # Filter low line.identity or low coverage hits.
        coverage: float = 1.0 * (line.alignment_length) / length
        if coverage < minimum_coverage or line.identity < minimum_identity:
            continue

        query_protein_id: str = line.query_id
        target_protein_id: str = line.subject_id
        identity_score: float = (line.identity + coverage) / 2
        target_data: tuple[str, float] = (target_protein_id, identity_score)

        # Save only the best match for each query,
        # more than one allowed, if best matches
        # has got the same score
        query_to_target.setdefault(query_protein_id, [target_data])
        if query_protein_id in query_to_target:
            # append matches if same identity_score
            if identity_score == query_to_target[query_protein_id][0][1]:
                query_to_target[query_protein_id].append(target_data)

            # update matches if better identity_score
            elif identity_score > query_to_target[query_protein_id][0][1]:
                query_to_target[query_protein_id] = [target_data]

    q2t_frame = pd.DataFrame.from_dict(query_to_target, orient="index")

    # If there's more than one record for query_to_target, re-arrange and re-format the dataframe.
    if len(q2t_frame.index) > 1:
        q2t_series: pd.Series[Any] | pd.DataFrame = (
            q2t_frame.unstack().dropna().sort_index(level=1)
        )
        q2t_series = q2t_series.reset_index(level=[0])
        q2t_series["extId"] = q2t_series.index
        q2t_series.reset_index(inplace=True)
        q2t_series.rename(
            columns={
                "index": "index",
                "level_0": "oldIndex",
                0: "value",
                "extId": "extId",
            },
            inplace=True,
        )

        ext_to_meta_diamond = pd.DataFrame(
            q2t_series["value"].tolist(), columns=["protid", "score"]
        )
        ext_to_meta_diamond["extid"] = q2t_series["extId"]
        ext_to_meta_diamond = ext_to_meta_diamond[["extid", "score", "protid"]]

        # Select only rows that were not matched by md5.
        ext_to_meta_diamond = ext_to_meta_diamond[
            ~ext_to_meta_diamond["extid"].isin(ext2meta["extid"])
        ]
        ext2meta = pd.concat([ext2meta, ext_to_meta_diamond], axis=0)
    return ext2meta


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
    file_name: str,
    reference_hashes: dict[str, list[str]],
    minimum_coverage: float,
    minimum_identity: float,
) -> str:
    db_id: int = 0
    version: int = 0
    file_prefix: str = os.path.basename(file_name)
    if VERBOSE:
        print(f"[INFO] Starting to process: {file_prefix}")

    target_hashes: dict[str, list[str]] = get_md5_fasta(file_name)
    target_to_reference_tmp: pd.DataFrame = crosslink_with_md5(
        target_hashes, reference_hashes
    )

    number_all: list[str] = run_command_with_return(
        f"zcat {file_name} | grep -c '>'"
    )
    number_all_result: int = int(number_all[0])
    number_matched: int = target_to_reference_tmp["extid"].nunique()

    leftover: int = number_all_result - number_matched
    if VERBOSE:
        print(f"Number of leftovers from crosslinking by MD5 hashes: {leftover}")

    has_not_matched_by_md5: bool = leftover > 0
    if has_not_matched_by_md5:
        target_to_reference_tmp = crosslink_with_diamond(
            target_to_reference_tmp,
            tmp_dir,
            reference_file,
            file_name,
            minimum_coverage,
            minimum_identity,
            1,
        )

    target_to_reference_tmp["dbid"] = db_id
    target_to_reference_tmp["version"] = version
    target_to_reference_tmp = target_to_reference_tmp[
        ["dbid", "extid", "version", "protid", "score"]
    ]

    ext2metaFn = os.path.join(
        output_directory, f"{file_prefix}.ext2meta.tbl.gz"
    )
    target_to_reference_tmp.drop_duplicates(inplace=True)
    target_to_reference_tmp.to_csv(
        ext2metaFn, compression="gzip", index=False, sep="\t", header=False
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


def run_shell_command(cmd: str, skip_error: bool) -> None:
    try:
        process = sp.Popen(cmd, shell=True)
    except Exception as e:
        if skip_error:
            print("Error ocurred, but you chose to ommit it", file=sys.stderr)
            pass
        else:
            raise e

    process.communicate(b"Y\n")
    process.wait()
    if process.returncode != 0 and skip_error:
        print("Error ocurred, but you chose to ommit it", file=sys.stderr)
        return
    elif process.returncode != 0 and not skip_error:
        raise ChildProcessError(f"ERROR: Execution of cmd [{cmd}] failed.")


def run_command_with_return(cmd: str) -> list[str]:
    process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    process_stdout: IO[bytes] | None = process.stdout
    if not process_stdout:
        return []
    process_stdout_lines: list[bytes] = process_stdout.readlines()
    decoded_process_stdout_lines: list[str] = [
        line.decode("utf-8").strip() for line in process_stdout_lines
    ]
    return decoded_process_stdout_lines


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
) -> None:
    reference_hashes: dict[str, list[str]] = get_md5_fasta(reference_file)

    make_diamond_db(tmp_directory, reference_file)

    process_taxid_args: list[
        tuple[str, str, str, str, dict[str, list[str]], float, float]
    ] = [
        (
            output_directory,
            tmp_directory,
            reference_file,
            fasta_file,
            reference_hashes,
            minimum_coverage,
            minimum_identity,
        )
        for fasta_file in target_fasta_files
    ]

    # with ProcessPoolExecutor(max_workers=max_threads) as executor:
    #     parallelized_processes_results: list[str] = list(
    #         executor.map(process_taxid, process_taxid_args)
    #     )
    #     print_parallelized_processes_logs(
    #         target_fasta_files, parallelized_processes_results
    #     )

    # FOR DEBUGGING:
    for args in process_taxid_args:
        process_taxid(args[0], args[1], args[2], args[3], args[4], args[5], args[6])



def main() -> int:
    t0 = datetime.now()
    args = CrosslinkBySequenceArgs.get_arguments()
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
    )

    dt = datetime.now() - t0
    print(f"#Time elapsed: {dt}")
    return 0
