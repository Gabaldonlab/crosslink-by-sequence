from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import Generator

import pandas as pd

from crosslink_by_sequence.utils import get_sequence_lengths
from crosslink_by_sequence.utils import run_shell_command


__DIAMOND_BIN_PATH: str = str(Path(__file__).parent / "bin" / "diamond")


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
        cls, output_file_path: str, verbose: bool
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
                    if verbose:
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


def make_diamond_db(tmp_dir: str, reference_file: str, verbose: bool) -> str:
    """
    Creates the cache database file from the reference proteome.

    @return str - Path of the cache database file.
    """
    quiet_placeholder: str = "--quiet"
    if verbose:
        quiet_placeholder = ""

    db_file_name: str = f"{os.path.basename(reference_file)}.db"
    output_file_path: str = os.path.join(tmp_dir, db_file_name)

    cmd: str = (
        f"{__DIAMOND_BIN_PATH} makedb -d {output_file_path} --in {reference_file} {quiet_placeholder}"
    )
    if not os.path.isfile(output_file_path + ".dmnd"):
        run_shell_command(cmd, False)
        print(cmd)
    else:
        print("[INFO] Diamond DB file exists for ", reference_file)
    return db_file_name

def _restack_target_to_reference_frame(
    target_to_reference: pd.DataFrame, query_to_target_frame: pd.DataFrame
) -> pd.DataFrame:
    query_to_target_series: pd.Series[Any] | pd.DataFrame = (
        query_to_target_frame.unstack().dropna().sort_index(level=1)
    )
    query_to_target_series = query_to_target_series.reset_index(level=[0])
    query_to_target_series["extId"] = query_to_target_series.index
    query_to_target_series.reset_index(inplace=True)
    query_to_target_series.rename(
        columns={
            "index": "index",
            "level_0": "oldIndex",
            0: "value",
            "extId": "extId",
        },
        inplace=True,
    )

    target_to_reference_diamond_frame = pd.DataFrame(
        query_to_target_series["value"].tolist(), columns=["protid", "score"]
    )
    target_to_reference_diamond_frame["extid"] = query_to_target_series[
        "extId"
    ]
    target_to_reference_diamond_frame = target_to_reference_diamond_frame[
        ["extid", "score", "protid"]
    ]

    # Select only rows that were not matched by md5.
    target_to_reference_diamond_frame = target_to_reference_diamond_frame[
        ~target_to_reference_diamond_frame["extid"].isin(
            target_to_reference["extid"]
        )
    ]
    target_to_reference = pd.concat(
        [target_to_reference, target_to_reference_diamond_frame], axis=0
    )
    return target_to_reference


def _run_diamond(
    tmp_dir: str, missing_file_name: str, threads: int, reference_file: str, db_file_path: str, verbose: bool
) -> str:
    """
    Runs Diamond as a process with the given parameters.

    @return str - The output file's path.
    """
    quiet_placeholder: str = "--quiet"
    if verbose:
        quiet_placeholder = ""

    database_file_path: str = os.path.join(tmp_dir, db_file_path)
    output_file_name: str = f"{os.path.basename(missing_file_name)}diamond"
    output_file_path: str = os.path.join(tmp_dir, output_file_name)

    # -tileSize=8 caused `Internal error genoFind.c 2225` error
    cmd: str = (
        f"{__DIAMOND_BIN_PATH} blastp "
        f" {quiet_placeholder}"
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
    target_to_reference: pd.DataFrame,
    tmp_dir: str,
    reference_file: str,
    missing_file_name: str,
    db_file_path: str,
    minimum_coverage: float,
    minimum_identity: float,
    threads: int,
    verbose: bool,
) -> pd.DataFrame:
    output_file_path: str = _run_diamond(
        tmp_dir, missing_file_name, threads, reference_file, db_file_path, verbose
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
        DiamondOutputFileLine.yield_output_file_lines(
            output_file_path, verbose
        )
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

    query_to_target_frame = pd.DataFrame.from_dict(
        query_to_target, orient="index"
    )

    # If there's more than one record for query_to_target, re-arrange and re-format the dataframe.
    if len(query_to_target_frame.index) > 1:
        target_to_reference = _restack_target_to_reference_frame(
            target_to_reference, query_to_target_frame
        )
    return target_to_reference
