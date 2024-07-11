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
import sys
import subprocess
import time
from typing import IO
from Bio import SeqIO
from typing import Callable, Any


def time_it(func: Any) -> Callable:
    def wrapper_function(*args: Any, **kwargs: Any) -> None:
        t0: float = time.monotonic()
        func(*args, **kwargs)
        print(f"# Time elapsed: {time.monotonic() - t0}")

    return wrapper_function


def run_shell_command(cmd: str, skip_error: bool) -> None:
    try:
        process = subprocess.Popen(cmd, shell=True)
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
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process_stdout: IO[bytes] | None = process.stdout
    if not process_stdout:
        return []
    process_stdout_lines: list[bytes] = process_stdout.readlines()
    decoded_process_stdout_lines: list[str] = [
        line.decode("utf-8").strip() for line in process_stdout_lines
    ]
    return decoded_process_stdout_lines


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
