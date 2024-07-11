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

import argparse
import sys

from dataclasses import dataclass
from pathlib import Path
import argparse, os, sys
import gzip, hashlib
import subprocess as sp
from typing import IO, Iterable

from Bio import SeqIO
from datetime import datetime

from multiprocessing import Pool


import pandas as pd
import numpy as np


@dataclass
class CrosslinkBySequenceArgs:
    verbose: bool
    target_fasta_gzip_files: list[str]
    target_reference_species_fasta_gzip_file: str
    output_directory: str
    tmp_directory: str
    minimum_coverage: float
    minimum_identity: float
    max_threads: int

    @classmethod
    def get_arguments(cls, args: list[str] =sys.argv[1:]) -> CrosslinkBySequenceArgs:
        desc: str = """
            The tool takes as input a reference proteome file with protein IDs in the fasta headers and proteome for cross linking.
            The output of the tool is a table with the following header fields:
                - "dbid"
                - "extid"
                - "version"
                - "protid"
                - "score"
            """
        epilog: str = "Author: majerdaniel93@gmail.com"
        usage: str = "%(prog)s [options] -f */fasta/*.faa.gz"
        parser = argparse.ArgumentParser(
            usage=usage, description=desc, epilog=epilog
        )

        parser.add_argument("--verbose", default=False, action="store_true")
        parser.add_argument("--version", action="version", version="1.0.0")
        parser.add_argument(
            "--target_fasta_gzip_files",  # old: fastas
            nargs="+",
            required=True,
            default=[],
            type=str,
            help="Path to the target proteome fasta files gzipped.",
        )
        parser.add_argument(
            "--target_reference_species_fasta_gzip_file",  # old: target_reference_species_fasta_gzip_file
            required=True,
            help="Reference species fasta file gzipped.",
        )
        parser.add_argument(
            "--output_directory",  # old: output_directory
            default="fasta.crosslinked",
            type=str,
            help="output directory  [%(default)s]",
        )
        parser.add_argument(
            "--tmp_directory",  # old: tmp_directory
            default="",
            type=str,
            help="TMP directory to store blat psl files [%(default)s]",
        )
        parser.add_argument(
            "--minimum_coverage",  # old: minimum_coverage
            default=0.95,
            type=float,
            help="min. coverage for blat hits  [%(default)s]",
        )
        parser.add_argument(
            "--minimum_identity",  # old: minimum_identity
            default=0.98,
            type=float,
            help="Min. identity for blat hits  [%(default)s]",
        )
        parser.add_argument(
            "--max_threads",  # old: threads
            default=4,
            type=int,
            help="number of cores [%(default)s]",
        )

        return CrosslinkBySequenceArgs(**vars(parser.parse_args(args)))
