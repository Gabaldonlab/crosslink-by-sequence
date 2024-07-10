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
import subprocess as sp
import sys
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path
from typing import IO

import pandas as pd
from Bio import SeqIO

from crosslink_by_sequence.argument_parser import CrosslinkBySequenceArgs


VERBOSE: bool = False


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
    # print(dict(list(hash2prot.items())[0:2]))
    return hash_to_protein


def make_diamond_db(tmp_dir: str, reference_file: str) -> None:
    db_file_name: str = f"{os.path.basename(reference_file)}db"
    output_file_path: str = os.path.join(tmp_dir, db_file_name)
    cmd: str = f"diamond makedb -d {output_file_path} --in {reference_file} "
    if not os.path.isfile(output_file_path + ".dmnd"):
        run_command(cmd, False)
        print(cmd)
    else:
        print("[INFO] Diamond DB file exists for ", reference_file)


def get_sequence_length(fasta_file_path: str) -> dict[str, int]:
    """
    @return dict[str, int]
        { external_id: <sequence_length: int> }
    """
    sequence_lengths: dict[str, int] = {}
    with gzip.open(fasta_file_path, mode="rt") as handle:
        for fasta_record in SeqIO.parse(handle, "fasta"):
            external_id: str = str(fasta_record.id)
            sequence_length: int = len(str(fasta_record.seq))
            sequence_lengths[external_id] = sequence_length
    return sequence_lengths


def crosslink_diamond(
    ext2meta: pd.DataFrame,
    tmp_dir: str,
    reference_file: str,
    missing_file_name: str,
    minimum_coverage: float,
    minimum_identity: float,
    threads: int,
):
    db_file_name: str = f"{os.path.basename(reference_file)}db"

    # run diamond
    database_file_path: str = os.path.join(tmp_dir, db_file_name)
    output_file_name: str = f"{os.path.basename(missing_file_name)}diamond"
    output_file_path: str = os.path.join(tmp_dir, output_file_name)

    # -tileSize=8 caused `Internal error genoFind.c 2225` error
    cmd: str = (
        f"diamond blastp "
        f" --db {database_file_path}"
        f" --query {missing_file_name}"
        f" --out {output_file_path}"
        " --outfmt 6"
        " --header"
        " --strand both"
        f" --threads {threads}"
    )

    # TODO - check existance of the file
    run_command(cmd, False)
    lenSequence: dict[str, int] = get_sequence_length(missing_file_name)
    lenSequenceSubj = get_sequence_length(reference_file)

    # parse diamond output
    q2t = {}
    Qset = set()
    Tset = set()
    for line in open(output_file_path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        # Query ID, Subject ID, Percentage of identical matches, Alignment length, Number of mismatches, Number of gap openings, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Expected value, Bit score
        lData: list[str] = line.split("\t")

        # skip lines with errors
        if len(lData) != 12:
            if VERBOSE:
                print(
                    f"  Error: blat2xlinks: Expected 12 elements in line, got {len(lData)} {str(lData)}!\n",
                    file=sys.stderr,
                )
            continue

        Q: str = lData[0]
        T: str = lData[1]
        identity: str = lData[2]
        aLength: str = lData[3]
        nMis: str = lData[4]
        nGaps: str = lData[5]
        sQ: str = lData[6]
        eQ: str = lData[7]
        sS: str = lData[8]
        eS: str = lData[9]
        expV: str = lData[10]
        score: str = lData[11]

        aLength: int = int(aLength)
        identity: float = float(identity) / 100

        length: int = lenSequenceSubj[T]
        if lenSequence[Q] > lenSequenceSubj[T]:
            length = lenSequence[Q]

        cov: float = 1.0 * (aLength) / length

        # store Q and T
        Qset.add(Q)
        Tset.add(T)

        # Filter low identity or low coverage hits.
        if cov < minimum_coverage or identity < minimum_identity:
            continue

        Qprotid: str = Q
        Tprotid: str = T
        iScore: float = (identity + cov) / 2
        Tdata: tuple[str, float] = (Tprotid, iScore)

        # Save only the best match for each query,
        # more than one allowed, if best matches
        # has got the same score
        q2t.setdefault(Qprotid, [Tdata])
        if Qprotid in q2t:
            # append matches if same iScore
            if iScore == q2t[Qprotid][0][1]:
                q2t[Qprotid].append(Tdata)

            # update matches if better iScore
            elif iScore > q2t[Qprotid][0][1]:
                q2t[Qprotid] = [Tdata]

    q2t_frame = pd.DataFrame.from_dict(q2t, orient="index")
    if len(q2t_frame.index) > 1:
        q2t_frame = q2t_frame.unstack().dropna().sort_index(level=1)
        q2t_frame = q2t_frame.reset_index(level=[0])
        q2t_frame["extId"] = q2t_frame.index
        q2t_frame = q2t_frame.reset_index()

        q2t_frame.columns = ["index", "oldIndex", "value", "extId"]

        ext2metaBlat = pd.DataFrame(q2t_frame["value"].tolist())
        ext2metaBlat.columns = ["protid", "score"]
        ext2metaBlat["extid"] = q2t_frame["extId"]
        ext2metaBlat = ext2metaBlat[["extid", "score", "protid"]]

        # Select only rows that were not matched by md5.
        ext2metaBlat = ext2metaBlat[
            ~ext2metaBlat["extid"].isin(ext2meta["extid"])
        ]
        ext2meta = pd.concat([ext2meta, ext2metaBlat], axis=0)
    return ext2meta


def crosslink_md5(
    hash_to_meta: dict[str, list[str]], hash_to_external: dict[str, list[str]]
) -> pd.DataFrame:
    """Link proteins using md5.
    Return ext2meta for given dbID
    """
    ext2meta = pd.DataFrame(columns=["extid", "score", "protid"])
    for md5 in hash_to_external:
        # if same seq already present
        if md5 in hash_to_meta:
            # metaID is a list, can be a few elements
            metaID = hash_to_meta[md5]
            for extID in hash_to_external[md5]:
                # add ext2meta entry
                for metIDEl in metaID:
                    ext2meta = ext2meta.append(
                        pd.DataFrame(
                            [[extID, "1", metIDEl]], columns=ext2meta.columns
                        )
                    )
    return ext2meta


# output_directory, tmp_directory, refFile, fnames, taxid, minimum_coverage, minimum_identity, protein_limit
def process_taxid(
    output_directory: str,
    tmp_dir: str,
    reference_file: str,
    file_name: str,
    hash_to_meta: dict[str, list[str]],
    minimum_coverage: float,
    minimum_identity: float,
):

    name: str = os.path.basename(file_name)
    print("[INFO] Xlink ", name)
    dbid: int = 0
    version: int = 0
    filePrefix: str = name

    hash_to_ext: dict[str, list[str]] = get_md5_fasta(file_name)
    ext_to_meta_tmp: pd.DataFrame = crosslink_md5(hash_to_meta, hash_to_ext)

    number_all: list[bytes] = run_command_with_return(
        f"zcat {file_name} | grep -c '>'"
    )
    number_all_result: int = int(number_all[0].decode("utf-8").strip())
    number_matched: int = ext_to_meta_tmp["extid"].nunique()

    leftover: int = number_all_result - number_matched
    print("Rest ", leftover)

    if leftover > 0:
        ext_to_meta_tmp = crosslink_diamond(
            ext_to_meta_tmp,
            tmp_dir,
            reference_file,
            file_name,
            minimum_coverage,
            minimum_identity,
            1,
        )

    ext_to_meta_tmp["dbid"] = dbid
    ext_to_meta_tmp["version"] = version

    ext_to_meta_tmp = ext_to_meta_tmp[
        ["dbid", "extid", "version", "protid", "score"]
    ]

    # Format ext2meta to the format compartible with the current ext2meta schema in the DB
    ext2metaFn = os.path.join(
        output_directory, "%s.ext2meta.tbl.gz" % filePrefix
    )
    ext_to_meta_tmp.to_csv(
        ext2metaFn, compression="gzip", index=False, sep="\t", header=False
    )

    # Calculate number of matched and non-matched proteins.
    # Print some stats with the percentage of orphans per file.
    number_matched = ext_to_meta_tmp["extid"].nunique()
    numberOrphans = number_all - number_matched
    percOrphans = numberOrphans * 100 / number_all

    log_file_path: str = os.path.join(output_directory, "{filePrefix}.log")
    with open(log_file_path, "wt", encoding="UTF-8") as opened_log_file:
        opened_log_file.write(
            f"{filePrefix}\t{number_all}\t{number_matched}\t{numberOrphans}\t{percOrphans:.2f}\n"
        )

    return filePrefix


def process(
    target_fasta_files: list[str],
    reference_file: str,
    output_directory: str,
    tmp_directory: str,
    max_threads: int,
    minimum_coverage: float,
    minimum_identity: float,
):
    """xlink proteomes"""
    hash2meta = get_md5_fasta(reference_file)

    make_diamond_db(tmp_directory, reference_file)

    multiprocessing_pool = Pool(max_threads)
    tdata = []
    for fastaFile in target_fasta_files:
        tdata.append(
            [
                output_directory,
                tmp_directory,
                reference_file,
                fastaFile,
                hash2meta,
                minimum_coverage,
                minimum_identity,
            ]
        )

    parallelized_processes: list[str] = multiprocessing_pool.starmap(
        process_taxid, tdata
    )

    for idx, data in enumerate(parallelized_processes, 1):
        if not data:
            continue
        fileNamePrefix: str = data
        print(f" {idx} / {len(target_fasta_files)}  {fileNamePrefix}    \r")
        print(f"[INFO] done {fileNamePrefix} ")


def run_command(cmd: str, skip_error: bool) -> None:
    if skip_error:
        try:
            process = sp.Popen(cmd, shell=True)
        except:
            pass
        process.communicate("Y\n")
        if process.wait() != 0:
            print("Error ocurred, but you chose to ommit it")
            return
    else:
        try:
            process = sp.Popen(cmd, shell=True)
        except OSError as e:
            print("Error: Execution cmd failed", file=sys.stderr)
            raise e

        process.communicate("Y\n")
        if process.wait() != 0:
            raise ChildProcessError(f"ERROR: Execution of cmd [{cmd}] failed.")


def run_command_with_return(cmd: str) -> list[bytes]:
    process_stdout: IO[bytes] | None = sp.Popen(
        cmd, shell=True, stdout=sp.PIPE
    ).stdout
    return process_stdout.readlines()


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
    process(
        args.target_fasta_gzip_files,
        args.target_reference_species_fasta_gzip_file,
        args.output_directory,
        args.tmp_directory,
        args.max_threads,
        args.minimum_coverage,
        args.minimum_identity,
    )

    dt = datetime.now() - t0
    sys.stderr.write()

    print(f"#Time elapsed: {dt}")
    return 0
