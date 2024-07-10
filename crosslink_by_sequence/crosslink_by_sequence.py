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
from typing import Optional

from dataclasses import dataclass

import argparse, os, sys
import gzip, hashlib
import subprocess as sp
from typing import IO, Iterable

from Bio import SeqIO
from datetime import datetime

from multiprocessing import Pool


import pandas as pd
import numpy as np


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
    ext2meta,
    tmp_dir,
    reference_file,
    missing_file_name: str,
    covTh,
    idTh,
    threads,
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
    Qset, Tset = set(), set()
    for line in open(output_file_path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # parse with error capturing
        # Query ID, Subject ID, Percentage of identical matches, Alignment length, Number of mismatches, Number of gap openings, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Expected value, Bit score
        lData = line.split("\t")
        if len(lData) == 12:
            (
                Q,
                T,
                identity,
                aLength,
                nMis,
                nGaps,
                sQ,
                eQ,
                sS,
                eS,
                expV,
                score,
            ) = lData
            aLength = int(aLength)
            identity = float(identity) / 100
            if lenSequence[Q] > lenSequenceSubj[T]:
                length = lenSequence[Q]
            else:
                length = lenSequenceSubj[T]
            cov = 1.0 * (aLength) / length
            # store Q and T
            Qset.add(Q)
            Tset.add(T)
        # skip lines with errors
        else:
            if verbose:
                sys.stderr.write(
                    "  Error: blat2xlinks: Expected 21 elements in line, got %s %s!\n"
                    % (len(lData), str(lData))
                )
            continue

        # filter low identity or low coverage hits
        if cov < covTh or identity < idTh:
            continue

        Qprotid = Q
        Tprotid = T
        iScore = (identity + cov) / 2
        Tdata = (Tprotid, iScore)

        # save only the best match for each query, more than one allowed, if best matches has got the same score
        if Qprotid in q2t:
            # append matches if same iScore
            if iScore == q2t[Qprotid][0][1]:
                # add entry in t2q
                #  if Tprotid not in t2q:
                #      t2q[Tprotid]=[]
                #  t2q[Tprotid].append((Qprotid, iScore))
                # add entry to q2t
                q2t[Qprotid].append(Tdata)
            # update matches if better iScore
            elif iScore > q2t[Qprotid][0][1]:
                q2t[Qprotid] = [Tdata]
        else:
            q2t[Qprotid] = [Tdata]
            # add entry in t2q
            # if Tprotid not in t2q:
            #    t2q[Tprotid]=[]
            # t2q[Tprotid].append((Qprotid, iScore))
    # print("T2Q:", len(t2q))
    # print(dict(list(t2q.items())[0:2]))
    # print("Q2T:")
    # print(dict(list(q2t.items())[0:2]))

    df = pd.DataFrame.from_dict(q2t, orient="index")
    # print(df)
    if len(df.index) > 1:
        df = df.unstack().dropna().sort_index(level=1)
        df = df.reset_index(level=[0])
        df["extId"] = df.index
        df = df.reset_index()

        df.columns = ["index", "oldIndex", "value", "extId"]

        ext2metaBlat = pd.DataFrame(df["value"].tolist())
        ext2metaBlat.columns = ["protid", "score"]
        ext2metaBlat["extid"] = df["extId"]
        # extid score  protid
        ext2metaBlat = ext2metaBlat[["extid", "score", "protid"]]
        # print(ext2metaBlat.head())
        # print(ext2meta.head())
        ##select only rows that were not matched by md5
        ext2metaBlat = ext2metaBlat[
            ~ext2metaBlat["extid"].isin(ext2meta["extid"])
        ]
        ext2meta = pd.concat([ext2meta, ext2metaBlat], axis=0)
    return ext2meta


def xLink_blat(ext2meta, tmpDir, refFile, missingFn, covTh, idTh):

    # run blat
    outputFileName = tmpDir + os.path.basename(missingFn)
    psl = "%s.psl" % outputFileName
    # -tileSize=8 caused `Internal error genoFind.c 2225` error
    cmd = (
        "blat -prot -out=psl -minIdentity=90 -minMatch=1 -out=psl -noHead %s %s %s"
        % (refFile, missingFn, psl)
    )
    # print(cmd)
    # TODO - check existance of the file
    run_command(cmd, False)
    # parse blat output
    q2t = {}
    Qset, Tset = set(), set()
    for line in open(psl):
        line = line.strip()
        if not line:
            continue
        # parse with error capturing
        lData = line.split("\t")
        if len(lData) == 21:
            (
                match,
                mMatch,
                rMatch,
                Ns,
                Qgapcount,
                Qgapbases,
                Tgapcount,
                Tgapbases,
                strand,
                Q,
                Qsize,
                Qstart,
                Qend,
                T,
                Tsize,
                Tstart,
                Tend,
                bcount,
                bsizes,
                qStarts,
                tStarts,
            ) = lData
            match, mMatch, Qsize, Tsize = (
                int(match),
                int(mMatch),
                int(Qsize),
                int(Tsize),
            )
            # take longer query or target size as general size of compared proteins
            if Qsize > Tsize:
                length = Qsize
            else:
                length = Tsize
            # length=Q_size
            cov = 1.0 * (match + mMatch) / length
            identity = 1.0 * match / length
            # store Q and T
            Qset.add(Q)
            Tset.add(T)
        # skip lines with errors
        else:
            if verbose:
                sys.stderr.write(
                    "  Error: blat2xlinks: Expected 21 elements in line, got %s %s!\n"
                    % (len(lData), str(lData))
                )
            continue

        # filter low identity or low coverage hits
        if cov < covTh or identity < idTh:
            continue

        Qprotid = Q
        Tprotid = T
        iScore = (identity + cov) / 2
        Tdata = (Tprotid, iScore)

        # save only the best match for each query, more than one allowed, if best matches has got the same score
        if Qprotid in q2t:
            # append matches if same iScore
            if iScore == q2t[Qprotid][0][1]:
                # add entry in t2q
                #  if Tprotid not in t2q:
                #      t2q[Tprotid]=[]
                #  t2q[Tprotid].append((Qprotid, iScore))
                # add entry to q2t
                q2t[Qprotid].append(Tdata)
            # update matches if better iScore
            elif iScore > q2t[Qprotid][0][1]:
                q2t[Qprotid] = [Tdata]
        else:
            q2t[Qprotid] = [Tdata]
            # add entry in t2q
            # if Tprotid not in t2q:
            #    t2q[Tprotid]=[]
            # t2q[Tprotid].append((Qprotid, iScore))
    # print("T2Q:", len(t2q))
    # print(dict(list(t2q.items())[0:2]))
    # print("Q2T:")
    # print(dict(list(q2t.items())[0:2]))

    df = pd.DataFrame.from_dict(q2t, orient="index")
    df = df.unstack().dropna().sort_index(level=1)
    df = df.reset_index(level=[0])
    df["extId"] = df.index
    df = df.reset_index()

    df.columns = ["index", "oldIndex", "value", "extId"]

    ext2metaBlat = pd.DataFrame(df["value"].tolist())
    ext2metaBlat.columns = ["protid", "score"]
    ext2metaBlat["extid"] = df["extId"]
    # select only rows that were not matched by md5
    ext2metaBlat = ext2metaBlat[~ext2metaBlat["extid"].isin(ext2meta["extid"])]
    ext2meta = pd.concat([ext2meta, ext2metaBlat], axis=0)
    # ext2meta2 = pd.DataFrame.from_dict(q2t)
    return ext2meta


def xLink_md5(hash2meta, hash2ext):
    """Link proteins using md5.
    Return ext2meta for given dbID
    """
    ext2meta = pd.DataFrame(columns=["extid", "score", "protid"])
    for md5 in hash2ext:
        # if same seq already present
        if md5 in hash2meta:
            # metaID is a list, can be a few elements
            metaID = hash2meta[md5]
            for extID in hash2ext[md5]:
                # add ext2meta entry
                for metIDEl in metaID:
                    ext2meta = ext2meta.append(
                        pd.DataFrame(
                            [[extID, "1", metIDEl]], columns=ext2meta.columns
                        )
                    )

    # print(dict(list(ext2meta.items())[0:2]))
    # print(ext2meta.head())
    return ext2meta


# outDir, tmpDir, refFile, fnames, taxid, minCoverage, minIdentity, protLimit
def process_taxid(
    outDir,
    tmp_dir,
    reference_file,
    file_name,
    hash2meta,
    taxid,
    minimum_coverage,
    minimum_identity,
    protLimit,
):

    name: str = os.path.basename(file_name)
    print("[INFO] Xlink ", name)
    dbid: int = 0
    version: int = 0
    filePrefix: str = name

    hash2ext = get_md5_fasta(file_name)

    ext_to_meta_tmp = xLink_md5(hash2meta, hash2ext)

    number_all: list[bytes] = run_command_with_return(
        f"zcat {file_name} | grep -c '>'"
    )
    number_all_result: int = int(number_all[0].decode("utf-8").strip())
    number_matched: int = ext_to_meta_tmp["extid"].nunique()

    leftover: int = number_all_result - number_matched
    print("Rest ", leftover)

    if leftover > 0:
        # run blat comparison between two files in the temporary folder
        # this function returns ext2meta dictionary with keys - protein names from fileName and values - protein names from refFile
        # ext2metaTmp = xLink_blat(ext2metaTmp, tmpDir, refFile, fileName, minCoverage, minIdentity)
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

    # format ext2meta to the format compartible with the current ext2meta schema in the DB
    ext2metaFn = os.path.join(outDir, "%s.ext2meta.tbl.gz" % filePrefix)
    ext_to_meta_tmp.to_csv(
        ext2metaFn, compression="gzip", index=False, sep="\t", header=False
    )
    # format fasta file with orphans
    # outOrphans = outDir+filePrefix+".orphans.faa.gz"
    # print("[INFO] writing orphans for ",name )
    # saveGzipFasta(fileName,ext2metaTmp["extid"].unique(), outOrphans, 'invert' )

    # at the same time calculate number of matched and non-matched proteins
    # print some stats with % of orphans per file
    number_matched = ext_to_meta_tmp["extid"].nunique()
    numberOrphans = number_all - number_matched
    percOrphans = numberOrphans * 100 / number_all

    logFn = os.path.join(outDir, "%s.log" % filePrefix)
    log = open(logFn, "wt")
    print(
        "%s\t%s\t%s\t%s\t%.2f"
        % (filePrefix, number_all, number_matched, numberOrphans, percOrphans),
        file=log,
    )
    log.close()

    return filePrefix


def get_taxids(fileName):
    """Return taxids and species2strains"""
    taxid = fileName.split("/")[-1].split(".")[1]
    return taxid


def fasta_generator(
    output_directory: str,
    tmp_dir: str,
    reference_file_path: str,
    file_paths: Iterable[str],
    taxid,
    minCoverage,
    minIdentity,
    protLimit,
):
    for fastaFile in file_paths:
        yield (
            output_directory,
            tmp_dir,
            reference_file_path,
            fastaFile,
            taxid,
            minCoverage,
            minIdentity,
            protLimit,
        )


def saveGzipFasta(file_path: str, idList, outName, param):
    handle_in = gzip.open(file_path, "rt")
    handle_out = gzip.open(outName, "wt")
    for record in SeqIO.parse(handle_in, "fasta"):
        if np.isin(record.id, idList):
            if param == "invert":
                continue
            else:
                handle_out.write(record.format("fasta"))
        else:
            if param == "invert":
                handle_out.write(record.format("fasta"))
            else:
                continue


def process(
    fnames,
    refFile,
    outDir,
    tmpDir,
    nprocs,
    minCoverage,
    minIdentity,
    protLimit,
):
    """xlink proteomes"""
    ###get taxids and species2strains
    # taxid = get_taxids(refFile)
    taxid = 99999999
    hash2meta = get_md5_fasta(refFile)
    ###process files
    # print("[INFO] Processing taxa ", taxid)
    make_diamond_db(tmpDir, refFile)

    multiprocessing_pool = Pool(nprocs)
    # tdata = fasta_generator(outDir, tmpDir, refFile, fnames, taxid, minCoverage, minIdentity, protLimit)
    tdata = []
    for fastaFile in fnames:
        tdata.append(
            [
                outDir,
                tmpDir,
                refFile,
                fastaFile,
                hash2meta,
                taxid,
                minCoverage,
                minIdentity,
                protLimit,
            ]
        )

    parallelized_processes: list[str] = multiprocessing_pool.starmap(
        process_taxid, tdata
    )

    for idx, data in enumerate(parallelized_processes, 1):
        if not data:
            continue
        fileNamePrefix = data
        sys.stderr.write(
            " %s / %s  %s    \r" % (idx, len(fnames), fileNamePrefix)
        )
        print("[INFO] done %s " % (fileNamePrefix))


def create_folder(name):
    if not os.path.exists(name):
        cmd = "mkdir " + name
        try:
            run_command(cmd, False)
        except:
            print("Unable to create directory ", name)


def run_command(cmd: str, skip_error: bool) -> None:
    if skip_error:
        try:
            process = sp.Popen(cmd, shell=True)
        except:
            pass
        process.communicate("Y\n")
        if process.wait() != 0:
            print("Error ocurred, but you chose to ommit it")
    else:
        try:
            process = sp.Popen(cmd, shell=True)
        except OSError as e:
            sys.exit("Error: Execution cmd failed")
        process.communicate("Y\n")
        if process.wait() != 0:
            sys.exit("ERROR: Execution cmd failed")


def run_command_with_return(cmd: str) -> list[bytes]:
    process_stdout: IO[bytes] | None = sp.Popen(
        cmd, shell=True, stdout=sp.PIPE
    ).stdout
    return process_stdout.readlines()


@dataclass
class CrosslinkBySequenceArgs:
    verbose: bool
    target_fasta_gzip_file: list[str]
    target_reference_species_fasta_gzip_file: str
    output_directory: str
    tmp_directory: str
    minimum_coverage: float
    minimum_identity: float
    protein_limit: int
    max_threads: int

    @classmethod
    def get_arguments(cls, args=sys.argv[1:]) -> CrosslinkBySequenceArgs:
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
            "--target_reference_species_fasta_gzip_file",  # old: refSpecie
            required=True,
            help="Reference species fasta file gzipped.",
        )
        parser.add_argument(
            "--output_directory",  # old: outDir
            default="fasta.crosslinked",
            type=str,
            help="output directory  [%(default)s]",
        )
        parser.add_argument(
            "--tmp_directory",  # old: tmpDir
            default="",
            type=str,
            help="TMP directory to store blat psl files [%(default)s]",
        )
        parser.add_argument(
            "--minimum_coverage",  # old: minCoverage
            default=0.95,
            type=float,
            help="min. coverage for blat hits  [%(default)s]",
        )
        parser.add_argument(
            "--minimum_identity",  # old: minIdentity
            default=0.98,
            type=float,
            help="Min. identity for blat hits  [%(default)s]",
        )
        parser.add_argument(
            "--protein_limit",  # old: protLimit
            default=0,
            type=int,
            help="only taxa having >l proteins [%(default)s]",
        )
        parser.add_argument(
            "--max_threads",  # old: threads
            default=4,
            type=int,
            help="number of cores [%(default)s]",
        )

        return CrosslinkBySequenceArgs(**vars(parser.parse_args(args)))


def main() -> int:
    t0 = datetime.now()
    args = CrosslinkBySequenceArgs.get_arguments()

    if args.verbose:
        print(f"Options: {str(args)}")
    # create out and tmp directory if they dont exists
    create_folder(o.outDir)
    create_folder(o.tmpDir)

    # xLink proteomes
    process(
        args.fastas,
        args.refSpecie,
        args.outDir,
        args.tmpDir,
        args.threads,
        args.minCoverage,
        args.minIdentity,
        args.protLimit,
    )

    dt = datetime.now() - t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)

    return 0
