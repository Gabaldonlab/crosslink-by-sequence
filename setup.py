from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from setuptools import setup

def remove_old_egg_info() -> None:
    tar_command = ["rm", "-rf", "crosslink_by_sequence.egg-info"]
    try:
        subprocess.run(tar_command, check=True)
        print("Old .egg-info directory clean up completed successfully.")
    except subprocess.CalledProcessError as e:
        print(
            f"Error occurred during old .egg-info directory clean up: {e}",
            file=sys.stderr,
        )
        raise e


@dataclass
class DiamondSetup:
    diamond_download_url: str = (
        "http://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz"
    )
    diamond_tar_file: str = diamond_download_url.split("/")[-1]

    @classmethod
    def _download_diamond(cls) -> None:
        # Execute the wget command to download the file
        wget_command: list[str] = ["wget", cls.diamond_download_url]
        try:
            subprocess.run(wget_command, check=True)
            print("Download completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error occurred during download: {e}", file=sys.stderr)
            raise e

    @classmethod
    def _extract_diamond(cls) -> None:
        output_dir = Path("./crosslink_by_sequence/bin")
        shutil.rmtree(str(output_dir), ignore_errors=True)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Execute the tar command to extract the file
        tar_command = ["tar", "-C", str(output_dir), "-xzf", cls.diamond_tar_file]
        try:
            subprocess.run(tar_command, check=True)
            print("Extraction completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error occurred during extraction: {e}", file=sys.stderr)
            raise e

    @classmethod
    def _remove_diamond_junks(cls) -> None:
        # Execute the tar command to extract the file
        tar_command = ["rm", "-rf", cls.diamond_tar_file]
        try:
            subprocess.run(tar_command, check=True)
            print("Diamond clean up completed successfully.")
        except subprocess.CalledProcessError as e:
            print(
                f"Error occurred during Diamond clean up: {e}", file=sys.stderr
            )
            raise e

    @classmethod
    def setup(cls) -> None:
        cls._download_diamond()
        cls._extract_diamond()
        cls._remove_diamond_junks()


remove_old_egg_info()
DiamondSetup.setup()
setup()
