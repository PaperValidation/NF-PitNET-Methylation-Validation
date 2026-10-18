"""Centralized dataset paths loaded from the project root ``.env``.

This is the single source of truth for the two data locations:

* ``DATASET_PATH``      -> raw / permanent datasets (E:/Medical Data/..., on disk E:)
* ``DATASET_PATH_TEMP`` -> intermediate / temp data   (F:/Data Temp/..., local, not synced)

Scripts should import from here instead of hardcoding the base path, so that a
path change only needs to happen in ``.env``.

Usage (from a script living inside ``scripts/``)::

    from _paths import (
        DATASET_PATH, DATASET_PATH_TEMP,
        DATASET_PATH_WSL, DATASET_PATH_TEMP_WSL,
    )
    # DATASET_PATH      = "E:/Medical Data/Neurosurgery/NF-PitNET"
    # DATASET_PATH_WSL  = "/mnt/e/Medical Data/Neurosurgery/NF-PitNET"
"""
import os
import re

_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_ENV_PATH = os.path.join(_PROJECT_ROOT, ".env")


def _load_env() -> dict:
    """Minimal .env parser (no external dependency)."""
    env: dict = {}
    try:
        with open(_ENV_PATH, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                m = re.match(r"^([A-Za-z0-9_]+)\s*=\s*'?(.*?)'?\.?\s*$", line)
                if m:
                    env[m.group(1)] = m.group(2)
    except FileNotFoundError:
        pass
    return env


_ENV = _load_env()


def _win_to_wsl(win: str) -> str:
    """Convert a Windows path like ``X:/a/b`` to ``/mnt/x/a/b``."""
    if not win or len(win) < 2:
        return win
    drive = win[0].lower()
    rest = win[2:].replace("\\", "/").lstrip("/")
    return f"/mnt/{drive}/{rest}"


DATASET_PATH = _ENV.get("DATASET_PATH", "E:/Medical Data/Neurosurgery/NF-PitNET")
DATASET_PATH_TEMP = _ENV.get("DATASET_PATH_TEMP", "F:/Data Temp/NF-PitNET")
PROJECT_PATH = _PROJECT_ROOT

DATASET_PATH_WSL = _win_to_wsl(DATASET_PATH)
DATASET_PATH_TEMP_WSL = _win_to_wsl(DATASET_PATH_TEMP)
PROJECT_PATH_WSL = _win_to_wsl(PROJECT_PATH)


if __name__ == "__main__":
    print("PROJECT_PATH         =", PROJECT_PATH)
    print("DATASET_PATH         =", DATASET_PATH)
    print("DATASET_PATH_TEMP    =", DATASET_PATH_TEMP)
    print("DATASET_PATH_WSL     =", DATASET_PATH_WSL)
    print("DATASET_PATH_TEMP_WSL=", DATASET_PATH_TEMP_WSL)
