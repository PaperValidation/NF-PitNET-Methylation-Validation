#!/usr/bin/env python3
"""Download GSE115783 series matrix from NCBI FTP (HTTP, no proxy needed).

GEO's ftp.ncbi.nlm.nih.gov is reachable over plain HTTP from this host, unlike
www.ncbi.nlm.nih.gov (HTTPS handshake times out). Resumable + progress.
"""
import os, sys, time, urllib.request

URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115783/matrix/GSE115783_series_matrix.txt.gz"
OUT = r"G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/GSE115783/GSE115783_series_matrix.txt.gz"

def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    headers = {"User-Agent": "Mozilla/5.0"}
    # resume if partial
    start = 0
    if os.path.exists(OUT):
        start = os.path.getsize(OUT)
    mode = "ab" if start else "wb"
    if start:
        headers["Range"] = f"bytes={start}-"
        print(f"[resume] {start} bytes already present", flush=True)
    req = urllib.request.Request(URL, headers=headers)
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=60) as resp, open(OUT, mode) as f:
        total = int(resp.headers.get("Content-Length", 0)) + start
        got = start
        chunk = 1024 * 1024
        while True:
            b = resp.read(chunk)
            if not b:
                break
            f.write(b)
            got += len(b)
            dt = time.time() - t0
            spd = got / dt / 1e6 if dt > 0 else 0
            pct = 100 * got / total if total else 0
            print(f"\r[dl] {got/1e6:.1f}/{total/1e6:.1f} MB ({pct:.0f}%) {spd:.1f} MB/s", end="", flush=True)
    print(f"\n[done] {OUT} ({os.path.getsize(OUT)/1e6:.1f} MB)", flush=True)

if __name__ == "__main__":
    main()
