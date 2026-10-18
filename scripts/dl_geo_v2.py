"""Robust GEO/ftp downloader via Windows-host SOCKS5 proxy (127.0.0.1:8428).

Fixes two issues seen with plain urllib:
 1. NCBI ftp-host requires a Range header on supplementary files (else 403).
 2. Use remote DNS resolution (rdns=True) to mirror curl's socks5h behavior.
Skips the >16KB write bug by streaming through Python (not curl -o).
"""
import sys, os, re, time, urllib.request
import socks

from _paths import DATASET_PATH as GDRIVE

def patch_proxy(rdns=True):
    orig = urllib.request.socket.create_connection
    def cc(*a, **k):
        host = a[0][0]; port = a[0][1]
        sx = socks.socksocket()
        sx.set_proxy(socks.SOCKS5, "127.0.0.1", 8428, rdns=rdns)
        sx.settimeout(120); sx.connect((host, port)); return sx
    urllib.request.socket.create_connection = cc
    return orig

def winpath(p):
    m = re.match(r"^/([a-zA-Z])/(.*)$", p)
    if m:
        return f"{m.group(1)}:/{m.group(2)}"
    return p

def dl(url, out, chunk=1024*1024):
    out = winpath(out)
    d = os.path.dirname(out)
    if d:
        os.makedirs(d, exist_ok=True)
    orig = patch_proxy(rdns=True)
    try:
        req = urllib.request.Request(url, headers={
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)",
            "Range": "bytes=0-",
        })
        t = time.time(); total = 0
        with urllib.request.urlopen(req, timeout=180) as r, open(out, "wb") as f:
            while True:
                b = r.read(chunk)
                if not b:
                    break
                f.write(b); total += len(b)
        print(f"OK {total:,} bytes -> {out} in {time.time()-t:.1f}s")
    except Exception as e:
        print("ERR", type(e).__name__, e)
        return False
    finally:
        urllib.request.socket.create_connection = orig
    return True

if __name__ == "__main__":
    url = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(GDRIVE, os.path.basename(url))
    dl(url, out)
