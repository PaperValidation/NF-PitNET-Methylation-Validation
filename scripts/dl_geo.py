import sys, os, time, re
import urllib.request
import socks

from _paths import DATASET_PATH as GDRIVE
os.makedirs(GDRIVE, exist_ok=True)

def make_opener():
    orig = urllib.request.socket.create_connection
    def cc(*a, **k):
        host = a[0][0]; port = a[0][1]
        sx = socks.socksocket(); sx.set_proxy(socks.SOCKS5, "127.0.0.1", 8428)
        sx.settimeout(60); sx.connect((host, port)); return sx
    urllib.request.socket.create_connection = cc
    return urllib.request.build_opener(), orig

def winpath(p):
    # accept /g/... or G:/... -> G:/...
    m = re.match(r"^/([a-zA-Z])/(.*)$", p)
    if m:
        return f"{m.group(1)}:/{m.group(2)}"
    return p

def dl(url, out):
    out = winpath(out)
    d = os.path.dirname(out); 
    if d: os.makedirs(d, exist_ok=True)
    op, orig = make_opener()
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    try:
        t = time.time(); data = op.open(req, timeout=120).read()
        with open(out, "wb") as f: f.write(data)
        print(f"OK {len(data):,} bytes -> {out} in {time.time()-t:.1f}s")
    except Exception as e:
        print("ERR", type(e).__name__, e)
    finally:
        urllib.request.socket.create_connection = orig

if __name__ == "__main__":
    # url out  (out can be relative -> placed under GDRIVE)
    url = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else None
    if not out:
        out = os.path.join(GDRIVE, os.path.basename(url))
    dl(url, out)
