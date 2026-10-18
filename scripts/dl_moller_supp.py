import re, hashlib, os, sys
import requests
from _paths import DATASET_PATH as BASE

PROXIES = None  # direct connection to NCBI works; proxy IP was blocked
HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Referer": "https://pmc.ncbi.nlm.nih.gov/articles/PMC13070495/",
}
ZIP_URL = "https://pmc.ncbi.nlm.nih.gov/articles/instance/13070495/bin/noaf269_supplementary_data.zip"
OUT = os.path.join(BASE, "Moller2025_supplement", "noaf269_supplementary_data.zip")
os.makedirs(os.path.dirname(OUT), exist_ok=True)

def get_challenge():
    r = requests.get(ZIP_URL, proxies=PROXIES, headers=HEADERS, timeout=60)
    html = r.text
    m = re.search(r'POW_CHALLENGE\s*=\s*"([^"]+)"', html)
    d = re.search(r'POW_DIFFICULTY\s*=\s*"(\d+)"', html)
    if not m:
        return None, None, html
    return m.group(1), int(d.group(1)), html

def solve(challenge, difficulty):
    prefix = "0" * difficulty
    nonce = 0
    while True:
        h = hashlib.sha256((challenge + str(nonce)).encode("utf-8")).hexdigest()
        if h.startswith(prefix):
            return nonce, h
        nonce += 1
        if nonce % 200000 == 0:
            print("  ...tried", nonce)

def main():
    challenge, difficulty, html = get_challenge()
    if challenge is None:
        print("No POW challenge found. Body[:300]:", html[:300])
        sys.exit(1)
    print(f"Challenge={challenge} difficulty={difficulty}")
    nonce, h = solve(challenge, difficulty)
    print(f"Solved nonce={nonce} hash={h}")
    cookie = f"cloudpmc-viewer-pow={challenge},{nonce}"
    # download with cookie
    sess = requests.Session()
    r = sess.get(ZIP_URL, proxies=PROXIES,
                 headers={**HEADERS, "Cookie": cookie,
                           "Accept": "application/zip,*/*"},
                 timeout=120)
    ct = r.headers.get("Content-Type", "")
    print(f"Download status={r.status_code} len={len(r.content)} ct={ct}")
    if r.status_code == 200 and len(r.content) > 50000 and "text/html" not in ct:
        with open(OUT, "wb") as f:
            f.write(r.content)
        print("SAVED", OUT)
        # list zip
        import zipfile, io
        with zipfile.ZipFile(io.BytesIO(r.content)) as z:
            print("=== ZIP entries ===")
            for n in z.namelist():
                print("  ", n)
    else:
        print("Download failed. Body[:300]:", r.content[:300])

if __name__ == "__main__":
    main()
