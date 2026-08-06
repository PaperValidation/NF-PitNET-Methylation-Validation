import urllib.request, re, os, sys
try:
    import socks, socket
    socket.socket = socks.socksocket
except Exception as e:
    print("no socks module:", e)

def open_via_proxy(url):
    # patch create_connection to route through SOCKS5 127.0.0.1:8428
    import urllib.request as ur, socket as s, socks
    s.socket = socks.socksocket
    orig = s.create_connection
    def cc(*a, **k):
        host, port = a[0][0], a[0][1]
        sk = socks.socksocket(); sk.set_proxy(socks.SOCKS5, "127.0.0.1", 8428)
        sk.connect((host, port)); return sk
    s.create_connection = cc
    req = ur.Request(url, headers={"User-Agent":"Mozilla/5.0"})
    return ur.build_opener().open(req, timeout=60).read()

outdir = r"G:\CMCC Sync\Medical Data\Neurosurgery\NF-PitNET"
os.makedirs(outdir, exist_ok=True)
url = "https://pmc.ncbi.nlm.nih.gov/articles/PMC13070495/?report=printable"
try:
    data = open_via_proxy(url).decode("utf-8", "ignore")
    print("HTML_LEN", len(data))
    html_path = os.path.join(outdir, "pmc13070495_printable.html")
    open(html_path, "w", encoding="utf-8").write(data)
    print("SAVED", html_path)
    # strip tags for text scanning
    text = re.sub(r"<[^>]+>", " ", data)
    text = re.sub(r"&[a-z]+;", " ", text)
    for kw in ["562","CCDC141","ATP1A3","ANK2","CASC15","HDC","NCKAP5","GUCA1C","Supplementary File S2","probes"]:
        print(f"--{kw}:", "FOUND" if kw in data else "no")
    gses = sorted(set(re.findall(r"GSE\d{5,7}", data)))
    print("GSEs in article:", gses)
    # look for supplementary file links
    supp = sorted(set(re.findall(r"https?://[^\"\'<\s]*supplementary[^\"\'<\s]*", data, re.I)))
    supp += sorted(set(re.findall(r"https?://[^\"\'<\s]*\.zip[^\"\'<\s]*", data, re.I)))
    print("SUPP_URLS:", supp[:30])
except Exception as e:
    print("ERR", type(e).__name__, e)
