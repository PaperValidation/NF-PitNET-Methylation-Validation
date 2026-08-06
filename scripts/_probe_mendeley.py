import urllib.request, json, socket

def fetch(url, proxy=None, timeout=90, accept="application/json"):
    handlers = []
    if proxy:
        import socks
        socks.set_default_proxy(socks.PROXY_TYPE_SOCKS5, proxy[0], proxy[1])
        socket.socket = socks.socksocket
    req = urllib.request.Request(url, headers={
        "User-Agent": "Mozilla/5.0",
        "Accept": accept,
    })
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = r.read()
            return r.status, r.headers.get("Content-Type",""), data
    except urllib.error.HTTPError as e:
        try: body = e.read()[:1200]
        except: body=b""
        return e.code, "HTTPError", body
    except Exception as e:
        return None, type(e).__name__, repr(e).encode()

proxies = [("127.0.0.1",8428), ("127.0.0.1",1987)]
endpoints = [
    "https://data.mendeley.com/api/v1/datasets/vr2wb38r6t/versions/1/files",
    "https://data.mendeley.com/api/v1/datasets/vr2wb38r6t/versions/1",
]

for label, px in [("DIRECT", None)] + [("PROXY%s:%d"%(h,p), (h,p)) for h,p in proxies]:
    print("="*60)
    print("MODE:", label)
    for u in endpoints:
        st, ct, data = fetch(u, proxy=px, timeout=90)
        print("  URL:", u)
        print("    status/err:", st, "| CT:", ct, "| len:", len(data))
        if isinstance(data, bytes) and len(data) < 2500:
            print("    body:", data.decode("utf-8","replace"))
        elif isinstance(data, bytes):
            # peek for json vs html
            head = data[:200].decode("utf-8","replace")
            print("    head:", head.replace("\n"," "))
    # reset socket for next mode
    import socks
    socks.set_default_proxy(None)
