import urllib.request, re, socket
socket.setdefaulttimeout(30)
HDR = {
 "User-Agent":"Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36",
 "Accept":"text/html,application/xhtml+xml,*/*;q=0.8",
 "Referer":"https://www.ncbi.nlm.nih.gov/geo/",
 "Connection":"keep-alive"}
# 1) fetch the series page HTML and grep for series_matrix / download links
url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115783"
try:
    req=urllib.request.Request(url, headers=HDR)
    html=urllib.request.urlopen(req, timeout=30).read().decode("utf-8","ignore")
    print("PAGE len:", len(html))
    links=re.findall(r'href="([^"]*series_matrix[^"]*)"', html)
    print("series_matrix links found:", len(links))
    for l in links[:10]:
        print("  ", l)
    # also any geo/download links
    dls=re.findall(r'(https://www\.ncbi\.nlm\.nih\.gov/geo/download/[^\"]+)', html)
    print("geo/download links:", len(dls))
    for l in dls[:10]:
        print("  ", l)
except Exception as e:
    print("ERR", type(e).__name__, str(e)[:120])
# 2) read esearch body
try:
    req=urllib.request.Request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSE115783", headers=HDR)
    body=urllib.request.urlopen(req, timeout=20).read().decode("utf-8","ignore")
    print("ESEARCH body:", body[:500])
except Exception as e:
    print("ESEARCH ERR", type(e).__name__, str(e)[:120])
