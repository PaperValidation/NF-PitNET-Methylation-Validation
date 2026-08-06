"""Merge 3 standalone PNGs into Figure 4 (A/B/C) multi-panel, 300 dpi."""
from PIL import Image, ImageDraw, ImageFont
import os

SRC = r'C:\Users\tidus\PycharmProjects\Neurosurgery-NF-PitNET\results\Moller2025\strict_projection'
OUT = r'C:\Users\tidus\PycharmProjects\Neurosurgery-NF-PitNET\results\Moller2025\figure4_merged.png'

files = {
    'A': 'km_strict_riskhi.png',
    'B': 'forest_strict.png',
    'C': 'fig5c_lineage_bar.png'
}

# Load images
imgs = {}
for label, fn in files.items():
    p = os.path.join(SRC, fn)
    imgs[label] = Image.open(p)
    print(f'Panel {label}: {imgs[label].size}')

# Determine layout: A + B row 1, C row 2 centered (since C is narrower)
# or A top, B+C bottom depending on aspect ratios
a_w, a_h = imgs['A'].size  # 900x675
b_w, b_h = imgs['B'].size  # 900x570
c_w, c_h = imgs['C'].size  # 825x570

# Row 1: A alone (full width KM plot)
# Row 2: B (forest) + C (lineage bar) side by side
row1_w = max(a_w, b_w + c_w + 40)  # A is 900, B+C is 900+825+40=1765 → need B+C to fit
# Actually let's do: Row 1 = A (KM), Row 2 = B+C side by side
# Panel widths: B=900, C=825, total=1725 with gap. A=900. Make width = max(900+20, 1725+20) = 1745
margin = 20
gap = 30

# Row heights
row1_h = a_h   # 675
row2_h = max(b_h, c_h)  # 570

total_w = max(a_w, b_w + c_w + gap) + 2 * margin  # max(900, 1725) = 1725 + 40 = 1765. Let's use 1800.
total_w = 1800
total_h = margin + row1_h + gap + row2_h + margin  # 20+675+30+570+20=1315

# Create canvas at 300 dpi scale
scale = 2  # upscale for 300dpi print
canvas = Image.new('RGB', (total_w * scale, total_h * scale), 'white')
draw = ImageDraw.Draw(canvas)

# Panel A - centered
a_scaled = imgs['A'].resize((a_w * scale, a_h * scale), Image.LANCZOS)
ax = (total_w * scale - a_scaled.width) // 2
ay = margin * scale
canvas.paste(a_scaled, (ax, ay))
# Add "A" label
draw.text((margin * scale + 5, ay + 5), 'A', fill='black')

# Panel B - left side of row 2
b_scaled = imgs['B'].resize((b_w * scale, b_h * scale), Image.LANCZOS)
bx = (total_w * scale - (b_w + c_w + gap) * scale) // 2
by = (margin + row1_h + gap) * scale
canvas.paste(b_scaled, (bx, by))
draw.text((bx + 5, by + 5), 'B', fill='black')

# Panel C - right side of row 2
c_scaled = imgs['C'].resize((c_w * scale, c_h * scale), Image.LANCZOS)
cx = bx + b_scaled.width + gap * scale
cy = by
canvas.paste(c_scaled, (cx, cy))
draw.text((cx + 5, cy + 5), 'C', fill='black')

# Save at 300 dpi
canvas.save(OUT, dpi=(300, 300))
print(f'Saved: {OUT} ({canvas.size[0]}x{canvas.size[1]})')
print(f'Print size at 300dpi: {canvas.size[0]/300:.1f} x {canvas.size[1]/300:.1f} inches')
