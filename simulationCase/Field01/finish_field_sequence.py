"""Embed original 600-dpi MATLAB field tiles in vector PDF/SVG exports.

MATLAB's painter embeds display-resolution versions by default. Replace only
the six existing image objects, without moving plots or touching data/text.
"""
from pathlib import Path
import base64
import json
import re
import sys
import fitz

root = Path(sys.argv[1])
pdf = root / 'Fig_field_sequence_strict_AB.pdf'
tiles = [root / 'source_notes' / f'field_{c}_600dpi.png' for c in 'abcdef']
doc = fitz.open(pdf)
assert len(doc) == 1
refs = []
for xref in range(1, doc.xref_length()):
    if doc.xref_get_key(xref, 'Subtype')[1] != '/Image':
        continue
    w = int(doc.xref_get_key(xref, 'Width')[1])
    h = int(doc.xref_get_key(xref, 'Height')[1])
    if .98 < w / h < 1.02:
        refs.append(xref)
assert len(refs) == 6, f'Expected six square field images, got {refs}'
text_before = doc[0].get_text()
for xref, tile in zip(refs, tiles):
    doc[0].replace_image(xref, filename=str(tile))
assert doc[0].get_text() == text_before
tmp = pdf.with_name(pdf.stem + '_export_tmp.pdf')
doc.save(tmp, garbage=4, deflate=True)
doc.close()
tmp.replace(pdf)

svg = pdf.with_suffix('.svg')
content = svg.read_text(encoding='utf-8')
matches = list(re.finditer(r'(xlink:href=")data:image/[^\"]+(\")', content))
assert len(matches) == 7, 'Expected six field tiles followed by the colourbar'
for m, tile in reversed(list(zip(matches[:6], tiles))):
    value = m[1] + 'data:image/png;base64,' + base64.b64encode(tile.read_bytes()).decode('ascii') + m[2]
    content = content[:m.start()] + value + content[m.end():]
svg.write_text(content, encoding='utf-8')
with fitz.open(pdf) as check:
    info = check[0].get_image_info()
    squares = [v for v in info if .98 < v['width'] / v['height'] < 1.02]
    assert len(squares) == 6 and min(v['width'] for v in squares) >= 1000
    quality = {'pdf_page_mm': [x * 25.4 / 72 for x in (check[0].rect.width, check[0].rect.height)],
               'field_tiles': [{'pixels': [v['width'], v['height']],
                                'effective_dpi_x': v['width'] / ((v['bbox'][2] - v['bbox'][0]) / 72)} for v in squares],
               'native_vector_fonts': sorted(set(v[3] for v in check[0].get_fonts()))}
    (root / 'export_quality.json').write_text(json.dumps(quality, indent=2), encoding='utf-8')
    check[0].get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False).save(root / 'source_notes' / 'pdf_check.png')
print('PDF/SVG: six full-resolution field tiles; original vector text and data charts retained.')
