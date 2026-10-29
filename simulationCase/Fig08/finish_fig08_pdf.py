"""Restore the 600 dpi geometry bitmap lost by MATLAB's PDF image downsampling.

Only replace the image object inside its existing PDF pattern. Its transform,
draw order, vector data axes, and vector text are unchanged. Requires PyMuPDF.
"""
from pathlib import Path
import base64
import re
import sys
import fitz

root = Path(sys.argv[1])
pdf = root / 'Fig08_finite_three_dimensional.pdf'
bitmap = root / 'source_notes' / 'geometry_view_600dpi.png'
doc = fitz.open(pdf)
image_refs = [x for x in range(1, doc.xref_length())
              if doc.xref_get_key(x, 'Subtype')[1] == '/Image']
assert len(doc) == 1 and len(image_refs) == 1, 'Expected one geometry image'
text_before = doc[0].get_text()
doc[0].replace_image(image_refs[0], filename=str(bitmap))
assert doc[0].get_text() == text_before, 'Vector text changed'
tmp = pdf.with_name(pdf.stem + '_export_tmp.pdf')
doc.save(tmp, garbage=4, deflate=True)
doc.close()
tmp.replace(pdf)
with fitz.open(pdf) as check:
    info = check[0].get_image_info()
    assert len(info) == 1 and info[0]['width'] >= 1700
svg = pdf.with_suffix('.svg')
if svg.exists():
    content = svg.read_text(encoding='utf-8')
    data_url = 'data:image/png;base64,' + base64.b64encode(bitmap.read_bytes()).decode('ascii')
    content, count = re.subn(r'(xlink:href=")data:image/[^\"]+(\")',
                            lambda m: m[1] + data_url + m[2], content)
    assert count == 1, 'Expected one SVG geometry image'
    svg.write_text(content, encoding='utf-8')
print('PDF/SVG: native vector data/text + full-resolution geometry bitmap')
