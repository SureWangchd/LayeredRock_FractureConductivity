"""Restore the two native 600 dpi geometry images in MATLAB's vector PDF."""
from pathlib import Path
import re
import sys
import fitz

root = Path(sys.argv[1])
pdf = root/'Fig09_CT_digital_rock.pdf'
doc = fitz.open(pdf)
assert len(doc) == 1
text_before = doc[0].get_text()
patterns = []
sections = []
for xref in range(1,doc.xref_length()):
    obj = doc.xref_object(xref)
    if '/PatternType 1' not in obj:
        continue
    matrix = re.search(r'/Matrix\s*\[([^]]+)\]',obj)
    image = re.search(r'/Im\d+\s+(\d+)\s+0\s+R',obj)
    if matrix and image:
        x = float(matrix[1].split()[4]); ref = int(image[1])
        if int(doc.xref_get_key(ref,'Width')[1])>128:
            patterns.append((x,ref))
        else:
            sections.append(ref)
assert len(patterns) == 2, 'Expected two 3D geometry panels'
assert len(sections) == 2, 'Expected two categorical S1 material sections'
for j,(_,xref) in enumerate(sorted(patterns),start=1):
    doc[0].replace_image(xref,filename=str(root/'source_notes'/f'geometry_{j}_600dpi.png'))
for xref in sections:
    doc[0].replace_image(xref,filename=str(root/'source_notes'/'ct_S1_section_native128.png'))
assert doc[0].get_text() == text_before
tmp = pdf.with_name(pdf.stem+'_export_tmp.pdf')
doc.save(tmp,garbage=4,deflate=True);doc.close();tmp.replace(pdf)
with fitz.open(pdf) as check:
    info=check[0].get_image_info()
    assert len([x for x in info if x['width']>1600])==2
    assert len([x for x in info if x['width']==128 and x['height']==128])==2
print('PDF: 600 dpi CT geometry, native 128x128 material images; vector data, fracture-section polygons and text retained')
