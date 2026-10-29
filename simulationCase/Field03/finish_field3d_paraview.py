"""Restore full-resolution tiles in the ParaView/MATLAB PDF and SVG.

Matches image pixels, not PDF/SVG appearance order. For each format, compare
all n embedded rasters with all n source tiles, then enumerate all n!
assignments. Existing vector labels and image placement remain unchanged.
No source-data scans, hashes, COMSOL calls, or new visualization calculations.

Run with the Python environment containing PyMuPDF and Pillow (no NumPy):
  python finish_field3d_paraview.py F:/.../Field3D/first_version_refined --panels abcdef
  # Defaults to Fig_3D_current_fields.pdf/svg; auto-discovers panel_a.png onward.
  # Other output filenames remain supported by --stem.
"""
from __future__ import annotations

import argparse
import base64
import io
import itertools
import json
import math
from pathlib import Path
import re
import shutil

import fitz
from PIL import Image, ImageChops, ImageDraw, ImageFilter, ImageStat


RESAMPLE = getattr(Image, "Resampling", Image)
TRANSPOSE = getattr(Image, "Transpose", Image)
ORIENTATIONS = {
    "identity": None,
    "flip_left_right": TRANSPOSE.FLIP_LEFT_RIGHT,
    "flip_top_bottom": TRANSPOSE.FLIP_TOP_BOTTOM,
    "rotate_180": TRANSPOSE.ROTATE_180,
    "rotate_90": TRANSPOSE.ROTATE_90,
    "rotate_270": TRANSPOSE.ROTATE_270,
    "transpose": TRANSPOSE.TRANSPOSE,
    "transverse": TRANSPOSE.TRANSVERSE,
}
FILTERS = {"lanczos": RESAMPLE.LANCZOS, "bicubic": RESAMPLE.BICUBIC, "bilinear": RESAMPLE.BILINEAR, "box": RESAMPLE.BOX, "nearest": RESAMPLE.NEAREST}
SCORE_BLUR_PX = 1.0
IMAGE_TAG = re.compile(r"<(?:[\w.-]+:)?image\b[^>]*>", re.I | re.S)
DATA_HREF = re.compile(r"\b(?:xlink:href|href)\s*=\s*(?P<q>['\"])(?P<uri>data:image/[^'\"]+)(?P=q)", re.I | re.S)


def rgb_image(data):
    im = Image.open(io.BytesIO(data)) if isinstance(data, bytes) else Image.open(data)
    im.load()
    if im.mode == "RGBA" or "transparency" in im.info:
        rgba = im.convert("RGBA")
        white = Image.new("RGBA", rgba.size, "white")
        white.alpha_composite(rgba)
        return white.convert("RGB")
    return im.convert("RGB")


def orient(im, name):
    op = ORIENTATIONS[name]
    return im if op is None else im.transpose(op)


def compatible_size(a, b):
    tolerance = 0.025 + 2.0 / min(a)
    return abs(math.log((a[0]/a[1])/(b[0]/b[1]))) < tolerance


def plausible_raster(size, tiles):
    return min(size) >= 24 and any(compatible_size(size, im.size) or compatible_size(size, im.size[::-1]) for im in tiles)


def foreground_mask(im):
    r,g,b = im.split()
    return ImageChops.darker(ImageChops.darker(r,g),b).point(lambda x: 255 if x < 248 else 0)


def raw_pixel_score(candidate, scaled_source):
    diff = ImageChops.difference(candidate, scaled_source)
    mask = ImageChops.lighter(foreground_mask(candidate), foreground_mask(scaled_source))
    if mask.getbbox() is None:
        return 0.0
    return sum(ImageStat.Stat(diff, mask).mean)/(3*255)


def pixel_score(candidate, scaled_source):
    # MATLAB's 96-dpi image embedding aliases subpixel lines, even when the
    # final source is 600 dpi. Compare one-pixel low-pass copies so sampling
    # phase does not dominate identity matching. Exported tiles are untouched.
    if SCORE_BLUR_PX > 0:
        kernel=ImageFilter.GaussianBlur(SCORE_BLUR_PX)
        candidate=candidate.filter(kernel);scaled_source=scaled_source.filter(kernel)
    return raw_pixel_score(candidate,scaled_source)


def best_pair(candidate, source):
    # Compare at the embedded raster size; cap only unexpectedly large or
    # already-restored images to keep re-runs inexpensive.
    target = candidate.copy()
    target.thumbnail((512,512), RESAMPLE.LANCZOS)
    best = None
    for name in ORIENTATIONS:
        shape = source.size[::-1] if name in ("rotate_90","rotate_270","transpose","transverse") else source.size
        if not compatible_size(candidate.size, shape):
            continue
        for kernel, resampling in FILTERS.items():
            # Resize before the orientation operation to avoid repeatedly
            # transposing a multi-megapixel source for small matching images.
            pre_shape = target.size[::-1] if name in ("rotate_90","rotate_270","transpose","transverse") else target.size
            scaled = orient(source.resize(pre_shape, resampling), name)
            score = pixel_score(target, scaled)
            if best is None or score < best["mae"]:
                best = {"mae":score,"raw_mae_without_comparison_blur":raw_pixel_score(target,scaled),
                        "orientation":name,"resampling":kernel,"comparison_pixels":list(target.size)}
    return best


def global_match(candidates, tiles, labels, kind, max_mae):
    n = len(tiles)
    if len(candidates) != n:
        raise RuntimeError(f"{kind}: expected {n} plausible panel images; found {len(candidates)}. No files changed. Candidates: {[v['id'] for v in candidates]}")
    matrix = [[best_pair(c["image"],tile) for tile in tiles] for c in candidates]
    rankings=[]
    for permutation in itertools.permutations(range(n)):
        chosen=[matrix[i][j] for i,j in enumerate(permutation)]
        if all(p is not None for p in chosen):
            rankings.append((sum(p["mae"] for p in chosen),permutation))
    if not rankings:
        raise RuntimeError(f"{kind}: no valid {n}-panel correspondence. No files changed.")
    rankings.sort(key=lambda item:item[0])
    best_cost, permutation = rankings[0]
    gap = rankings[1][0]-best_cost if len(rankings)>1 else None
    if gap is not None and gap <= 1e-10:
        raise RuntimeError(f"{kind}: pixel matching is tied between different panel assignments. Do not guess the order; inspect duplicate/indistinguishable tiles. No files changed.")
    assignment=[]
    for i,j in enumerate(permutation):
        detail=matrix[i][j]
        if detail["mae"] > max_mae:
            raise RuntimeError(f"{kind}: panel {labels[j]} has foreground MAE {detail['mae']:.6f} > {max_mae}. The tile may have a different crop. No files changed.")
        assignment.append({"candidate_index":i,"candidate_id":candidates[i]["id"],"source_index":j,"panel":labels[j],
                           "embedded_pixels_before":list(candidates[i]["image"].size),"source_pixels":list(tiles[j].size),**detail})
    report={"method":f"foreground RGB mean absolute pixel difference; global enumeration of all {math.factorial(n)} permutations; independent PDF/SVG matching",
            "comparison_blur_pixels":SCORE_BLUR_PX,"comparison_blur_scope":"Only low-resolution matching copies; final embedded source pixels are not blurred or resampled.",
            "number_of_panels":n,"permutations_enumerated":math.factorial(n),"valid_permutations":len(rankings),"best_total_mae":best_cost,"second_best_total_mae":rankings[1][0] if len(rankings)>1 else None,
            "second_best_gap":gap,"candidate_ids":[c["id"] for c in candidates],"source_panels":labels,
            "mae_matrix":[[None if p is None else p["mae"] for p in row] for row in matrix],"assignment":assignment,
            "top_assignments":[{"total_mae":score,"panels_by_candidate":[labels[j] for j in perm]} for score,perm in rankings[:5]]}
    return assignment,report


def pdf_image_references(doc):
    """Find executed Image XObjects, including MATLAB's tiling Patterns.

    Page.get_images() omits image resources inside tiling Patterns. These
    are still ordinary Image XObjects, not BI/ID/EI inline images. Follow
    named Do and scn/SCN resource references without changing content streams.
    """
    found=[];seen_images=set();visited=set()

    def resource_map(owner, kind):
        typ,value=doc.xref_get_key(owner,"Resources/"+kind)
        if typ=="xref":value=doc.xref_object(int(value.split()[0]))
        return {m.group(1):int(m.group(2)) for m in
                re.finditer(r"/([^\s<>\[\](){}/%]+)\s+(\d+)\s+\d+\s+R",value)}

    def walk(owner, content, path, inherited_owner=None):
        if owner in visited:return
        visited.add(owner)
        for kind,operator in [("XObject",rb"Do"),("Pattern",rb"(?:scn|SCN)")]:
            refs=resource_map(owner,kind)
            if not refs and inherited_owner is not None:refs=resource_map(inherited_owner,kind)
            used={m.group(1).decode("latin1") for m in re.finditer(
                rb"/([^\s<>\[\](){}/%]+)\s+"+operator+rb"\b",content)}
            for name,xref in refs.items():
                if name not in used:continue
                next_path=path+[f"{kind}/{name}"]
                if doc.xref_is_image(xref):
                    if xref not in seen_images:
                        seen_images.add(xref)
                        found.append({"xref":xref,"resource_path":" -> ".join(next_path)})
                elif doc.xref_is_stream(xref):
                    walk(xref,doc.xref_stream(xref),next_path,owner)

    content=b"\n".join(doc.xref_stream(xref) for xref in doc[0].get_contents())
    walk(doc[0].xref,content,["page"])
    return found


def pdf_candidates(doc, tiles):
    out=[]
    for reference in pdf_image_references(doc):
        xref=reference["xref"]
        width=int(doc.xref_get_key(xref,"Width")[1]);height=int(doc.xref_get_key(xref,"Height")[1])
        if not plausible_raster((width,height),tiles):continue
        image=rgb_image(doc.extract_image(xref)["image"])
        mask_type,mask_ref=doc.xref_get_key(xref,"SMask")
        smask=int(mask_ref.split()[0]) if mask_type=="xref" else 0
        if smask:
            alpha=Image.open(io.BytesIO(doc.extract_image(smask)["image"])).convert("L")
            rgba=image.convert("RGBA");rgba.putalpha(alpha)
            base=Image.new("RGBA",rgba.size,"white");base.alpha_composite(rgba);image=base.convert("RGB")
        out.append({"id":f"xref_{xref}","xref":xref,"image":image,"resource_path":reference["resource_path"]})
    return out


def svg_candidates(content, tiles):
    out=[]
    for tag_index,match in enumerate(IMAGE_TAG.finditer(content)):
        href=DATA_HREF.search(match.group())
        if not href:
            continue
        uri=href.group("uri")
        if ";base64," not in uri:
            raise RuntimeError("SVG image URI is not base64 encoded; no files changed.")
        image=rgb_image(base64.b64decode(uri.split(",",1)[1]))
        if plausible_raster(image.size,tiles):
            begin=match.start()+href.start("uri");end=match.start()+href.end("uri")
            out.append({"id":f"svg_image_{tag_index}","image":image,"begin":begin,"end":end,"tag_index":tag_index})
    return out


def png_bytes(im):
    buffer=io.BytesIO();im.save(buffer,format="PNG");return buffer.getvalue()


def make_montage(path, candidates, assignment, tiles):
    canvas=Image.new("RGB",(860,len(assignment)*190),"white");draw=ImageDraw.Draw(canvas)
    for row,entry in enumerate(assignment):
        low=candidates[entry["candidate_index"]]["image"].copy()
        matched=orient(tiles[entry["source_index"]],entry["orientation"])
        matched=matched.resize(low.size,FILTERS[entry["resampling"]])
        draw.text((8,row*190+4),f"{entry['candidate_id']} -> panel {entry['panel']} | MAE {entry['mae']:.7f} | {entry['orientation']}",fill="black")
        for column,image in enumerate([low,matched]):
            image.thumbnail((410,160),RESAMPLE.LANCZOS)
            canvas.paste(image,(10+430*column,row*190+25))
    canvas.save(path)


def main():
    global SCORE_BLUR_PX
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("--stem",default="Fig_3D_current_fields")
    parser.add_argument("--panels",default=None,help="Explicit consecutive panel letters (e.g. abcd). Otherwise discover source_tiles/panel_[a-z].png.")
    parser.add_argument("--max-mae",type=float,default=0.04)
    parser.add_argument("--comparison-blur-px",type=float,default=1.0,help="Anti-alias comparison copies only; source/export pixels remain unchanged. Set 0 for raw matching.")
    args=parser.parse_args();root=args.root.resolve()
    if not 0 <= args.comparison_blur_px <= 2:
        raise ValueError("Comparison blur must be 0 to 2 embedded-image pixels.")
    SCORE_BLUR_PX=args.comparison_blur_px
    pdf=root/(args.stem+".pdf");svg=root/(args.stem+".svg")
    labels=list(args.panels) if args.panels is not None else sorted(
        p.stem[-1] for p in (root/"source_tiles").glob("panel_*.png")
        if re.fullmatch(r"panel_[a-z]",p.stem))
    if not 1 <= len(labels) <= 8 or labels != list("abcdefgh"[:len(labels)]):
        raise RuntimeError("Expected 1 to 8 consecutive source panels starting at a; specify --panels if other tiles remain. No files changed.")
    paths=[root/"source_tiles"/f"panel_{letter}.png" for letter in labels]
    missing=[str(p) for p in [pdf,svg,*paths] if not p.is_file()]
    if missing:
        raise FileNotFoundError("Required exports are not ready; no files changed:\n"+"\n".join(missing))
    tiles=[rgb_image(p) for p in paths]
    doc=fitz.open(pdf)
    if len(doc)!=1:
        raise RuntimeError("Expected one-page figure PDF; no files changed.")
    original_text=doc[0].get_text();original_fonts=sorted({v[3] for v in doc[0].get_fonts()})
    pc=pdf_candidates(doc,tiles);original_svg=svg.read_text(encoding="utf-8");sc=svg_candidates(original_svg,tiles)
    # BOTH correspondences must pass before either deliverable is changed.
    pa,pr=global_match(pc,tiles,labels,"PDF",args.max_mae)
    sa,sr=global_match(sc,tiles,labels,"SVG",args.max_mae)
    verification=root/"source_tiles"/"embedding_verification";verification.mkdir(parents=True,exist_ok=True)
    make_montage(verification/"pdf_match_montage.png",pc,pa,tiles)
    make_montage(verification/"svg_match_montage.png",sc,sa,tiles)
    for entry in pa:
        im=orient(tiles[entry["source_index"]],entry["orientation"])
        xref=pc[entry["candidate_index"]]["xref"]
        # Replace the Image object directly. Its existing Pattern/Form/page
        # transform is retained; no new page resource or content is inserted.
        # RGB source tiles are opaque (any source alpha was composited white).
        doc.update_object(xref,f"<< /Type /XObject /Subtype /Image /Width {im.width} /Height {im.height} /BitsPerComponent 8 /ColorSpace /DeviceRGB >>")
        doc.update_stream(xref,im.tobytes(),compress=True)
        restored=rgb_image(doc.extract_image(xref)["image"])
        if restored.size!=im.size or ImageChops.difference(restored,im).getbbox() is not None:
            raise RuntimeError(f"Lossless PDF replacement check failed for panel {entry['panel']}; originals unchanged.")
        entry["embedded_pixels_after"]=list(restored.size);entry["replacement_pixels_exact"]=True
        entry["resource_path"]=pc[entry["candidate_index"]]["resource_path"]
    if doc[0].get_text()!=original_text or sorted({v[3] for v in doc[0].get_fonts()})!=original_fonts:
        raise RuntimeError("PDF vector text/font preservation check failed; originals unchanged.")
    new_svg=original_svg
    for entry in sorted(sa,key=lambda v:sc[v["candidate_index"]]["begin"],reverse=True):
        candidate=sc[entry["candidate_index"]];im=orient(tiles[entry["source_index"]],entry["orientation"])
        uri="data:image/png;base64,"+base64.b64encode(png_bytes(im)).decode("ascii")
        new_svg=new_svg[:candidate["begin"]]+uri+new_svg[candidate["end"]:]
        entry["embedded_pixels_after"]=list(im.size);entry["replacement_pixels_exact"]=True
    strip=lambda s:DATA_HREF.sub(lambda m:m.group(0).replace(m.group("uri"),"IMAGE_PAYLOAD"),s)
    if strip(original_svg)!=strip(new_svg):
        raise RuntimeError("SVG content outside image payloads changed; originals unchanged.")
    # Save to new temporary files, preserving a reviewable first MATLAB export.
    tmp_pdf=pdf.with_name(pdf.stem+"_hires_tmp.pdf");tmp_svg=svg.with_name(svg.stem+"_hires_tmp.svg")
    # Keep xref numbers for direct post-save verification; no image digest scans.
    doc.save(tmp_pdf,garbage=0,deflate=True);doc.close()
    tmp_svg.write_text(new_svg,encoding="utf-8")
    with fitz.open(tmp_pdf) as check:
        assert check[0].get_text()==original_text
        for entry in pa:
            xref=pc[entry["candidate_index"]]["xref"]
            im=rgb_image(check.extract_image(xref)["image"])
            expected=orient(tiles[entry["source_index"]],entry["orientation"])
            assert im.size==expected.size and ImageChops.difference(im,expected).getbbox() is None
        page_mm=[check[0].rect.width*25.4/72,check[0].rect.height*25.4/72]
        check[0].get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False).save(verification/"pdf_after_hires.png")
    for p in [pdf,svg]:
        backup=p.with_name(p.stem+"_matlab_lowdpi_original"+p.suffix)
        if not backup.exists():shutil.copy2(p,backup)
    tmp_pdf.replace(pdf);tmp_svg.replace(svg)
    report={"source_tiles":[str(p) for p in paths],"pdf":pr,"svg":sr,"pdf_page_mm":page_mm,
            "vector_text_unchanged":True,"vector_fonts_unchanged":True,"native_vector_fonts":original_fonts,
            "svg_nonimage_content_unchanged":True,"matching_uses_appearance_order":False,"hashes_computed":False,
            "note":"Only image payloads replaced. Any detected source orientation is preserved; source tiles and simulation data are untouched.",
            "visual_correspondence_checks":[str(verification/"pdf_match_montage.png"),str(verification/"svg_match_montage.png")]}
    (root/"field3d_embedding_validation.json").write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("PDF mapping:",", ".join(f"{e['candidate_id']}->{e['panel']}" for e in pa))
    print("SVG mapping:",", ".join(f"{e['candidate_id']}->{e['panel']}" for e in sa))
    print(f"Restored {len(labels)} full-resolution panels in PDF/SVG; vector text and placement retained.")


if __name__=="__main__":
    main()
