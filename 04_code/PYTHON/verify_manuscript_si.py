"""Check the current manuscript/SI package and write an explicit audit report."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path

from lxml import etree as ET
from PIL import Image
from io import BytesIO

NS = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
      'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
      'v':'urn:schemas-microsoft-com:vml',
      'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
W = '{'+NS['w']+'}'


def load(path):
    with zipfile.ZipFile(path) as z:
        doc = ET.fromstring(z.read('word/document.xml'))
        parts = {n:z.read(n) for n in z.namelist() if n.startswith('word/media/')}
        rels = ET.fromstring(z.read('word/_rels/document.xml.rels'))
        crc_ok = z.testzip() is None
    paras = [''.join(p.xpath('.//w:t/text()',namespaces=NS)) for p in doc.xpath('//w:p',namespaces=NS)]
    return {'path':str(path),'root':doc,'paras':paras,'text':'\n'.join(paras),'media':parts,
            'rels':{r.get('Id'):r.get('Target') for r in rels},'crc_ok':crc_ok}


def normalize_caption(s):
    return re.sub(r'^Figure (\d+)[ |.]+',r'Figure \1. ',s).replace('–','-')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path.cwd())
    p.add_argument('--manuscript',type=Path)
    p.add_argument('--si',type=Path)
    p.add_argument('--baseline-main',type=Path)
    p.add_argument('--output',type=Path,required=True)
    args = p.parse_args()
    root = args.root.resolve()
    ms = load(args.manuscript or root/'01_manuscript/Zhang_main_manuscript_polished.docx')
    si = load(args.si or root/'01_manuscript/Zhang_supplementary_information.docx')
    checks = {'valid_docx_zip_crc':ms['crc_ok'] and si['crc_ok']}
    title = 'Joint associations of residential greenness and ozone with biological ageing'
    checks['matching_title'] = title in ms['paras'] and title in si['paras']
    checks['main_table_count_one'] = len(ms['root'].xpath('//w:tbl',namespaces=NS))==1
    checks['si_table_count_34'] = len(si['root'].xpath('//w:tbl',namespaces=NS))==34
    checks['si_figure_labels_S1_S8'] = all(re.search(rf'Supplementary Figure S{i}\.',si['text']) for i in range(1,9))
    checks['si_table_labels_S1_S30'] = all(re.search(rf'Supplementary Table S{i}\.',si['text']) for i in range(1,31))
    tables = set(int(x) for x in re.findall(r'Tables? S(\d+)',ms['text']))
    figures = set(int(x) for x in re.findall(r'(?:Figures?|Figs?\.) S(\d+)',ms['text']))
    checks['main_supplement_calls_resolve'] = tables <= set(range(1,31)) and figures <= set(range(1,9))
    checks['core_samples_both_documents'] = all(v in d['text'] for d in (ms,si) for v in ['6,472','3,281','5,914','3,154','6,471'])
    checks['primary_p_values_both_documents'] = all(v in d['text'] for d in (ms,si) for v in ['0.016','0.003'])
    checks['joint_did_p_values_both_documents'] = all(v in d['text'] for d in (ms,si) for v in ['0.605','0.779'])
    checks['no_mde_exclusion_claim'] = not re.search(r'(?:exclude|rule out)\s+interaction(?:s| contrasts)?\s+(?:larger|greater)',si['text'],flags=re.I)
    checks['no_unsupported_leakage_free_claim'] = 'leakage-free' not in si['text'].lower()
    checks['no_selection_bias_exoneration'] = 'indicating limited exposure-related selection' not in si['text']
    checks['no_word_field_errors'] = not re.search(r'Error!|Reference source not found|\?\?\?\?',ms['text']+'\n'+si['text'])
    captions = {}
    for i in range(1,5):
        caps = [s for s in ms['paras'] if re.match(rf'^Figure {i}\s*[|.]',s)]
        captions[str(i)] = caps
        checks[f'figure_{i}_inline_and_end_legend_match'] = len(caps)==2 and normalize_caption(caps[0])==normalize_caption(caps[1])
    names = ['Figure1_study_design.png','Figure2_primary_associations.png','Figure3_joint_patterns.png','Figure4_exposure_response.png']
    for i,name in enumerate(names,1):
        checks[f'main_figure_{i}_media_sync'] = ms['media'].get(f'word/media/image{i}.png') == (root/'02_figures/main'/name).read_bytes()
    si_names = ['FigS1_participant_flow.png','FigS2_quartile_scores.png','FigS3_iqr_models.png','FigS4_counterfactual_age.png',
                'FigS5_lightgbm_scenarios.png','FigS6_dose_response_qgcomp.png','FigS7_phenoage_calibration.png','FigS8_lightgbm_nested.png']
    for i,name in enumerate(si_names,1):
        checks[f'si_figure_{i}_media_sync'] = si['media'].get(f'word/media/image{i}.png') == (root/'02_figures/supplementary'/name).read_bytes()
    ratios = []
    for node in si['root'].xpath('//v:shape[v:imagedata]',namespaces=NS):
        rid = node.find('v:imagedata',NS).get('{'+NS['r']+'}id')
        target = 'word/'+si['rels'][rid].lstrip('/')
        style = dict(re.findall(r'([a-z-]+):([^;]+)',node.get('style','')))
        if target not in si['media'] or 'width' not in style or 'height' not in style:
            continue
        im = Image.open(BytesIO(si['media'][target]))
        def points(value):
            m = re.fullmatch(r'\s*([0-9.]+)\s*(pt|in|cm|mm)?\s*', value)
            if not m:
                return None
            number, unit = float(m.group(1)), m.group(2) or 'pt'
            return number * {'pt': 1.0, 'in': 72.0, 'cm': 72.0/2.54, 'mm': 72.0/25.4}[unit]
        width = points(style['width'])
        height = points(style['height'])
        if width is None or height is None:
            continue
        ratios.append({'media':target,'display_ratio':width/height,'pixel_ratio':im.width/im.height,
                       'ok':abs((width/height)/(im.width/im.height)-1)<0.015})
    checks['si_image_aspect_ratios'] = len(ratios)==8 and all(r['ok'] for r in ratios)
    if args.baseline_main:
        old = load(args.baseline_main)
        old_table = old['root'].xpath('//w:tbl',namespaces=NS)[0]
        new_table = ms['root'].xpath('//w:tbl',namespaces=NS)[0]
        checks['main_table_numerical_cells_unchanged'] = old_table.xpath('.//w:t/text()',namespaces=NS)==new_table.xpath('.//w:t/text()',namespaces=NS)
        checks['main_figures_unchanged_from_baseline'] = old['media']==ms['media']
    result = {'checks':checks,'all_pass':all(checks.values()),'captions':captions,'si_aspect_ratios':ratios,
              'main':ms['path'],'si':si['path'],
              'scope':'Structural and selected cross-document checks, not proof of raw-data reconstruction or visual quality.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'all_pass':result['all_pass'],'checks':checks},indent=2))
    return 0 if result['all_pass'] else 1


if __name__=='__main__':
    raise SystemExit(main())
