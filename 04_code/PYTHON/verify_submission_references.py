"""Validate numbered DOCX citations, bibliography order and bookmark linkage.

Run with a manuscript path and an explicit report path. This tool does not
verify source authenticity; retain a separate DOI/source evidence ledger.
"""
from __future__ import annotations

import argparse
import json
import re
import zipfile
from pathlib import Path

from lxml import etree as ET

W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
NS = {'w': W}
Q = lambda name: '{' + W + '}' + name


def text(node):
    return ''.join(node.xpath('.//w:t/text()', namespaces=NS))


def compress(numbers):
    numbers = sorted(set(numbers))
    groups = []
    for n in numbers:
        if groups and n == groups[-1][-1] + 1:
            groups[-1].append(n)
        else:
            groups.append([n])
    return ','.join(str(g[0]) if len(g) == 1 else f'{g[0]}-{g[-1]}' for g in groups)


def audit(path):
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read('word/document.xml'))
        numbering = ET.fromstring(z.read('word/numbering.xml'))
        zip_error = z.testzip()
    body = root.find(Q('body'))
    paragraphs = body.findall(Q('p'))
    ref_idx = next(i for i,p in enumerate(paragraphs) if text(p).strip() == 'References')
    refs = []
    for p in paragraphs[ref_idx+1:]:
        if text(p).strip() == 'Figure Legends':
            break
        anchors = p.xpath('./w:bookmarkStart', namespaces=NS)
        anchors = [a for a in anchors if re.fullmatch(r'R\d+', a.get(Q('name'),''))]
        if anchors:
            refs.append((int(anchors[0].get(Q('name'))[1:]), p))
        elif p.find('./'+Q('pPr')+'/'+Q('numPr')) is not None:
            refs.append((-1, p))
    expected = set(range(1,len(refs)+1))
    starts = root.xpath('//w:bookmarkStart', namespaces=NS)
    ends = root.xpath('//w:bookmarkEnd', namespaces=NS)
    names = [s.get(Q('name')) for s in starts]
    start_ids = [s.get(Q('id')) for s in starts]
    end_ids = [s.get(Q('id')) for s in ends]
    groups = []
    issues = []
    all_super = True
    link_match = True
    citation_grammar = True
    ordered_groups = True
    merged = True
    all_cited = set()
    first_seen = []
    for pidx, p in enumerate(paragraphs[:ref_idx],1):
        elems = list(p)
        i = 0
        while i < len(elems):
            node = elems[i]
            if node.tag != Q('hyperlink') or not re.fullmatch(r'R\d+',node.get(Q('anchor'),'')):
                i += 1
                continue
            members = []
            while i < len(elems):
                e = elems[i]
                is_link = e.tag == Q('hyperlink') and re.fullmatch(r'R\d+', e.get(Q('anchor'),''))
                is_sep = e.tag == Q('r') and re.fullmatch(r'[,\-–]+',text(e))
                if not (is_link or is_sep):
                    break
                members.append(e); i += 1
            raw = ''.join(text(e) for e in members)
            valid = bool(re.fullmatch(r'\d+(?:[-–]\d+)?(?:,\d+(?:[-–]\d+)?)*',raw))
            citation_grammar &= valid
            nums = []
            if valid:
                for chunk in raw.split(','):
                    bounds = list(map(int,re.split('[-–]',chunk)))
                    if len(bounds) == 2:
                        if bounds[0] >= bounds[1]:
                            ordered_groups = False
                        nums.extend(range(bounds[0],bounds[1]+1))
                    else:
                        nums += bounds
                ordered_groups &= nums == sorted(set(nums))
                merged &= raw == compress(nums)
                for n in nums:
                    if n not in all_cited:
                        first_seen.append(n)
                    all_cited.add(n)
            for e in members:
                for r in ([e] if e.tag == Q('r') else e.findall(Q('r'))):
                    va = r.find('./'+Q('rPr')+'/'+Q('vertAlign'))
                    all_super &= va is not None and va.get(Q('val')) == 'superscript'
                if e.tag == Q('hyperlink'):
                    link_match &= e.get(Q('anchor')) == 'R'+text(e) and int(text(e)) in expected
            groups.append({'paragraph':pidx,'display':raw,'expanded':nums,'context':text(p)})
    body_text = '\n'.join(text(p) for p in paragraphs[:ref_idx])
    author_year_left = re.findall(r'(?:et al\.,?\s*\(?|VanderWeele and Knol,\s*)20\d\d', body_text)
    nums = []
    list_format = True
    for _,p in refs:
        n = p.find('./'+Q('pPr')+'/'+Q('numPr')+'/'+Q('numId'))
        list_format &= n is not None
        if n is not None:
            nums.append(n.get(Q('val')))
    for num_id in set(nums):
        actual = numbering.xpath(f'./w:num[@w:numId="{num_id}"]/w:abstractNumId/@w:val',namespaces=NS)
        list_format &= len(actual)==1
        if actual:
            fmt = numbering.xpath(f'./w:abstractNum[@w:abstractNumId="{actual[0]}"]/w:lvl[@w:ilvl="0"]/w:numFmt/@w:val',namespaces=NS)
            list_format &= fmt == ['decimal']
    checks = {
        'bibliography_present_and_all_cited': bool(refs) and all_cited == expected,
        'all_citation_numbers_and_separators_superscript': bool(groups) and all_super,
        'no_missing_or_trailing_separators': citation_grammar,
        'no_descending_ranges_or_unsorted_groups': ordered_groups,
        'first_appearance_strictly_increasing': first_seen == list(range(1,len(refs)+1)),
        'bookmark_anchor_and_visible_number_match': link_match and len(names)==len(set(names)) and sorted(start_ids)==sorted(end_ids),
        'bibliography_order_and_automatic_numbering': [n for n,_ in refs] == list(range(1,len(refs)+1)) and list_format,
        'all_consecutive_numbers_compressed_including_pairs': merged,
        'no_residual_author_year_citation_syntax': not author_year_left,
        'docx_zip_valid': zip_error is None,
    }
    return {'document':str(path),'reference_count':len(refs),'citation_groups':len(groups),
            'first_appearance':first_seen,'unused':sorted(expected-all_cited),
            'missing':sorted(all_cited-expected),'checks':checks,'all_pass':all(checks.values()),
            'source_authenticity':'Not tested by this structural verifier; see separate live source audit.',
            'occurrences':groups}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manuscript',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    result = audit(args.manuscript)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='occurrences'},ensure_ascii=True,indent=2))
    return 0 if result['all_pass'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
