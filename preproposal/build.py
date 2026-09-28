# -*- coding: utf-8 -*-
"""Rebuild word/document.xml from content.py, reusing the original's four
flowchart drawings, section properties and header references verbatim."""
import re, os, shutil, subprocess, sys
from content import C

SRC = 'unpacked/word/document.xml'
raw = open(SRC, encoding='utf-8').read()

# --- salvage the pieces we must not regenerate -------------------------------
head = raw[:raw.index('<w:body>') + len('<w:body>')]
sectpr = re.search(r'<w:sectPr.*?</w:sectPr>', raw, re.S).group(0)

# the four <w:drawing> blocks live in one paragraph, in document order
drawings = re.findall(r'<w:drawing>.*?</w:drawing>', raw, re.S)
assert len(drawings) == 4, f'expected 4 drawings, found {len(drawings)}'

ESC = {'&': '&amp;', '<': '&lt;', '>': '&gt;'}
def esc(s):
    return ''.join(ESC.get(ch, ch) for ch in s)

def runs(text):
    """Split **bold** / *italic* markup into formatted <w:r> runs."""
    out = []
    for tok in re.split(r'(\*\*[^*]+\*\*|\*[^*]+\*)', text):
        if not tok:
            continue
        if tok.startswith('**') and tok.endswith('**'):
            rpr, body = '<w:b w:val="1"/>', tok[2:-2]
        elif tok.startswith('*') and tok.endswith('*'):
            rpr, body = '<w:i w:val="1"/>', tok[1:-1]
        else:
            rpr, body = '', tok
        out.append(
            f'<w:r><w:rPr>{rpr}<w:rtl w:val="0"/></w:rPr>'
            f'<w:t xml:space="preserve">{esc(body)}</w:t></w:r>')
    return ''.join(out)

DBL = '<w:spacing w:line="480" w:lineRule="auto"/>'
SGL = '<w:spacing w:line="240" w:lineRule="auto" w:after="120"/>'

def para(text, ind='', spacing=DBL, jc='', extra='', pre=''):
    # CT_PPr enforces child order: keepNext/keepLines -> spacing -> ind -> jc -> rPr.
    # Emitting these out of sequence makes strict readers reject the part.
    return (f'<w:p><w:pPr>{extra}{spacing}{ind}{jc}<w:rPr/></w:pPr>'
            f'{pre}{runs(text)}</w:p>')

# --- table scaffolding (dual widths; DXA everywhere) -------------------------
TW, C0, C1 = 9360, 2900, 6460
def cell(text, w, bold=False, shade=False):
    sh = ('<w:shd w:val="clear" w:color="auto" w:fill="EDEDED"/>' if shade else '')
    body = f'**{text}**' if bold else text
    return (f'<w:tc><w:tcPr><w:tcW w:w="{w}" w:type="dxa"/>{sh}</w:tcPr>'
            f'<w:p><w:pPr>{SGL}<w:ind w:left="72" w:right="72" w:firstLine="0"/>'
            f'<w:rPr/></w:pPr>{runs(body)}</w:p></w:tc>')

BORDER = ('<w:tblBorders>' + ''.join(
    f'<w:{s} w:val="single" w:sz="4" w:space="0" w:color="808080"/>'
    for s in ('top', 'left', 'bottom', 'right', 'insideH', 'insideV')) +
    '</w:tblBorders>')

# --- walk the content --------------------------------------------------------
body, tbl, fignum = [], [], 0
for kind, text in C:
    if kind != 'TBLROW' and tbl:
        rows = ''.join(
            '<w:tr>' + cell(a, C0, bold=(i == 0), shade=(i == 0)) +
            cell(b, C1, bold=(i == 0), shade=(i == 0)) + '</w:tr>'
            for i, (a, b) in enumerate(tbl))
        body.append(
            f'<w:tbl><w:tblPr><w:tblW w:w="{TW}" w:type="dxa"/>{BORDER}'
            f'<w:tblLayout w:type="fixed"/></w:tblPr>'
            f'<w:tblGrid><w:gridCol w:w="{C0}"/><w:gridCol w:w="{C1}"/></w:tblGrid>'
            f'{rows}</w:tbl>'
            f'<w:p><w:pPr>{SGL}<w:rPr/></w:pPr></w:p>')
        tbl = []

    if kind == 'TITLE':
        body.append(para(f'**{text}**',
                         jc='<w:jc w:val="center"/>',
                         ind='<w:ind w:left="0" w:firstLine="0"/>'))
    elif kind in ('H1', 'H2'):
        lvl = '0' if kind == 'H1' else '1'
        sz = '28' if kind == 'H1' else '24'
        body.append(
            f'<w:p><w:pPr><w:pStyle w:val="Heading{1 if kind=="H1" else 2}"/>'
            f'<w:keepNext w:val="1"/>'
            f'<w:spacing w:before="280" w:after="120" w:line="240" w:lineRule="auto"/>'
            f'<w:ind w:left="0" w:firstLine="0"/>'
            f'<w:outlineLvl w:val="{lvl}"/>'
            f'<w:rPr><w:b w:val="1"/><w:color w:val="000000"/><w:sz w:val="{sz}"/>'
            f'<w:szCs w:val="{sz}"/></w:rPr></w:pPr>'
            f'<w:r><w:rPr><w:b w:val="1"/><w:color w:val="000000"/>'
            f'<w:sz w:val="{sz}"/><w:szCs w:val="{sz}"/><w:rtl w:val="0"/></w:rPr>'
            f'<w:t xml:space="preserve">{esc(text)}</w:t></w:r></w:p>')
    elif kind == 'P':
        body.append(para(text, ind='<w:ind w:left="0" w:firstLine="720"/>'))
    elif kind == 'PL':
        body.append(para(text, ind='<w:ind w:left="720" w:firstLine="0"/>'))
    elif kind == 'FIG':
        i = int(text)
        body.append(
            f'<w:p><w:pPr>{SGL}<w:ind w:left="0" w:firstLine="0"/>'
            f'<w:jc w:val="center"/><w:rPr/></w:pPr>'
            f'<w:r><w:rPr><w:rtl w:val="0"/></w:rPr>{drawings[i]}</w:r></w:p>')
    elif kind == 'CAP':
        body.append(para(text, spacing=SGL,
                         ind='<w:ind w:left="360" w:right="360" w:firstLine="0"/>',
                         extra='<w:keepLines w:val="1"/>'))
        body.append(f'<w:p><w:pPr>{SGL}<w:rPr/></w:pPr></w:p>')
    elif kind == 'TBLROW':
        tbl.append(tuple(text.split('|', 1)))
    elif kind == 'REF':
        body.append(para(text, spacing=DBL,
                         ind='<w:ind w:left="720" w:hanging="720"/>'))

out = head + ''.join(body) + sectpr + '</w:body></w:document>'
open(SRC, 'w', encoding='utf-8').write(out)
print(f'document.xml rewritten: {len(out):,} bytes, {out.count("<w:p>"):,} paragraphs')
