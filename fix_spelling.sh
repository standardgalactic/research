#!/usr/bin/env bash
set -euo pipefail
DIR="${1:-.}"
python3 - "$DIR" <<'PY'
from pathlib import Path
import re, sys, shutil, difflib
folder=Path(sys.argv[1]); stem='Why_AI_Safety_Branding_Accelerates_Risk'
# Only transcription/spelling/obvious grammar repairs. No fact-checking or claim rewrites.
fixes=[
 ('Evidenciary Authority','Evidentiary Authority'),
 ('Flaksian','Flyxion'),
 ('can sport 100 tons','can support 100 tons'),
 ('claimgraphcheck.py','claimgraph_check.py'),
 ('und, unresult','und, unresolved'),
 ('to higher political representation','to hire political representation'),
 ('antipomorphic mediation','anthropomorphic mediation'),
 ('paired extraction protocol','paired abstraction protocol'),
 ('the cultural reach of political resonance','the cultural reach or political resonance'),
 ('some consensus views that\'s scraped','some consensus views it scraped'),
 ('my gmails','my Gmail messages'),
 ('S1 prospectus','S-1 prospectus'),
 ('an S1','an S-1'),
 ('the S1','the S-1'),
 ('S1 risk factors','S-1 risk factors'),
]
# "misanthropic statement" is ambiguous: could be a deliberate coined name.
# Leave it unchanged rather than guessing. Similarly, leave all reported
# figures, dates, litigation, and IPO assertions unchanged.
for ext in ('.txt','.vtt'):
    path=folder/(stem+ext)
    if not path.is_file(): raise SystemExit(f'Missing file: {path}')
    original=path.read_text(encoding='utf-8-sig')
    updated=original
    counts=[]
    for old,new in fixes:
        # Allow a phrase to wrap over lines in subtitle cues.
        pat=r'(?<!\w)'+r'\s+'.join(map(re.escape,old.split()))+r'(?!\w)'
        updated,n=re.subn(pat,new,updated,flags=re.IGNORECASE)
        if n: counts.append((old,new,n))
    if ext=='.vtt':
        times=lambda s: re.findall(r'^\s*(?:\d\d:)?\d\d:\d\d\.\d{3}\s*-->[^\n]*$',s,re.M)
        if times(original)!=times(updated): raise SystemExit('VTT timing changed: refusing write')
    if updated!=original:
        backup=path.with_name(path.name+'.bak')
        if not backup.exists(): shutil.copy2(path,backup)
        path.write_text(updated,encoding='utf-8')
    print(f'{path.name}: {sum(n for _,_,n in counts)} replacements')
    for old,new,n in counts: print(f'  {n} × {old!r} → {new!r}')
print('Done. Original files backed up as .bak (if changed).')
PY
