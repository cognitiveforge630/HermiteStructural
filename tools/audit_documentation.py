"""Check packaged Markdown links and literal PowerShell Python entry points.

This is a static integrity check, not a promise that expensive GUI/solver
commands have been run. Output paths under outputs/ are generated, not shipped.
Run: python -m tools.audit_documentation
"""
from __future__ import annotations

import importlib.util
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MARKDOWN_LINK = re.compile(r'!?\[[^\]]*\]\(([^)]+)\)')
FENCE = re.compile(r'^```(?:powershell|pwsh|shell|bash)\s*$')
PYTHON_COMMAND = re.compile(r'^\s*(?:python|py)(?:\.exe)?\s+(.*)$', re.IGNORECASE)


def markdown_files():
    return sorted(path for path in ROOT.rglob('*.md')
                  if 'outputs' not in path.parts)


def check_docs():
    failures=[]
    link_total=0
    command_total=0
    for file in markdown_files():
        txt=file.read_text(encoding='utf-8')
        for match in MARKDOWN_LINK.finditer(txt):
            dest=match.group(1).split(' ',1)[0].split('#',1)[0]
            if not dest or re.match(r'^\w+://|^mailto:',dest):
                continue
            link_total+=1
            if not (file.parent / dest.replace('%20',' ')).exists():
                failures.append(f'{file.relative_to(ROOT)}: broken link: {dest}')
        in_code=False
        for line_num, line in enumerate(txt.splitlines(),start=1):
            if FENCE.match(line):
                in_code=True
                continue
            if in_code and line.startswith('```'):
                in_code=False
                continue
            if not in_code or line.lstrip().startswith('#'):
                continue
            match=PYTHON_COMMAND.match(line)
            if not match:continue
            args=match.group(1).strip().split()
            if not args:continue
            if args[0]=='-m' and len(args)>=2:
                mod=args[1]
                if mod=='pip' or mod=='unittest':continue
                command_total+=1
                if '<' in mod or '{' in mod:continue
                if not (ROOT / (mod.replace('.','/')+'.py')).is_file() and not (ROOT/mod.replace('.','/')/'__main__.py').is_file():
                    failures.append(f'{file.relative_to(ROOT)}:{line_num}: nonexistent -m module: {mod}')
            elif args[0].endswith('.py'):
                command_total+=1
                target=args[0].replace('\\','/').removeprefix('./')
                if not (ROOT/target).is_file():
                    failures.append(f'{file.relative_to(ROOT)}:{line_num}: nonexistent Python script: {target}')
    return {'markdown_files':len(markdown_files()),'checked_internal_links':link_total,
            'checked_python_commands':command_total,'failures':failures}


def main():
    result=check_docs()
    for key,value in result.items():
        if key!='failures':print(f'{key}: {value}')
    for failure in result['failures']:
        print('FAIL:',failure)
    if result['failures']:
        raise SystemExit(f"FAIL: {len(result['failures'])} documentation errors")
    print('PASS: all resolvable Markdown links and literal Python entrypoints exist')


if __name__ == '__main__':
    main()
