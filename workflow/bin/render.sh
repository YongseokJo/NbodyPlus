#!/usr/bin/env bash
set -euo pipefail

# Robust template renderer using Python string replacement.
# Placeholders look like: {{NAME}}
# Usage:
#   render.sh template.in out NAME value NAME2 value2 ...

TEMPLATE=${1:?template}
OUT=${2:?out}
shift 2

if [[ $(( $# % 2 )) -ne 0 ]]; then
  echo "render.sh: expected key/value pairs" >&2
  exit 2
fi

# Build a small Python one-liner to perform literal replacements.
python3 - "$@" <<PY
import sys
from pathlib import Path

tmpl = Path(r"$TEMPLATE").read_text()
args = sys.argv[1:]
if len(args) % 2 != 0:
  print('render: expected key/value pairs', file=sys.stderr)
  sys.exit(2)
for i in range(0, len(args), 2):
  k = args[i]
  v = args[i+1]
  key = '{{' + k + '}}'
  tmpl = tmpl.replace(key, v)

Path(r"$OUT").write_text(tmpl)
PY
