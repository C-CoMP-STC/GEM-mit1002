#!/usr/bin/env bash
# Build the continuous-curation manifesto as a PDF supplement.
#
#   bash docs/build.sh
#
# Citations in continuous-curation.md are pandoc keys like [@wilson2017good].
# They render literally on GitHub and are resolved only here. The bibliography
# is docs/references.bib; the citation style is docs/asm.csl (ASM numbered).
#
# Requires pandoc and a LaTeX engine. On macOS:
#   brew install pandoc
#   brew install --cask basictex     # then: sudo tlmgr install lmodern
#
# Output (docs/continuous-curation.pdf) is a build artifact -- add it to
# .gitignore unless you decide to commit it alongside the other generated files.

set -euo pipefail

cd "$(dirname "$0")"

SRC=continuous-curation.md
OUT=continuous-curation.pdf

command -v pandoc >/dev/null || { echo "error: pandoc not found" >&2; exit 1; }

# --citeproc is built in from pandoc 2.11; older builds need a separate filter.
PANDOC_MAJOR=$(pandoc --version | head -1 | sed -E 's/^pandoc ([0-9]+)\.([0-9]+).*/\1\2/')
if [ "${PANDOC_MAJOR:-0}" -lt 211 ] 2>/dev/null; then
    echo "error: pandoc >= 2.11 required for --citeproc; found $(pandoc --version | head -1)" >&2
    echo "       macOS: brew install pandoc" >&2
    exit 1
fi

# Warn about placeholder author lists rather than silently shipping them.
if grep -vn '^%' references.bib | grep -q 'TODO-AUTHORS'; then
    echo "warning: references.bib still has TODO-AUTHORS entries:" >&2
    grep -n 'TODO-AUTHORS' references.bib | grep -v ':%' | sed 's/^/  /' >&2
fi

pandoc "$SRC" \
    --citeproc \
    --bibliography=references.bib \
    --csl=asm.csl \
    --metadata-file=metadata.yaml \
    --lua-filter=promote-headings.lua \
    --include-in-header=latex-header.tex \
    --toc \
    --number-sections \
    -V geometry:margin=1in \
    -V linkcolor:blue \
    -o "$OUT"

echo "wrote docs/$OUT"
