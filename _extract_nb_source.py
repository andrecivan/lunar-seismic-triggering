"""Extractor de source code de notebooks (sin outputs JSON).
Uso: python _extract_nb_source.py nb1.ipynb nb2.ipynb ...
Imprime a stdout para auditoría rápida.
Script de uso único.
"""
import nbformat
import sys
from pathlib import Path

for path in sys.argv[1:]:
    p = Path(path)
    nb = nbformat.read(p, as_version=4)
    n_code = sum(1 for c in nb.cells if c.cell_type == 'code')
    n_md   = sum(1 for c in nb.cells if c.cell_type == 'markdown')
    print(f"\n{'='*100}")
    print(f"# FILE: {p.name}   |   total_cells={len(nb.cells)}  code={n_code}  md={n_md}")
    print(f"{'='*100}")
    for i, cell in enumerate(nb.cells, 1):
        tag = "MD" if cell.cell_type == 'markdown' else "PY"
        # Conteo de líneas para ubicar funciones largas
        nlines = len(cell.source.splitlines())
        print(f"\n--- Cell {i} [{tag}] ({nlines} lines) ---")
        print(cell.source)
