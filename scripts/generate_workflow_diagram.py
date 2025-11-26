#!/usr/bin/env python3
"""Generate workflow diagram"""
import subprocess

# Generate DOT
result = subprocess.run(['snakemake', '--dag'], capture_output=True, text=True)

# Save text version
with open('docs/workflow_dag.txt', 'w') as f:
    f.write(result.stdout)
print("✓ Saved: docs/workflow_dag.txt")

# Try PNG
try:
    subprocess.run(['dot', '-Tpng'], input=result.stdout.encode(),
                   stdout=open('docs/workflow_diagram.png', 'wb'), check=True)
    print("✓ Saved: docs/workflow_diagram.png")
except FileNotFoundError:
    print("⚠ Install graphviz: sudo apt install graphviz")