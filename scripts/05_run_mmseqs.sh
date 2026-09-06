#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 4 ]]; then
  echo "Usage: $0 QUERY_FASTA COMBINED_PWM_REFERENCE_FASTA [OUTPUT_PREFIX] [THREADS]" >&2
  exit 2
fi

query_fasta=$1
reference_fasta=$2
output_prefix=${3:-TF_sequences}
threads=${4:-8}
tmp_dir="${output_prefix}_tmp"
target_prefix="${output_prefix}_reference"

command -v mmseqs >/dev/null 2>&1 || {
  echo "Error: mmseqs is not available on PATH." >&2
  exit 1
}
[[ -s "$query_fasta" ]] || { echo "Error: missing query FASTA: $query_fasta" >&2; exit 1; }
[[ -s "$reference_fasta" ]] || { echo "Error: missing reference FASTA: $reference_fasta" >&2; exit 1; }

mkdir -p "$tmp_dir"
mmseqs createdb "$query_fasta" "${output_prefix}.db"
mmseqs createdb "$reference_fasta" "${target_prefix}.db"
mmseqs search "${output_prefix}.db" "${target_prefix}.db" "${output_prefix}.ali" "$tmp_dir" \
  --threads "$threads" -s 7.5 --max-seq-id 1.0 --num-iterations 4 --alignment-mode 3
mmseqs convertalis "${output_prefix}.db" "${target_prefix}.db" "${output_prefix}.ali" \
  "${output_prefix}.comparison" \
  --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

echo "MMseqs2 alignment complete: ${output_prefix}.comparison"

