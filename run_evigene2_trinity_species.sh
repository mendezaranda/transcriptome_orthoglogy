#!/bin/bash
#
#SBATCH --job-name=evigene_trinity
#SBATCH -N1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --gres=localtmp:80G
#SBATCH -t 12:00:00
#SBATCH -o /home/dmendez/outdir/%j_%x
#SBATCH -e /home/dmendez/errdir/%j_%x
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.mendezaranda@mdc-berlin.de



# ================================================================
# De novo transcriptome cleanup pipeline for Orthofinder
# Daniel Mendez Aranda, 2025-12
# ================================================================

set -eo pipefail

########################
# USAGE
########################
# sbatch run_evigene_trinity_species.sh CryNat /path/to/CryNat_merged.prefixed.fa

species="${1:-}"
input_fasta="${2:-}"

BASE_OUTDIR="/fast/AG_Lewin/dmendez/transcriptomes/AMR/new"
ORTHOFINDER_DIR="/fast/AG_Lewin/dmendez/orthofinder/20251219/"
NCPU="${NCPU:-16}"
CDHIT_ID="${CDHIT_ID:-0.97}"

# Optional knobs
MIN_AA="${MIN_AA:-100}"                 # select protein size minimum
COMPLETE_ONLY="${COMPLETE_ONLY:-1}"   # 1 = keep only ORF type:complete

export EVIGENEHOME="${EVIGENEHOME:-/fast/AG_Lewin/dmendez/tools/evidentialgene}"

if [[ -z "$species" || -z "$input_fasta" ]]; then
  echo "Usage: $0 <species> <merged_trinity_fasta>"
  exit 1
fi
if [[ ! -f "$input_fasta" ]]; then
  echo "ERROR: input fasta not found: $input_fasta"
  exit 1
fi

source /home/dmendez/.bashrc
conda activate /fast/AG_Lewin/dmendez/.conda/envs/Trinity
module load mpi || true

if [[ -z "${EVIGENEHOME:-}" || ! -x "$EVIGENEHOME/scripts/prot/tr2aacds.pl" ]]; then
  echo "ERROR: Could not locate Evigene (tr2aacds.pl). Set EVIGENEHOME explicitly."
  exit 1
fi

outdir="${BASE_OUTDIR}/${species}"
mkdir -p "${outdir}"
cd "${outdir}"

merged_fa="${outdir}/${species}_merged.fa"
cdhit_fa="${outdir}/${species}_cdhit.uniq.fa"
okay_mrna="${outdir}/${species}_okay.mrna"
final_tx="${outdir}/${species}_okay.mrna"
final_pep="${final_tx}.transdecoder.pep"

# NOTE: we keep an unprefixed collapsed file, then prefix once for OrthoFinder
collapsed_pep="${outdir}/${species}.primary_by_locus.pep.fa"
orthofinder_pep="${ORTHOFINDER_DIR}/${species}.pep.fasta"

count_fa () { grep -c '^>' "$1" 2>/dev/null || true; }

echo ">>> Starting Evigene+cleanup pipeline for ${species}"
date
echo "Input Trinity fasta: ${input_fasta}"
echo "Output dir:          ${outdir}"
echo

########################
# STEP 1: LINK MERGED FASTA
########################
if [[ ! -s "${merged_fa}" ]]; then
  ln -sf "${input_fasta}" "${merged_fa}"
fi

########################
# STEP 2: CD-HIT-EST
########################
if [[ ! -s "${cdhit_fa}" ]]; then
  echo ">>> Running cd-hit-est on ${merged_fa}"
  cd-hit-est -i "${merged_fa}" -o "${cdhit_fa}" \
    -c "${CDHIT_ID}" -aS 0.90 -g 1 -T "${NCPU}" -M 0 2>&1 | tee cdhit.log
else
  echo ">>> Skipping cd-hit-est, ${cdhit_fa} already exists"
fi

echo ">>> Checking for duplicate FASTA IDs (first token) in ${cdhit_fa}"
dups=$(grep '^>' "${cdhit_fa}" | awk '{print $1}' | sort | uniq -d | head -n 5 || true)
if [[ -n "${dups}" ]]; then
  echo "ERROR: duplicate IDs detected after cd-hit-est. Examples:"
  echo "${dups}"
  echo "Fix header collisions (prefix each Trinity run BEFORE merge) and rerun."
  exit 2
fi

########################
# STEP 3: Evigene tr2aacds
########################
base="${species}_cdhit"
evigene_work="${outdir}/${base}_evigene"
mkdir -p "${evigene_work}"
cd "${evigene_work}"

if [[ ! -s "okayset/${base}.okay.mrna" && ! -s "${okay_mrna}" ]]; then
  echo ">>> Running Evigene tr2aacds on ${cdhit_fa}"
  ln -sf "${cdhit_fa}" "${base}.mrna"
  "$EVIGENEHOME/scripts/prot/tr2aacds.pl" \
    -mrna "${base}.mrna" \
    -NCPU "${NCPU}" \
    -MINCDS=90 2>&1 | tee "${outdir}/tr2aacds.log"

  ok_mrna_path="okayset/${base}.okay.mrna"
  if [[ -s "${ok_mrna_path}" ]]; then
    cp "${ok_mrna_path}" "${okay_mrna}"
  else
    echo "WARNING: Evigene produced no okay.mrna; falling back to cd-hit transcripts."
    cp "${cdhit_fa}" "${okay_mrna}"
  fi
else
  echo ">>> Skipping Evigene, ${okay_mrna} already exists or okayset present"
fi

cd "${outdir}"

########################
# STEP 4: TransDecoder
########################
if [[ ! -s "${final_pep}" ]]; then
  echo ">>> Predicting proteins with TransDecoder..."
  TransDecoder.LongOrfs -t "${final_tx}" 2>&1 | tee final_transdecoder_longorfs.log
  TransDecoder.Predict  -t "${final_tx}" 2>&1 | tee final_transdecoder_predict.log
else
  echo ">>> Skipping TransDecoder, ${final_pep} already exists"
fi

########################
# STEP 5: Collapse to ONE protein per Evigene locus (strip t1/t2/...)
########################
if [[ ! -s "${collapsed_pep}" ]]; then
  echo ">>> Collapsing proteins per Evigene locus (t1/t2/... collapse)..."
  python3 - <<PY
import re
from pathlib import Path

inp = Path("${final_pep}")
out = Path("${collapsed_pep}")
min_aa = int("${MIN_AA}")
complete_only = int("${COMPLETE_ONLY}")

hdr_re = re.compile(
    r'^>(\\S+)\\s+GENE\\.([^~\\s]+)~~([^\\s]+)\\s+ORF.*?score=([0-9.]+)\\s+len:([0-9]+)'
)

def fasta_iter(path):
    h = None
    seq = []
    with path.open() as fh:
        for line in fh:
            line=line.rstrip("\\n")
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(seq)
                h=line
                seq=[]
            else:
                seq.append(line.strip())
        if h is not None:
            yield h, "".join(seq)

best = {}  # locus -> (score, length, locus, seq)
n_total = 0
n_used = 0
n_complete = 0
n_matched = 0

for header, seq in fasta_iter(inp):
    n_total += 1
    seq = seq.replace("*","")
    if min_aa and len(seq) < min_aa:
        continue

    m = hdr_re.match(header)
    if not m:
        # Without rich headers we cannot do t-collapse reliably -> keep unique first token
        locus = header[1:].split()[0]
        score = 0.0
        L = len(seq)
    else:
        n_matched += 1
        full_id, gene_id, orf_id, score_str, len_str = m.groups()
        score = float(score_str)
        L = int(len_str)

        if complete_only:
            mt = re.search(r'ORF type:([^ ]+)', header)
            if not mt:
                continue
            if mt.group(1) == "complete":
                n_complete += 1
            else:
                continue

        # collapse Evigene transcript variants: ...t1, ...t2, ... -> locus
        locus = re.sub(r't\\d+$', '', gene_id)

    n_used += 1
    prev = best.get(locus)
    if prev is None or (score, L) > (prev[0], prev[1]):
        best[locus] = (score, L, locus, seq)

with out.open("w") as oh:
    for locus, (score, L, locus, seq) in best.items():
        oh.write(f">{locus}\\n")
        for i in range(0, len(seq), 60):
            oh.write(seq[i:i+60] + "\\n")

print("Input proteins:", n_total)
print("Used (after min_aa/complete filters):", n_used)
print("Matched rich TD headers:", n_matched)
if complete_only:
    print("Kept complete ORFs:", n_complete)
print("Primary loci written:", len(best))
print("Wrote:", out)
PY
else
  echo ">>> Skipping collapse, ${collapsed_pep} already exists"
fi

########################
# STEP 6: Prefix + copy to OrthoFinder (single-token headers)
########################
tmp_prefixed="${outdir}/${species}.orthofinder.faa"
awk -v sp="${species}" '
  /^>/ {
    id = substr($0,2);
    sub(/[[:space:]].*$/, "", id);
    print ">" sp "_" id;
    next
  }
  { print }
' "${collapsed_pep}" > "${tmp_prefixed}"

cp "${tmp_prefixed}" "${orthofinder_pep}"

########################
# REPORT
########################
echo ">>> Counts:"
echo "merged      : $(count_fa "${merged_fa}")"
echo "cdhit       : $(count_fa "${cdhit_fa}")"
echo "okay_mrna   : $(count_fa "${okay_mrna}")"
echo "transdecoder: $(count_fa "${final_pep}")"
echo "collapsed   : $(count_fa "${collapsed_pep}")"
echo "orthofinder : $(count_fa "${orthofinder_pep}")"

if [[ "$(count_fa "${collapsed_pep}")" -ge "$(count_fa "${final_pep}")" ]]; then
  echo "WARNING: collapsed is not smaller than TransDecoder pep; check that headers contain 'GENE.' and Evigene IDs."
fi

echo ">>> OrthoFinder proteome written to: ${orthofinder_pep}"

conda deactivate || true
echo ">>> Completed Evigene+cleanup pipeline for ${species}"
date
