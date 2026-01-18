#!/bin/bash
#
#SBATCH --job-name=bpp_prep
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --gres=localtmp:40G
#SBATCH --mem=32G
#SBATCH -t 6:00:00
#SBATCH -o /home/dmendez/outdir/%j_%x
#SBATCH -e /home/dmendez/errdir/%j_%x

set -eo pipefail

# ================================================================
# Prepare input for BPP from Orthofinder single-copy orthologues
# Daniel Mendez Aranda, 2026.01
# ================================================================

# -----------------------
# Inputs
# -----------------------

MODE="${MODE:-array}"  # array | merge

RESULTS_DIR="${RESULTS_DIR:-/fast/AG_Lewin/dmendez/orthofinder/20251219/primary/input_clean/OrthoFinder/Results_Jan09}"
IN_DIR="${IN_DIR:-$RESULTS_DIR/bpp_from_SCO/aln}"
WORK="${WORK:-$RESULTS_DIR/bpp_from_SCO}"

TRIMAL="${TRIMAL:-/fast/AG_Lewin/dmendez/.conda/envs/funannotate/bin/trimal}"
MAFFT_THREADS="${MAFFT_THREADS:-1}"
CPUS="${SLURM_CPUS_PER_TASK:-16}"

ALN_DIR="$WORK/aln"
TRIM_DIR="$WORK/aln_trim"
PHY_DIR="$WORK/phy_blocks"
FINAL_PHY="$WORK/bpp_loci.phy"
IMAP="$WORK/imap.txt"
SKIPLOG="$WORK/skipped_loci.txt"

mkdir -p "$ALN_DIR" "$TRIM_DIR" "$PHY_DIR"
: > "$SKIPLOG" 2>/dev/null || true

# ----------------- merge mode -----------------
if [[ "$MODE" == "merge" ]]; then
  ls "$PHY_DIR"/*.phy >/dev/null 2>&1 || { echo "No PHYLIP blocks in $PHY_DIR"; exit 2; }

  cat "$PHY_DIR"/*.phy > "$FINAL_PHY"

  python3 - <<PY
import re
seqfile = "${FINAL_PHY}"
imap_out = "${IMAP}"

pairs=set()
with open(seqfile) as f:
    for line in f:
        line=line.strip()
        if not line:
            continue
        if re.match(r'^\\d+\\s+\\d+$', line):
            continue
        name=line.split()[0]
        sp=name.split("_",1)[0]
        pairs.add((name,sp))

with open(imap_out,"w") as out:
    for name,sp in sorted(pairs):
        out.write(f"{name} {sp}\\n")

nloci = 0
with open(seqfile) as f:
    for line in f:
        if re.match(r'^\\d+\\s+\\d+\\s*$', line):
            nloci += 1

print("Wrote:", seqfile)
print("Loci blocks:", nloci)
print("Wrote:", imap_out, "mappings:", len(pairs))
PY
  exit 0
fi

# ----------------- array mode -----------------

# Absolute path check for trimal
if [[ ! -x "$TRIMAL" ]]; then
  echo "trimal not executable: $TRIMAL"
  exit 2
fi

mapfile -t FILES < <(ls -1 "$IN_DIR"/OG*.fa 2>/dev/null | sort)
nfiles="${#FILES[@]}"
if [[ "$nfiles" -eq 0 ]]; then
  echo "No OG*.fa files found in: $IN_DIR"
  exit 2
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
TASK_COUNT="${SLURM_ARRAY_TASK_COUNT:-1}"

SHARD_LIST="$WORK/shard_${TASK_ID}.list"
: > "$SHARD_LIST"
for i in "${!FILES[@]}"; do
  t=$(( (i % TASK_COUNT) + 1 ))
  if [[ "$t" -eq "$TASK_ID" ]]; then
    echo "${FILES[$i]}" >> "$SHARD_LIST"
  fi
done

process_one () {
  local f="$1"
  local aln_dir="$2"
  local trim_dir="$3"
  local phy_dir="$4"
  local trimal="$5"
  local mafft_threads="$6"
  local skiplog="$7"

  local base
  base=$(basename "$f" .fa)

  local aln="$aln_dir/${base}.aln.fa"
  local trim="$trim_dir/${base}.trim.fa"
  local phy="$phy_dir/${base}.phy"

  # MAFFT: skip if alignment exists
  if [[ ! -s "$aln" ]]; then
    mafft --auto --thread "$mafft_threads" "$f" > "$aln"
  fi

  # trimAl: always regenerate trim
  "$trimal" -automated1 -in "$aln" -out "$trim" >/dev/null 2>&1 || cp "$aln" "$trim"

  # if trim not produced, log + skip
  if [[ ! -s "$trim" ]]; then
    echo -e "${base}\tNO_TRIM_FILE" >> "$skiplog"
    return 0
  fi

  # FASTA -> PHYLIP block
  python3 - <<PY
import sys

inp="${trim}"
outp="${phy}"
base="${base}"
skiplog="${skiplog}"

def read_fasta(path):
    seqs={}
    name=None
    parts=[]
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name]="".join(parts)
                name=line[1:].split()[0]
                parts=[]
            else:
                parts.append(line.replace(" ", ""))
        if name is not None:
            seqs[name]="".join(parts)
    return seqs

seqs=read_fasta(inp)
names=list(seqs.keys())
if not names:
    with open(skiplog,"a") as s:
        s.write(base + "\\tEMPTY\\n")
    sys.exit(0)

lens={len(seqs[n]) for n in names}
if len(lens)!=1:
    with open(skiplog,"a") as s:
        s.write(base + "\\tINCONSISTENT_LENGTHS\\n")
    sys.exit(0)

nseq=len(names)
nsite=lens.pop()

with open(outp,"w") as out:
    out.write(f"{nseq} {nsite}\\n")
    for n in names:
        out.write(f"{n} {seqs[n]}\\n")
    out.write("\\n")
PY
}

export -f process_one

# Pass everything as args (no exported env → avoids "environment: line 5" issues)
while IFS= read -r f; do
  printf '%s\0' "$f"
done < "$SHARD_LIST" | xargs -0 -P "$CPUS" -I{} bash -lc \
  'process_one "$1" "$2" "$3" "$4" "$5" "$6" "$7"' _ {} \
  "$ALN_DIR" "$TRIM_DIR" "$PHY_DIR" "$TRIMAL" "$MAFFT_THREADS" "$SKIPLOG"
