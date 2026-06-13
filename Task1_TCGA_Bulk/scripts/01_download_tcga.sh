#!/usr/bin/env bash
# =============================================================================
# 01_download_tcga.sh — Download TCGA-LUAD data via GDC API + gdc-client
#
# Downloads:
#   1. RNA-seq STAR Counts  (~601 files, ~2.5 GB)
#   2. Somatic Mutation MAF (~618 files, ~60 MB)
#   3. Clinical metadata    (GDC API → local TSV, no bulk file download)
#
# Usage (from Task1_TCGA_Bulk/):
#   bash scripts/01_download_tcga.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${PROJECT_DIR}/data"
MANIFEST_DIR="${DATA_DIR}/manifests"
LOG_FILE="${DATA_DIR}/download.log"

export PATH="${HOME}/.local/bin:${HOME}/.aspera/connect/bin:${PATH}"

mkdir -p "${MANIFEST_DIR}" \
         "${DATA_DIR}/star_counts" \
         "${DATA_DIR}/somatic_maf" \
         "${DATA_DIR}/clinical"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_FILE}"
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || { log "ERROR: '$1' not found in PATH"; exit 1; }
}

require_cmd python3
require_cmd gdc-client

log "=== TCGA-LUAD download started ==="
log "Project dir: ${PROJECT_DIR}"
cd "${PROJECT_DIR}"
export DATA_DIR

# ---------------------------------------------------------------------------
# Step 1: Query GDC API → generate manifests + clinical TSV
# ---------------------------------------------------------------------------
python3 - <<'PY'
import csv
import json
import os
import sys
import urllib.error
import urllib.request

PROJECT = "TCGA-LUAD"
BASE = "https://api.gdc.cancer.gov"
DATA_DIR = os.environ.get("DATA_DIR", "data")
MANIFEST_DIR = os.path.join(DATA_DIR, "manifests")
os.makedirs(MANIFEST_DIR, exist_ok=True)


def post(endpoint, payload):
    url = f"{BASE}/{endpoint}"
    data = json.dumps(payload).encode()
    req = urllib.request.Request(
        url, data=data, headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.load(resp)


def project_filter():
    return {"op": "=", "content": {"field": "cases.project.project_id", "value": PROJECT}}


def fetch_all_files(extra_filters, fields):
    filters = {"op": "and", "content": [project_filter()] + extra_filters}
    hits = []
    offset = 0
    page = 500
    while True:
        payload = {
            "filters": filters,
            "fields": fields,
            "format": "JSON",
            "size": page,
            "from": offset,
        }
        resp = post("files", payload)
        batch = resp["data"]["hits"]
        if not batch:
            break
        hits.extend(batch)
        offset += len(batch)
        total = resp["data"]["pagination"]["total"]
        print(f"  fetched {len(hits)}/{total} files...", flush=True)
        if len(hits) >= total:
            break
    return hits


def write_manifest(hits, out_path):
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["id", "filename", "md5", "size", "state"])
        for h in hits:
            writer.writerow([
                h["file_id"],
                h["file_name"],
                h.get("md5sum", ""),
                h.get("file_size", 0),
                "live",
            ])
    total_mb = sum(h.get("file_size", 0) or 0 for h in hits) / 1e6
    print(f"  manifest: {out_path}  ({len(hits)} files, {total_mb:.0f} MB)")


def first_record(obj):
    """GDC may return nested fields as dict or list."""
    if not obj:
        return {}
    if isinstance(obj, list):
        return obj[0] if obj else {}
    if isinstance(obj, dict):
        return obj
    return {}


# --- Manifests ---
print("\n[1/3] STAR Counts manifest")
star_hits = fetch_all_files(
    [
        {"op": "=", "content": {"field": "files.data_category", "value": "Transcriptome Profiling"}},
        {"op": "=", "content": {"field": "files.data_type", "value": "Gene Expression Quantification"}},
        {"op": "=", "content": {"field": "files.analysis.workflow_type", "value": "STAR - Counts"}},
    ],
    "file_id,file_name,md5sum,file_size",
)
write_manifest(star_hits, os.path.join(MANIFEST_DIR, "luad_star_counts.manifest.tsv"))

print("\n[2/3] Somatic MAF manifest")
maf_hits = fetch_all_files(
    [
        {"op": "=", "content": {"field": "files.data_category", "value": "Simple Nucleotide Variation"}},
        {"op": "=", "content": {"field": "files.data_type", "value": "Masked Somatic Mutation"}},
    ],
    "file_id,file_name,md5sum,file_size",
)
write_manifest(maf_hits, os.path.join(MANIFEST_DIR, "luad_somatic_maf.manifest.tsv"))

# --- Clinical via API ---
print("\n[3/3] Clinical metadata via GDC cases API")
clinical_out = os.path.join(DATA_DIR, "clinical", "tcga_luad_clinical.tsv")
fields = (
    "submitter_id,"
    "samples.submitter_id,samples.sample_type,samples.sample_type_id,"
    "demographic.vital_status,demographic.days_to_birth,"
    "diagnoses.vital_status,diagnoses.days_to_death,diagnoses.days_to_last_follow_up,"
    "diagnoses.primary_diagnosis,diagnoses.tumor_stage"
)
case_filter = {"op": "=", "content": {"field": "project.project_id", "value": PROJECT}}
rows = []
offset = 0
page = 500
while True:
    payload = {
        "filters": case_filter,
        "fields": fields,
        "format": "JSON",
        "size": page,
        "from": offset,
    }
    resp = post("cases", payload)
    batch = resp["data"]["hits"]
    if not batch:
        break
    for case in batch:
        case_id = case.get("submitter_id", "")
        demo = first_record(case.get("demographic"))
        diag = first_record(case.get("diagnoses"))
        vital = demo.get("vital_status") or diag.get("vital_status") or ""
        days_death = diag.get("days_to_death")
        days_follow = diag.get("days_to_last_follow_up")
        for sample in case.get("samples") or []:
            rows.append({
                "case_submitter_id": case_id,
                "sample_submitter_id": sample.get("submitter_id", ""),
                "sample_type": sample.get("sample_type", ""),
                "sample_type_id": sample.get("sample_type_id", ""),
                "vital_status": vital,
                "days_to_death": days_death if days_death is not None else "",
                "days_to_last_follow_up": days_follow if days_follow is not None else "",
                "primary_diagnosis": diag.get("primary_diagnosis", ""),
                "tumor_stage": diag.get("tumor_stage", ""),
            })
    offset += len(batch)
    total = resp["data"]["pagination"]["total"]
    print(f"  fetched {offset}/{total} cases...", flush=True)
    if offset >= total:
        break

if not rows:
    print("ERROR: no clinical rows retrieved", file=sys.stderr)
    sys.exit(1)

fieldnames = list(rows[0].keys())
with open(clinical_out, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
print(f"  clinical TSV: {clinical_out}  ({len(rows)} sample rows)")
PY

log "Manifests and clinical TSV generated."

# ---------------------------------------------------------------------------
# Step 2: gdc-client download — STAR Counts
# ---------------------------------------------------------------------------
STAR_MANIFEST="${MANIFEST_DIR}/luad_star_counts.manifest.tsv"
STAR_DIR="${DATA_DIR}/star_counts"
N_PROC="${GDC_N_PROCESSES:-8}"

if [[ -f "${STAR_MANIFEST}" ]]; then
  STAR_N=$(tail -n +2 "${STAR_MANIFEST}" | wc -l)
  log "Downloading STAR Counts (${STAR_N} files) → ${STAR_DIR}"
  gdc-client download \
    -m "${STAR_MANIFEST}" \
    -d "${STAR_DIR}" \
    -n "${N_PROC}" \
    --no-verify \
    2>&1 | tee -a "${LOG_FILE}"
  log "STAR Counts download finished."
else
  log "ERROR: missing ${STAR_MANIFEST}"
  exit 1
fi

# ---------------------------------------------------------------------------
# Step 3: gdc-client download — Somatic MAF
# ---------------------------------------------------------------------------
MAF_MANIFEST="${MANIFEST_DIR}/luad_somatic_maf.manifest.tsv"
MAF_DIR="${DATA_DIR}/somatic_maf"

if [[ -f "${MAF_MANIFEST}" ]]; then
  MAF_N=$(tail -n +2 "${MAF_MANIFEST}" | wc -l)
  log "Downloading Somatic MAF (${MAF_N} files) → ${MAF_DIR}"
  gdc-client download \
    -m "${MAF_MANIFEST}" \
    -d "${MAF_DIR}" \
    -n "${N_PROC}" \
    --no-verify \
    2>&1 | tee -a "${LOG_FILE}"
  log "Somatic MAF download finished."
else
  log "ERROR: missing ${MAF_MANIFEST}"
  exit 1
fi

# ---------------------------------------------------------------------------
# Step 4: Summary
# ---------------------------------------------------------------------------
log "=== Download summary ==="
log "STAR Counts files : $(find "${STAR_DIR}" -type f | wc -l)"
log "MAF files         : $(find "${MAF_DIR}" -type f | wc -l)"
log "Clinical TSV      : ${DATA_DIR}/clinical/tcga_luad_clinical.tsv"
log "Total data size   : $(du -sh "${DATA_DIR}" | cut -f1)"
log "=== TCGA-LUAD download complete ==="
