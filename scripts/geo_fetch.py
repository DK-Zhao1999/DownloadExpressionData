#!/usr/bin/env python3
"""Download and inventory processed gene-expression data from GEO.

This helper intentionally uses only the Python standard library plus pandas and
numpy, so it can run in a minimal analysis environment.  It does NOT run raw
read quantification; its job is to acquire the author-processed matrix and the
sample metadata, then make the available matrix files easy to inspect.

Commands
--------
    geo_fetch.py fetch GSE50760 --outdir raw/GSE50760
    geo_fetch.py inspect raw/GSE50760
    geo_fetch.py assemble --input raw/GSE50760/per_sample --out raw/GSE50760/matrix.tsv.gz

The fetch command writes:
    <outdir>/<GSE>_series_matrix.txt.gz      raw GEO series matrix (if present)
    <outdir>/series_matrix.tsv.gz            clean genes x samples table (if present)
    <outdir>/samples.tsv                     one row per sample
    <outdir>/series_meta.json                series-level metadata
    <outdir>/suppl/                          series-level supplementary files
    <outdir>/per_sample/<GSM>/               per-sample supplementary files
    <outdir>/fetch_manifest.json             URLs, sizes, checksums, warnings
    <outdir>/matrix_candidates.tsv           ranked candidate matrices

Only processed data are downloaded by default.  FASTQ/BAM/SRA/CEL/RAW.tar are
skipped unless --suppl all is used; this is deliberate: the skill prepares a
matrix, not a raw sequencing project.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import shutil
import sys
import tarfile
import time
import urllib.parse
import urllib.request
import zipfile
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo"
GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
USER_AGENT = "dsh-geo-expression-matrix/1.0 (research data preparation)"
TABLE_SUFFIXES = {
    ".txt",
    ".tsv",
    ".csv",
    ".tab",
    ".xls",
    ".xlsx",
    ".txt.gz",
    ".tsv.gz",
    ".csv.gz",
    ".tab.gz",
    ".xls.gz",
    ".xlsx.gz",
}
ARCHIVE_SUFFIXES = {".tar", ".tar.gz", ".tgz", ".zip"}
RAW_SUFFIXES = {
    ".fastq",
    ".fq",
    ".fastq.gz",
    ".fq.gz",
    ".bam",
    ".sam",
    ".cram",
    ".sra",
    ".cel",
    ".idat",
    ".raw",
    ".mzml",
    ".wiff",
}
NA_VALUES = ["", "NA", "N/A", "null", "NULL", "NaN", "nan", "None"]


# ---------------------------------------------------------------------------
# Generic HTTP / file helpers
# ---------------------------------------------------------------------------


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def json_safe(value):
    """Replace NaN/Inf so manifests are valid, portable JSON."""
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def https_url(url: str) -> str:
    """GEO metadata often contains ftp:// links; prefer HTTPS."""
    if url.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
        return "https://ftp.ncbi.nlm.nih.gov/" + url[len("ftp://ftp.ncbi.nlm.nih.gov/") :]
    return url


def http_request(url: str, method: str = "GET", timeout: int = 60):
    req = urllib.request.Request(
        https_url(url),
        method=method,
        headers={"User-Agent": USER_AGENT, "Accept": "*/*"},
    )
    return urllib.request.urlopen(req, timeout=timeout)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def remote_size(url: str, timeout: int = 30) -> int | None:
    try:
        with http_request(url, method="HEAD", timeout=timeout) as response:
            value = response.headers.get("Content-Length")
            return int(value) if value and value.isdigit() else None
    except Exception:
        return None


def download_file(
    url: str,
    destination: Path,
    timeout: int = 60,
    force: bool = False,
) -> dict:
    """Stream a URL to disk with retries, a .part file, and a checksum."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    url = https_url(url)
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return {
            "url": url,
            "path": str(destination),
            "size": destination.stat().st_size,
            "sha256": sha256_file(destination),
            "skipped": True,
        }

    last_error = None
    for attempt in range(3):
        temp = destination.with_name(destination.name + ".part")
        try:
            with http_request(url, timeout=timeout) as response, temp.open("wb") as out:
                expected = response.headers.get("Content-Length")
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            actual = temp.stat().st_size
            if expected and expected.isdigit() and actual != int(expected):
                raise IOError(f"size mismatch: expected {expected}, got {actual}")
            os.replace(temp, destination)
            return {
                "url": url,
                "path": str(destination),
                "size": actual,
                "sha256": sha256_file(destination),
                "skipped": False,
            }
        except Exception as exc:  # network errors vary by platform
            last_error = exc
            temp.unlink(missing_ok=True)
            if attempt < 2:
                time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"failed to download {url}: {last_error}")


def list_links(url: str, timeout: int = 60) -> list[dict]:
    """Parse a simple Apache-style directory index."""
    with http_request(url, timeout=timeout) as response:
        html = response.read().decode("utf-8", "replace")
    links = []
    seen = set()
    for href, label in re.findall(
        r'<a\s+href="([^"]+)"[^>]*>([^<]+)</a>', html, flags=re.IGNORECASE
    ):
        href = href.strip()
        label = label.strip()
        if not href or href.startswith(("../", "/", "?", "#")):
            continue
        if href.lower().startswith("http"):
            continue
        name = urllib.parse.unquote(href)
        if name in seen or label in {"Parent Directory", ""}:
            continue
        seen.add(name)
        links.append(
            {
                "name": name,
                "url": urllib.parse.urljoin(url, href),
            }
        )
    return links


def safe_member_path(destination: Path, member: str) -> Path | None:
    target = (destination / member).resolve()
    try:
        target.relative_to(destination.resolve())
    except ValueError:
        return None
    return target


def extract_archive(path: Path, destination: Path, max_members: int = 2000) -> list[Path]:
    """Extract only table-like members, refusing traversal and huge archives."""
    destination.mkdir(parents=True, exist_ok=True)
    extracted: list[Path] = []
    lower = path.name.lower()

    def is_table_member(name: str) -> bool:
        low = name.lower()
        if low.endswith(tuple(RAW_SUFFIXES)):
            return False
        return any(low.endswith(suffix) for suffix in TABLE_SUFFIXES)

    if lower.endswith((".tar", ".tar.gz", ".tgz")):
        with tarfile.open(path, "r:*") as archive:
            for member in archive.getmembers()[:max_members]:
                if not member.isfile() or not is_table_member(member.name):
                    continue
                target = safe_member_path(destination, member.name)
                if target is None:
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                with archive.extractfile(member) as src, target.open("wb") as out:
                    if src is not None:
                        shutil.copyfileobj(src, out)
                extracted.append(target)
    elif lower.endswith(".zip"):
        with zipfile.ZipFile(path) as archive:
            for info in archive.infolist()[:max_members]:
                if info.is_dir() or not is_table_member(info.filename):
                    continue
                target = safe_member_path(destination, info.filename)
                if target is None:
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(info) as src, target.open("wb") as out:
                    shutil.copyfileobj(src, out)
                extracted.append(target)
    return extracted


# ---------------------------------------------------------------------------
# GEO accession metadata and series matrix parsing
# ---------------------------------------------------------------------------


def gse_stub(accession: str) -> str:
    digits = accession[3:]
    return "GSE" + (digits[:-3] if len(digits) > 3 else "") + "nnn"


def gpl_stub(platform: str) -> str:
    digits = platform[3:]
    return "GPL" + (digits[:-3] if len(digits) > 3 else "") + "nnn"


def download_platform_annotation(
    platform: str, outdir: Path, timeout: int, force: bool
) -> dict | None:
    """Download the optional GEO platform annotation table (arrays only)."""
    stub = gpl_stub(platform)
    for filename in (f"{platform}.annot.gz", f"{platform}.annot.txt.gz"):
        url = f"{GEO_FTP}/platforms/{stub}/{platform}/annot/{filename}"
        if remote_size(url, timeout=min(timeout, 30)) is None:
            continue
        destination = outdir / "platform_annotation" / filename
        info = download_file(url, destination, timeout=timeout, force=force)
        info["platform"] = platform
        print(f"[fetch] annotation: {platform} -> {destination}")
        return info
    print(f"[fetch] no annotation table found for {platform}", file=sys.stderr)
    return None


def split_tsv_values(line: str) -> list[str]:
    if "\t" not in line:
        return []
    rest = line.split("\t", 1)[1]
    try:
        row = next(csv.reader([rest], delimiter="\t", quotechar='"'))
    except csv.Error:
        row = rest.split("\t")
    return [value.strip() for value in row]


def fetch_series_metadata(accession: str, timeout: int = 60) -> dict:
    url = f"{GEO_ACC}?acc={accession}&targ=self&form=text&view=brief"
    with http_request(url, timeout=timeout) as response:
        text = response.read().decode("utf-8", "replace")
    fields: dict[str, list[str]] = OrderedDict()
    for raw in text.splitlines():
        line = raw.strip()
        if not line.startswith("!Series_"):
            continue
        body = line[len("!Series_") :]
        if "\t" in body:
            key, value = body.split("\t", 1)
        elif " = " in body:
            key, value = body.split(" = ", 1)
        else:
            key, value = body, ""
        fields.setdefault(key.strip(), []).append(value.strip().strip('"'))

    def one(key: str, default: str = "") -> str:
        values = fields.get(key)
        return values[0] if values else default

    sample_ids: list[str] = []
    for value in fields.get("sample_id", []):
        sample_ids.extend(re.findall(r"GSM\d+", value))
    platform_ids = re.findall(r"GPL\d+", " ".join(fields.get("platform_id", [])))
    return {
        "accession": accession,
        "url": url,
        "title": one("title"),
        "status": one("status"),
        "submission_date": one("submission_date"),
        "last_update_date": one("last_update_date"),
        "type": one("type"),
        "summary": one("summary"),
        "overall_design": one("overall_design"),
        "platform_ids": sorted(set(platform_ids)),
        "sample_ids": sample_ids,
        "pubmed_ids": fields.get("pubmed_id", []),
        "relations": fields.get("relation", []),
    }


def parse_series_matrix(path: Path) -> dict:
    """Parse GEO series-matrix metadata and the optional table."""
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    sample_fields: dict[str, list[list[str]]] = OrderedDict()
    series_fields: dict[str, list[str]] = OrderedDict()
    table_lines: list[str] = []
    in_table = False

    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                if line.strip():
                    table_lines.append(line)
                continue
            if line.startswith("!Sample_"):
                key, values = line[len("!Sample_") :].split("\t", 1)[0], split_tsv_values(line)
                sample_fields.setdefault(key, []).append(values)
            elif line.startswith("!Series_"):
                key = line[len("!Series_") :].split("\t", 1)[0]
                values = split_tsv_values(line)
                value = values[0] if values else ""
                series_fields.setdefault(key, []).append(value)

    n_samples = 0
    for chunks in sample_fields.values():
        for chunk in chunks:
            n_samples = max(n_samples, len(chunk))
    sample_ids = []
    if "geo_accession" in sample_fields and sample_fields["geo_accession"]:
        sample_ids = [str(v) for v in sample_fields["geo_accession"][0]]
    for value in series_fields.get("sample_id", []):
        for gsm in re.findall(r"GSM\d+", value):
            if gsm not in sample_ids:
                sample_ids.append(gsm)
    if not sample_ids:
        sample_ids = [f"S{i + 1}" for i in range(n_samples)]
    n_samples = max(n_samples, len(sample_ids))

    columns: dict[str, list[str]] = OrderedDict()
    for key, chunks in sample_fields.items():
        for index, values in enumerate(chunks):
            column = key if index == 0 else f"{key}_{index + 1}"
            padded = [str(v) for v in values] + [""] * (n_samples - len(values))
            columns[column] = padded

    # Parse key:value characteristics into useful columns such as tissue.
    for key in list(columns):
        if not key.startswith("characteristics"):
            continue
        parsed: dict[str, list[str]] = OrderedDict()
        for value in columns[key]:
            if ":" not in value:
                continue
            name, content = value.split(":", 1)
            name = re.sub(r"[^0-9A-Za-z]+", "_", name.strip()).strip("_").lower()
            if not name:
                continue
            name = re.sub(r"_+", "_", name)
            parsed.setdefault(name, []).append(content.strip())
        for name, values in parsed.items():
            if len(values) != n_samples:
                continue
            target = name
            suffix = 2
            while target in columns:
                target = f"{name}_{suffix}"
                suffix += 1
            columns[target] = values

    # The canonical sample ID column must be first and unique.
    sample_ids = [str(v) for v in columns.get("geo_accession", sample_ids)]
    if len(sample_ids) < n_samples:
        sample_ids.extend(f"S{i + 1}" for i in range(len(sample_ids), n_samples))
    sample_ids = sample_ids[:n_samples]
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("series matrix contains duplicate sample IDs")
    ordered = OrderedDict()
    ordered["sample_id"] = sample_ids
    for key, values in columns.items():
        if key == "geo_accession":
            continue
        ordered[key] = values[:n_samples]
    samples = pd.DataFrame(ordered)
    samples["sample_id"] = samples["sample_id"].astype(str)

    matrix = None
    matrix_stats = None
    if table_lines:
        text = "\n".join(table_lines)
        frame = pd.read_csv(
            io.StringIO(text),
            sep="\t",
            quotechar='"',
            dtype=str,
            keep_default_na=True,
            na_values=NA_VALUES,
        )
        if frame.shape[0] > 0 and frame.shape[1] > 1:
            first = frame.columns[0]
            gene_ids = frame[first].astype(str)
            raw_values = frame.iloc[:, 1:]
            values = raw_values.apply(pd.to_numeric, errors="coerce")
            nonnumeric = int(
                ((raw_values.notna()) & values.isna()).to_numpy().sum()
            )
            matrix = pd.concat(
                [gene_ids.rename("gene_id"), values], axis=1
            )
            matrix = matrix[matrix["gene_id"].notna() & (matrix["gene_id"] != "")]
            array = matrix.iloc[:, 1:].to_numpy(dtype=float)
            matrix_stats = {
                "rows": int(matrix.shape[0]),
                "samples": int(matrix.shape[1] - 1),
                "na_cells": int(np.isnan(array).sum()),
                "nonnumeric_cells": nonnumeric,
                "min": float(np.nanmin(array)) if np.isfinite(array).any() else None,
                "max": float(np.nanmax(array)) if np.isfinite(array).any() else None,
            }

    series_meta = {}
    for key, values in series_fields.items():
        series_meta[key] = values[0] if len(values) == 1 else values
    return {
        "samples": samples,
        "matrix": matrix,
        "matrix_stats": matrix_stats,
        "series_meta": series_meta,
    }


# ---------------------------------------------------------------------------
# Table inventory / candidate ranking
# ---------------------------------------------------------------------------


def detect_separator(path: Path) -> tuple[str, int]:
    """Return (separator, skiprows) for a delimited text file.

    Leading blank/comment/annotation lines are skipped.  This matters for GEO
    platform annotation files and for tables with a short free-text preamble.
    """
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    lines: list[str] = []
    skiprows = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for index, raw in enumerate(fh):
            line = raw.rstrip("\r\n")
            if not line.strip() or line.startswith(("!", "#")):
                skiprows = index + 1
                continue
            lines.append(line)
            if len(lines) >= 12:
                break
    if not lines:
        return "\t", skiprows
    counts = {sep: lines[0].count(sep) for sep in ("\t", ",", ";")}
    separator = max(counts, key=counts.get)
    if counts[separator] == 0:
        separator = r"\s+"
    return separator, skiprows


def read_table_preview(path: Path, nrows: int = 30) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path, nrows=nrows, dtype=str)
    sep, skiprows = detect_separator(path)
    return pd.read_csv(
        path,
        sep=sep,
        nrows=nrows,
        skiprows=skiprows,
        quotechar='"',
        engine="python",
        dtype=str,
        keep_default_na=True,
        na_values=NA_VALUES,
    )


def infer_value_type(name: str, columns: list[str], minimum: float, maximum: float) -> str:
    haystack = " ".join([name] + [str(c) for c in columns]).lower()
    if "tpm" in haystack:
        return "TPM"
    if "fpkm" in haystack or "rpkm" in haystack:
        return "FPKM/RPKM"
    if "count" in haystack or "read count" in haystack:
        return "counts"
    if "log2" in haystack or "log-2" in haystack:
        return "log2"
    if "normalized" in haystack or "normalised" in haystack:
        return "normalized"
    if minimum < -0.01:
        return "log2/centered (negative values)"
    if maximum > 1000:
        return "linear (TPM/counts-like range)"
    if maximum <= 30:
        return "log2 or normalized"
    return "unknown linear"


def inspect_table(path: Path) -> dict | None:
    name = path.name
    low = name.lower()
    if low.endswith(tuple(RAW_SUFFIXES)) or low.endswith(tuple(ARCHIVE_SUFFIXES)):
        return None
    if not any(low.endswith(suffix) for suffix in TABLE_SUFFIXES):
        return None
    try:
        frame = read_table_preview(path)
    except Exception as exc:
        return {
            "file": str(path),
            "name": name,
            "error": str(exc),
            "score": -999,
            "recommended": "no",
        }
    if frame.empty or frame.shape[1] < 2:
        return None
    first = frame.columns[0]
    numeric_columns = []
    for column in frame.columns[1:]:
        values = pd.to_numeric(frame[column], errors="coerce")
        if values.notna().mean() >= 0.8:
            numeric_columns.append(column)
    if not numeric_columns:
        return None
    numeric = frame[numeric_columns].apply(pd.to_numeric, errors="coerce")
    minimum = float(np.nanmin(numeric.to_numpy(dtype=float)))
    maximum = float(np.nanmax(numeric.to_numpy(dtype=float)))
    value_type = infer_value_type(name, list(numeric_columns), minimum, maximum)
    duplicate_first = bool(frame[first].duplicated().any())
    layout = "long/duplicated-ID" if duplicate_first and len(numeric_columns) <= 2 else "wide"
    score = 20
    if value_type == "TPM":
        score = 100
    elif value_type == "FPKM/RPKM":
        score = 65
    elif value_type == "normalized":
        score = 45
    elif value_type == "log2":
        score = 35
    elif value_type == "counts":
        score = 15
    elif value_type.startswith("log2/centered"):
        score = 25
    elif value_type.startswith("linear"):
        score = 30
    if len(numeric_columns) >= 10:
        score += 40
    elif len(numeric_columns) >= 3:
        score += 15
    elif len(numeric_columns) == 1:
        score -= 20
    if frame.shape[0] >= 1000:
        score += 10
    if layout.startswith("long"):
        score -= 25
    if any("gene_assignment" in str(c).lower() for c in frame.columns):
        score -= 20
    return {
        "file": str(path),
        "name": name,
        "rows_previewed": int(frame.shape[0]),
        "columns": int(frame.shape[1]),
        "first_column": str(first),
        "numeric_columns": len(numeric_columns),
        "layout": layout,
        "value_type_guess": value_type,
        "min_preview": minimum,
        "max_preview": maximum,
        "score": score,
        "recommended": "no",
    }


def inspect_directory(path: Path, max_files: int = 80) -> list[dict]:
    files: list[Path] = []
    if path.is_file():
        files = [path]
    else:
        for candidate in sorted(path.rglob("*")):
            if not candidate.is_file():
                continue
            if candidate.name.endswith(".part") or candidate.name.startswith("."):
                continue
            if candidate.name in {
                "fetch_manifest.json",
                "assemble_manifest.json",
                "samples.tsv",
                "matrix_candidates.tsv",
                "filelist.txt",
                "clinical_clean.csv",
                "survival_analysis.csv",
                "sample_map.tsv",
            }:
                continue
            if re.search(
                r"(_sample_qc|_qc_summary|_metadata|_processing_log|_sessionInfo"
                r"|_removed_samples|_subset_summary)\.(tsv|csv|txt)$",
                candidate.name,
                re.I,
            ):
                continue
            files.append(candidate)
            if len(files) >= max_files:
                break
    results = []
    for file in files:
        result = inspect_table(file)
        if result:
            results.append(result)
    results.sort(key=lambda item: item.get("score", -999), reverse=True)
    if results:
        results[0]["recommended"] = "yes"
    return results


def write_candidates(results: list[dict], outdir: Path) -> Path:
    path = outdir / "matrix_candidates.tsv"
    columns = [
        "file",
        "name",
        "rows_previewed",
        "columns",
        "first_column",
        "numeric_columns",
        "layout",
        "value_type_guess",
        "min_preview",
        "max_preview",
        "score",
        "recommended",
        "error",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for result in results:
            writer.writerow(result)
    return path


# ---------------------------------------------------------------------------
# Assemble per-sample tables
# ---------------------------------------------------------------------------


def resolve_column(frame: pd.DataFrame, spec: str | None, default_index: int) -> str:
    if spec is None:
        if frame.shape[1] <= default_index:
            raise ValueError(f"table has only {frame.shape[1]} columns")
        return str(frame.columns[default_index])
    if spec in frame.columns:
        return str(spec)
    if spec.isdigit():
        index = int(spec)
        if index >= frame.shape[1]:
            raise ValueError(f"column index {index} out of range")
        return str(frame.columns[index])
    raise ValueError(f"column not found: {spec}")


def assemble_files(
    files: list[Path],
    output: Path,
    gene_column: str | None = None,
    value_column: str | None = None,
    sample_regex: str | None = r"(GSM\d+)",
    duplicate: str = "error",
) -> dict:
    series = {}
    used = []
    for path in files:
        try:
            frame = read_table_preview(path, nrows=10_000_000)
        except Exception as exc:
            print(f"[assemble] skip {path.name}: {exc}", file=sys.stderr)
            continue
        if frame.shape[1] < 2:
            continue
        try:
            gene_col = resolve_column(frame, gene_column, 0)
            val_col = resolve_column(frame, value_column, 1)
        except ValueError as exc:
            print(f"[assemble] skip {path.name}: {exc}", file=sys.stderr)
            continue
        genes = frame[gene_col].astype(str).str.strip()
        values = pd.to_numeric(frame[val_col], errors="coerce")
        keep = genes.notna() & (genes != "") & (genes.str.lower() != "nan")
        genes = genes[keep]
        values = values[keep]
        if genes.empty:
            continue
        if genes.duplicated().any():
            if duplicate == "error":
                raise ValueError(f"{path.name}: duplicate gene IDs; use --duplicate mean|max|first")
            aggregated = pd.DataFrame({"gene": genes, "value": values})
            if duplicate == "first":
                aggregated = aggregated.drop_duplicates("gene", keep="first")
            elif duplicate == "mean":
                aggregated = aggregated.groupby("gene", as_index=False)["value"].mean()
            elif duplicate == "max":
                aggregated = aggregated.groupby("gene", as_index=False)["value"].max()
            genes = aggregated["gene"]
            values = aggregated["value"]
        sample = path.name
        if sample_regex:
            match = re.search(sample_regex, path.name)
            if match:
                sample = match.group(1) if match.groups() else match.group(0)
        if not sample:
            sample = path.name
        if sample in series:
            sample = f"{sample}_{len(series) + 1}"
        series[sample] = pd.Series(values.to_numpy(), index=genes.to_numpy(), name=sample)
        used.append(str(path))

    if len(series) < 1:
        raise ValueError("no usable per-sample tables found")
    matrix = pd.concat(series.values(), axis=1)
    matrix.index.name = "gene_id"
    matrix = matrix.sort_index()
    output.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if str(output).lower().endswith(".gz") else None
    matrix.to_csv(output, sep="\t", compression=compression, na_rep="NA")
    return {
        "output": str(output),
        "samples": int(matrix.shape[1]),
        "genes": int(matrix.shape[0]),
        "sample_names": list(matrix.columns),
        "input_files": used,
        "value_type_hint": value_type_from_names(used),
    }


def value_type_from_names(names: list[str]) -> str | None:
    haystack = " ".join(names).lower()
    if "tpm" in haystack:
        return "TPM"
    if "fpkm" in haystack or "rpkm" in haystack:
        return "FPKM"
    if "count" in haystack:
        return "counts"
    return None


# ---------------------------------------------------------------------------
# fetch command
# ---------------------------------------------------------------------------


def download_series_matrix(accession: str, outdir: Path, timeout: int, force: bool) -> dict | None:
    matrix_dir = f"{GEO_FTP}/series/{gse_stub(accession)}/{accession}/matrix/"
    try:
        links = list_links(matrix_dir, timeout=timeout)
    except Exception as exc:
        print(f"[fetch] cannot list matrix directory: {exc}", file=sys.stderr)
        return None
    target = None
    for link in links:
        if re.search(r"_series_matrix\.txt(\.gz)?$", link["name"], re.I):
            target = link
            break
    if target is None:
        print("[fetch] no series-matrix file found", file=sys.stderr)
        return None
    destination = outdir / target["name"]
    info = download_file(target["url"], destination, timeout=timeout, force=force)
    info["name"] = target["name"]
    return info


def collect_supplementary(
    samples: pd.DataFrame,
    series_meta: dict,
    accession: str,
) -> list[dict]:
    entries: list[dict] = []
    seen = set()

    def add(url: str, sample_id: str | None, origin: str) -> None:
        for token in re.findall(r"(?:https?|ftp)://[^\s\"']+", str(url)):
            token = https_url(token.strip())
            if not token or token in seen:
                continue
            seen.add(token)
            name = Path(urllib.parse.urlparse(token).path).name
            if not name:
                continue
            entries.append(
                {
                    "url": token,
                    "name": name,
                    "sample_id": sample_id,
                    "origin": origin,
                }
            )

    for column in samples.columns:
        if not str(column).startswith("supplementary_file"):
            continue
        for sample_id, value in zip(samples["sample_id"], samples[column]):
            if isinstance(value, str) and value.strip():
                add(value, str(sample_id), "sample")
    series_files = series_meta.get("supplementary_file", [])
    if isinstance(series_files, str):
        series_files = [series_files]
    for value in series_files:
        if isinstance(value, list):
            for item in value:
                add(item, None, "series")
        elif value:
            add(value, None, "series")

    try:
        for link in list_links(
            f"{GEO_FTP}/series/{gse_stub(accession)}/{accession}/suppl/"
        ):
            add(link["url"], None, "series-dir")
    except Exception:
        pass
    return entries


def should_skip_supplement(
    entry: dict,
    mode: str,
    has_matrix: bool,
    has_per_sample_tables: bool,
) -> tuple[bool, str]:
    name = entry["name"].lower()
    if name == "filelist.txt":
        return False, ""
    if mode == "all":
        return False, ""
    if name.endswith(tuple(RAW_SUFFIXES)):
        return True, "raw data type"
    if "raw" in name and name.endswith(tuple(ARCHIVE_SUFFIXES)):
        if not has_matrix and not has_per_sample_tables and entry.get("origin") == "series":
            return False, ""
        return True, "RAW archive"
    if name.endswith(tuple(ARCHIVE_SUFFIXES)):
        if has_matrix or has_per_sample_tables:
            return True, "archive while processed tables are available"
    return False, ""


def fetch_command(args: argparse.Namespace) -> int:
    accession = args.accession.upper()
    if not re.fullmatch(r"GSE\d+", accession):
        print("error: accession must look like GSE12345", file=sys.stderr)
        return 2
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[fetch] accession: {accession}")
    series_meta = fetch_series_metadata(accession, timeout=args.timeout)
    print(f"[fetch] title: {series_meta.get('title', '')[:160]}")
    if series_meta.get("platform_ids"):
        print(f"[fetch] platform(s): {', '.join(series_meta['platform_ids'])}")
        if len(series_meta["platform_ids"]) > 1:
            print(
                "[fetch] WARNING: multiple platforms; do not combine without a batch plan",
                file=sys.stderr,
            )

    manifest = {
        "accession": accession,
        "fetched_at_utc": now_utc(),
        "geo_metadata": series_meta,
        "series_matrix": None,
        "samples_tsv": None,
        "supplementary": [],
        "annotations": [],
        "extracted": [],
        "candidates": [],
        "warnings": [],
    }

    if args.annot and series_meta.get("platform_ids"):
        for platform in series_meta["platform_ids"]:
            try:
                info = download_platform_annotation(
                    platform, outdir, args.timeout, args.force
                )
                if info:
                    manifest["annotations"].append(info)
            except Exception as exc:
                manifest["warnings"].append(
                    f"failed to download annotation for {platform}: {exc}"
                )

    series_info = download_series_matrix(accession, outdir, args.timeout, args.force)
    samples = None
    matrix_info = None
    if not series_info:
        manifest["warnings"].append("no GEO series matrix file was downloaded")
    if series_info:
        manifest["series_matrix"] = series_info
        try:
            parsed = parse_series_matrix(Path(series_info["path"]))
            samples = parsed["samples"]
            samples.to_csv(outdir / "samples.tsv", sep="\t", index=False, lineterminator="\n")
            manifest["samples_tsv"] = str(outdir / "samples.tsv")
            if parsed["matrix"] is not None:
                matrix_path = outdir / "series_matrix.tsv.gz"
                parsed["matrix"].to_csv(
                    matrix_path, sep="\t", index=False, compression="gzip", na_rep="NA"
                )
                matrix_info = {"path": str(matrix_path), **parsed["matrix_stats"]}
                manifest["series_matrix"]["parsed"] = matrix_info
            else:
                manifest["warnings"].append(
                    "series matrix has no expression table; looking for supplementary matrices"
                )
            merged_meta = dict(series_meta)
            merged_meta.update(parsed["series_meta"])
            merged_meta["accession"] = accession
            merged_meta["platform_ids"] = sorted(
                set(
                    series_meta.get("platform_ids", [])
                    + re.findall(r"GPL\d+", json.dumps(parsed["series_meta"]))
                )
            )
            (outdir / "series_meta.json").write_text(
                json.dumps(json_safe(merged_meta), ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
            series_meta = merged_meta
        except Exception as exc:
            manifest["warnings"].append(f"failed to parse series matrix: {exc}")
            print(f"[fetch] WARNING: failed to parse series matrix: {exc}", file=sys.stderr)

    if samples is None:
        # No matrix file: still provide the series sample list as a minimal table.
        ids = series_meta.get("sample_ids", [])
        samples = pd.DataFrame({"sample_id": ids if ids else []})
        samples.to_csv(outdir / "samples.tsv", sep="\t", index=False, lineterminator="\n")
        manifest["samples_tsv"] = str(outdir / "samples.tsv")

    suppl_entries = collect_supplementary(samples, series_meta, accession)
    has_per_sample_tables = any(
        "/samples/" in entry["url"]
        and any(entry["name"].lower().endswith(suffix) for suffix in TABLE_SUFFIXES)
        for entry in suppl_entries
    )
    print(f"[fetch] supplementary candidates: {len(suppl_entries)}")
    for entry in suppl_entries:
        skip, reason = should_skip_supplement(
            entry, args.suppl, matrix_info is not None, has_per_sample_tables
        )
        if skip:
            entry["skipped"] = reason
            manifest["supplementary"].append(entry)
            continue
        if args.suppl == "none":
            entry["skipped"] = "--suppl none"
            manifest["supplementary"].append(entry)
            continue
        size = None
        if "/samples/" not in entry["url"]:
            size = remote_size(entry["url"], timeout=min(args.timeout, 30))
        entry["remote_size"] = size
        if (
            args.suppl != "all"
            and size is not None
            and size > args.max_suppl_mb * 1024 * 1024
        ):
            entry["skipped"] = f"larger than --max-suppl-mb={args.max_suppl_mb}"
            manifest["supplementary"].append(entry)
            print(f"[fetch] skip large file: {entry['name']} ({size / 1e6:.1f} MB)")
            continue
        if entry.get("sample_id") and "/samples/" in entry["url"]:
            destination = outdir / "per_sample" / entry["sample_id"] / entry["name"]
        else:
            destination = outdir / "suppl" / entry["name"]
        try:
            print(f"[fetch] download: {entry['name']}")
            info = download_file(entry["url"], destination, timeout=args.timeout, force=args.force)
            entry.update(info)
            entry["destination"] = str(destination)
            manifest["supplementary"].append(entry)
            if args.delay > 0:
                time.sleep(args.delay)
        except Exception as exc:
            entry["error"] = str(exc)
            manifest["supplementary"].append(entry)
            manifest["warnings"].append(f"failed to download {entry['name']}: {exc}")
            print(f"[fetch] WARNING: {exc}", file=sys.stderr)

    # Extract table-like members from any downloaded archives.
    for archive in outdir.rglob("*"):
        if not archive.is_file() or not archive.name.lower().endswith(
            (".tar", ".tar.gz", ".tgz", ".zip")
        ):
            continue
        try:
            extracted = extract_archive(archive, outdir / "extracted" / archive.stem)
            manifest["extracted"].extend(str(path) for path in extracted)
            if extracted:
                print(f"[fetch] extracted {len(extracted)} table(s) from {archive.name}")
        except Exception as exc:
            manifest["warnings"].append(f"failed to extract {archive.name}: {exc}")

    # If the series matrix is empty and we have per-sample tables, assemble them.
    per_sample_files = sorted((outdir / "per_sample").rglob("*.txt.gz")) + sorted(
        (outdir / "per_sample").rglob("*.tsv.gz")
    )
    if matrix_info is None and len(per_sample_files) >= 2:
        try:
            hint = value_type_from_names([path.name for path in per_sample_files])
            suffix = f"_{hint}" if hint else ""
            assembled_path = outdir / f"per_sample_matrix{suffix}.tsv.gz"
            result = assemble_files(
                per_sample_files,
                assembled_path,
                gene_column=None,
                value_column=None,
                sample_regex=r"(GSM\d+)",
                duplicate="mean",
            )
            manifest["assembled_per_sample"] = result
            manifest["warnings"].append(
                "per-sample tables were assembled with duplicate=mean; "
                "verify duplicate gene IDs and sample names before publication"
            )
            print(
                f"[fetch] assembled {result['samples']} samples x {result['genes']} genes: "
                f"{assembled_path}"
            )
        except Exception as exc:
            manifest["warnings"].append(f"could not assemble per-sample tables: {exc}")
            print(f"[fetch] WARNING: could not assemble per-sample tables: {exc}", file=sys.stderr)

    candidates = inspect_directory(outdir)
    manifest["candidates"] = candidates
    write_candidates(candidates, outdir)
    (outdir / "fetch_manifest.json").write_text(
        json.dumps(json_safe(manifest), ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )

    print("\n[fetch] candidate matrices (ranked):")
    for item in candidates[:10]:
        print(
            f"  {item.get('score', ''):>4}  {item.get('value_type_guess', ''):<28} "
            f"{item.get('name', '')}"
        )
    print(f"\n[fetch] done: {outdir}")
    print("[fetch] next: inspect matrix_candidates.tsv and choose the TPM/processed matrix")
    return 0


# ---------------------------------------------------------------------------
# inspect / assemble commands
# ---------------------------------------------------------------------------


def inspect_command(args: argparse.Namespace) -> int:
    path = Path(args.path).expanduser().resolve()
    if not path.exists():
        print(f"error: not found: {path}", file=sys.stderr)
        return 2
    outdir = Path(args.outdir).expanduser().resolve() if args.outdir else (
        path if path.is_dir() else path.parent
    )
    outdir.mkdir(parents=True, exist_ok=True)
    results = inspect_directory(path, max_files=args.max_files)
    write_candidates(results, outdir)
    if not results:
        print("no numeric table candidates found")
        return 0
    print(f"{'score':>5}  {'guess':<30} {'rows':>6} {'cols':>5}  file")
    for item in results:
        print(
            f"{item.get('score', ''):>5}  {item.get('value_type_guess', ''):<30} "
            f"{item.get('rows_previewed', ''):>6} {item.get('columns', ''):>5}  "
            f"{item.get('name', '')}"
        )
    print(f"\nwrote: {outdir / 'matrix_candidates.tsv'}")
    return 0


def assemble_command(args: argparse.Namespace) -> int:
    input_dir = Path(args.input).expanduser().resolve()
    output = Path(args.out).expanduser().resolve()
    files = sorted(
        path
        for path in input_dir.rglob(args.pattern)
        if path.is_file()
        and "series_matrix" not in path.name.lower()
        and not path.name.endswith(".part")
    )
    if not files:
        print(f"error: no files matching {args.pattern} under {input_dir}", file=sys.stderr)
        return 2
    result = assemble_files(
        files,
        output,
        gene_column=args.gene_column,
        value_column=args.value_column,
        sample_regex=args.sample_regex,
        duplicate=args.duplicate,
    )
    result["created_at_utc"] = now_utc()
    manifest_path = output.parent / "assemble_manifest.json"
    manifest_path.write_text(
        json.dumps(json_safe(result), ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(
        f"[assemble] {result['samples']} samples x {result['genes']} genes -> {output}"
    )
    print(f"[assemble] manifest: {manifest_path}")
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download processed GEO expression matrices and sample metadata."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    fetch = sub.add_parser("fetch", help="download a GEO series matrix and processed supplements")
    fetch.add_argument("accession", help="GEO series accession, e.g. GSE50760")
    fetch.add_argument("--outdir", required=True, help="dataset output directory")
    fetch.add_argument(
        "--suppl",
        choices=["auto", "all", "none"],
        default="auto",
        help="auto: processed/table-like files; all: every supplementary file",
    )
    fetch.add_argument("--max-suppl-mb", type=int, default=500, help="skip larger files in auto mode")
    fetch.add_argument(
        "--annot",
        action="store_true",
        help="also download the GPL platform annotation table (arrays)",
    )
    fetch.add_argument("--timeout", type=int, default=60, help="HTTP timeout in seconds")
    fetch.add_argument(
        "--delay",
        type=float,
        default=0.2,
        help="seconds between supplementary downloads (be polite to NCBI)",
    )
    fetch.add_argument("--force", action="store_true", help="re-download existing files")
    fetch.set_defaults(func=fetch_command)

    inspect = sub.add_parser("inspect", help="rank local matrix-like tables")
    inspect.add_argument("path", help="file or directory to inspect")
    inspect.add_argument("--outdir", default=None, help="where to write matrix_candidates.tsv")
    inspect.add_argument("--max-files", type=int, default=80)
    inspect.set_defaults(func=inspect_command)

    assemble = sub.add_parser("assemble", help="merge per-sample two-column tables into a matrix")
    assemble.add_argument("--input", required=True, help="directory containing per-sample tables")
    assemble.add_argument("--out", required=True, help="output .tsv.gz matrix")
    assemble.add_argument("--pattern", default="*.txt.gz", help="glob pattern under --input")
    assemble.add_argument("--gene-column", default=None, help="gene column name or 0-based index")
    assemble.add_argument("--value-column", default=None, help="value column name or 0-based index")
    assemble.add_argument(
        "--sample-regex",
        default=r"(GSM\d+)",
        help="regex used to derive sample names from file names; first group wins",
    )
    assemble.add_argument(
        "--duplicate",
        choices=["error", "mean", "max", "first"],
        default="error",
        help="how to collapse duplicate gene IDs within a per-sample file",
    )
    assemble.set_defaults(func=assemble_command)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        return 130
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())


