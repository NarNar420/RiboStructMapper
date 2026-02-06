import os
from typing import List, Tuple, Optional

from Bio import SeqIO
from Bio.Seq import Seq


def _parse_gtf_attributes(attr_str: str) -> dict:
    attrs = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " in part:
            key, val = part.split(" ", 1)
            val = val.strip().strip('"')
            attrs[key] = val
        else:
            attrs[part] = ""
    return attrs


def _find_seq_record(key: str, seq_records: dict):
    if key in seq_records:
        return seq_records[key]
    if len(seq_records) == 1:
        # fallback: return the only sequence
        return next(iter(seq_records.values()))
    # try matching against description or id startswith
    for rec in seq_records.values():
        if key == rec.id or key == rec.name:
            return rec
        if key in rec.description:
            return rec
        if rec.id.split()[0] == key:
            return rec
    return None


def parse_genomic_data(fasta_path: str, gtf_path: str, gene_id: Optional[str] = None, transcript_id: Optional[str] = None) -> Tuple[str, List[Tuple[str, int]]]:
    """
    Parse FASTA and GTF to extract the concatenated CDS nucleotide sequence and a coordinate map.

    Args:
        fasta_path: Path to FASTA file containing genomic sequences (chromosomes/contigs).
        gtf_path: Path to GTF/GFF file containing CDS features.
        gene_id: Optional gene_id to select. If provided, only transcripts belonging to this gene are considered.
        transcript_id: Optional transcript_id to select. If provided, that transcript's CDS is used.

    Returns:
        tuple: (nucleotide_sequence: str, coordinate_map: list of (seqname, genomic_pos) for each nucleotide in transcript order)

    Raises:
        FileNotFoundError: if either input file is missing.
        ValueError: if no CDS features are found or fasta sequence cannot be matched.
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")

    seq_records = {rec.id: rec for rec in SeqIO.parse(fasta_path, "fasta")}
    if not seq_records:
        raise ValueError(f"No sequences found in FASTA: {fasta_path}")

    # Parse GTF and collect CDS features
    cds_features = []
    with open(gtf_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            seqname, source, feature, start, end, score, strand, frame, attrs = parts[:9]
            if feature != "CDS":
                continue
            attrd = _parse_gtf_attributes(attrs)
            cds_features.append({
                "seqname": seqname,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "attributes": attrd,
            })

    if not cds_features:
        raise ValueError(f"No CDS features found in GTF: {gtf_path}")

    # Group by transcript_id
    transcripts = {}
    for f in cds_features:
        tid = f["attributes"].get("transcript_id") or f["attributes"].get("transcript") or ""
        gid = f["attributes"].get("gene_id") or f["attributes"].get("gene") or ""
        if gene_id and gid and gid != gene_id:
            continue
        transcripts.setdefault(tid, {"records": [], "gene_id": gid, "seqname": f["seqname"], "strand": f["strand"]})
        transcripts[tid]["records"].append((f["start"], f["end"]))

    if transcript_id:
        if transcript_id not in transcripts:
            raise ValueError(f"transcript_id '{transcript_id}' not found in GTF (after filtering by gene_id if provided).")
        chosen_tid = transcript_id
    else:
        # Choose the transcript with the longest total CDS length
        best = None
        best_len = -1
        for tid, v in transcripts.items():
            total = sum(e - s + 1 for s, e in v["records"])
            if total > best_len:
                best = tid
                best_len = total
        if best is None:
            raise ValueError("No transcript found after filtering GTF (maybe gene_id filter removed all transcripts).")
        chosen_tid = best

    chosen = transcripts[chosen_tid]
    seqname = chosen["seqname"]
    strand = chosen["strand"]
    seq_rec = _find_seq_record(seqname, seq_records)
    if seq_rec is None:
        raise ValueError(f"Sequence '{seqname}' referenced in GTF not found in FASTA.")

    # Normalize and sort intervals in transcript order
    intervals = chosen["records"]
    if strand == "+":
        intervals = sorted(intervals, key=lambda x: x[0])
    else:
        intervals = sorted(intervals, key=lambda x: x[0], reverse=True)

    seq_parts: List[str] = []
    coordinate_map: List[Tuple[str, int]] = []

    for s, e in intervals:
        if s < 1 or e > len(seq_rec.seq):
            # allow but warn — here we raise to keep behavior strict
            raise ValueError(f"CDS coordinates {s}-{e} are out of bounds for sequence '{seqname}' (length {len(seq_rec.seq)}).")
        seg = seq_rec.seq[s - 1 : e]
        if strand == "+":
            seq_parts.append(str(seg))
            coordinate_map.extend([(seqname, pos) for pos in range(s, e + 1)])
        else:
            # for negative strand, take reverse-complement of the exon and map coordinates from end->start
            seq_parts.append(str(seg.reverse_complement()))
            coordinate_map.extend([(seqname, pos) for pos in range(e, s - 1, -1)])

    nucleotide_seq = "".join(seq_parts)

    if len(nucleotide_seq) != len(coordinate_map):
        raise RuntimeError("Internal error: nucleotide sequence length disagrees with coordinate map length.")

    return nucleotide_seq, coordinate_map


__all__ = ["parse_genomic_data"]
