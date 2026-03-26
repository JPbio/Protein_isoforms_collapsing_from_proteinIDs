#!/usr/bin/env python3

"""
Protein isoform collapsing from protein IDs or BLAST/DIAMOND tables.

This script supports two modes:

1) List mode (-l)
   Input:
     - a text file with one protein accession per line
     - an NCBI gene mapping file (gene2accession or gene2refseq)

   Output:
     - protein-to-gene mapping table
     - isoform groups inferred from shared GeneID
     - grouping report
     - unmapped proteins
     - summary statistics

2) Table mode (-t)
   Input:
     - a BLAST/DIAMOND-like tabular file with the expected header
     - an NCBI gene mapping file (gene2accession or gene2refseq)

   Output:
     - annotated full table
     - filtered representative table
     - grouping report showing mapped and unmapped relationships
     - summary statistics

Core logic:
- First try exact accession -> GeneID mapping
- If exact mapping fails, try version-stripped mapping
  Example:
    XP_015836279.1 -> XP_015836279
    XP_015836279.2 -> XP_015836279

Why this is needed:
- RefSeq proteins are often mappable to GeneIDs
- Many non-RefSeq / GenBank / nr accessions are not
- Therefore, in table mode, the script also applies heuristics to collapse
  likely redundant unmapped proteins based on alignment behavior

Author intent:
- Provide a practical preprocessing step for phylogeny, domain analysis,
  homology-based discovery, and redundancy reduction in large protein hit sets
"""

import sys
import gzip
import argparse
import os
import csv
from collections import defaultdict


# Expected header for table mode (-t).
# The script checks that all these columns are present.
EXPECTED_HEADER = [
    "qseqid", "qlen", "qstart", "qend", "length", "qcovhsp", "pident",
    "sseqid", "slen", "sstart", "send", "scovhsp", "evalue", "bitscore",
    "qstrand", "qframe", "stitle", "staxids", "slineages", "sscinames",
    "skingdoms", "sphylums", "db_label"
]


def status(message):
    """
    Print progress and status messages to stderr.

    Why stderr?
    - keeps progress separate from data files
    - easier to redirect stdout if needed in future
    """
    sys.stderr.write(f"[INFO] {message}\n")


def open_maybe_gzip(path):
    """
    Open a plain text or .gz file transparently.

    Parameters
    ----------
    path : str
        File path.

    Returns
    -------
    file handle
        Open text handle.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def infer_prefix(path):
    """
    Infer output prefix from the input filename.

    Example
    -------
    'yakuba.tsv' -> 'yakuba'
    'hits.final.txt' -> 'hits'
    """
    base = os.path.basename(path)
    return base.split(".")[0]


def strip_version(accession):
    """
    Remove trailing accession version.

    Examples
    --------
    XP_015836279.1 -> XP_015836279
    NP_123456      -> NP_123456
    """
    return accession.rsplit(".", 1)[0] if "." in accession else accession


def safe_int(x, default=0):
    """
    Convert a value to int safely.

    Useful for sorting numeric columns that may occasionally be empty
    or malformed.
    """
    try:
        return int(str(x).strip())
    except Exception:
        return default


def gene_sort_value(gene_id):
    """
    Sorting helper for GeneID columns.

    Behavior
    --------
    - numeric GeneIDs sort numerically
    - NA goes last
    - non-numeric strings sort after numeric GeneIDs but before NA
    """
    gene_id = str(gene_id).strip()
    if gene_id == "NA":
        return (1, float("inf"))
    if gene_id.isdigit():
        return (0, int(gene_id))
    return (0, gene_id)


# =========================================================
# Shared mapping engine
# =========================================================

def build_gene_mapping(accessions, gene_file):
    """
    Build protein accession -> GeneID mapping.

    Strategy
    --------
    1. Exact accession match is preferred
    2. If exact match is unavailable, use version-stripped match

    Parameters
    ----------
    accessions : list[str]
        Protein accessions to map.
    gene_file : str
        gene2accession or gene2refseq file (.gz or plain text)

    Returns
    -------
    mapping : dict
        accession -> GeneID
    match_type : dict
        accession -> one of:
        - exact
        - version_stripped
        - unmapped
    """
    wanted_exact = set(accessions)
    wanted_base = {strip_version(x) for x in accessions}

    # Databases built from the NCBI mapping file.
    exact_db = {}
    base_db = {}

    status(f"Reading gene mapping file: {gene_file}")

    with open_maybe_gzip(gene_file) as f:
        for i, line in enumerate(f, 1):
            if i % 1_000_000 == 0:
                status(f"Processed {i:,} lines from gene mapping file")

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue

            # NCBI columns:
            # field[1] = GeneID
            # field[5] = protein_accession.version
            gene_id = fields[1].strip()
            protein_acc = fields[5].strip()

            if not protein_acc or protein_acc == "-":
                continue

            base_acc = strip_version(protein_acc)

            # Save only the first observed mapping per accession/base.
            # This is consistent with previous script behavior.
            if protein_acc in wanted_exact and protein_acc not in exact_db:
                exact_db[protein_acc] = gene_id

            if base_acc in wanted_base and base_acc not in base_db:
                base_db[base_acc] = gene_id

    mapping = {}
    match_type = {}

    for acc in accessions:
        base = strip_version(acc)

        if acc in exact_db:
            mapping[acc] = exact_db[acc]
            match_type[acc] = "exact"
        elif base in base_db:
            mapping[acc] = base_db[base]
            match_type[acc] = "version_stripped"
        else:
            match_type[acc] = "unmapped"

    exact_n = sum(1 for v in match_type.values() if v == "exact")
    stripped_n = sum(1 for v in match_type.values() if v == "version_stripped")
    unmapped_n = sum(1 for v in match_type.values() if v == "unmapped")

    status("Finished scanning gene mapping file")
    status(f"Mapped accessions: {len(mapping)}")
    status(f"  exact matches: {exact_n}")
    status(f"  version-stripped matches: {stripped_n}")
    status(f"  unmapped: {unmapped_n}")

    return mapping, match_type


# =========================================================
# List mode (-l)
# =========================================================

def read_protein_list(path):
    """
    Read a plain text list of protein accessions.

    Rules
    -----
    - one accession per line
    - empty lines ignored
    - lines beginning with '#' ignored
    - duplicate accessions removed, preserving first occurrence order
    """
    proteins = []
    seen = set()

    status(f"Reading protein list from: {path}")

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            acc = line.split()[0]
            if acc not in seen:
                proteins.append(acc)
                seen.add(acc)

    status(f"Loaded {len(proteins)} unique protein accessions")
    return proteins


def build_gene_to_proteins(mapping):
    """
    Reverse accession -> GeneID mapping into GeneID -> [proteins].
    """
    gene_to_proteins = defaultdict(list)
    for protein, gene_id in mapping.items():
        gene_to_proteins[gene_id].append(protein)
    return gene_to_proteins


def write_list_mapping_table(protein_list, mapping, match_type, out_tsv):
    """
    Write protein-to-gene table for list mode.

    Includes all input proteins:
    - mapped proteins get their GeneID
    - unmapped proteins get GeneID = NA

    Sorted by GeneID to make isoform inspection easier.
    """
    status(f"Writing (sorted by GeneID): {out_tsv}")

    rows = []
    for protein in protein_list:
        gene_id = mapping.get(protein, "NA")
        mtype = match_type.get(protein, "unmapped")
        rows.append((protein, gene_id, mtype))

    rows = sorted(rows, key=lambda x: (gene_sort_value(x[1]), x[0]))

    with open(out_tsv, "w") as out:
        out.write("Protein_accession\tGeneID\tMatch_type\n")
        for protein, gene_id, mtype in rows:
            out.write(f"{protein}\t{gene_id}\t{mtype}\n")


def write_list_isoform_groups(gene_to_proteins, out_isoforms):
    """
    Write only GeneIDs associated with more than one input protein.

    This is the simplest GeneID-based isoform grouping report.
    """
    status(f"Writing: {out_isoforms}")

    isoform_rows = []
    for gene_id, proteins in gene_to_proteins.items():
        if len(proteins) > 1:
            isoform_rows.append((str(gene_id).strip(), len(proteins), sorted(proteins)))

    isoform_rows = sorted(isoform_rows, key=lambda x: (gene_sort_value(x[0]), x[2][0]))

    with open(out_isoforms, "w") as out:
        out.write("GeneID\tN_proteins\tProteins\n")
        for gene_id, n_proteins, proteins in isoform_rows:
            out.write(f"{gene_id}\t{n_proteins}\t{','.join(proteins)}\n")

    status(f"Isoform groups written: {len(isoform_rows)}")


def write_list_groups_report(protein_list, mapping, out_groups):
    """
    Write a complete grouping report for list mode.

    Mapped proteins are grouped by GeneID.
    Unmapped proteins are reported one by one as unmapped_single.
    """
    status(f"Writing group report: {out_groups}")

    mapped_groups = defaultdict(list)
    unmapped = []

    for protein in protein_list:
        gene_id = mapping.get(protein, "NA")
        if gene_id == "NA":
            unmapped.append(protein)
        else:
            mapped_groups[gene_id].append(protein)

    with open(out_groups, "w") as out:
        out.write("Group_type\tGroup_id\tN_proteins\tProteins\n")

        for gene_id in sorted(mapped_groups, key=gene_sort_value):
            proteins = sorted(mapped_groups[gene_id])
            out.write(f"mapped_gene\t{gene_id}\t{len(proteins)}\t{','.join(proteins)}\n")

        for protein in sorted(unmapped):
            out.write(f"unmapped_single\tNA\t1\t{protein}\n")


def write_list_unmapped(protein_list, mapping, out_unmapped):
    """
    Write unmapped proteins as a simple text list.
    """
    status(f"Writing: {out_unmapped}")

    unmapped = [p for p in protein_list if p not in mapping]

    with open(out_unmapped, "w") as out:
        for protein in unmapped:
            out.write(protein + "\n")

    status(f"Unmapped proteins written: {len(unmapped)}")


def write_list_summary(protein_list, mapping, match_type, gene_to_proteins, out_summary):
    """
    Write summary statistics for list mode.
    """
    status(f"Writing: {out_summary}")

    total = len(protein_list)
    mapped = len(mapping)
    unmapped = total - mapped
    genes = len(gene_to_proteins)
    iso_genes = sum(1 for prots in gene_to_proteins.values() if len(prots) > 1)
    iso_prots = sum(len(prots) for prots in gene_to_proteins.values() if len(prots) > 1)
    exact_n = sum(1 for v in match_type.values() if v == "exact")
    stripped_n = sum(1 for v in match_type.values() if v == "version_stripped")

    with open(out_summary, "w") as out:
        out.write(f"Total proteins:\t{total}\n")
        out.write(f"Mapped proteins:\t{mapped}\n")
        out.write(f"Unmapped proteins:\t{unmapped}\n")
        out.write(f"Exact matches:\t{exact_n}\n")
        out.write(f"Version-stripped matches:\t{stripped_n}\n")
        out.write(f"Unique GeneIDs:\t{genes}\n")
        out.write(f"Genes with isoforms:\t{iso_genes}\n")
        out.write(f"Proteins in isoform genes:\t{iso_prots}\n")


def run_list_mode(protein_file, gene_file, prefix):
    """
    Execute list mode workflow.
    """
    proteins = read_protein_list(protein_file)
    mapping, match_type = build_gene_mapping(proteins, gene_file)
    gene_to_proteins = build_gene_to_proteins(mapping)

    write_list_mapping_table(proteins, mapping, match_type, f"{prefix}.protein_to_gene.tsv")
    write_list_isoform_groups(gene_to_proteins, f"{prefix}.isoform_groups.tsv")
    write_list_groups_report(proteins, mapping, f"{prefix}.groups.tsv")
    write_list_unmapped(proteins, mapping, f"{prefix}.unmapped.txt")
    write_list_summary(proteins, mapping, match_type, gene_to_proteins, f"{prefix}.summary.txt")


# =========================================================
# Table mode (-t)
# =========================================================

def read_subject_ids_from_table(path):
    """
    Read unique sseqid values from a BLAST/DIAMOND-style table.

    Also validates that the expected header columns are present.
    """
    status(f"Reading subject accessions from table: {path}")
    sseqids = []
    seen = set()

    with open_maybe_gzip(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input table appears empty or has no header.")
        missing = [c for c in EXPECTED_HEADER if c not in reader.fieldnames]
        if missing:
            raise ValueError(
                "Input table is missing expected columns:\n  " + ", ".join(missing)
            )

        for row in reader:
            sseqid = row["sseqid"].strip()
            if not sseqid:
                continue
            if sseqid not in seen:
                seen.add(sseqid)
                sseqids.append(sseqid)

    status(f"Loaded {len(sseqids)} unique sseqid values from table")
    return sseqids


def read_full_table(path, mapping, match_type):
    """
    Read full table and annotate each row with mapping metadata.

    Added columns
    -------------
    GeneID
    Match_type
    Selected
    Selection_reason
    Representative_of
    Isoform_signature

    Isoform_signature is used for unmapped redundancy heuristics:
    staxids|slen|sstart|send|scovhsp
    """
    status(f"Reading full table and annotating rows: {path}")
    rows = []

    with open_maybe_gzip(path) as f:
        reader = csv.DictReader(f, delimiter="\t")

        for idx, row in enumerate(reader, 1):
            sseqid = row["sseqid"].strip()
            staxids = row["staxids"].strip()
            gene_id = mapping.get(sseqid, "NA")
            mtype = match_type.get(sseqid, "unmapped")

            row["_row_index"] = idx
            row["GeneID"] = gene_id
            row["Match_type"] = mtype
            row["Selected"] = "no"
            row["Selection_reason"] = "not_selected"
            row["Representative_of"] = ""
            row["Isoform_signature"] = (
                f"{staxids}|{row['slen']}|{row['sstart']}|{row['send']}|{row['scovhsp']}"
            )
            rows.append(row)

    status(f"Annotated {len(rows)} rows")
    return rows


def representative_key(row):
    """
    Ranking for choosing a representative inside a mapped gene group.

    Priority
    --------
    1. larger slen
    2. exact mapping over version_stripped
    3. first occurrence in file
    """
    slen = safe_int(row["slen"], 0)
    match_rank = 0 if row["Match_type"] == "exact" else 1
    return (-slen, match_rank, row["_row_index"])


def select_representatives(rows):
    """
    Select representative rows in table mode.

    Logic
    -----
    Mapped rows:
      - grouped by (staxids, GeneID)
      - if only one row -> selected as unique
      - if multiple rows -> longest slen wins (ties broken as above)

    Unmapped rows:
      - if only one unmapped row for a taxid -> selected as unmatched,
        unless an identical signature is already represented
      - remaining unmapped rows are grouped by Isoform_signature
      - first occurrence becomes unmatched_signature representative,
        unless identical signature already represented

    Returns
    -------
    all_rows : list[dict]
        Annotated rows, including selected/not-selected labels
    selected_rows : list[dict]
        Only rows with Selected == yes
    """
    mapped_by_tax_gene = defaultdict(list)
    unmapped_by_taxid = defaultdict(list)

    for row in rows:
        if row["GeneID"] != "NA":
            key = (row["staxids"].strip(), row["GeneID"])
            mapped_by_tax_gene[key].append(row)
        else:
            key = row["staxids"].strip()
            unmapped_by_taxid[key].append(row)

    # First process mapped rows.
    for key, group in mapped_by_tax_gene.items():
        group_sorted = sorted(group, key=representative_key)
        winner = group_sorted[0]

        if len(group) == 1:
            winner["Selected"] = "yes"
            winner["Selection_reason"] = "unique"
        else:
            winner["Selected"] = "yes"
            winner["Selection_reason"] = "biggest"

        winner["Representative_of"] = f"{key[0]}|{key[1]}"

        for loser in group_sorted[1:]:
            loser["Selection_reason"] = f"same_gene_as_selected:{winner['sseqid']}"

    # Track already selected signatures to avoid duplicate unmapped reps.
    selected_signature_to_row = {}
    for row in rows:
        if row["Selected"] == "yes":
            sig = row["Isoform_signature"]
            if sig not in selected_signature_to_row:
                selected_signature_to_row[sig] = row

    # Unmapped rows that are unique for their taxid.
    for staxid, group in unmapped_by_taxid.items():
        if len(group) == 1:
            row = group[0]
            sig = row["Isoform_signature"]

            if sig in selected_signature_to_row:
                row["Selection_reason"] = (
                    f"same_isoform_signature_as_selected:{selected_signature_to_row[sig]['sseqid']}"
                )
            else:
                row["Selected"] = "yes"
                row["Selection_reason"] = "unmatched"
                row["Representative_of"] = f"{staxid}|UNMATCHED"
                selected_signature_to_row[sig] = row

    # Remaining unmapped rows grouped by signature.
    unmapped_by_signature = defaultdict(list)
    for row in rows:
        if row["GeneID"] == "NA" and row["Selected"] != "yes":
            unmapped_by_signature[row["Isoform_signature"]].append(row)

    for sig, group in unmapped_by_signature.items():
        if sig in selected_signature_to_row:
            rep = selected_signature_to_row[sig]
            for row in group:
                row["Selection_reason"] = f"same_isoform_signature_as_selected:{rep['sseqid']}"
            continue

        winner = min(group, key=lambda r: r["_row_index"])
        winner["Selected"] = "yes"
        winner["Selection_reason"] = "unmatched_signature"
        winner["Representative_of"] = f"{winner['staxids']}|SIGNATURE"
        selected_signature_to_row[sig] = winner

        for row in group:
            if row["_row_index"] != winner["_row_index"]:
                row["Selection_reason"] = f"same_isoform_signature_as_selected:{winner['sseqid']}"

    selected_rows = [r for r in rows if r["Selected"] == "yes"]
    return rows, selected_rows


def sort_table_rows_for_output(rows):
    """
    Sort table-mode rows for cleaner inspection.

    Primary keys:
    - taxid
    - GeneID
    - sseqid
    - original row index
    """
    return sorted(
        rows,
        key=lambda r: (
            r["staxids"].strip(),
            gene_sort_value(r["GeneID"]),
            r["sseqid"],
            r["_row_index"]
        )
    )


def write_annotated_table(rows, out_path):
    """
    Write full annotated table in table mode.
    """
    status(f"Writing annotated table: {out_path}")

    out_fields = EXPECTED_HEADER + [
        "GeneID", "Match_type", "Isoform_signature",
        "Selected", "Selection_reason", "Representative_of"
    ]

    rows = sort_table_rows_for_output(rows)

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=out_fields,
            delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_filtered_table(rows, out_path):
    """
    Write only selected representative rows in table mode.
    """
    status(f"Writing filtered representatives: {out_path}")

    out_fields = EXPECTED_HEADER + [
        "GeneID", "Match_type", "Isoform_signature",
        "Selected", "Selection_reason", "Representative_of"
    ]

    rows = sort_table_rows_for_output(rows)

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=out_fields,
            delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_table_groups_report(all_rows, selected_rows, out_path):
    """
    Write grouping report for table mode.

    Includes:
    - mapped_gene groups: grouped by (TaxID, GeneID)
    - unmapped_signature groups: clustered unmapped proteins
    - unmapped_single groups: isolated unmapped proteins

    This report is especially useful to inspect which proteins were collapsed
    together and which representative was retained.
    """
    status(f"Writing group report: {out_path}")

    mapped_groups = defaultdict(list)
    unmapped_signature_groups = defaultdict(list)

    for row in all_rows:
        sseqid = row["sseqid"]
        staxid = row["staxids"].strip()
        gene_id = row["GeneID"]
        sig = row["Isoform_signature"]

        if gene_id != "NA":
            mapped_groups[(staxid, gene_id)].append(row)
        else:
            unmapped_signature_groups[(staxid, sig)].append(row)

    with open(out_path, "w") as out:
        out.write(
            "Group_type\tTaxID\tGroup_id\tN_proteins\tProteins\t"
            "Selected_representative\tSelection_reason\n"
        )

        mapped_keys = sorted(
            mapped_groups.keys(),
            key=lambda x: (safe_int(x[0], 0), gene_sort_value(x[1]))
        )

        for staxid, gene_id in mapped_keys:
            group_rows = mapped_groups[(staxid, gene_id)]
            proteins = sorted({r["sseqid"] for r in group_rows})

            selected_rep = "NA"
            selected_reason = "NA"
            for row in group_rows:
                if row["Selected"] == "yes":
                    selected_rep = row["sseqid"]
                    selected_reason = row["Selection_reason"]
                    break

            out.write(
                f"mapped_gene\t{staxid}\t{gene_id}\t{len(proteins)}\t"
                f"{','.join(proteins)}\t{selected_rep}\t{selected_reason}\n"
            )

        unmapped_keys = sorted(
            unmapped_signature_groups.keys(),
            key=lambda x: (safe_int(x[0], 0), x[1])
        )

        for staxid, sig in unmapped_keys:
            group_rows = unmapped_signature_groups[(staxid, sig)]
            proteins = sorted({r["sseqid"] for r in group_rows})

            selected_rep = "NA"
            selected_reason = "NA"
            for row in group_rows:
                if row["Selected"] == "yes":
                    selected_rep = row["sseqid"]
                    selected_reason = row["Selection_reason"]
                    break

            if len(proteins) == 1:
                group_type = "unmapped_single"
            else:
                group_type = "unmapped_signature"

            out.write(
                f"{group_type}\t{staxid}\t{sig}\t{len(proteins)}\t"
                f"{','.join(proteins)}\t{selected_rep}\t{selected_reason}\n"
            )


def write_table_summary(all_rows, selected_rows, out_path):
    """
    Write summary statistics for table mode.
    """
    status(f"Writing: {out_path}")

    total_rows = len(all_rows)
    total_selected = len(selected_rows)

    mapped_rows = sum(1 for r in all_rows if r["GeneID"] != "NA")
    unmapped_rows = total_rows - mapped_rows

    selected_mapped = sum(1 for r in selected_rows if r["GeneID"] != "NA")
    selected_unmapped = total_selected - selected_mapped

    exact_matches = sum(1 for r in all_rows if r["Match_type"] == "exact")
    stripped_matches = sum(1 for r in all_rows if r["Match_type"] == "version_stripped")

    reason_counts = defaultdict(int)
    for r in selected_rows:
        reason_counts[r["Selection_reason"]] += 1

    unique_genes_selected = len({
        (r["staxids"].strip(), r["GeneID"])
        for r in selected_rows
        if r["GeneID"] != "NA"
    })

    with open(out_path, "w") as f:
        f.write(f"Total rows:\t{total_rows}\n")
        f.write(f"Mapped rows:\t{mapped_rows}\n")
        f.write(f"Unmapped rows:\t{unmapped_rows}\n")
        f.write(f"Exact matches:\t{exact_matches}\n")
        f.write(f"Version-stripped matches:\t{stripped_matches}\n")
        f.write(f"Selected representatives:\t{total_selected}\n")
        f.write(f"Selected mapped representatives:\t{selected_mapped}\n")
        f.write(f"Selected unmapped representatives:\t{selected_unmapped}\n")
        f.write(f"Unique mapped gene groups selected:\t{unique_genes_selected}\n")
        f.write("\nSelected representatives by reason:\n")
        for reason in sorted(reason_counts):
            f.write(f"{reason}\t{reason_counts[reason]}\n")


def run_table_mode(table_file, gene_file, prefix):
    """
    Execute table mode workflow.
    """
    sseqids = read_subject_ids_from_table(table_file)
    mapping, match_type = build_gene_mapping(sseqids, gene_file)
    rows = read_full_table(table_file, mapping, match_type)
    all_rows, selected_rows = select_representatives(rows)

    write_annotated_table(all_rows, f"{prefix}.annotated.tsv")
    write_filtered_table(selected_rows, f"{prefix}.filtered.tsv")
    write_table_groups_report(all_rows, selected_rows, f"{prefix}.groups.tsv")
    write_table_summary(all_rows, selected_rows, f"{prefix}.summary.txt")


# =========================================================
# CLI
# =========================================================

def parse_args():
    """
    Parse command-line arguments.

    The script requires exactly one input mode:
    - -l for protein list mode
    - -t for table mode

    And always requires:
    - -g2a for the NCBI mapping file
    """
    parser = argparse.ArgumentParser(
        prog="Protein_isoforms_collapsing_from_proteinIDs.py",
        description=(
            "Map protein accessions to GeneIDs and collapse likely redundant isoforms.\n\n"
            "This tool supports two modes:\n"
            "  1) List mode (-l): map one accession per line to GeneIDs and report groups\n"
            "  2) Table mode (-t): annotate and collapse redundancy in BLAST/DIAMOND-style results\n"
        ),
        epilog=(
            "Examples:\n"
            "  List mode:\n"
            "    python Protein_isoforms_collapsing_from_proteinIDs.py \\\n"
            "      -l proteins.txt \\\n"
            "      -g2a gene2accession.gz\n\n"
            "  List mode with custom prefix:\n"
            "    python Protein_isoforms_collapsing_from_proteinIDs.py \\\n"
            "      -l proteins.txt \\\n"
            "      -g2a gene2refseq.gz \\\n"
            "      -prefix my_run\n\n"
            "  Table mode:\n"
            "    python Protein_isoforms_collapsing_from_proteinIDs.py \\\n"
            "      -t hits.tsv \\\n"
            "      -g2a gene2accession.gz\n\n"
            "  Table mode with custom prefix:\n"
            "    python Protein_isoforms_collapsing_from_proteinIDs.py \\\n"
            "      -t hits.tsv \\\n"
            "      -g2a gene2accession.gz \\\n"
            "      -prefix charon_DMD_nr_refseq_FINAL\n\n"
            "Input requirements:\n"
            "  -l : plain text file with one protein accession per line\n"
            "  -t : tab-separated table with this header:\n"
            "       qseqid qlen qstart qend length qcovhsp pident sseqid slen sstart send\n"
            "       scovhsp evalue bitscore qstrand qframe stitle staxids slineages\n"
            "       sscinames skingdoms sphylums db_label\n\n"
            "  -g2a : NCBI gene2accession or gene2refseq file (.gz or plain text)\n"
            "         Examples:\n"
            "         https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz\n"
            "         https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz\n\n"
            "Matching behavior:\n"
            "  - exact accession -> GeneID match is attempted first\n"
            "  - if exact match fails, version-stripped matching is attempted\n"
            "    Example: XP_015836279.1 and XP_015836279.2 both fallback to XP_015836279\n\n"
            "Outputs in list mode (-l):\n"
            "  <prefix>.protein_to_gene.tsv\n"
            "  <prefix>.isoform_groups.tsv\n"
            "  <prefix>.groups.tsv\n"
            "  <prefix>.unmapped.txt\n"
            "  <prefix>.summary.txt\n\n"
            "Outputs in table mode (-t):\n"
            "  <prefix>.annotated.tsv\n"
            "  <prefix>.filtered.tsv\n"
            "  <prefix>.groups.tsv\n"
            "  <prefix>.summary.txt\n\n"
            "Notes:\n"
            "  - Outside RefSeq, many proteins do not have GeneID mappings.\n"
            "  - In table mode, unmapped proteins are heuristically grouped using:\n"
            "      staxids + slen + sstart + send + scovhsp\n"
            "  - This is intended to reduce redundancy, not to define true transcript isoforms.\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-l",
        help="Protein list file, one accession per line"
    )
    group.add_argument(
        "-t",
        help="BLAST/DIAMOND-style input table"
    )

    parser.add_argument(
        "-g2a",
        required=True,
        help="NCBI gene2accession or gene2refseq file (.gz or plain text)"
    )
    parser.add_argument(
        "-prefix",
        required=False,
        help="Output prefix; default is inferred from the input filename"
    )

    return parser.parse_args()


def main():
    """
    Main entry point.
    """
    args = parse_args()

    input_path = args.l if args.l else args.t
    prefix = args.prefix if args.prefix else infer_prefix(input_path)

    status(f"Using output prefix: {prefix}")
    status("Starting program")

    if args.l:
        status("Running in list mode (-l)")
        run_list_mode(args.l, args.g2a, prefix)
    else:
        status("Running in table mode (-t)")
        run_table_mode(args.t, args.g2a, prefix)

    status("Done")


if __name__ == "__main__":
    main()
