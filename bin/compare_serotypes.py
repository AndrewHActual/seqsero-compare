#!/usr/bin/env python3
"""
Parse SeqSero2 SeqSero_result.tsv files from three sources and build a
per-sample comparison of O group, phase-1 H (H1), and phase-2 H (H2) calls.

SeqSero2 TSV columns of interest:
  - Sample name
  - O antigen prediction
  - H1 antigen prediction(fliC)
  - H2 antigen prediction(fljB)
"""
import argparse
import sys
from pathlib import Path
import pandas as pd

# Map the antigen fields we care about to likely SeqSero2 column names.
COLUMN_ALIASES = {
    "sample": ["Sample name", "Sample", "sample"],
    "O":      ["O antigen prediction", "O antigen"],
    "H1":     ["H1 antigen prediction(fliC)", "H1 antigen prediction", "H1 antigen"],
    "H2":     ["H2 antigen prediction(fljB)", "H2 antigen prediction", "H2 antigen"],
}


def pick_column(df, aliases):
    for name in aliases:
        if name in df.columns:
            return name
    return None


def parse_result(tsv_path):
    """Return dict {sample, O, H1, H2} from a single SeqSero2 result TSV."""
    df = pd.read_csv(tsv_path, sep="\t", dtype=str).fillna("")
    df.columns = [c.strip() for c in df.columns]
    if df.empty:
        return None

    row = df.iloc[0]
    result = {}
    for key, aliases in COLUMN_ALIASES.items():
        col = pick_column(df, aliases)
        val = row[col].strip() if col else ""
        # Normalize empty / "N/A" style values so comparisons are fair.
        if val in ("", "N/A", "NA", "-", "* No prediction"):
            val = "NA"
        result[key] = val

    # Fall back to filename stem if sample name column was missing/blank.
    if result.get("sample", "NA") in ("", "NA"):
        result["sample"] = Path(tsv_path).stem.split("_")[0]
    return result


def load_dir(directory):
    """Load all *.tsv in a directory keyed by sample id."""
    out = {}
    for tsv in sorted(Path(directory).glob("*.tsv")):
        parsed = parse_result(tsv)
        if parsed:
            out[parsed["sample"]] = parsed
    return out


def match_label(values):
    """Given a list of calls, return MATCH / MISMATCH ignoring NA-only sets."""
    real = [v for v in values if v != "NA"]
    if len(set(real)) <= 1 and len(real) == len(values):
        return "MATCH"
    if len(set(real)) <= 1 and real:
        return "MATCH_WITH_MISSING"
    return "MISMATCH"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amplicon-dir",  required=True)
    ap.add_argument("--assembly-dir",  required=True)
    ap.add_argument("--simulated-dir", required=True)
    ap.add_argument("--out-tsv",       required=True)
    ap.add_argument("--out-csv",       required=True)
    ap.add_argument("--out-mqc",       required=True,
                    help="MultiQC custom-content table (*_mqc.tsv)")
    args = ap.parse_args()

    amplicon  = load_dir(args.amplicon_dir)
    assembly  = load_dir(args.assembly_dir)
    simulated = load_dir(args.simulated_dir)

    samples = sorted(set(amplicon) | set(assembly) | set(simulated))
    if not samples:
        sys.exit("No SeqSero2 results found in any source directory.")

    rows = []
    for s in samples:
        a  = amplicon.get(s,  {})
        b  = assembly.get(s,  {})
        c  = simulated.get(s, {})

        row = {"sample": s}
        for antigen in ("O", "H1", "H2"):
            av = a.get(antigen, "NA")
            bv = b.get(antigen, "NA")
            cv = c.get(antigen, "NA")
            row[f"{antigen}_amplicon"]  = av
            row[f"{antigen}_assembly"]  = bv
            row[f"{antigen}_simulated"] = cv
            row[f"{antigen}_match"]     = match_label([av, bv, cv])

        # Overall concordance across all three antigens.
        overall = [row[f"{ag}_match"] for ag in ("O", "H1", "H2")]
        if all(m == "MATCH" for m in overall):
            row["overall"] = "ALL_MATCH"
        elif any(m == "MISMATCH" for m in overall):
            row["overall"] = "MISMATCH"
        else:
            row["overall"] = "MATCH_WITH_MISSING"
        rows.append(row)

    cols = ["sample"]
    for ag in ("O", "H1", "H2"):
        cols += [f"{ag}_amplicon", f"{ag}_assembly", f"{ag}_simulated", f"{ag}_match"]
    cols += ["overall"]

    out = pd.DataFrame(rows)[cols]
    out.to_csv(args.out_tsv, sep="\t", index=False)
    out.to_csv(args.out_csv, index=False)
    print(out.to_string(index=False))

    # ---- MultiQC custom-content table ----
    # A leading YAML header comment tells MultiQC how to render this as a
    # custom table in its own section. The '# plot_type: table' + config
    # block is parsed by MultiQC's custom-content module.
    mqc_header = [
        "# id: 'seqsero2_comparison'",
        "# section_name: 'SeqSero2 Serotype Comparison'",
        "# description: 'O group, phase-1 H (fliC) and phase-2 H (fljB) calls "
        "from amplicon reads, assembly, and simulated amplicons, with per-antigen "
        "concordance flags.'",
        "# format: 'tsv'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'seqsero2_comparison_table'",
        "#     title: 'SeqSero2 Serotype Concordance'",
        "#     col1_header: 'Sample'",
    ]
    with open(args.out_mqc, "w") as fh:
        fh.write("\n".join(mqc_header) + "\n")
        out.to_csv(fh, sep="\t", index=False)


if __name__ == "__main__":
    main()