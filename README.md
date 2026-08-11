# seqsero-compare

An nf-core-style Nextflow pipeline that compares SeqSero2 *Salmonella*
serotype calls across three data sources for the same isolate and reports
where they agree and disagree.

## Overview

For each sample the pipeline runs SeqSero2 on:

1. **Amplicon reads** — real amplicon sequencing generated with a large
   primer panel (paired-end FASTQ).
2. **Assembled genome** — the isolate assembly FASTA.
3. **Simulated amplicons** — in-silico amplicons cut from the assembly with
   cutadapt using the same primer panel, then serotyped as reads.

All three converge on SeqSero2's `SeqSero_result.tsv`, and a comparison step
extracts the **O group**, **phase-1 H antigen (fliC)**, and
**phase-2 H antigen (fljB)** calls from each source, flags per-antigen
concordance, and rolls the results into a MultiQC report.

## Workflow

```
                 ┌───────────────────────┐
 amplicon reads ─► SEQSERO2_READS         ─┐
                 └───────────────────────┘ │
                 ┌───────────────────────┐ │
 assembly fasta ─► SEQSERO2_ASSEMBLY      ─┼─► COMPARE_SEROTYPES ─► MULTIQC
        │        └───────────────────────┘ │        ▲
        │        ┌───────────────────────┐ │        │
        └───────► CUTADAPT_SIMULATE ──────►SEQSERO2_SIMULATED
                 └───────────────────────┘   (simulated amplicons)
```

## Usage

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --primer_fasta primers.tsv \
    -profile docker
```

### Samplesheet (`--input`)

```csv
sample,amplicon_fastq_1,amplicon_fastq_2,assembly_fasta
SAL001,/data/SAL001_R1.fastq.gz,/data/SAL001_R2.fastq.gz,/data/SAL001.fasta
```

### Primer panel (`--primer_fasta`)

A tab-separated file, one amplicon per line: `name<TAB>forward<TAB>reverse`.
The reverse primer is given 5'→3' as ordered; the pipeline reverse-complements
it internally. Lines beginning with `#` are ignored.

## Outputs

Results are written to `results/run_<YYYY-MM-DD_HH-MM-SS>/`:

| Path | Contents |
|------|----------|
| `seqsero2/amplicon_reads/`     | SeqSero2 results, real amplicon reads |
| `seqsero2/assembly/`           | SeqSero2 results, assembly |
| `simulated_amplicons/`         | cutadapt-derived amplicon FASTAs + logs |
| `seqsero2/simulated_amplicon/` | SeqSero2 results, simulated amplicons |
| `comparison/`                  | `serotype_comparison.{tsv,csv}` |
| `multiqc/`                     | `multiqc_report.html` |

The comparison table has, per sample, the O/H1/H2 call from each source, a
per-antigen flag (`MATCH` / `MISMATCH` / `MATCH_WITH_MISSING`), and an
overall verdict.

## Output directory timestamping

Each run writes to a timestamped subdirectory so repeat runs never overwrite
each other. Control this with:

- `--outdir_base <dir>` — base directory; a run subfolder is appended
  (default `./results`).
- `--run_name <name>` — use a **stable** subfolder name instead of a
  timestamp. Pass the same `--run_name` with `-resume` to reuse the same
  output directory (see below).
- `--outdir <dir>` — pin an exact directory with no base/run_name logic.

### Resuming a run

By default the outdir contains a launch timestamp, which changes every run —
so a bare `-resume` would write to a *new* folder. To make `-resume` reuse the
same output directory, give the run a stable name:

```bash
# First run
nextflow run main.nf --input samplesheet.csv --primer_fasta primers.tsv \
    --run_name batch2026_08 -profile docker

# Resume the same run into the same folder
nextflow run main.nf --input samplesheet.csv --primer_fasta primers.tsv \
    --run_name batch2026_08 -profile docker -resume
```

Nextflow's work cache still lives in `./work`, so `-resume` skips completed
tasks regardless; `--run_name` just keeps the published outputs in one place.

## Testing

A bundled `test` profile runs the whole pipeline end-to-end on a tiny dataset:

```bash
nextflow run main.nf -profile test,docker
```

This caps resources for a laptop or CI runner and uses a stable `run_name` of
`test`. You must supply three small real *Salmonella* test files under
`assets/test_data/` (`TEST01.fasta`, `TEST01_R1.fastq.gz`,
`TEST01_R2.fastq.gz`) — see the test-data note in the repo. Random sequence
won't work, because SeqSero2 needs biologically real input to make a call.

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir_base` | `./results` | Parent directory for run outputs |
| `--run_name` | `null` | Stable run subfolder name (enables `-resume` reuse) |
| `--seqsero2_reads_mode` | `a` | SeqSero2 mode for reads (`a`=microassembly, `k`=k-mer) |
| `--cutadapt_error_rate` | `0.1` | Allowed primer mismatch rate |
| `--cutadapt_min_overlap`| `10` | Minimum primer overlap |
| `--cutadapt_min_len`    | `50` | Minimum amplicon length to keep |
| `--multiqc_title`       | `null` | Custom MultiQC report title |

## Design decisions & notes

**Three sources, one format.** All three code paths converge on SeqSero2's
`SeqSero_result.tsv`, so the comparison step parses a single consistent
format. The read-based module is reused for both real and simulated amplicons
via an aliased include (`SEQSERO2_READS as SEQSERO2_SIMULATED`), the idiomatic
nf-core way to avoid duplicating a process.

**SeqSero2 modes.** Reads default to `-m a` (allele microassembly), the
recommended workflow for read data; the assembly uses `-m k` (k-mer), which is
what SeqSero2 uses on FASTA input. Both are exposed as params so you can
override them — e.g. switch amplicons to `-m k` if microassembly struggles
with panel data.

**Amplicon simulation with cutadapt.** Simulation uses cutadapt's
linked-adapter mode (`-g FWD...REVCOMP(REV)`) with `--discard-untrimmed`, so
only sequence flanked by a matched primer pair survives — that is the
in-silico amplicon. A shell loop reverse-complements each reverse primer
(including IUPAC ambiguity codes) and builds one `-g` argument per pair.

> **Caveat:** plain cutadapt scans linearly and reports the first match per
> pair per sequence. If your panel has many primers, or you need *every* hit
> on a contig, consider running cutadapt iteratively or using a dedicated
> in-silico PCR tool. `seqkit amplicon` is purpose-built for this and worth
> considering; cutadapt is wired in here as requested.

**Concordance with missing-data awareness.** The comparison uses a
`MATCH` / `MISMATCH` / `MATCH_WITH_MISSING` scheme rather than a plain
match/no-match. Simulated or real amplicons may not recover an antigen the
assembly does, and you probably don't want that silently scored as a hard
mismatch — `MATCH_WITH_MISSING` marks agreement among the sources that made a
call while flagging that one source returned nothing.

**MultiQC.** The comparison script emits a custom-content table
(`*_mqc.tsv`) with an embedded YAML header, and cutadapt logs are picked up by
MultiQC's built-in cutadapt module, so the report shows both the serotype
concordance table and amplicon-recovery stats in one place.

## Things to verify for your setup

- **SeqSero2 column headers** have shifted slightly between versions. The
  parser tries several aliases, but confirm against your installed 1.3.x
  output.
- **Primers per amplicon.** The TSV format assumes one forward/reverse pair
  per amplicon; if your panel design differs, adjust the parser and the
  cutadapt module.

## Tools

SeqSero2 · cutadapt · MultiQC · Nextflow (DSL2)