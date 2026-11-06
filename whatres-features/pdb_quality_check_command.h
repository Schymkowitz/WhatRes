// pdb_quality_check_command.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// Command: "pdb_quality_check"
//
// Validates one or more PDB structures against the structural reference
// database built by the "structure_survey" command.
//
// For each PDB:
//   1. NativeQCInputs        → dihedrals and FoldX-style occupancy burial,
//                              from the coordinates
//   2. BondGeometry          → bond lengths and angles
//   3. CavityDetector        → cavity summary
//   4. SurfacePatcher        → surface patches
//   5. Look up each residue in the reference DB (AA × side-chain burial bin)
//   6. Flag outlier metrics (outside p1–p99 or user-defined thresholds)
//   7. Score residues and whole structure
//
// Output:
//   - <prefix>quality_report.tsv    per-residue report with percentiles
//   - <prefix>quality_summary.txt   per-PDB scores and overview
//   - <prefix>quality_outliers.tsv  only the flagged residues
//
// Scoring:
//   Per metric:   -log2(min(percentile, 1-percentile) / 50)
//                 → 0 at median, ~1 at p25/p75, ~4.3 at p5/p95, ~5.6 at p1/p99
//   Per residue:  sum of per-metric scores (worst metrics dominate)
//   Per PDB:      mean residue score + worst-5 penalty
//
// Usage:
//   whatres-features -do pdb_quality_check -in pdb_list.txt -opt Options.txt
//
// Requires:
//   <qc_reference_db>/path/to/structure_reference.bin; in the options file.

#pragma once
#include "commands.h"
#include "optionclass.h"

class pdb_quality_check_command : public sequence_command
{
public:
    pdb_quality_check_command(Option_class& O)
        : sequence_command("pdb_quality_check",
                           "Validate PDB structures against reference distributions",
                           O)
    {}

    void Execute() override;

    std::string HelpText() override
    {
        return R"(Score one or more PDB structures against the real-structure reference
database built by `structure_survey`, flagging residues whose local
geometry falls outside the normal range.

Use it for
  - benchmarking predicted models (AlphaFold and similar) against the
    distributions of real structures, surfacing their systematic
    geometric deviations;
  - quality-control of experimental structures before relying on them for
    modelling.

Per structure it computes the backbone dihedrals and a FoldX-style
side-chain burial (the FoldX occupancy envelope, NativeQCInputs.h) directly
from the coordinates, plus BondGeometry (bond lengths/angles), CavityDetector
and SurfacePatcher, then looks up each residue in the reference DB by
amino-acid x side-chain burial bin and flags metrics outside the p1-p99
range. The per-residue energy metrics are not scored.

Scoring
  per metric     -log2(min(pct, 1-pct) / 50)
                 -> 0 at the median, ~1 at p25/p75, ~4.3 at p5/p95,
                    ~5.6 at p1/p99
  per residue    sum of per-metric scores (worst metrics dominate)
  per structure  mean residue score + a worst-5 penalty

Input
  -in   a .txt file listing one PDB path per line.

Output
  <prefix>quality_report.tsv    per-residue report with percentiles
  <prefix>quality_summary.txt   per-PDB scores and overview
  <prefix>quality_outliers.tsv  only the flagged residues

Requires
  <qc_reference_db>/path/to/structure_reference.bin; in the options file
  (the database produced by `structure_survey`).)";
    }
};
