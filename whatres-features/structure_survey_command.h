// structure_survey_command.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// Command: "structure_survey"
//
// Processes a large set of PDB structures (e.g. 25k high-resolution set) and
// builds a binary reference database of per-residue structural metric
// distributions, binned by amino acid type × side-chain burial.
//
// For each PDB:
//   1. NativeQCInputs        → phi/psi/omega and FoldX-style side-chain
//                              burial, from the coordinates
//   2. BondGeometry          → covalent bond lengths and angles
//   3. SurfacePatcher        → surface patches, charge clusters (per-chain)
//   4. CavityDetector        → buried cavities (per-structure)
//
// Output:
//   - structure_reference.bin      binary reference database
//   - structure_survey.tsv         combined per-residue TSV (for inspection)
//   - structure_survey_summary.txt   human-readable cell counts
//
// Usage:
//   whatres-features -do structure_survey -in pdb_list.txt -opt Options.txt

#pragma once
#include "commands.h"
#include "optionclass.h"

class structure_survey_command : public sequence_command
{
public:
    structure_survey_command(Option_class& O)
        : sequence_command("structure_survey",
                           "Build a structural-metric reference database from real PDB structures",
                           O)
    {}

    void Execute() override;

    std::string HelpText() override
    {
        return R"(Process a large set of PDB structures and build a binary reference
database of per-residue structural-metric distributions, binned by
amino-acid type x side-chain burial.

Why it exists
  Structure-prediction methods (AlphaFold and similar) appear to make
  *systematic* geometric errors. structure_survey was built to quantify
  that: collect distributions of structural metrics from real
  X-ray / cryo-EM structures across resolution thresholds, and
  from predicted models, so the systematic deviation of predicted
  structures from real ones can be measured. The reference database it
  produces is also consumed by `pdb_quality_check`, which scores how far
  a given structure departs from the real-structure norm.

Per structure it runs:
  NativeQCInputs        - phi/psi/omega and FoldX-style side-chain burial
                          (FoldX occupancy envelope), from the coordinates
  BondGeometry          - covalent bond lengths and angles
  SurfacePatcher        - surface patches and charge clusters (per chain)
  CavityDetector        - buried cavities (per structure)

Input
  -in   a .txt file listing one PDB path per line.

Output
  <prefix>structure_reference.bin       binary reference database
  <prefix>structure_survey.tsv          combined per-residue metrics
  <prefix>structure_survey_summary.txt  human-readable cell counts
  <prefix>structure_survey_cavities.tsv buried-cavity table
  <prefix>structure_survey_patches.tsv  surface-patch / charge-cluster table

Note
  The energy columns of the per-residue table are left empty.)";
    }
};
