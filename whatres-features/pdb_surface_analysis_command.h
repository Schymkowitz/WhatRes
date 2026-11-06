// pdb_surface_analysis_command.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// Command: "pdb surface analysis"
//
// Reads a list of PDB files and runs SurfacePatcher on every chain in each
// structure. Outputs one TSV row per chain containing:
//   - global surface statistics (exposure fraction, mean burial, net charge)
//   - metrics for the largest hydrophobic patch
//   - metrics for the largest positive charge cluster
//   - metrics for the largest negative charge cluster
//
// When surface_verbose is true an additional per-patch detail file is also
// written.
//
// Options (Options.txt prefix: surface_*):
//   surface_burial_radius      Cβ-neighbour sphere radius in Å       (default 10.0)
//   surface_burial_threshold   exposure cutoff (neighbours < this)    (default 16)
//   surface_adj_cutoff         Cβ–Cβ adjacency distance in Å         (default  8.0)
//   surface_min_patch_size     minimum residues to report a patch     (default  2)
//   surface_hphob_set          hydrophobic residue alphabet           (default "AVILMFWYPCG")
//   surface_verbose            also write a per-patch detail file     (default false)
//
// Cavity detection (prefix: cavity_*):
//   cavity_grid_spacing        voxel edge length in Å                 (default  1.0)
//   cavity_probe_radius        solvent probe radius in Å              (default  1.4)
//   cavity_min_volume          minimum cavity volume in ų            (default 10.0)
//   cavity_lining_dist         centroid–atom cutoff for lining list   (default  3.5)
//
// Interface void analysis (prefix: interface_*):
//   interface_grid_spacing     voxel edge for interface box in Å      (default  0.5)
//   interface_contact_dist     inter-chain cutoff for interface res    (default  5.0)
//   interface_shell_dist       void voxel needs both chains within    (default  4.0)
//   interface_box_padding      padding around interface box in Å      (default  3.0)

#pragma once
#include "commands.h"
#include "optionclass.h"

class pdb_surface_analysis_command : public sequence_command
{
public:
    pdb_surface_analysis_command(Option_class& O)
        : sequence_command("pdb_surface_analysis",
                           "Surface patch and charge-cluster analysis of PDB structures "
                           "for aggregation / developability scoring", O)
    {}

    void Execute() override;

    std::string HelpText() override
    {
        return R"(Analyse the molecular surface of every chain in a set of PDB
structures - hydrophobic patches and charge clusters - for aggregation
and antibody-developability scoring.

Context
  Large exposed hydrophobic patches and strong surface charge clusters
  are classic aggregation and developability liabilities.

Per chain it runs SurfacePatcher and reports one TSV row with:
  - global surface statistics: exposure fraction, mean burial, net charge
  - the largest hydrophobic patch
  - the largest positive charge cluster
  - the largest negative charge cluster
It can additionally run buried-cavity detection and inter-chain
interface-void analysis.

Input
  -in   a .txt file listing one PDB path per line.

Output
  <prefix>_surface_analysis.txt  one row per chain
  <prefix>_cavities.txt          buried cavities, one row per structure
  <prefix>_interface_voids.txt   inter-chain interface voids
  with surface_verbose true, also <prefix>_surface_patches.txt (per patch)
  and <prefix>_cavity_details.txt (per cavity)

Key Options.txt families
  surface_*    patch detection - burial radius / threshold, Cb-Cb
               adjacency cutoff, minimum patch size, hydrophobic
               alphabet, verbose flag
  cavity_*     buried-cavity detection - grid spacing, probe radius,
               minimum volume, lining distance
  interface_*  interface-void analysis - grid spacing, contact and shell
               distances, box padding
  (defaults in optionclass.cpp; example values in Options_example.txt))";
    }
};
