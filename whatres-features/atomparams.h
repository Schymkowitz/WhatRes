// atomparams.h
// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz
//
// Interface to the compiled atom-parameter library (libatomparams.a): per-atom
// volume, minimal occupancy (Occ) and maximal occupancy (Occmax) of the FoldX
// 2014 parameter set, used for burial (NativeQCInputs.h via FxLookup in
// FxAtoms.h). The parameter values are FoldX's and are distributed only in
// compiled form, under the FoldX licence. A plain C interface, so one library
// links with any C++ compiler on its platform.
#pragma once

extern "C" {
// Number of parameter rows.
int  ap_n_atoms(void);
// Row i (0-based): residue (3-letter), atom name, volume (A^3), Occ, Occmax.
// Returns 0 if i is out of range.
int  ap_atom(int i, const char** res3, const char** atom,
              double* volume, double* occ, double* occmax);
// Parameters given to atoms that are not in the table.
void ap_fallback(double* volume, double* occ, double* occmax);
}
