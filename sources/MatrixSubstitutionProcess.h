
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#ifndef MATRIXSUB_H
#define MATRIXSUB_H

#include <unordered_map>

#include "SubstitutionProcess.h"
#include "MatrixProfileProcess.h"

class MatrixSubstitutionProcess : public virtual SubstitutionProcess, public virtual MatrixProfileProcess	{

	public:

	MatrixSubstitutionProcess() : propagate_aux(0), propagate_aux_size(0),
	                              sitepropagate_aux(0), sitepropagate_aux_size(0),
	                              expdiag_aux(0), expdiag_aux_size(0) {}
	virtual ~MatrixSubstitutionProcess() {
		delete[] propagate_aux;
		delete[] sitepropagate_aux;
		delete[] expdiag_aux;
	}

	virtual int GetNstate(int site) {return GetMatrix(site)->GetNstate();}
	// virtual int GetNstate() {return GetMatrix(0)->GetNstate();}

	virtual const double* GetStationary(int site)	{
		return GetMatrix(site)->GetStationary();
	}

	protected:

	// CPU Level 3: implementations of likelihood propagation and substitution mapping methods
	void Propagate(double*** from, double*** to, double time, bool condalloc = false);

	void SitePropagate(int site, double** from, double** to, double time, bool condalloc = false);

	BranchSitePath** SamplePaths(int* stateup, int* statedown, double time);
	BranchSitePath** SampleRootPaths(int* rootstate);
	BranchSitePath* ResampleAcceptReject(int maxtrial, int stateup, int statedown, double rate, double totaltime, SubMatrix* matrix);
	BranchSitePath* ResampleUniformized(int stateup, int statedown, double rate, double totaltime, SubMatrix* matrix);

	void SimuPropagate(int* stateup, int* statedown, double time);

	// Persistent scratch buffers for Propagate/SitePropagate.
	// Replaces a per-call new[]/delete[] pair (one O(Nsite*Nrate*Nstate) and
	// one O(Nstate)) that was the dominant heap-allocator pressure during
	// likelihood evaluation. Lazily sized on first use, never shrunk.
	double* propagate_aux;
	size_t  propagate_aux_size;
	double* sitepropagate_aux;
	int     sitepropagate_aux_size;

	// Per-call cache of exp(length * eigenval[k]) keyed by (matrix, rate cat).
	// In SumOverRateAllocations mode, length = time * rate[j] is invariant
	// across sites that share a matrix, so each site re-evaluates the same
	// nstate exp() calls. We dedupe by SubMatrix* identity within a single
	// Propagate call. Reused across calls (cleared at the top of Propagate),
	// never freed mid-run. matrix_to_slot is also persistent: clear() reuses
	// bucket storage across calls.
	double* expdiag_aux;
	size_t  expdiag_aux_size;
	std::unordered_map<SubMatrix*, int> matrix_to_slot;
};

#endif

