
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "MatrixSubstitutionProcess.h"
#include "Random.h"

#include <cmath>
#include <iostream>
#include <vector>
using namespace std;


void MatrixSubstitutionProcess::SimuPropagate(int* stateup, int* statedown, double time)	{

	const int nstate = GetMatrix(sitemin)->GetNstate();
	double cumul[nstate];
	double expdiag[nstate];
	for(int i=sitemin; i<sitemax; i++)	{

		int up = stateup[i];

		SubMatrix* matrix = GetMatrix(i);
		double** eigenvect = matrix->GetEigenVect();
		double** inveigenvect = matrix->GetInvEigenVect();
		double* eigenval = matrix->GetEigenVal();

		int j = ratealloc[i];
		double length = time * GetRate(i,j);
		for (int k=0; k<nstate; k++)	{
			expdiag[k] = exp(length * eigenval[k]);
		}

		double totprob = 0;
		for (int k=0; k<nstate; k++)	{
			double tot = 0;
			for (int l=0; l<nstate; l++)	{
				tot += eigenvect[up][l] * expdiag[l] * inveigenvect[l][k];
			}
			totprob += tot;
			cumul[k] = totprob;
		}
		if (fabs(totprob - 1) > 1e-6)	{
			cerr << "error in MatrixSubstitutionProcess::SimuPropagate: tot prob is not 1\n";
			cerr << totprob << '\n';
			exit(1);
		}

		double u = rnd::GetRandom().Uniform();
		int k = 0;
		while ((k<nstate) && (u>cumul[k]))	{
			k++;
		}
		if (k == nstate)	{
			cerr << "error in MatrixSubstitutionProcess::SimuPropagate: overflow\n";
			exit(1);
		}
		statedown[i] = k;
	}
}

//-------------------------------------------------------------------------
//	(CPU level 3)
//
//	* conditional likelihood propagation
//
//	(CPU level 3)
//-------------------------------------------------------------------------

void MatrixSubstitutionProcess::Propagate(double*** from, double*** to, double time, bool condalloc)	{

	// propchrono.Start();
	int i,j,k,l,offset;
	double length,max,maxup;
	const int nstate = GetMatrix(sitemin)->GetNstate();
	const int nrate = GetNrate(0);
	// Lazily (re)allocate the persistent scratch buffer. Same size and
	// indexing as the previous per-call new[]; this is purely a storage
	// lift, no FP arithmetic changes.
	const size_t needed = (size_t) GetNsite() * nrate * nstate;
	if (needed > propagate_aux_size)	{
		delete[] propagate_aux;
		propagate_aux = new double[needed];
		propagate_aux_size = needed;
	}
	double* aux = propagate_aux;

	// Dedupe exp(length * eigenval[k]) across sites that share a matrix.
	// Only valid when GetRate(i,j) is independent of i for fixed j -- i.e.
	// when SumOverRateAllocations is true. Propagate is called only in that
	// regime by current callers (the collapsed/data-augmentation phase uses
	// SimuPropagate, not Propagate), so we expect can_dedupe to be true on
	// the hot path. If it ever isn't, we fall through to the inline form.
	const bool can_dedupe = SumOverRateAllocations();
	if (can_dedupe)	{
		matrix_to_slot.clear();
		int nmat = 0;
		for (int s=sitemin; s<sitemax; s++)	{
			if (ActiveSite(s))	{
				SubMatrix* m = GetMatrix(s);
				if (matrix_to_slot.find(m) == matrix_to_slot.end())	{
					matrix_to_slot[m] = nmat++;
				}
			}
		}
		const size_t exp_needed = (size_t) nmat * nrate * nstate;
		if (exp_needed > expdiag_aux_size)	{
			delete[] expdiag_aux;
			expdiag_aux = new double[exp_needed];
			expdiag_aux_size = exp_needed;
		}
		// Same expression as the inline form (length * eigenval[k]).
		// exp() is deterministic per platform libm, so the cached value
		// is bit-identical to the value the inline form would compute.
		for (auto& p : matrix_to_slot)	{
			SubMatrix* m = p.first;
			int mi = p.second;
			double* eigenval_m = m->GetEigenVal();
			double* base = expdiag_aux + (size_t) mi * nrate * nstate;
			for (int jj=0; jj<nrate; jj++)	{
				double lj = time * GetRate(0, jj);
				double* slot = base + jj * nstate;
				for (int kk=0; kk<nstate; kk++)	{
					slot[kk] = exp(lj * eigenval_m[kk]);
				}
			}
		}
	}

	for(i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            SubMatrix* matrix = GetMatrix(i);
            double** eigenvect = matrix->GetEigenVect();
            double** inveigenvect = matrix->GetInvEigenVect();
            double* eigenval = matrix->GetEigenVal();
            const double* expdiag_for_matrix = 0;
            if (can_dedupe)	{
                int mi = matrix_to_slot[matrix];
                expdiag_for_matrix = expdiag_aux + (size_t) mi * nrate * nstate;
            }
            for(j=0; j<GetNrate(i); j++)	{
                if ((!condalloc) || (ratealloc[i] == j))	{
                    double* up = from[i][j];
                    double* down = to[i][j];
                    //SubMatrix* matrix = GetMatrix(i);
                    length = time * GetRate(i,j);

                    //double** eigenvect = matrix->GetEigenVect();
                    //double** inveigenvect= matrix->GetInvEigenVect();
                    //double* eigenval = matrix->GetEigenVal();

                    // substitution matrix Q = P L P^{-1} where L is diagonal (eigenvalues) and P is the eigenvector matrix
                    // we need to compute 
                    // down = exp(length * Q) . up
                    // which we express as 
                    // down = P ( exp(length * L) . (P^{-1} . up) )  
            
                    // thus we successively do the following matrix.vector products

                    // P^{-1} . up  -> aux
                    // exp(length * L) . aux  -> aux 	(where exp(length*L) is diagonal, so this is linear)
                    // P . aux -> down

                    /*
                    int nstate = GetNstate();
                    double* aux = new double[nstate];
                    */
                    offset = nstate*(i*GetNrate(0) + j);
                    double* const aux_site = aux + offset;

                    // P^{-1} . up  -> aux_site
                    // Local accumulator + hoisted row pointer lets the
                    // compiler keep the running sum in an FP register and
                    // emit a single store per output element. Bit-identical
                    // to the original (same FP add order: ((0 + a0) + a1) + ...).
                    for(k=0; k<nstate; k++)	{
                        double sum = 0.0;
                        const double* row = inveigenvect[k];
                        for(l=0; l<nstate; l++)	{
                            sum += row[l] * up[l];
                        }
                        aux_site[k] = sum;
                    }

                    // exp(length * L) . aux_site  -> aux_site
                    // Read from cached expdiag table when available; the
                    // cached entry is exp(time * rate[j] * eigenval[k]),
                    // which equals exp(length * eigenval[k]) bit-for-bit
                    // since length = time * GetRate(i,j) = time * rate[j]
                    // when SumOverRateAllocations is true.
                    if (can_dedupe)	{
                        const double* expdiag_jk = expdiag_for_matrix + j * nstate;
                        for(k=0; k<nstate; k++)	{
                            aux_site[k] *= expdiag_jk[k];
                        }
                    }
                    else	{
                        for(k=0; k<nstate; k++)	{
                            aux_site[k] *= exp(length * eigenval[k]);
                        }
                    }

                    // P . aux_site -> down
                    for(k=0; k<nstate; k++)	{
                        double sum = 0.0;
                        const double* row = eigenvect[k];
                        for(l=0; l<nstate; l++)	{
                            sum += row[l] * aux_site[l];
                        }
                        down[k] = sum;
                    }

                    // exit in case of numerical errors
                    for(k=0; k<nstate; k++)	{
                        if (std::isnan(down[k]))	{
                            cerr << "error in back prop\n";
                            for(l=0; l<nstate; l++)	{
                                cerr << up[l] << '\t' << down[l] << '\t' << matrix->Stationary(l) << '\n';
                            }
                            exit(1);
                        }
                    }
                    maxup = 0.0;
                    for(k=0; k<nstate; k++)	{
                        if (up[k] < 0.0)	{
                            cerr << "error in backward propagate: negative prob : " << up[k] << "\n";
                            exit(1);
                        }
                        if (maxup < up[k])	{
                            maxup = up[k];
                        }
                    }
                    max = 0.0;
                    for(k=0; k<nstate; k++)	{
                        if (down[k] < 0.0)	{
                            infprobcount++;
                            down[k] = 0.0;
                        }
                        if (max < down[k])	{
                            max = down[k];
                        }
                    }
                    /*
                    if (maxup == 0.0)	{
                        cerr << "error in backward propagate: null up array\n";
                        cerr << "site : " << i << '\n';
                        for(l=0; l<nstate; l++)	{
                            cerr << matrix->Stationary(l) << '\n';
                        }
                        cerr << time << '\t' << length << '\n';
                        cerr << GetDim() << '\n';
                        cerr << GetMinStat(i) << '\n';
                        exit(1);
                    }
                    if (max == 0.0)	{
                        cerr << "error in backward propagate: null array\n";
                        for(k=0; k<nstate; k++)	{
                            cerr << up[k] << '\t' << down[k] << '\n';
                        }
                        cerr << length << '\n';
                        cerr << '\n';
                        exit(1);
                    }
                    */

                    // this is the offset (in log)
                    down[nstate] = up[nstate];
                }
            }
        }
    }
    // aux is the persistent member buffer; do not delete here.
}

// Fused Initialize+Propagate for leaf branches.  When leafstates[i] is a
// known state s (>= 0), up[] is the indicator at s, so P^{-1}*up reduces
// to column s of P^{-1} — skipping the first O(nstate^2) matvec entirely.
// When leafstates[i] == -1 (missing data), up is all-1s and we fall back
// to the full row-sum path (same cost as Propagate).  The result is
// bit-identical to Initialize(aux)+Propagate(aux,to) because:
//   x*0.0 = 0.0, x*1.0 = x, and 0.0+x = x  in IEEE 754,
// so the reduced sum over the indicator equals the exact column entry.
void MatrixSubstitutionProcess::PropagateTip(const int* leafstates, double*** to, double time, double*** /*aux*/, bool condalloc)	{

	int i,j,k,l;
	double length,max;
	const int nstate = GetMatrix(sitemin)->GetNstate();
	const int nrate = GetNrate(0);

	if ((size_t)nstate > propagate_aux_size)	{
		delete[] propagate_aux;
		propagate_aux = new double[nstate];
		propagate_aux_size = nstate;
	}

	const bool can_dedupe = SumOverRateAllocations();
	if (can_dedupe)	{
		matrix_to_slot.clear();
		int nmat = 0;
		for (int s=sitemin; s<sitemax; s++)	{
			if (ActiveSite(s))	{
				SubMatrix* m = GetMatrix(s);
				if (matrix_to_slot.find(m) == matrix_to_slot.end())	{
					matrix_to_slot[m] = nmat++;
				}
			}
		}
		const size_t exp_needed = (size_t) nmat * nrate * nstate;
		if (exp_needed > expdiag_aux_size)	{
			delete[] expdiag_aux;
			expdiag_aux = new double[exp_needed];
			expdiag_aux_size = exp_needed;
		}
		for (auto& p : matrix_to_slot)	{
			SubMatrix* m = p.first;
			int mi = p.second;
			double* eigenval_m = m->GetEigenVal();
			double* base = expdiag_aux + (size_t) mi * nrate * nstate;
			for (int jj=0; jj<nrate; jj++)	{
				double lj = time * GetRate(0, jj);
				double* slot = base + jj * nstate;
				for (int kk=0; kk<nstate; kk++)	{
					slot[kk] = exp(lj * eigenval_m[kk]);
				}
			}
		}
	}

	for(i=sitemin; i<sitemax; i++)	{
        if (ActiveSite(i))  {
            SubMatrix* matrix = GetMatrix(i);
            double** eigenvect = matrix->GetEigenVect();
            double** inveigenvect = matrix->GetInvEigenVect();
            double* eigenval = matrix->GetEigenVal();
            int state = leafstates[i];
            const double* expdiag_for_matrix = 0;
            if (can_dedupe)	{
                int mi = matrix_to_slot[matrix];
                expdiag_for_matrix = expdiag_aux + (size_t) mi * nrate * nstate;
            }
            for(j=0; j<GetNrate(i); j++)	{
                if ((!condalloc) || (ratealloc[i] == j))	{
                    double* down = to[i][j];
                    length = time * GetRate(i,j);

                    double* scratch = propagate_aux;

                    if (state == -1)	{
                        // Missing data: up is all-1s. Row sums of P^{-1}.
                        for(k=0; k<nstate; k++)	{
                            double sum = 0.0;
                            const double* row = inveigenvect[k];
                            for(l=0; l<nstate; l++)	{
                                sum += row[l];
                            }
                            scratch[k] = sum;
                        }
                    }
                    else	{
                        // Known state: column read replaces matvec.
                        for(l=0; l<nstate; l++)	{
                            scratch[l] = inveigenvect[l][state];
                        }
                    }

                    // exp(length * L) . scratch
                    if (can_dedupe)	{
                        const double* expdiag_jk = expdiag_for_matrix + j * nstate;
                        for(k=0; k<nstate; k++)	{
                            scratch[k] *= expdiag_jk[k];
                        }
                    }
                    else	{
                        for(k=0; k<nstate; k++)	{
                            scratch[k] *= exp(length * eigenval[k]);
                        }
                    }

                    // P . scratch -> down
                    for(k=0; k<nstate; k++)	{
                        double sum = 0.0;
                        const double* row = eigenvect[k];
                        for(l=0; l<nstate; l++)	{
                            sum += row[l] * scratch[l];
                        }
                        down[k] = sum;
                    }

                    for(k=0; k<nstate; k++)	{
                        if (std::isnan(down[k]))	{
                            cerr << "error in PropagateTip\n";
                            for(l=0; l<nstate; l++)	{
                                cerr << down[l] << '\t' << matrix->Stationary(l) << '\n';
                            }
                            exit(1);
                        }
                    }
                    max = 0.0;
                    for(k=0; k<nstate; k++)	{
                        if (down[k] < 0.0)	{
                            infprobcount++;
                            down[k] = 0.0;
                        }
                        if (max < down[k])	{
                            max = down[k];
                        }
                    }
                    down[nstate] = 0;
                }
            }
        }
    }
}

void MatrixSubstitutionProcess::SitePropagate(int i, double** from, double** to, double time, bool condalloc)	{

	// propchrono.Start();
	int j,k,l;
	double length,max,maxup;

	// should be dependent on site
	// const int nstate = GetMatrix(GetSiteMin())->GetNstate();

	const int needed_small = GetNstate(i);
	if (needed_small > sitepropagate_aux_size)	{
		delete[] sitepropagate_aux;
		sitepropagate_aux = new double[needed_small];
		sitepropagate_aux_size = needed_small;
	}
	double* aux = sitepropagate_aux;

	SubMatrix* matrix = GetMatrix(i);
	double** eigenvect = matrix->GetEigenVect();
	double** inveigenvect = matrix->GetInvEigenVect();
	double* eigenval = matrix->GetEigenVal();

	int nstate = matrix->GetNstate();

	for(j=0; j<GetNrate(i); j++)	{

		if ((!condalloc) || (ratealloc[i] == j))	{

			double* up = from[j];
			double* down = to[j];

			length = time * GetRate(i,j);

			// substitution matrix Q = P L P^{-1} where L is diagonal (eigenvalues) and P is the eigenvector matrix
			// we need to compute 
			// down = exp(length * Q) . up
			// which we express as 
			// down = P ( exp(length * L) . (P^{-1} . up) )  
	
			// thus we successively do the following matrix.vector products

			// P^{-1} . up  -> aux
			// exp(length * L) . aux  -> aux 	(where exp(length*L) is diagonal, so this is linear)
			// P . aux -> down

			// P^{-1} . up  -> aux
			for(k=0; k<nstate; k++)	{
				double sum = 0.0;
				const double* row = inveigenvect[k];
				for(l=0; l<nstate; l++)	{
					sum += row[l] * up[l];
				}
				aux[k] = sum;
			}

			// exp(length * L) . aux  -> aux
			for(k=0; k<nstate; k++)	{
				aux[k] *= exp(length * eigenval[k]);
			}

			// P . aux -> down
			for(k=0; k<nstate; k++)	{
				double sum = 0.0;
				const double* row = eigenvect[k];
				for(l=0; l<nstate; l++)	{
					sum += row[l] * aux[l];
				}
				down[k] = sum;
			}

			// exit in case of numerical errors
			for(k=0; k<nstate; k++)	{
				if (isnan(down[k]))	{
					cerr << "error in back prop\n";
					for(l=0; l<nstate; l++)	{
						cerr << up[l] << '\t' << down[l] << '\t' << matrix->Stationary(l) << '\n';
					}
					exit(1);
				}
			}
			maxup = 0.0;
			for(k=0; k<nstate; k++)	{
				if (up[k] < 0.0)	{
					cerr << "error in backward propagate: negative prob : " << up[k] << "\n";
					exit(1);
				}
				if (maxup < up[k])	{
					maxup = up[k];
				}
			}
			max = 0.0;
			for(k=0; k<nstate; k++)	{
				if (down[k] < 0.0)	{
					infprobcount++;
					down[k] = 0.0;
				}
				if (max < down[k])	{
					max = down[k];
				}
			}
			/*
			if (maxup == 0.0)	{
				cerr << "error in backward propagate: null up array\n";
				cerr << "site : " << i << '\n';
				for(l=0; l<nstate; l++)	{
					cerr << matrix->Stationary(l) << '\n';
				}
				cerr << time << '\t' << length << '\n';
				exit(1);
			}
			if (max == 0.0)	{
				cerr << "error in backward propagate: null array\n";
				for(k=0; k<nstate; k++)	{
					cerr << up[k] << '\t' << down[k] << '\n';
				}
				cerr << length << '\n';
				cerr << '\n';
				exit(1);
			}
			*/
			// this is the offset (in log)
			down[nstate] = up[nstate];
		}
	}
	// aux is the persistent member buffer; do not delete here.
}

