
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/

#include <cassert>
#include "ExpoConjugateGTRPhyloProcess.h"
#include "Parallel.h"
#include <string.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* ExpoConjugateGTRPhyloProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void ExpoConjugateGTRPhyloProcess::CreateSuffStat()	{

	PhyloProcess::CreateSuffStat();
	if (siteprofilesuffstatcount)	{
		cerr << "error in ExpoConjugateGTRPhyloProcess::CreateSuffStat\n";
		cerr << myid << '\n';
		exit(1);
	}
	allocsiteprofilesuffstatcount = new int[GetNsite()*GetDim()];
	allocsiteprofilesuffstatbeta = new double[GetNsite()*GetDim()];
	tmpcount = new int[GetNsite()*GetDim()];
	tmpbeta = new double[GetNsite()*GetDim()];
	siteprofilesuffstatcount = new int*[GetNsite()];
	siteprofilesuffstatbeta = new double*[GetNsite()];
	// for (int i=sitemin; i<sitemax; i++)	{
	for (int i=0; i<GetNsite(); i++)	{
		siteprofilesuffstatcount[i] = allocsiteprofilesuffstatcount + i*GetDim();
		siteprofilesuffstatbeta[i] = allocsiteprofilesuffstatbeta + i*GetDim();
	}
}

void ExpoConjugateGTRPhyloProcess::DeleteSuffStat()	{

	if (siteprofilesuffstatcount)	{
		// for (int i=sitemin; i<sitemax; i++)	{
		/*
		for (int i=0; i<GetNsite(); i++)	{
			delete[] siteprofilesuffstatcount[i];
			delete[] siteprofilesuffstatbeta[i];
		}
		*/
		delete[] siteprofilesuffstatcount;
		delete[] siteprofilesuffstatbeta;
		siteprofilesuffstatcount = 0;
		siteprofilesuffstatbeta = 0;
		delete[] allocsiteprofilesuffstatcount;
		delete[] allocsiteprofilesuffstatbeta;
		allocsiteprofilesuffstatcount = 0;
		allocsiteprofilesuffstatbeta = 0;
		delete[] tmpcount;
		delete[] tmpbeta;
		tmpcount = 0;
		tmpbeta = 0;
	}
	PhyloProcess::DeleteSuffStat();
}

void ExpoConjugateGTRPhyloProcess::UpdateRRSuffStat()	{

	for (int k=0; k<GetNrr(); k++)	{
		rrsuffstatcount[k] = 0;
		rrsuffstatbeta[k] = 0;
	}
	for (int j=1; j<GetNbranch(); j++)	{
		// AddRRSuffStat(rrsuffstatcount,rrsuffstatbeta,submap[j],blarray[j]);
		AddRRSuffStat(rrsuffstatcount,rrsuffstatbeta,submap[j],blarray[j],missingmap[j]);
	}
}

void ExpoConjugateGTRPhyloProcess::UpdateSiteRateSuffStat()	{

	// cerr << "in update site rate : " << GetTotalLength() << '\n';
	for (int i=sitemin; i<sitemax; i++)	{
		siteratesuffstatcount[i] = 0;
		siteratesuffstatbeta[i] = 0;
	}
	for (int j=1; j<GetNbranch(); j++)	{
		// AddSiteRateSuffStat(siteratesuffstatcount,siteratesuffstatbeta,submap[j],blarray[j]);
		AddSiteRateSuffStat(siteratesuffstatcount,siteratesuffstatbeta,submap[j],blarray[j],missingmap[j]);
	}
}

void ExpoConjugateGTRPhyloProcess::UpdateBranchLengthSuffStat()	{

	branchlengthsuffstatcount[0] = 0;
	branchlengthsuffstatbeta[0] = 0;
	for (int j=1; j<GetNbranch(); j++)	{
		int& count = branchlengthsuffstatcount[j];
		double& beta = branchlengthsuffstatbeta[j];
		count = 0;
		beta = 0;
		// AddBranchLengthSuffStat(count,beta,submap[j]);
		AddBranchLengthSuffStat(count,beta,submap[j],missingmap[j]);
	}
}

void ExpoConjugateGTRPhyloProcess::UpdateSiteProfileSuffStat()	{

	for (int i=sitemin; i<sitemax; i++)	{
		for (int k=0; k<GetDim(); k++)	{
			siteprofilesuffstatcount[i][k] = 0;
			siteprofilesuffstatbeta[i][k] = 0;
		}
	}
	for (int j=0; j<GetNbranch(); j++)	{
		// AddSiteProfileSuffStat(siteprofilesuffstatcount,siteprofilesuffstatbeta,submap[j],blarray[j], (j == 0));
		AddSiteProfileSuffStat(siteprofilesuffstatcount,siteprofilesuffstatbeta,submap[j],blarray[j],missingmap[j]);
	}
}

// Compute the (recvcounts, displs) layout for Allgatherv over the
// site-partition shared by master and slaves.  Master (rank 0) owns
// no sites and contributes 0 elements.  Each slave i (1..nprocs-1)
// owns [(i-1)*width, i*width) in flat units of GlobalNstate, with the
// last slave receiving any remainder.  Layout depends only on nprocs
// and GetNsite(), so it is identical on every rank.
static void compute_siteprofile_layout(int nprocs, int nsite, int nstate,
                                        int* recvcounts, int* displs) {
	int width = nsite / (nprocs - 1);
	recvcounts[0] = 0;
	displs[0] = 0;
	for (int i = 1; i < nprocs; ++i) {
		int smin = (i - 1) * width;
		int smax = (i == nprocs - 1) ? nsite : i * width;
		recvcounts[i] = (smax - smin) * nstate;
		displs[i] = smin * nstate;
	}
}

void ExpoConjugateGTRPhyloProcess::GlobalUpdateSiteProfileSuffStat()	{

	// MPI2
	// ask slaves to update siteprofilesuffstats then Allgatherv into the
	// master and back to all slaves in a single collective.  The previous
	// implementation used a serial Recv-per-slave gather + Bcast pair,
	// which forced O(nprocs) sequential copies through stack ivector/dvector
	// on the master and an extra round-trip.  Bit-identical because there
	// is no summation -- each slave owns a unique site slice and the data
	// is just moved into place.
	assert(myid == 0);
	MESSAGE signal = UPDATE_SPROFILE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	int recvcounts[nprocs], displs[nprocs];
	compute_siteprofile_layout(nprocs, GetNsite(), GetGlobalNstate(), recvcounts, displs);

	// Master contributes nothing; receives full gather into allocsiteprofile*.
	MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_INT,
	               allocsiteprofilesuffstatcount, recvcounts, displs, MPI_INT,
	               MPI_COMM_WORLD);
	MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE,
	               allocsiteprofilesuffstatbeta, recvcounts, displs, MPI_DOUBLE,
	               MPI_COMM_WORLD);
}

void ExpoConjugateGTRPhyloProcess::SlaveUpdateSiteProfileSuffStat()	{

	UpdateSiteProfileSuffStat();

	// siteprofilesuffstatcount[i] is a pointer-view into
	// allocsiteprofilesuffstatcount + i*GetDim() (see Create()), so the
	// slave's just-computed slice already lives at the correct offset
	// (displs[myid]) in the flat buffer.  Allgatherv with MPI_IN_PLACE
	// reads from that slot and gathers all slices into every rank.
	int recvcounts[nprocs], displs[nprocs];
	compute_siteprofile_layout(nprocs, GetNsite(), GetGlobalNstate(), recvcounts, displs);

	MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_INT,
	               allocsiteprofilesuffstatcount, recvcounts, displs, MPI_INT,
	               MPI_COMM_WORLD);
	MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE,
	               allocsiteprofilesuffstatbeta, recvcounts, displs, MPI_DOUBLE,
	               MPI_COMM_WORLD);
}

void ExpoConjugateGTRPhyloProcess::GlobalUpdateRRSuffStat()	{

	// MPI2
	// should send message to slaves for updating their rrsuffstats
	// by calling UpdateRRSuffStat();
	// then collect all suff stats
	//
	// suff stats are contained in 2 arrays
	// int* rrsuffstatcount
	// double* rrsuffstatbeta

	// should be summed over all slaves (reduced)
	assert(myid == 0);
	int i,j,workload = Nrr;
	MPI_Status stat;
	MESSAGE signal = UPDATE_RRATE;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for(i=0; i<workload; ++i) {
		rrsuffstatcount[i] = 0;
		rrsuffstatbeta[i] = 0.0;
	}

	int ivector[workload];
	double dvector[workload];
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(ivector,workload,MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
		// MPI_Recv(ivector,workload,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<workload; ++j) {
			rrsuffstatcount[j] += ivector[j];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(dvector,workload,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		// MPI_Recv(dvector,workload,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<workload; ++j) {
			rrsuffstatbeta[j] += dvector[j];
		}
	}

	MPI_Bcast(rrsuffstatcount,Nrr,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(rrsuffstatbeta,Nrr,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void ExpoConjugateGTRPhyloProcess::SlaveUpdateRRSuffStat()	{

	UpdateRRSuffStat();
	int workload = Nrr;

	MPI_Send(rrsuffstatcount,workload,MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(rrsuffstatbeta,workload,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	MPI_Bcast(rrsuffstatcount,Nrr,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(rrsuffstatbeta,Nrr,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

int ExpoConjugateGTRPhyloProcess::GlobalCountMapping()	{

	GlobalUpdateSiteProfileSuffStat();
	return PhyloProcess::GlobalCountMapping();
}

int ExpoConjugateGTRPhyloProcess::CountMapping()	{

	int total = 0;	
	for(int i = sitemin; i < sitemax; i++){
		total += CountMapping(i);
	}
	return total;
}

int ExpoConjugateGTRPhyloProcess::CountMapping(int i)	{

	const int* tmp = GetSiteProfileSuffStatCount(i);
	int total = 0;
	for (int k=0; k<GetNstate(i); k++)	{
		total += tmp[k];
	}
	total--;
	return total;
}

