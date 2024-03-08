#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "pm18MpiHeader.h"

#include "typesMeshHeader.h"
#include "cfdHeader.h"
#include "constantsCFD.h"

#include "cfdFieldValuesInitialization.h"
#include "mpiCfdCore.h"

#include "cfdDataInitialization.h"

int cfdDataInitializationHost(gtypeHm4MeshAttributes meshAttributes, gtypeHm4CfdHostCore *cfdCore, pm18MpiTopologyType mpiTopology)
{
	// опнбепйю йнппейрмнярх йнмярюмр
	if(mpiTopology.rankTASK == 0)
	{
		double diff = hm4CheckCFDDefinedConstants();
		printf("hm4CheckCFDDefinedConstants => %E\n", diff); printf("\n"); fflush(stdout);
		if(diff > 1.0e-8) { fprintf(stderr, "[%d]: ERROR INIT CFD <000>\n", mpiTopology.rankTASK); exit(0); }
	}
	MPI_Barrier(mpiTopology.commTASK);
	// опнбепйю йнппейрмнярх йнмярюмр

	// хмхжхюкхгюжхъ гмювемхи мебнглсыеммнцн онрнйю
	{
		double estimatedDT;

		cfdCore->cellQ   = (double *)malloc(cfdCore->nE * 5 * sizeof(double)); if(cfdCore->cellQ   == NULL) exit(0);
		cfdCore->cellQph = (double *)malloc(cfdCore->nE * 5 * sizeof(double)); if(cfdCore->cellQph == NULL) exit(0);
		hm4CfdInitCellsStartValuesForUndisturbedFlow(cfdCore->nE, cfdCore->cellQ, cfdCore->cellQph);

		mpiCoreCalculateTimeStep(0, cfdCore->nEinterior,
			                     cfdCore->cellQph, cfdCore->cellH, &estimatedDT,
								 0, cfdCore->cellLevels,
								 mpiTopology.commTASK);
		if(mpiTopology.rankTASK == 0) { printf(">>>>>>>>>>>>>>>>> ESTIMATED DT: %g\n", estimatedDT); fflush(stdout); }
		MPI_Barrier(mpiTopology.commTASK);

		mpiCoreCheckMinMaxPhysAndPrimVariables(cfdCore->nE, cfdCore->cellQ, cfdCore->cellQph, stdout, mpiTopology);		
	}
	// хмхжхюкхгюжхъ гмювемхи мебнглсыеммнцн онрнйю

	// онярюмнбйю цпюмхвмшу сякнбхи
	if(mpiTopology.rankTASK == 0)
	hm4SetBoundaryConditionsForBoundaryFaces(meshAttributes.boundarySurfaceNames, 81, stdout,
			                                 cfdCore->nBCF, cfdCore->BCF,
											 nBSF, (int *)BoundaryConditions);
	else
	hm4SetBoundaryConditionsForBoundaryFaces(meshAttributes.boundarySurfaceNames, 81, NULL,
			                                 cfdCore->nBCF, cfdCore->BCF,
											 nBSF, (int *)BoundaryConditions);
	MPI_Barrier(mpiTopology.commTASK);
	// онярюмнбйю цпюмхвмшу сякнбхи

	// онрнйх он цпюмъл
	cfdCore->faceCE = (double *)malloc((cfdCore->nICF + cfdCore->nBCF) * 5 * sizeof(double));
	if(cfdCore->faceCE == NULL) exit(0);
	memset((void *)cfdCore->faceCE, 0, (cfdCore->nICF + cfdCore->nBCF) * 5 * sizeof(double));
	MPI_Barrier(mpiTopology.commTASK);
	// онрнйх он цпюмъл

	// цпюдхемрш тхгхвеяйху оепелеммшу
	cfdCore->celldQph = (double *)malloc(cfdCore->nE * 15 * sizeof(double));
	if(cfdCore->celldQph == NULL) exit(0);
	memset((void *)cfdCore->celldQph, 0, cfdCore->nE * 15 * sizeof(double));
	MPI_Barrier(mpiTopology.commTASK);
	// цпюдхемрш тхгхвеяйху оепелеммшу

	return 0;
}
