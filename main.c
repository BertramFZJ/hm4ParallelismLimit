#define _SET_AFFINITY_ 1

#define _CRT_SECURE_NO_WARNINGS

#define _ACCURACY_ 1
#define _LEVELS_NUMBER_ 1
#define _GAUSS_GRADIENT_TRANSFER_ 0
#define _ITERATIONS_NUMBER_ 500

#if( (_ACCURACY_ == 2) && (_GAUSS_GRADIENT_TRANSFER_ == 1) )
#define _LEVELS_NUMBER_ 1
#endif

#if( (_ACCURACY_ == 2) && (_GAUSS_GRADIENT_TRANSFER_ == 0) && (_LEVELS_NUMBER_ < 2) )
#define _LEVELS_NUMBER_ 2
#endif

#define _MIXED_MESH_NAME_ "./INPUT_DATA/sphereAMR_1.hmsb"

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include "pm18MpiHeader.h"

#if(_SET_AFFINITY_ == 1)

#ifndef __USE_GNU
#define __USE_GNU
#include <sched.h>
#include <sys/types.h>
#undef __USE_GNU
#else
#include <sched.h>
#include <sys/types.h>
#endif

#endif

#include "ioHM3D.h"

#include "typesMeshHeader.h"
#include "cfdHeader.h"

#include "hm4MeshPreprocessorCfdMPI.h"
#include "cfdGeometryInitialization.h"
#include "cfdDataInitialization.h"
#include "mpiTransferHost.h"
#include "mpiCfdCoreMain.h"

int main(int argc, char *argv[])
{
	pm18MpiTopologyType mpiTopology;
    gtypeHm4MeshAttributes meshAttributes;    
    gtypeHm4MeshTopology globalMesh, hostMesh;	
	gtypeHm4CfdHostCore cfdCoreHost;

	double taskTime = - 100.0;

	// »Õ»÷»¿À»«¿÷»ﬂ “ŒœŒÀŒ√»» MPI
	{
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &mpiTopology.sizeWORLD);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpiTopology.rankWORLD);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(mpiTopology.rankWORLD == 0)
		{ printf("MPI_COMM_WORLD: CHECK TOPOLOGY >>> %d MPI_PROCESSES\n", mpiTopology.sizeWORLD); fflush(stdout); }
		MPI_Barrier(MPI_COMM_WORLD);		
	}
	// »Õ»÷»¿À»«¿÷»ﬂ “ŒœŒÀŒ√»» MPI

	// œ–»¬ﬂ« ¿   ﬂƒ–¿Ã
#if(_SET_AFFINITY_ == 1)
	{
		int corePerNode = 128; // ◊»—ÀŒ CPU-ﬂƒ≈– Õ¿ ŒƒÕŒÃ ”«À≈
		int mpiProcessPerNode = 128; // ◊»—ÀŒ MPI-œ–Œ÷≈——Œ¬ Õ¿ ŒƒÕŒÃ ”«À≈

		int coresPerMpiProcess = 1;
		int myNodeId = -1;
		int myCoreId = -1;

		cpu_set_t mask; int get, idProcess;

		if(mpiTopology.rankWORLD == 0) { printf("\n"); printf("SET AFFINITY FOR MPI PROCESSES ..... \n"); fflush(stdout); }
		MPI_Barrier(MPI_COMM_WORLD);

		coresPerMpiProcess = corePerNode / mpiProcessPerNode;
		myNodeId = mpiTopology.rankWORLD / mpiProcessPerNode;
		myCoreId = mpiTopology.rankWORLD % mpiProcessPerNode;
		myCoreId = myCoreId * coresPerMpiProcess;

		// SET AFFINITY
		get = sched_getaffinity(0, sizeof(cpu_set_t), &mask); if(get != 0) exit(0);
        CPU_ZERO(&mask);
        CPU_SET(myCoreId, &mask);
		get = sched_setaffinity(0, sizeof(cpu_set_t), &mask); if(get != 0) exit(0);
		// SET AFFINITY

		// CHECK AFFINITY
		if(mpiTopology.rankWORLD != 0) MPI_Recv((void *)&get, 1, MPI_INT, mpiTopology.rankWORLD - 1, 187, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		{ int i;
			for(i=0; i<CPU_SETSIZE; i++) if(CPU_ISSET(i, &mask) != 0) { printf("[mpi %3d] [core %3d]\n", mpiTopology.rankWORLD, i); fflush(stdout); } }
		if(mpiTopology.rankWORLD != mpiTopology.sizeWORLD - 1) MPI_Send(&get, 1, MPI_INT, mpiTopology.rankWORLD + 1, 187, MPI_COMM_WORLD);
		// CHECK AFFINITY

		MPI_Barrier(MPI_COMM_WORLD);
		if(mpiTopology.rankWORLD == 0) { printf("\n"); fflush(stdout); }
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif
	// œ–»¬ﬂ« ¿   ﬂƒ–¿Ã

	// »Õ»÷»¿À»«¿÷»ﬂ –¿¡Œ◊≈… √–”œœ€
	{
		mpiTopology.rankTASK = mpiTopology.rankWORLD;
		mpiTopology.sizeTASK = mpiTopology.sizeWORLD;
		mpiTopology.commTASK = MPI_COMM_WORLD;
	}
	// »Õ»÷»¿À»«¿÷»ﬂ –¿¡Œ◊≈… √–”œœ€

    hm4SerialMeshPreprocessorCfdMPI(_MIXED_MESH_NAME_, _LEVELS_NUMBER_, &meshAttributes, &globalMesh, &hostMesh, mpiTopology.commTASK);

    hm4RotateHostMeshCfdMPI(&hostMesh, mpiTopology.commTASK);

#if 0
    {
        char fname[300];
        int neProc = 0;
        int i;

        for(i=0; i<rotateMesh.meshElementsNumber; i++)
            if(rotateMesh.partitionUP[i] != mpiTopology.rankTASK)
            { neProc = i; break; }

        sprintf(fname, "./OUTPUT_DATA/decoCenter_%02d.plt", mpiTopology.rankTASK);
        IOHM3D_PlotHm4Mesh(rotateMesh.meshNodesNumber, rotateMesh.meshNodesCoordinates, neProc, rotateMesh.meshEX, rotateMesh.meshEA, fname);

        sprintf(fname, "./OUTPUT_DATA/decoAll_%02d.plt", mpiTopology.rankTASK);
        IOHM3D_PlotHm4Mesh(rotateMesh.meshNodesNumber, rotateMesh.meshNodesCoordinates, rotateMesh.meshElementsNumber, rotateMesh.meshEX, rotateMesh.meshEA, fname);

    }
#endif

	// »Õ»÷»¿À»«¿÷»ﬂ —“–” “”– œ–»≈Ã¿-œ≈–≈ƒ¿◊» ƒ¿ÕÕ€’ —Œ —œ≈÷»¿À‹Õ€Ã» ¬≈Ÿ≈—“¬≈ÕÕ€Ã» ¡”‘≈–¿Ã» ƒÀﬂ «Õ¿◊≈Õ»… ‘”Õ ÷»…
	InterconnectInitializationHostMPI(&hostMesh, &mpiTopology.hostInterconnect, mpiTopology.commTASK);

	cfdGeometryInitializationHost(&hostMesh, &cfdCoreHost, mpiTopology.commTASK);

	cfdDataInitializationHost(meshAttributes, &cfdCoreHost, mpiTopology);

	// «¿Ã≈Õ¿ —œ≈÷»¿À‹ÕŒ√Œ ¡”‘≈–¿ œ–»≈Ã¿ ƒ¿ÕÕ€’ Õ¿ —–≈« Ã¿——»¬¿ –¿—◊≈“Õ€’ œ¿–¿Ã≈“–Œ¬, ◊“Œ »— Àﬁ◊¿≈“ œ–ŒÃ≈∆”“Œ◊Õ”ﬁ ¡”‘≈–»«¿÷»ﬁ œ–»Õ»Ã¿≈Ã€’ ƒ¿ÕÕ€’
	// !!!!! Õ¿ƒŒ —ƒ≈À¿“‹ ›“Œ Œœ÷»≈… Ë œ–Œ¬≈—“» “≈—“»–Œ¬¿Õ»≈ !!!!!
	RebootHostsInterconnectRecvRequestsMPI(&mpiTopology.hostInterconnect, cfdCoreHost.cellQph + cfdCoreHost.nEinterior*5, mpiTopology.commTASK);	

	if(_ACCURACY_ == 1)
	{
		if(_LEVELS_NUMBER_ >  1)
			taskTime = cfdCoreModelMainLevelAc1(_ITERATIONS_NUMBER_, &cfdCoreHost, mpiTopology);

		if(_LEVELS_NUMBER_ == 1)
			taskTime = cfdCoreModelMainAc1     (_ITERATIONS_NUMBER_, &cfdCoreHost, mpiTopology);
	}


	if(_ACCURACY_ == 2)
	{
		if(_LEVELS_NUMBER_ >  2)
		{
			if(_LEVELS_NUMBER_ % 2 != 0)
			{
				MPI_Barrier(mpiTopology.commTASK);
				if(mpiTopology.rankTASK == 0)
				{ fprintf(stderr, "_LEVELS_NUMBER_ %% 2 != 0\n"); fflush(stderr); exit(0); }
				MPI_Barrier(mpiTopology.commTASK);				
			}

			taskTime = cfdCoreModelMainLevelAc2(_ITERATIONS_NUMBER_, &cfdCoreHost, mpiTopology);
		}
		
		if(_LEVELS_NUMBER_ == 2)
			taskTime = cfdCoreModelMainAc2     (_ITERATIONS_NUMBER_, &cfdCoreHost, mpiTopology);

		if(_LEVELS_NUMBER_ == 1)
		{
			InitHostsInterconnectGradientRequestsMPI(&mpiTopology.hostInterconnect, cfdCoreHost.celldQph + cfdCoreHost.nEinterior * 15, mpiTopology.commTASK);
			taskTime = cfdCoreModelMainTransfer2Ac2(_ITERATIONS_NUMBER_, &cfdCoreHost, mpiTopology);
		}
	}

#if 1
	fflush(stdout); fflush(stderr);
	MPI_Barrier(mpiTopology.commTASK);
	if(mpiTopology.rankTASK == 0)
	{
		FILE *file = fopen("./OUTPUT_DATA/report.txt", "a");

		fprintf(file, "MESH %d NODES %d CELLS\n", meshAttributes.meshNodesNumber, meshAttributes.meshElementsNumber);
		fprintf(file, "MPI COMM SIZE %d\n", mpiTopology.sizeTASK);
		fprintf(file, "ACCURACY %d\n", _ACCURACY_);
		fprintf(file, "LEVELS NUMBER %d\n", _LEVELS_NUMBER_);		
		if(_ACCURACY_ == 2) fprintf(file, "GAUSS GRADIENT TRANSFER %d\n", _GAUSS_GRADIENT_TRANSFER_);
		fprintf(file, "ITERATIONS NUMBER %d\n", _ITERATIONS_NUMBER_);
		fprintf(file, "TIME %.4lf SPEED %.4lf iter per sec\n", taskTime, (double)_ITERATIONS_NUMBER_ / taskTime);
		fprintf(file, "\n");
		
		fflush(file); fclose(file);
	}
	MPI_Barrier(mpiTopology.commTASK);
#endif

	MPI_Barrier(MPI_COMM_WORLD);
	if(mpiTopology.rankWORLD == 0) { fprintf(stderr, "\nBY\n"); fflush(stderr); }
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();	

	return 0;
}
