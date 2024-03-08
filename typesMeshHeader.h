#ifndef _gtypeHm4MeshAttributes_
typedef struct
{
	char *meshName;

    int meshNodesNumber;
    int meshElementsNumber;

    int meshBoundaryFacesNumber;
	int meshBoundarySurfacesNumber;
	char *boundarySurfaceNames;
} gtypeHm4MeshAttributes;
#define _gtypeHm4MeshAttributes_
#endif

#ifndef _gtypeHm4MeshTopology_
typedef struct
{
	int meshElementsNumber;
	int *meshEX;
	int *meshEA;

	int meshNodesNumber;
	double *meshNodesCoordinates;
	
	int meshBoundaryFacesNumber;
	int *meshBoundaryFacesTopology;

	int *partitionUP;
	int *indexesUP;

	int *partitionDOWN;
	int *indexesDOWN;
} gtypeHm4MeshTopology;
#define _gtypeHm4MeshTopology_
#endif
