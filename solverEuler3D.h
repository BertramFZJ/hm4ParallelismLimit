int solverRusanov3D(double *QphR, double *QphL,
					double *flux,
                    double    nX, double    nY, double    nZ,
                    double    sX, double    sY, double    sZ,
                    double    tX, double    tY, double    tZ);

int solverCUSP3D(double *QphR, double *QphL,
				 double *flux,
				 double    nX, double    nY, double    nZ,
				 double    sX, double    sY, double    sZ,
				 double    tX, double    tY, double    tZ);


#if 0
int solverRoePike3D(double *QphR, int shockFixR,
					double *QphL, int shockFixL,
					double *flux,
                    double    nX, double    nY, double    nZ,
                    double    sX, double    sY, double    sZ,
                    double    tX, double    tY, double    tZ);

int solverHLL3D(double *QphR, int shockFixR,
				double *QphL, int shockFixL,
				double *flux,
                double    nX, double    nY, double    nZ,
                double    sX, double    sY, double    sZ,
                double    tX, double    tY, double    tZ);

int solverAUSM3D(double *QphR, int shockFixR,
				 double *QphL, int shockFixL,
				 double *flux,
                 double    nX, double    nY, double    nZ,
                 double    sX, double    sY, double    sZ,
                 double    tX, double    tY, double    tZ);

int solverTurb3D(double *QphR, int shockFixR,
				 double *QphL, int shockFixL,
				 double *flux,
				 double    nX, double    nY, double    nZ,
				 double    sX, double    sY, double    sZ,
				 double    tX, double    tY, double    tZ);

#endif
