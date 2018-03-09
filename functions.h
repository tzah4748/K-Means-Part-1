
typedef struct 
{
	double x = 0;
	double y = 0;
	double diameter = 0;
	double xSum = 0;
	double ySum = 0;
	int pointsInCluster = 0;
} Cluster;

typedef struct 
{
	double x = 0;
	double y = 0;
	double vx = 0;
	double vy = 0;
	Cluster* center = NULL;
	int centerChanged;
} Point;

void getInitialsFromFile(FILE* f, int* numOfPoints, int* numOfClusters, int* timeInterval, double* deltaTime, int* maxIterations, double* qualityMeasure);
void getPointsFromFile(FILE *f, Point* allPoints);
void initializeClusters(Cluster* clusterCenters, Point* allPoints, int numOfClusters);

double xPosition(Point p, double time);
double yPosition(Point p, double time);

double distClusterToPoint(Cluster c, Point p);
double distPointToPoint(Point p1, Point p2);

double calculateDiameter(Cluster* someCluster, Point* allPoints, int numOfPoints);
double calculateQM(Cluster* clusterCenters, int numOfClusters);
void writeResultsToFile(FILE* f, const char* outputPath, double time, int iter, double terminationQuality, Cluster* clusterCenters, int numOfClusters);

void printPoints(Point* p, int numOfPoints);
void createRandomPointsFile(FILE *f);