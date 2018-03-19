#include "functions.h"

int main(int argc, char *argv[])
{
	clock_t t1, t2;
	t1 = clock();
	const char* inputFile = "input.txt";
	const char* outputFile = "output.txt";
	/*Initials*/
	int numOfPoints, numOfClusters, maxIterations, iter, lastIndexGiven, terminationCondition = false;
	double timeInterval, deltaTime, qualityMeasure, time, terminationQuality, x, y, distFromCluster, oldDistFromCluster;
	Point *allPoints;
	Cluster *clusterCenters;
	
	//Generate new random points every excution of the code
	//createRandomPointsFile(inputFile);

	allPoints = getPointsFromFile(inputFile, &numOfPoints, &numOfClusters, &timeInterval, &deltaTime, &maxIterations, &qualityMeasure);
	clusterCenters = initializeClusters(allPoints, numOfClusters);
	for (time = 0; time < timeInterval; time += deltaTime)
	{
		for (iter = 0; iter < maxIterations; iter++)
		{
			for (int i = 0; i < numOfPoints; i++)
			{	/*Get point[i]'s Position*/
				x = xPosition(allPoints[i], time);
				y = yPosition(allPoints[i], time);
				/*At the beginning of the iteration(iter) the point's center is it's old center, meaning it hasn't changed yet*/
				allPoints[i].centerChanged = false;
				/*Define a center to a point if it doesn't exists*/
				if (allPoints[i].clusterIndex == -1)
				{
					allPoints[i].clusterIndex = 0;
					lastIndexGiven = -1;
				}
				else
					lastIndexGiven = allPoints[i].clusterIndex;
				/*Calculate the distance from point[i] to it's own cluster center*/
				oldDistFromCluster = sqrt(pow((clusterCenters[allPoints[i].clusterIndex].x - x), 2) + pow((clusterCenters[allPoints[i].clusterIndex].y - y), 2));
				/*For each point[i] calculate a new cluster center*/
				for (int j = 0; j < numOfClusters; j++)
				{
					/*The new distance from point[i] to cluster[j]*/
					distFromCluster = sqrt(pow((clusterCenters[j].x - x), 2) + pow((clusterCenters[j].y - y), 2));
					/*If true:
					- Define new distance as "old"
					- Define new center to the point[i]
					- Point[i]'s center has now changed.*/
					if (distFromCluster < oldDistFromCluster)
					{
						oldDistFromCluster = distFromCluster;
						lastIndexGiven = j;
					}
				}
				if (lastIndexGiven != allPoints[i].clusterIndex)
				{
					allPoints[i].centerChanged = true;
					allPoints[i].clusterIndex = lastIndexGiven;
				}
				/*After defining a point to its center, add the point to the center
				- This simple logic works due to the fact a cluster's new center is an average of all points in the cluster*/
				clusterCenters[allPoints[i].clusterIndex].xSum += x;
				clusterCenters[allPoints[i].clusterIndex].ySum += y;
				clusterCenters[allPoints[i].clusterIndex].pointsInCluster += 1;
				/*Termination condition/flag will stay true as long as at least one point[i] defined a new center*/
				terminationCondition |= allPoints[i].centerChanged;
			}
			/*False Positive - the flag is false if all points stayed in their clusters*/
			if (!terminationCondition) //False Positive
				break;
			terminationCondition = false;
			/*Calculate cluster[j]'s new center and reset cluster's info for the points he "own"*/
			for (int j = 0; j < numOfClusters; j++)
			{
				clusterCenters[j].x = clusterCenters[j].xSum / clusterCenters[j].pointsInCluster;
				clusterCenters[j].y = clusterCenters[j].ySum / clusterCenters[j].pointsInCluster;
				clusterCenters[j].xSum = 0;
				clusterCenters[j].ySum = 0;
				clusterCenters[j].pointsInCluster = 0;
			}
		}
		/*After the "iter" loop has stopped, calculate each cluster's diameter*/
		for (int j = 0; j < numOfClusters; j++)
			clusterCenters[j].diameter = calculateDiameter(allPoints, numOfPoints, j);
		/*Check if the quality of current clusters is enough to stop time iterations loop*/
		terminationQuality = calculateQM(clusterCenters, numOfClusters);
		if (terminationQuality <= qualityMeasure)
			break;
	}
	writeResultsToFile(outputFile, time, iter, terminationQuality, clusterCenters, numOfClusters);

	/*Printed results (for debugging)*/
	
	printf("First occurence at t = %f, iter = %d with q = %f\n\nCenter of the Clusters:\n\n", time, iter, terminationQuality);
	for (int i = 0; i < numOfClusters; i++)
	{
	printf("%f %f\n\n", clusterCenters[i].x, clusterCenters[i].y);
	}
	
	printf("Done!\nFile Location:\n%s\n", outputFile);
	free(allPoints);
	free(clusterCenters);
	t2 = clock();
	printf("K-Means Time: %f\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
	
	return 0;
}

void checkFile(FILE* f)
{
	if (!f)
	{
		printf("Failed opening the file. Exiting!\n");
		exit(1);
	}
}

Point* getPointsFromFile(const char* inputFile, int* numOfPoints, int* numOfClusters, double* timeInterval, double* deltaTime, int* maxIterations, double* qualityMeasure)
{
	FILE* f;
	fopen_s(&f, inputFile, "r");
	checkFile(f);
	fscanf_s(f, "%d %d %lf %lf %d %lf", numOfPoints, numOfClusters, timeInterval, deltaTime, maxIterations, qualityMeasure);
	Point* allPoints = (Point*)calloc(*numOfPoints, sizeof(Point));
	for (int i = 0; i < *numOfPoints; i++)
		fscanf_s(f, "%lf %lf %lf %lf", &allPoints[i].x, &allPoints[i].y, &allPoints[i].vx, &allPoints[i].vy);
	fclose(f);
	return allPoints;
}

Cluster* initializeClusters(Point* allPoints, int numOfClusters)
{
	Cluster* clusterCenters = (Cluster*)calloc(numOfClusters, sizeof(Cluster));
	for (int i = 0; i < numOfClusters; i++)
	{
		clusterCenters[i].x = allPoints[i].x;
		clusterCenters[i].y = allPoints[i].y;
	}
	return clusterCenters;
}

double xPosition(Point p, double time)
{
	return p.x + time*p.vx;
}

double yPosition(Point p, double time)
{
	return p.y + time*p.vy;
}

double calculateDiameter(Point* allPoints, int numOfPoints, int clusterIndex)
{
	double diameter = 0;
	double ret = 0;
	for (int i = 0; i < numOfPoints; i++)
	{
		Point a = allPoints[i];

		if (a.clusterIndex == clusterIndex)
		{
			for (int j = i + 1; j < numOfPoints; j++)
			{
				Point b = allPoints[j];
				if (b.clusterIndex == clusterIndex)
				{
					diameter = sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2));
					if (diameter > ret)
						ret = diameter;
				}
			}
		}
	}
	return ret;
}

double calculateQM(Cluster* clusterCenters, int numOfClusters)
{
	double q = 0;
	double count = 0;
	for (int i = 0; i < numOfClusters; i++)
	{
		for (int j = 0; j < numOfClusters; j++)
		{
			if (i != j)
			{
				q += clusterCenters[i].diameter / sqrt(pow((clusterCenters[i].x - clusterCenters[j].x), 2) + pow((clusterCenters[i].y - clusterCenters[j].y), 2));
				count += 1;
			}
		}
	}
	q = q / count;
	return q;
}

void writeResultsToFile(const char* outputFile, double time, int iter, double terminationQuality, Cluster* clusterCenters, int numOfClusters)
{
	FILE* f;
	fopen_s(&f, outputFile, "w");
	if (f == NULL)
	{
		printf("Failed opening the file. Exiting!\n");
		exit(1);
	}
	fprintf(f, "First occurence at t = %f, iter = %d with q = %f\n\nCenter of the Clusters:\n\n", time, iter, terminationQuality);
	for (int i = 0; i < numOfClusters; i++)
	{
		fprintf(f, "%f %f\n\n", clusterCenters[i].x, clusterCenters[i].y);
	}
	fclose(f);
}

void printPoints(Point* p, int numOfPoints)
{
	for (int i = 0; i < numOfPoints; i++)
		printf("%0.2f %0.2f\n", p[i].x, p[i].y);
}

void createRandomPointsFile(const char* inputFile)
{
	FILE* f;
	fopen_s(&f, inputFile, "w");
	checkFile(f);
	srand((unsigned int)time(NULL));
	fputs("10000 4 30 0.1 2000 7.3\n", f);
	for (int i = 0; i < 10000; i++)
	{
		fprintf(f, "%.2f %.2f %.2f %.2f\n", double((rand() % 20000) - 10000) / 100, double((rand() % 20000) - 10000) / 100, double((rand() % 40000) - 20000) / 100, double((rand() % 40000) - 20000) / 100);
	}
	fclose(f);
}
