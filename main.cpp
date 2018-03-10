#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include "functions.h"

int main(int argc, char *argv[])
{	
	FILE *f;
	const char* inputPath = "input.txt";
	const char* outputPath = "output.txt";
	/*Initials*/
	int numOfPoints, numOfClusters, timeInterval;
	double deltaTime;
	int maxIterations;
	double qualityMeasure;

	Point* allPoints;
	Cluster* clusterCenters;
	f = fopen(inputPath, "r");
	if (f == NULL)
	{
		printf("Failed opening the file. Exiting!\n");
		return 0;
	}
	getInitialsFromFile(f, &numOfPoints, &numOfClusters, &timeInterval, &deltaTime, &maxIterations, &qualityMeasure);
	allPoints = (Point*)calloc(numOfPoints, sizeof(Point));
	getPointsFromFile(f, allPoints);
	clusterCenters = (Cluster*)calloc(numOfClusters, sizeof(Cluster));
	initializeClusters(clusterCenters, allPoints, numOfClusters);

	double x, y, distFromCluster, oldDistFromCluster;
	int terminationCondition = false;
	double time;
	int iter;
	double terminationQuality;
	
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
				if(allPoints[i].center == NULL)
					allPoints[i].center = &clusterCenters[0];
				/*Calculate the distance from point[i] to it's own cluster center*/
				oldDistFromCluster = sqrt(pow((allPoints[i].center->x - x), 2) + pow((allPoints[i].center->y - y), 2));
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
						allPoints[i].center = &clusterCenters[j];
						allPoints[i].centerChanged = true;
					}
				}
				/*After defining a point to its center, add the point to the center
				- This simple logic works due to the fact a cluster's new center is an average of all points in the cluster*/
				allPoints[i].center->pointsInCluster += 1;
				allPoints[i].center->xSum += x;
				allPoints[i].center->ySum += y;
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
		{
			clusterCenters[j].diameter = calculateDiameter(&clusterCenters[j], allPoints, numOfPoints);
		}
		/*Check if the quality of current clusters is enough to stop time iterations loop*/
		terminationQuality = calculateQM(clusterCenters, numOfClusters);
		if (terminationQuality <= qualityMeasure)
			break;
	}

	writeResultsToFile(f, outputPath, time, iter, terminationQuality, clusterCenters, numOfClusters);
	
	/*Printed results (for debugging)*/
	/*
	printf("First occurence at t = %f, iter = %d with q = %f\n\nCenter of the Clusters:\n\n", time, iter, terminationQuality);
	for (int i = 0; i < numOfClusters; i++)
	{
		printf("%f %f\n\n", clusterCenters[i].x, clusterCenters[i].y);
	}
	*/
	free(allPoints);
	free(clusterCenters);
	return 0;
}

void getInitialsFromFile(FILE* f, int* numOfPoints, int* numOfClusters, int* timeInterval, double* deltaTime, int* maxIterations, double* qualityMeasure)
{	
	fscanf(f, "%d %d %d %lf %d %lf", numOfPoints, numOfClusters, timeInterval, deltaTime, maxIterations, qualityMeasure);
}
void getPointsFromFile(FILE *f, Point* allPoints)
{
	int i = 0;
	while (!feof(f))
	{
		fscanf(f, "%lf %lf %lf %lf", &allPoints[i].x, &allPoints[i].y, &allPoints[i].vx, &allPoints[i].vy);
		i += 1;
	}
	fclose(f);
}
void initializeClusters(Cluster* clusterCenters, Point* allPoints, int numOfClusters)
{
	for (int i = 0; i < numOfClusters; i++)
	{
		clusterCenters[i].x = allPoints[i].x;
		clusterCenters[i].y = allPoints[i].y;
	}
}

double xPosition(Point p, double time)
{
	return p.x + time*p.vx;
}
double yPosition(Point p, double time)
{
	return p.y + time*p.vy;
}

double distClusterToPoint(Cluster c, Point p)
{
	double ret = sqrt((c.x - p.x)*(c.x - p.x) + (c.y - p.y)*(c.y - p.y));
	return ret;
}
double distPointToPoint(Point p1, Point p2)
{
	double ret = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
	return ret;
}

double calculateDiameter(Cluster* someCluster, Point* allPoints, int numOfPoints)
{
	double diameter = 0;
	double ret = 0;
	for (int i = 0; i < numOfPoints; i++)
	{
		Point a = allPoints[i];
		
		if (a.center == someCluster)
		{
			for (int j = 0; j < numOfPoints; j++)
			{
				Point b = allPoints[j];
				if (i != j & b.center == someCluster)
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
void writeResultsToFile(FILE* f, const char* outputPath, double time, int iter, double terminationQuality, Cluster* clusterCenters, int numOfClusters)
{
	f = fopen(outputPath, "w");
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
void createRandomPointsFile(FILE *f)
{
	srand(time(NULL));
	fputs("5000 4 30 0.1 2000 7.3\n", f);
	for (int i = 0; i < 5000; i++)
	{
		fprintf(f, "%.2f %.2f %.2f %.2f\n", double((rand() % 20000) - 10000)/100, double((rand() % 20000) - 10000)/100, double((rand() % 40000) - 20000)/100, double((rand() % 40000) - 20000)/100);
	}
	fclose(f);
}
