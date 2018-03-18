#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

typedef struct cluster_t
{
	double x = 0;
	double y = 0;
	double xSum = 0;
	double ySum = 0;
	double diameter = 0;
	int pointsInCluster = 0;
} Cluster;

typedef struct point_t
{
	double x = 0;
	double y = 0;
	double vx = 0;
	double vy = 0;
	int centerChanged;
	int clusterIndex = -1;
} Point;

void checkFile(FILE* f);
Point* getPointsFromFile(const char* inputFile, int* numOfPoints, int* numOfClusters, double* timeInterval, double* deltaTime, int* maxIterations, double* qualityMeasure);
Cluster* initializeClusters(Point* allPoints, int numOfClusters);

double xPosition(Point p, double time);
double yPosition(Point p, double time);

double calculateDiameter(Point* allPoints, int numOfPoints, int clusterIndex);
double calculateQM(Cluster* clusterCenters, int numOfClusters);
void writeResultsToFile(const char* outputFile, double time, int iter, double terminationQuality, Cluster* clusterCenters, int numOfClusters);

void printPoints(Point* p, int numOfPoints);
void createRandomPointsFile(const char* inputFile);