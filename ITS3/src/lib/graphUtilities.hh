#ifndef GRAPH_UTILITIES_HH
#define GRAPH_UTILITIES_HH

#include <TGraph.h>

std::pair<double, double> FindEdge(TGraph & gr, int nIgnorePoints = 0, int nChecks = 10, int nSample = 200, const int nDerivativePoints = 40, const int nSmoothingPoints = 10);
std::pair<int, int> FindEdgeIndex(TGraph & gr, int nIgnorePoints = 0, int nChecks = 10, int nSample = 200, const int nDerivativePoints = 40, const int nSmoothingPoints = 10);
int FindMinimumIndex(TGraph & gr);
int FindMaximumIndex(TGraph & gr);  
int FindChangingDerivative(TGraph & grDer, const double& meanDer, const double& RMSDer, const int ignorePoints, const int nChecks, const bool backwards = false);
std::pair<double, double> GetMeanAndRMS(TGraph & gr, const int& begin, const int& end);
TGraph* SmoothGraph(TGraph& gr, const int nSmoothingPoints = 2);
TGraph* DerivativeGraph(TGraph& gr, const int nDerivativePoints = 40);
int FindPointIndex(TGraph& gr, const double& x, const int nIgnorePoints = 0);
int FindClosestPointIndex(TGraph& gr, const double y, const int nIgnorePoints = 0);

#endif