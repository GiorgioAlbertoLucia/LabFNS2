#include <numeric>
#include <algorithm>

#include <TGraph.h>
#include <TString.h>

#include "graphUtilities.hh"

/**
 * @brief Function to find the edge of the signal.
 * 
 * @param gr 
 * @param nIgnorePoints number of points to ignore at the beginning and at the end of the TGraph (due to the smoothing and/or 
 * peculiarity of the waveform)
 * @return std::tuple<double, double, double> edgeLeft, edgeRight
 */
std::pair<double, double> FindEdge(TGraph & gr, int nIgnorePoints, int nChecks, int nSample, const int nDerivativePoints, const int nSmoothingPoints)
{
    auto edgesIdx = FindEdgeIndex(gr, nIgnorePoints, nChecks, nSample, nDerivativePoints, nSmoothingPoints);
    return std::make_pair(gr.GetPointX(edgesIdx.first), gr.GetPointX(edgesIdx.second));
}

/**
 * @brief Function to find the edge of the signal.
 * 
 * @param gr 
 * @return std::tuple<double, double, double> edgeLeft, edgeRight
 */
std::pair<int, int> FindEdgeIndex(TGraph & gr, int nIgnorePoints, int nChecks, int nSample, const int nDerivativePoints, const int nSmoothingPoints)
{
    auto grSmooth = SmoothGraph(gr, nSmoothingPoints);                  // smooth the graph to find the edge more easily
    auto grDerivative = DerivativeGraph(*grSmooth, nDerivativePoints);  // take the derivative of the smoothed graph

    double meanDer{0.}, RMSDer{0.};
    if (nSample > gr.GetN())  nSample = int(gr.GetN()/2);
    if (nIgnorePoints == 0)  nIgnorePoints = nSmoothingPoints;
    auto resultDer = GetMeanAndRMS(*grDerivative, nIgnorePoints, nIgnorePoints+nSample);
    meanDer = resultDer.first;
    RMSDer = resultDer.second;

    double edgeLeft{0.}, edgeRight{0.};    
    int edgeLeftIndex = FindChangingDerivative(*grDerivative, meanDer, RMSDer, nIgnorePoints, nChecks);          // correct for smoothing
    int edgeRightIndex = FindChangingDerivative(*grDerivative, meanDer, RMSDer, nIgnorePoints, nChecks, true);   // correct for smoothing

    delete grSmooth;
    delete grDerivative;
    
    return std::make_pair(edgeLeftIndex, edgeRightIndex);
}

/**
 * @brief Function to find the minimum value in a graph.
 * 
 * @param gr 
 * @return double 
 */
int FindMinimumIndex(TGraph & gr)
{
    const int nSamples = gr.GetN();
    double * x = gr.GetX();
    double * y = gr.GetY();

    double minimum{999999.};
    int minimumIndex{0};
    for (int i = 0; i < nSamples; ++i)
    {
        if (y[i] < minimum)     
        {
            minimum = y[i];
            minimumIndex = i;
        }   
    }

    return minimumIndex;
}

/**
 * @brief Function to find the maximum value in a graph.
 * 
 * @param gr 
 * @return double 
 */
int FindMaximumIndex(TGraph & gr)
{
    const int nSamples = gr.GetN();
    double * x = gr.GetX();
    double * y = gr.GetY();

    double maximum{-9999999.};
    int maximumIndex{0};
    for (int i = 0; i < nSamples; ++i)
    {
        if (y[i] > maximum)     
        {
            maximum = y[i];
            maximumIndex = i;
        }   
    }

    return maximumIndex;
}

/**
 * @brief Find the point at which the derivative changes for n consecutive points with a value higher or lower than meanDer ± 3 * RMSDer.
 * 
 * @param grDer 
 * @param meanDer 
 * @param RMSDer 
 * @param ignorePoints Number of points to ignore at the beginning and at the end of the TGraph (due to smoothing).
 * @param nChecks Number of consecutive points with derivative beyond meanDer ± 3 * RMSDer to consider a change.
 * @param backwards If true, the function loops through the TGraph from the last point to the first one.
 * @return const int 
 */
int FindChangingDerivative(TGraph & grDer, const double & meanDer, const double & RMSDer, const int ignorePoints, const int nChecks, const bool backwards)
{
    const int nSamples = grDer.GetN();
    double * x = grDer.GetX();
    double * y = grDer.GetY();

    int changingPoint{0};
    int consecutiveCount{0};

    if (!backwards)
    {
        for (int i = ignorePoints; i < nSamples - nChecks + 1; ++i)
        {
            // Check if the derivative of n consecutive points is beyond meanDer ± 3 * RMSDer
            bool isChange = true;
            for (int j = 0; j < nChecks; ++j)
            {
                if (y[i + j] <= meanDer + 5 * RMSDer && y[i + j] >= meanDer - 5 * RMSDer)
                {
                    isChange = false;
                    break;
                }
            }

            if (isChange)
            {
                changingPoint = i;
                consecutiveCount = nChecks;
                break;
            }
        }
    }
    else
    {
        for (int i = nSamples - ignorePoints - 1; i > nChecks - 2; --i)
        {
            // Check if the derivative of n consecutive points is beyond meanDer ± 3 * RMSDer
            bool isChange = true;
            for (int j = 0; j < nChecks; ++j)
            {
                if (y[i - j] <= meanDer + 5 * RMSDer && y[i - j] >= meanDer - 5 * RMSDer)
                {
                    isChange = false;
                    break;
                }
            }

            if (isChange)
            {
                changingPoint = i;
                consecutiveCount = nChecks;
                break;
            }
        }
    }

    if (consecutiveCount < nChecks)
    {
        // If no change is found, return 0
        changingPoint = 0;
    }

    return changingPoint;
}

/**
 * @brief Get the Mean And RMS of a graph in a given range (points).
 * 
 * @param gr 
 * @param begin 
 * @param end 
 * @return std::pair<double, double> 
 */
std::pair<double, double> GetMeanAndRMS(TGraph & gr, const int & begin, const int & end)
{
    const int nSamples = gr.GetN();
    double * x = gr.GetX();
    double * y = gr.GetY();

    double mean{0.};
    for (int i = begin; i < end; ++i)
    {
        mean += y[i];
    }
    mean /= (end - begin);

    double RMS{0.};
    for (int i = begin; i < end; ++i)
    {
        RMS += (y[i] - mean) * (y[i] - mean);
    }
    RMS /= (end - begin);
    RMS = std::sqrt(RMS);

    return std::make_pair(mean, RMS);
}

/**
 * @brief Function to obtain a smoothed version of a TGraph, by doing a running average over nSmoothingPoints points.
 * 
 * @param gr 
 * @param nSmoothingPoints 
 * @return TGraph* 
 */
TGraph* SmoothGraph(TGraph& gr, const int nSmoothingPoints)
{
    const int nSamples = gr.GetN();
    double * x = gr.GetX();
    double * y = gr.GetY();

    auto grSmooth = new TGraph(0);
    grSmooth->SetName(Form("%s_smoothed", gr.GetName()));
    for (int i = 0; i < nSamples; ++i)
    {
        double sum{0.};
        for (int j = i - nSmoothingPoints; j < i + nSmoothingPoints; ++j)
        {
            if (j < 0 || j > nSamples) continue;
            sum += y[j];
        }
        grSmooth->SetPoint(i, x[i], sum / (2 * nSmoothingPoints + 1));
    }

    return grSmooth;
}

/**
 * @brief Function to obtain the derivative of a TGraph point by point. The derivative is computed as the average of the derivative done on 
 * nDerivativePoints points before and after the point in which the derivative is evaluated.
 * 
 * @param gr 
 * @param nDerivativePoints
 * @return TGraph* 
 */
TGraph* DerivativeGraph(TGraph& gr, const int nDerivativePoints)
{
    // compute the derivative of a TGraph point by point considering nDerivativePoints points before and after the point in which the derivative is evaluated
    const int nSamples = gr.GetN();
    double * x = gr.GetX();
    double * y = gr.GetY();

    auto grDerivative = new TGraph(0);
    grDerivative->SetName(Form("%s_derivative", gr.GetName()));
    
    for (int i = 1; i < nSamples - nDerivativePoints - 1; ++i)
    {
        double der{0.};
        for (int j = i; j < i + nDerivativePoints; j++)   der += (y[j + 1] - y[j - 1]) / (x[j + 1] - x[j - 1]);
        der /= nDerivativePoints;
        grDerivative->SetPoint(i, x[i], der);
    }
    grDerivative->SetPoint(0, x[0], grDerivative->GetPointY(1));
    grDerivative->SetPoint(nSamples - 1, x[nSamples - 1], grDerivative->GetPointY(nSamples - 2));

    // set the first and last point of the derivative graph to the second and second to last point of the original graph
    grDerivative->SetPoint(0, x[0], grDerivative->GetPointY(1));
    grDerivative->SetPoint(nSamples - 1, x[nSamples - 1], grDerivative->GetPointY(nSamples - 2));

    return grDerivative;
}

/**
 * @brief Find the first point in the TGraph with given value. The first skipPoints points are ignored.
 * 
 * @param gr 
 * @param x 
 * @param skipPoints 
 * @return int 
 */
int FindPointIndex(TGraph & gr, const double & x, const int nIgnorePoints)
{
    const int nSamples = gr.GetN();
    double * xGraph = gr.GetX();
    double * yGraph = gr.GetY();

    for (int i = nIgnorePoints; i < nSamples; ++i)  if (xGraph[i] == x) return i;
    return 0;
}

/**
 * @brief Find the point in TGraph with value closest to the target value.
 * 
 * @param gr
 * @param target
 * @return int
*/
int FindClosestPointIndex(TGraph & gr, const double target, const int nIgnorePoints)
{
    const int size = gr.GetN();
    double * array = gr.GetY();

    size_t effectiveSize = (size > nIgnorePoints) ? size - nIgnorePoints : 0;
    array += nIgnorePoints;

    // Create a vector of indices and sort it based on the values in the array
    size_t * indices = new size_t[effectiveSize];
    std::iota(indices, indices + effectiveSize, 0);
    std::sort(indices, indices + effectiveSize,
              [&array](size_t i1, size_t i2) { return array[i1] < array[i2]; });

    // Use binary search on the sorted indices
    auto it = std::lower_bound(indices, indices + effectiveSize, target,
                               [&array](size_t i, double value) { return array[i] < value; });

    if (it == indices) {
        // Target is less than or equal to the first element
        size_t result = *it;
        delete[] indices;
        return result;
    } else if (it == indices + effectiveSize) {
        // Target is greater than or equal to the last element
        size_t result = *(--it);
        delete[] indices;
        return result;
    }

    // Check the closest value between the current iterator and the one before
    size_t index1 = *it;
    size_t index2 = *(--it);

    size_t result = std::abs(array[index1] - target) < std::abs(array[index2] - target) ? index1 : index2;
    result += nIgnorePoints;

    delete[] indices;
    return result;
}
