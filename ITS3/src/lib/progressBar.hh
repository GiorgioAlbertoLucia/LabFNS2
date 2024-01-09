#ifndef PROGRESSBAR_HH
#define PROGRESSBAR_HH

#include <chrono>
#include <thread>

void updateProgressBar(int progress, int total, const std::chrono::steady_clock::time_point& startTime);

#endif // PROGRESSBAR_HH