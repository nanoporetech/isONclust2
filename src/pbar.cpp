#include <iostream>

// From: https://stackoverflow.com/a/14539953
void Pbar(double progress)
{
    int barWidth = 70;

    std::cerr << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
	if (i < pos)
	    std::cerr << "=";
	else if (i == pos)
	    std::cerr << ">";
	else
	    std::cerr << " ";
    }
    std::cerr << "] " << int(progress * 100.0) << " %\r";
    std::cerr.flush();
}
