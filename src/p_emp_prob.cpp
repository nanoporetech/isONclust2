#include "p_emp_prob.h"
#include <cstdlib>
#include <exception>
#include <iostream>
#include <sstream>
#include "p_emp_prob_data.h"
#include "util.h"

bool operator<(const ErrTarget& a, const ErrTarget& b)
{
    if (a.Error1 < b.Error1) {
	return true;
    }
    else if (a.Error1 == b.Error1) {
	if (a.Error2 < b.Error2) {
	    return true;
	}
    }
    return false;
}

MinSharedMap InitMinSharedMap(int kmerSize, int windowSize)
{
    MinSharedMap res;

    std::istringstream tokenStream(pMinText);
    std::string ts;

    while (getline(tokenStream, ts, '\n')) {
	auto tmp = parseMinTextLine(ts);
	auto k = std::get<0>(tmp);
	auto w = std::get<1>(tmp);
	auto p = std::get<2>(tmp);
	auto e1 = std::get<3>(tmp);
	auto e2 = std::get<4>(tmp);

	if (k == kmerSize && std::abs(w - windowSize) <= 2) {
	    ErrTarget keyOne{e1, e2};
	    ErrTarget keyTwo{e2, e1};

	    res[keyOne] = p;
	    res[keyTwo] = p;
	}
    }

    return res;
}

std::tuple<int, int, double, double, double> parseMinTextLine(
    const std::string& line)
{
    std::istringstream tokenStream(line);

    int k, w;
    double p, e1, e2;

    tokenStream >> k;
    tokenStream >> w;
    tokenStream >> p;
    tokenStream >> e1;
    tokenStream >> e2;

    return std::make_tuple(k, w, p, e1, e2);
}

double GetPMinShared(double e1, double e2, const MinSharedMap& msm)
{
    e1 = round(e1, 2);
    e2 = round(e2, 2);

    if (e1 > 0.15) {
	e1 = 0.15;
    }
    if (e1 < 0.01) {
	e1 = 0.01;
    }

    if (e2 > 0.15) {
	e2 = 0.15;
    }
    if (e2 < 0.01) {
	e2 = 0.01;
    }

    ErrTarget key{e1, e2};
    auto it = msm.find(key);
    if (it == msm.end()) {
	throw std::invalid_argument("Empirical probability lookup failure!");
    }
    else {
	return it->second;
    }
    return -1;
}
