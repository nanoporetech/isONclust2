#ifndef P_EMP_PROB_H_INCLUDED
#define P_EMP_PROB_H_INCLUDED

#include<vector>
#include<map>
#include<tuple>

struct ErrTarget {
    double Error1;
    double Error2;
};

typedef std::map<ErrTarget, double> MinSharedMap;

bool operator<(const ErrTarget& a, const ErrTarget& b);
MinSharedMap InitMinSharedMap(int kmerSize, int windowSize);
std::tuple<int, int, double, double, double> parseMinTextLine(const std::string& line);
double GetPMinShared(double e1, double e2, const MinSharedMap& msm);

#endif
