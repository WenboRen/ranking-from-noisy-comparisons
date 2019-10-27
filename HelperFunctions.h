#pragma once

#include<vector>
#include <algorithm> 
using namespace std;

bool operator == (const vector<double>& v1, const vector<double>& v2) {
	int len = min(v1.size(), v2.size());
	for (int i = 0; i < len; i++)
		if (v1[i] != v2[i])
			return false;
	return true;
}

double doubleRand(double lo, double hi) {
	return rand() / (double)RAND_MAX * (hi - lo) + lo;
}

bool intervalOverlapped(int l1, int r1, int l2, int r2) {
	return r1 >= l2 && r2 >= l1;
}

void randomReorder(vector<int> &S) {
	for (int t = 0, len = S.size(); t < len * len; t++) {
		int i = rand() % len;
		int j = rand() % len;
		swap(S[i], S[j]);
	}
}

void randomReorder(vector<double> &S, int r) {
	for (int t = 0; t <r * r; t++) {
		int i = rand() % r;
		int j = rand() % r;
		swap(S[i], S[j]);
	}
}

vector<double> stationaryDistribution(vector<vector<double>> &P, int iterations) {
	int n = P.size();
	vector<double> pi1(n, 0.0), pi2(n, 0.0);
	double init_value = 1.0 / n;
	for (int i = 0; i < n; i++) pi1[i] = init_value;
	while (iterations > 0) {
		for (int i = 0; i < n; i++) {
			pi2[i] = 0.0;
			for (int j = 0; j < n; j++) {
				pi2[i] += pi1[j] * P[j][i];
			}
		}
		pi1 = pi2;
		iterations--;
	}
	return pi1;
}

bool fullyConnected(vector<vector<bool>> &edges) {
	int n = edges.size();
	queue<int> q;
	vector<bool> visited(n, false);
	q.push(0);
	visited[0] = true;
	int numVisited = 1;
	while (q.size() > 0) {
		int top = q.front();
		for (int i = 0; i < n; i++) {
			if (edges[top][i] && !visited[i]) {
				q.push(i);
				visited[i] = true;
				numVisited++;
			}
		}
		q.pop();
	}
	return numVisited >= n;
}


bool Compare(const double pij) {
	if (rand() < pij * (double)RAND_MAX)
		return true;
	else return false;
}

double dlog(const double p, const double q) {
	if (p < 0.0000001) return 0;
	if (p > 10000000.0 * q) return INFINITY;
	return log(p / q);
}

double Dkl(const double p, const double q) {
	double x = p*dlog(p, q) + (1 - p)*dlog(1 - p, 1 - q);
	return(x);
}
