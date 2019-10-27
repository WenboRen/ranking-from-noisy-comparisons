#include <stdio.h>
#include <tchar.h>

#include<iostream>
#include "RankingAlgorithms.h"
#include "HelperFunctions.h"
#include<vector>
#include<algorithm>
#include<math.h>
#include<time.h>
using namespace std;


int main()
{
	//Choose an algorithm you want to test
	string Alg = "IterativeInsertionRanking";
	//string Alg = "LILUCB_BinarySearch";
	//string Alg = "ActiveRanking";
	//string Alg = "PLPAC_AMPR";
	//string Alg = "Hoeffding_LUCB_BinarySearch";
	//string Alg = "KL_LUCB_BinarySearch";
	//Choose an algorithm you want to test	

	//You may change the way to generate the data, i.e., p_{i,j}'s
	//string mode = "MNL";
	string mode = "Random";
	//string mode = "Homo";
	//string mode = "Easy";
	//You may change the way to generate the data, i.e., p_{i,j}'s

	int trials = 100;//Number of repetitions
	double delta = 0.01;//confidence

	srand(20191208);//Fix a seed so that the instances remain the same for tests of different algorithms

	//The number of n values you want to test
	vector<int> Ns = { 4,6,8,10,15,20,25,30,40,50,60,70,80,90,100 };

	//We do not provide datasets but provide the codes to generate the datasets
	vector<double> thetas;
	for (int i = 0; i < 100; i++) {
		double li = 0.9 * pow(1.2 / 0.8, i), ri = 1.1 * pow(1.2 / 0.8, i);
		thetas.push_back(li + (ri - li) *(double(rand()) / (double)RAND_MAX));
		//thetas.push_back(pow(1.2 / 0.8, i));
	}

	vector<double> randomNumberList;
	int randomNumberIterator = 0.0;
	if (mode == "Random") {
		const double Delta = 0.1;
		randomNumberList = vector<double>(40000, 0.0);
		for (int i = 0; i < 40000; i++)
			randomNumberList[i] = 0.8 * Delta + 0.7 * Delta * ((double)rand() / (double)RAND_MAX);
	}

	srand((unsigned int)(time(0)));

	for (int n : Ns) {
		randomReorder(thetas, n);
		vector<vector<double>> P(n, vector<double>(n, 0.0));
		vector<int> correctRanking;
		if (mode == "Homo") {
			for (int i = 0; i < n; i++)
				thetas.push_back((double(rand()) / (double)RAND_MAX));
			vector<pair<double, int>> fs;
			for (int i = 0; i < n; i++)
				fs.push_back(make_pair(thetas[i], i + 1));
			sort(fs.begin(), fs.end());
			for (int i = 0; i < n; i++)
				correctRanking.push_back(fs[i].second);//ascending
			double Delta = 0.1;
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					if (thetas[i] > thetas[j]) P[i][j] = 0.5 + Delta; else P[i][j] = 0.5 - Delta;
				}
			}
		}
		else if (mode == "MNL") {
			randomReorder(thetas, n);
			vector<pair<double, int>> fs;
			for (int i = 0; i < n; i++)
				fs.push_back(make_pair(thetas[i], i + 1));
			sort(fs.begin(), fs.end());
			for (int i = 0; i < n; i++)
				correctRanking.push_back(fs[i].second);//ascending
			double Delta = 0.1;
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					P[i][j] = thetas[i] / (thetas[i] + thetas[j]);
				}
			}
		}

		else if (mode == "Easy") {
			vector<int> A;
			for (int i = 0; i < n; i++)
				A.push_back(i + 1);
			randomReorder(A);//ascending
			correctRanking = A;
			double Delta = 0.1;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (i + 1 == j)
						P[A[i] - 1][A[j] - 1] = 0.5 - Delta;
					else if (i == j + 1)
						P[A[i] - 1][A[j] - 1] = 0.5 + Delta;
					else if (i < j)
						P[A[i] - 1][A[j] - 1] = 0.0;
					else
						P[A[i] - 1][A[j] - 1] = 1.0;
				}
			}
		}
		else if (mode == "Random") {
			vector<int> A;
			for (int i = 0; i < n; i++)
				A.push_back(i + 1);
			randomReorder(A);//ascending
			correctRanking = A;
			double Delta = 0.1;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (i < j)
						P[A[i] - 1][A[j] - 1] = 0.5 - randomNumberList[randomNumberIterator++];
					else
						P[A[i] - 1][A[j] - 1] = 0.5 + randomNumberList[randomNumberIterator++];
				}
			}
		}

		vector<double> Deltas(n);

		for (int i = 0; i < n; i++) {
			double minDelta = 1.0;
			for (int j = 0; j < n; ++j) {
				if (i != j && fabs(P[i][j] - 0.5) < minDelta) {
					minDelta = fabs(P[i][j] - 0.5);
				}
			}
			Deltas[i] = minDelta;
		}
		double SumDelta_2 = 0;
		for (int i = 0; i < n; i++) SumDelta_2 += 1.0 / (Deltas[i] * Deltas[i]);
		/////////////////////////////////////////////////////////////////////

		double Delta_min = 1.0;
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				Delta_min = min(Delta_min, fabs(P[i][j] - 0.5));
			}
		}

		double complexity = 0.0;
		double accumulated = 0.0;
		int correct = trials;
		cout << "n = " << n << ", delta = " << delta << " trials = " << trials << endl;

		vector<int> R;
		for (int t = 0; t < trials; t++) {
			complexity = 0.0;

			if (Alg == "IterativeInsertionRanking") {
				R = IterativeInsertionRanking(P, delta, &complexity);
			}
			else if (Alg == "LILUCB_BinarySearch") {
				double epsilon = 0.01, beta = 1.0, lambda = ((2.0 + beta) / beta) * ((2.0 + beta) / beta);
				R = LILUCB_BinarySearch(P, delta, epsilon, lambda, beta, &complexity);
			}
			else if (Alg == "ActiveRanking") {
				R = ActiveRanking(P, delta, &complexity); for (int i = 0; i < n; ++i) ++R[i];
			}
			else if (Alg == "PLPAC_AMPR") {
				vector<double> thetasprime; for (int k = 0; k < n; ++k) thetasprime.push_back(thetas[k]);
				R = PLPAC_AMPR(thetasprime, delta, 0.0, &complexity); for (int i = 0; i < n; ++i) ++R[i]; for (int i = 0; i < n / 2; i++) swap(R[i], R[n - 1 - i]);
			}
			else if (Alg == "Hoeffding_LUCB_BinarySearch") {
				R = Hoeffding_LUCB_BinarySearch(P, delta, &complexity);
			}
			else if (Alg == "KL_LUCB_BinarySearch") {
				R = KL_LUCB_BinarySearch(P, delta, &complexity);
			}

			//for (int i = 0; i < n; i++) cout << R[i] << "\t";
			accumulated += complexity;
			//cout << complexity << endl;
			for (int i = 0; i < n; i++) {
				if (R[i] != correctRanking[i]) {
					correct--;
					break;
				}
			}
		}
		accumulated /= (double)trials;

		cout << "aveDelta^{-2}=" << SumDelta_2 / n << "\n" << correct << " corrects out of " << trials << " trials\n" << "average complexity = " << accumulated << endl << "Complexity / H = " << (accumulated / SumDelta_2) << endl << endl;
	}

    return 0;
}



