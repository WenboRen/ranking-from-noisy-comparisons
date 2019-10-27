#pragma once

#include<vector>
#include<algorithm>
#include<math.h>
#include<random>
#include<stdio.h>
#include<queue>
using namespace std;

#include "HelperFunctions.h"

#define PI 3.1415927

bool AttemptingComparison(const double pij, const double epsilon, const double delta, double *complexity) {
	double bmax = ceil((log(2.0 / delta) / (2.0 * epsilon * epsilon)));
	int wi = 0;
	double pe;
	for (int t = 1; t <= bmax; t++) {
		if (Compare(pij))
			wi++;
		*complexity += 1;
		pe = (double)wi / (double)t;
		double bt = sqrt(log(PI*PI*(double)t*(double)t / (3.0 * delta)) / (2.0 * (double)t));
		if (pe > 0.5 + bt) 
			return true;
		if (pe < 0.5 - bt) 
			return false;
	}
	if (pe > 0.5) return true;
	else if (pe < 0.5) return false;
	else if (rand() < RAND_MAX / 2) return true;
	else return false;
}

struct BIT {
	BIT* lchild;
	BIT *rchild;
	BIT* parent;
	int left;
	int right;
	int counter;
};

//Inserted = true; Unsure = false
bool AttemptingInsertion(const int i, vector<int> &S, const vector<vector<double>> &P, const double epsilon, const double delta, double *complexity) {
	int n = S.size(); //S is ascending
	int h = 1 + ceil(log2(1.0000 + (double)n));
	if (n == 2)
		int agsadertge = 1 + 1;
	int tmax = ceil(max(4.0 * h, 512.0 * log(2.0 / delta) / 25.0));
	double q = 15.0 / 16.0;
	double q2 = sqrt(q), q3 = pow(q, 1.0 / 3.0);
	vector<double> Pis; //-infty, S[0], S[1],..., S[n-1], +infty
	Pis.push_back(1.0);
	for (int j = 1; j <= n; j++)
		Pis.push_back(P[i-1][S[j - 1]-1]);
	Pis.push_back(0.0);

	BIT* root = new BIT(); root->parent = NULL; root->counter = 0; 
	root->left = 0; root->right = 1 + n;

	queue<BIT*> bq;
	bq.push(root);
	while (!bq.empty()) {
		BIT* current = bq.front();
		bq.pop();
		if (current->right - current->left > 1) {
			int mid = (current->left + current->right) / 2;
			BIT *ltemp = new BIT(), *rtemp = new BIT();
			ltemp->parent = current;
			ltemp->left = current->left;
			ltemp->right = mid;
			ltemp->counter = 0;
			current->lchild = ltemp;
			bq.push(ltemp);

			rtemp->parent = current;
			rtemp->left = mid;
			rtemp->right = current->right;
			rtemp->counter = 0;
			current->rchild = rtemp;
			bq.push(rtemp);
		}
		else {
			current->lchild = NULL;
			current->rchild = NULL;
		}
	}

	BIT* X = root;
	for (int t = 1; t <= tmax; t++) {
		int Xmid = (X->left + X->right) / 2;
		if (X->left == 0 && X->right == n + 1) {
			double pij = Pis[Xmid];
			if (AttemptingComparison(pij, epsilon, 1 - q, complexity))
				X = X->rchild;
			else
				X = X->lchild;
		}
		else if (X->right - X->left == 1) { //i.e., leaf node
			if (AttemptingComparison(Pis[X->left], epsilon, 1 - q2, complexity) && !AttemptingComparison(Pis[X->right], epsilon, 1 - q2, complexity)){
				X->counter++;
				if (X->counter > 1.0 + 0.5 * (double)t + sqrt(0.5 * (double)t * log(PI*PI*(double)t*(double)t) / (3.0*delta))) {
					S.insert(S.begin() + X->left, i);
					return true;
				}
			}
			else {
				if (X->counter > 0)
					X->counter--;
				else
					X = X->parent;
			}
		}
		else {
			if (!AttemptingComparison(Pis[X->left], epsilon, 1 - q3, complexity) || AttemptingComparison(Pis[X->right], epsilon, 1 - q3, complexity))
				X = X->parent;
			else if (AttemptingComparison(Pis[Xmid], epsilon, 1 - q3, complexity))
				X = X->rchild;
			else
				X = X->lchild;
		}
	}

	bool inserted = false;
	while (!bq.empty()) bq.pop();
	bq.push(root);
	while (!bq.empty()) {
		BIT *current = bq.front();
		if (current->right - current->left == 1) {
			if (current->counter >= 1 + 5.0*tmax / 16.0) {
				S.insert(S.begin() + current->left, i);
				inserted = true;
			}
		}
		else {
			bq.push(current->lchild);
			bq.push(current->rchild);
		}
		delete current;
		bq.pop();
	}
	return inserted;
}

void IterativeAttemptingInsertion(const int i, vector<int> &S, const vector<vector<double>> &P, const double delta, double *complexity) {
	int t = 1;
	bool flag = false;
	while (!flag) {
		double epsilont = pow(0.5, (double)(1 + t));
		double deltat = 6.0 * delta / (PI * PI * (double)t * (double)t);
		flag = AttemptingInsertion(i, S, P, epsilont, deltat, complexity);
		t++;
	}
}

vector<int> IterativeInsertionRanking(const vector<vector<double>> &P, const double delta, double *complexity) {
	int n = P.size();
	vector<int> S;
	S.push_back(1);
	if (n <= 1) {
		return S;
	}
	double delta1 = delta / (double) n;

	for (int i = 2; i <= n; i++) {
		IterativeAttemptingInsertion(i, S, P, delta1, complexity);
	}
	return S;
}

bool LILCompare(const double p, const double delta, const double epsilon, const double lambda, const double beta, double *complexity) {
	double ratio = (1.0 + beta) * (1.0 + sqrt(epsilon));//To avoid repreated computation
	int T1 = 0, T2 = 0;//T2 is the 1/2-arm with deterministic rewards
	double S1 = 0, S2 = 0;
	const double sigma2 = 0.25;
	T1++; T2++;
	if (Compare(p)) S1 += 1.0;
	S2 += 0.5;

	*complexity += 1.0;

	int t = 1;
	while (true) {
		if (T1 >= 1 + lambda * T2 || T2 >= 1 + lambda * T1) {
			if (T1 > T2) return true;
			else return false;
		}
		double mu1 = S1 / (double)T1;
		double mu2 = S2 / (double)T2;
		double a1 = mu1 + ratio * sqrt(2.0*sigma2*(1.0 + epsilon) * log((fabs(log((1.0 + epsilon)*(double)T1) + 2.0) / delta)) / (double)T1);
		double a2 = mu2 + ratio * sqrt(2.0*sigma2*(1.0 + epsilon) * log(fabs((log((1.0 + epsilon)*(double)T2) + 2.0) / delta)) / (double)T2);
		if (a1 > a2) {
			T1++;
			if (Compare(p)) S1 += 1.0;
		}
		else {
			T2++;
			S2 += 0.5;
		}
		t++;
		*complexity += 1.0;
	}
}

bool LILInsertion(const int i, vector<int> &S, const vector<vector<double>> &P, const double delta, const double epsilon, const double alpha, const double beta, double *complexity) {
	int n = S.size(); //S is ascending
	double q = 2.0 / 3.0;
	double q2 = sqrt(q), q3 = pow(q, 1.0 / 3.0);
	vector<double> Pis; //-infty, S[0], S[1],..., S[n-1], +infty
	Pis.push_back(1.0);
	for (int j = 1; j <= n; j++)
		Pis.push_back(P[i - 1][S[j - 1] - 1]);
	Pis.push_back(0.0);

	const double a = 1.0;
	const double b = 6.0 * log((double)n) + 18.0 * log(1 / delta);
	const double c = log((double)n) * log((double)n);
	const double m = (b + sqrt(b * b - 4 * a * c)) / (2.0 * a);

	BIT* root = new BIT(); root->parent = NULL; root->counter = 0;
	root->left = 0; root->right = 1 + n;

	queue<BIT*> bq;
	bq.push(root);
	while (!bq.empty()) {
		BIT* current = bq.front();
		bq.pop();
		if (current->right - current->left > 1) {
			int mid = (current->left + current->right) / 2;
			BIT *ltemp = new BIT(), *rtemp = new BIT();
			ltemp->parent = current;
			ltemp->left = current->left;
			ltemp->right = mid;
			ltemp->counter = 0;
			current->lchild = ltemp;
			bq.push(ltemp);

			rtemp->parent = current;
			rtemp->left = mid;
			rtemp->right = current->right;
			rtemp->counter = 0;
			current->rchild = rtemp;
			bq.push(rtemp);
		}
		else {
			current->lchild = NULL;
			current->rchild = NULL;
		}
	}

	BIT* X = root;
	int t = 1;
	while(t < m){
		int Xmid = (X->left + X->right) / 2;
		if (X->left == 0 && X->right == n + 1) {
			double pij = Pis[Xmid];
			if (LILCompare(pij, 1 - q, epsilon, alpha, beta, complexity))
				X = X->rchild;
			else
				X = X->lchild;
		}
		else if (X->right - X->left == 1) { //i.e., leaf node
			if (LILCompare(Pis[X->left], 1 - q2, epsilon, alpha, beta, complexity) && !LILCompare(Pis[X->right], 1 - q2, epsilon, alpha, beta, complexity)) {
				X->counter++;
				//if (X->counter > log(((double)n))) {
					//S.insert(S.begin() + X->left, i);
					//break;
				//}
			}
			else {
				if (X->counter > 0)
					X->counter--;
				else
					X = X->parent;
			}
		}
		else {
			if (!LILCompare(Pis[X->left], 1 - q3, epsilon, alpha, beta, complexity) || LILCompare(Pis[X->right], 1 - q3, epsilon, alpha, beta, complexity))
				X = X->parent;
			else if (LILCompare(Pis[Xmid],  1 - q3, epsilon, alpha, beta, complexity))
				X = X->rchild;
			else
				X = X->lchild;
		}
		t++;
	}

	bool inserted = false;
	while (!bq.empty()) bq.pop();
	bq.push(root);
	while (!bq.empty()) {
		BIT *current = bq.front();
		if (current->right - current->left > 1) {
			bq.push(current->lchild);
			bq.push(current->rchild);
		}
		else if (current->right - current->left == 1 && current->counter > 0) {
			S.insert(S.begin() + X->left, i);
			return true;
		}
		delete current;
		bq.pop();
	}
}

vector<int> LILUCB_BinarySearch(const vector<vector<double>> &P, const double delta, const double epsilon, const double alpha, const double beta, double *complexity) {
	int n = P.size();
	vector<int> S;
	S.push_back(1);
	if (n <= 1) {
		return S;
	}
	double delta1 = delta / (double)n;

	for (int i = 2; i <= n; i++) {
		bool inserted = LILInsertion(i, S, P, delta1, epsilon, alpha, beta, complexity);
		if (!inserted) {
			vector<int> S2(n, 0);
			return S2;
		}
	}
	return S;
}

vector<int> ActiveRanking(const vector<vector<double>> &P, const double delta, double *complexity) {
	int n = P.size();
	vector<int> R(n,0), T, remain(n, 0);
	vector<double> estimated(n, 0.0);
	int t = 0;
	vector<bool> assured(n, false);
	int numAssured = 0;
	vector<pair<double, int>> fs;
	while (numAssured < n) {
		t++;
		for (int i = 0; i < n; i++) {
			if (assured[i]) continue;
			int j;
			while ((j = rand() % n) == i);
			int winner = Compare(P[i][j]) ? i : j;
			if (winner == i)
				estimated[i] = (1.0 + estimated[i] * (double)(t - 1)) / (double)t;
			else estimated[i] = estimated[i] * (double)(t - 1) / (double)t;
			*complexity += 1.0;
		}
		//arrange remaining items in the ascending order of empirical means
		fs.clear();
		for (int i = 0; i < n; i++) fs.push_back(make_pair(estimated[i], i));
		sort(fs.begin(), fs.end());
		double alpha_t = sqrt(log(125.0 * n * log(1.12 * (double)t) / delta) / (double)t);

		for (int i = 0; i < n ; ++i) {
			if (assured[i]) continue;

			bool left, right;

			//check right
			int j = i + 1;
			while (j < n && assured[j]) ++j;
			if (j >= n) 
				right = true;
			else if (estimated[fs[j].second] - estimated[fs[i].second] > 4.0 * alpha_t) 
				right = true;
			else 
				right = false;

			//check left
			j = i - 1;
			while (j >= 0 && assured[j]) --j;
			if (j < 0) 
				left = true;
			else if (estimated[fs[i].second] - estimated[fs[j].second] > 4.0 * alpha_t) 
				left = true;
			else 
				left = false;

			if (left && right) {
				assured[i] = true;
				++numAssured;
				R[i] = fs[i].second;
			}
		}
	}
	return R;
}

vector<int> BQS(const vector<double> &thetas, vector<int> A, const int B, int *remain) {
	int n = A.size();
	//randomReorder(A);
	if (B <= 0)
		return A;
	if (n <= 1) return A;
	int mid = rand() % n;
	vector<int> A0;
	vector<int> A1;
	vector<int> R;
	vector<bool> assigned(n, false);
	assigned[mid] = true;
	int numAssigned = 1;
	while (numAssigned < n) {
		int j = rand() % n;
		if (assigned[j] || j == mid) continue; //A[j] has been assigned to mid, A0, or A1
											   //for(int j = 0; j < n; j++){
		double r = ((double)rand() / (double)RAND_MAX) * (thetas[A[mid]] + thetas[A[j]]);
		if (r <= thetas[A[mid]]) A1.push_back(A[j]);
		else A0.push_back(A[j]);
		numAssigned++;
		assigned[j] = true;
	}
	*remain -= n - 1;
	BQS(thetas, A0, *remain, remain);
	BQS(thetas, A1, *remain, remain);
	for (int i = 0; i < A1.size(); i++) {
		R.push_back(A1[i]);
	}
	R.push_back(A[mid]);
	for (int i = 0; i < A0.size(); i++) {
		R.push_back(A0[i]);
	}
	return R;
}

vector<int> PLPAC_AMPR(const vector<double>& thetas, const double delta, const double epsilon, double* complexity) {
	int n = thetas.size();
	vector<vector<double>> P(n, vector<double>(n, 0.0));
	vector<vector<double>> N(n, vector<double>(n, 0.0));
	vector<vector<double>> W(n, vector<double>(n, 0.0));
	vector<bool> removed(n, false);
	vector<int> blow(n, 0);
	vector<int> bhigh(n, 0);
	vector<vector<double>> c(n, vector<double>(n, 0.0));
	vector<vector<int>> G(n, vector<int>(n, 0));

	for (int i = 0; i < n; i++) {
		blow[i] = 0.0;
		bhigh[i] = double(n);
	}
	bool completed = false;

	while (!completed) {
		for (int i = 0; i < n; i++) {
			for (int j = i+1; j < n; j++) {
				if (intervalOverlapped(blow[i], bhigh[i], blow[j], bhigh[j]))
					G[i][j] = G[j][i] = 1;
				else G[i][j] = G[j][i] = 0;
			}
		}
		queue<int> q;
		vector<bool> visited(n, false);
		int numVisited = 0;
		vector<int> component;
		while (numVisited < n) {
			component.clear();
			for (int i = 0; i < n; i++) {
				if (!visited[i]) {
					q.push(i);
					numVisited++;
					visited[i] = true;
				}
			}
			while (q.size() > 0) {
				int i = q.front();
				q.pop();
				component.push_back(i);
				for (int j = 0; j < n; j++) {
					if (!visited[j] && G[i][j] == 1) {
						q.push(j);
						numVisited++;
						visited[j] = true;
					}
				}
			}
			//Now we have found a component
			if (component.size() > 1) {
				int budget = 1 + (int)(3.0 * (component.size() + 1.0) * log(component.size()));
				int remain = budget;
				vector<int> R = BQS(thetas, component, budget, &remain);
				*complexity += budget - remain;
				for (int i = 0, len = R.size(); i < len; i++) {
					for (int j = i + 1; j < len; j++) {
						W[R[i]][R[j]] += 1.0;
						N[R[i]][R[j]] += 1.0;
						N[R[j]][R[i]] += 1.0;
					}
				}
			}
		}

/*		int budget = 1 + (int)(3.0 * (n + 1.0) * log(n));
		int remain = budget;
		vector<int> nset;
		for (int i = 0; i < n; i++) nset.push_back(i);
		vector<int> R = BQS(thetas, nset, budget, &remain);
		*complexity += budget - remain;
		for (int i = 0, len = R.size(); i < len; i++) {
			for (int j = i + 1; j < len; j++) {
				W[R[i]][R[j]] += 1.0;
				N[R[i]][R[j]] += 1.0;
				N[R[j]][R[i]] += 1.0;
			}
		}*/

		//All components have been BQSed
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					c[i][j] = sqrt(log(4.0 * (double)n * (double)n * N[i][j] * N[i][j] / delta) / (2.0 * N[i][j]));
					P[i][j] = W[i][j] / N[i][j];
				}
			}
		}
		for (int i = 0; i < n; i++) {
			blow[i] = bhigh[i] = 0;
			for (int j = 0; j < n; j++) {
				if (j != i && P[i][j] - c[i][j] > 0.5)
					blow[i]++;
			}
			bhigh[i] = blow[i];
			for (int j = 0; j < n; j++) {
				if (j != i && 0.5 >= P[i][j] - c[i][j] && 0.5 <= P[i][j] + c[i][j])
					bhigh[i]++;
			}
		}
		//Now determine whether ranking has been completed
		completed = true;
		for (int i = 0; i < n && completed; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j && intervalOverlapped(blow[i], bhigh[i], blow[j], bhigh[j])) {
					if (!(0.5 - epsilon < P[i][j] - c[i][j] && 0.5 + epsilon > P[i][j] + c[i][j])) {
						completed = false;
						break;
					}
				}
			}
		}
	}
	vector<int> R;
	for (int b = 0; b < n; b++) {
		for (int i = 0; i < n; i++) {
			if (blow[i] == b)
				R.push_back(i);
		}
	}
	return R;
}

bool HoeffdingCompare(const double p, const double delta, const double gamma, double *complexity) {
	int T1 = 0, T2 = 0;//T2 is the 1/2-arm with deterministic rewards
	double S1 = 0, S2 = 0;
	T1++; T2++;
	if (Compare(p)) S1 += 1.0;
	S2 += 0.5;
	*complexity += 1.0;
	double hat_mu2 = 0.5;
	const double k1 = 2.0 * (1.0 + 1.0 / (gamma - 1.0));

	int t = 2;
	while (true) {
		double hat_mu1 = S1 / (double)T1;
		double alpha = sqrt(log(k1 * pow((double)t, gamma) / delta) / (2.0 * (double)t));
		if (hat_mu1 - hat_mu2 > alpha) 
			return true;
		else if (hat_mu2 - hat_mu1 > alpha)
			return false;

		++T1; ++T2;
		if (Compare(p)) S1 += 1.0;
		S2 += 0.5;

		*complexity += 1.0;
		t += 2;
	}
}

bool HoeffdingInsertion(const int i, vector<int> &S, const vector<vector<double>> &P, const double delta, const double gamma, double *complexity) {
	int n = S.size(); //S is ascending
	double q = 2.0 / 3.0;
	double q2 = sqrt(q), q3 = pow(q, 1.0 / 3.0);
	vector<double> Pis; //-infty, S[0], S[1],..., S[n-1], +infty
	Pis.push_back(1.0);
	for (int j = 1; j <= n; j++)
		Pis.push_back(P[i - 1][S[j - 1] - 1]);
	Pis.push_back(0.0);

	BIT* root = new BIT(); root->parent = NULL; root->counter = 0;
	root->left = 0; root->right = 1 + n;

	const double a = 1.0;
	const double b = 6.0 * log((double)n) + 18.0 * log(1 / delta);
	const double c = log((double)n) * log((double)n);
	const double m = (b + sqrt(b * b - 4 * a * c)) / (2.0 * a);

	queue<BIT*> bq;
	bq.push(root);
	while (!bq.empty()) {
		BIT* current = bq.front();
		bq.pop();
		if (current->right - current->left > 1) {
			int mid = (current->left + current->right) / 2;
			BIT *ltemp = new BIT(), *rtemp = new BIT();
			ltemp->parent = current;
			ltemp->left = current->left;
			ltemp->right = mid;
			ltemp->counter = 0;
			current->lchild = ltemp;
			bq.push(ltemp);

			rtemp->parent = current;
			rtemp->left = mid;
			rtemp->right = current->right;
			rtemp->counter = 0;
			current->rchild = rtemp;
			bq.push(rtemp);
		}
		else {
			current->lchild = NULL;
			current->rchild = NULL;
		}
	}

	BIT* X = root;
	int t = 1;
	while (t < m) {
		int Xmid = (X->left + X->right) / 2;
		if (X->left == 0 && X->right == n + 1) {
			double pij = Pis[Xmid];
			if (HoeffdingCompare(pij, 1 - q, gamma, complexity))
				X = X->rchild;
			else
				X = X->lchild;
		}
		else if (X->right - X->left == 1) { //i.e., leaf node
			if (HoeffdingCompare(Pis[X->left], 1 - q2, gamma, complexity) && !HoeffdingCompare(Pis[X->right], 1 - q2,gamma, complexity)) {
				X->counter++;
			}
			else {
				if (X->counter > 0)
					X->counter--;
				else
					X = X->parent;
			}
		}
		else {
			if (!HoeffdingCompare(Pis[X->left], 1 - q3, gamma, complexity) || HoeffdingCompare(Pis[X->right], 1 - q3, gamma, complexity))
				X = X->parent;
			else if (HoeffdingCompare(Pis[Xmid], 1 - q3, gamma, complexity))
				X = X->rchild;
			else
				X = X->lchild;
		}
		t++;		
	}

	bool inserted = false;
	while (!bq.empty()) bq.pop();
	bq.push(root);
	while (!bq.empty()) {
		BIT *current = bq.front();
		if (current->right - current->left > 1) {
			bq.push(current->lchild);
			bq.push(current->rchild);
		}
		else if (current->right - current->left == 1 && current->counter > 0) {
			S.insert(S.begin() + X->left, i);
			return true;
		}
		delete current;
		bq.pop();
	}
	return false;
}

vector<int> Hoeffding_LUCB_BinarySearch(const vector<vector<double>> &P, const double delta, double *complexity) {
	const double gamma = 2.0;
	int n = P.size();
	vector<int> S;
	S.push_back(1);
	if (n <= 1) {
		return S;
	}
	double delta1 = delta / (double)n;

	for (int i = 2; i <= n; i++) {
		bool inserted = HoeffdingInsertion(i, S, P, delta1, gamma, complexity);
		if (!inserted) {
			vector<int> S2(n, 0);
			return S2;
		}
	}
	return S;
}

double KLUpperBound(const double p, const double beta, const double epsilon) {
	double theta = p;
	double dis;
	while (Dkl(p, theta) < beta) {
		theta += epsilon;
		if (theta >= 1.0)
			return 1.0;
	}
	return theta;
}

double KLLowerBound(const double p, const double beta, const double epsilon) {
	double theta = p;
	while (Dkl(p, theta) < beta) {
		theta -= epsilon;
		if (theta <= 0.0)
			return 0.0;
	}
	return theta;
}

bool KLCompare(const double p,  const double delta, double* complexity) {
	double wins1 = 0.0;
	double nums1 = 0.0;
	double ucb = 0.0;
	double lcb = 0.0;
	double t = 0;
	double _epsilon = 0.01;

	double gamma = 2.0;
	double k1 = 0.01 + 2.0 * (1.0 + 1.0 / (gamma - 1.0));

	t = 1;
	if (Compare(p)) 
		wins1 += 1.0;
	nums1 += 1.0;
	*complexity += 1.0;
	ucb = KLUpperBound(wins1 / nums1, log(k1 * pow((double)t, gamma) * 1.0 / delta) / 1.0, _epsilon);
	lcb = KLLowerBound(wins1 / nums1, log(k1 * pow((double)t, gamma) * 1.0 / delta) / 1.0, _epsilon);

	double B = INFINITY;

	while (true) {
		if (ucb - _epsilon < 0.5)
			return false;
		else if (lcb + _epsilon > 0.5)
			return true;

		nums1 += 1.0;
		if (Compare(p)) wins1 += 1.0;
		++t; *complexity += 1.0;
		ucb = KLUpperBound(wins1 / nums1, log(k1 * pow(t, gamma) * 1.0 / delta) / nums1, _epsilon);
		lcb = KLLowerBound(wins1 / nums1, log(k1 * pow(t, gamma) * 1.0 / delta) / nums1, _epsilon);
	}
}

bool KLLUCBInsertion(const int i, vector<int> &S, const vector<vector<double>> &P, const double delta, double *complexity) {
	int n = S.size(); //S is ascending
	double q = 2.0 / 3.0;
	double q2 = sqrt(q), q3 = pow(q, 1.0 / 3.0);
	vector<double> Pis; //-infty, S[0], S[1],..., S[n-1], +infty
	Pis.push_back(1.0);
	for (int j = 1; j <= n; j++)
		Pis.push_back(P[i - 1][S[j - 1] - 1]);
	Pis.push_back(0.0);

	BIT* root = new BIT(); root->parent = NULL; root->counter = 0;
	root->left = 0; root->right = 1 + n;

	const double a = 1.0;
	const double b = 6.0 * log((double)n) + 18.0 * log(1 / delta);
	const double c = log((double)n) * log((double)n);
	const double m = (b + sqrt(b * b - 4 * a * c)) / (2.0 * a);

	queue<BIT*> bq;
	bq.push(root);
	while (!bq.empty()) {
		BIT* current = bq.front();
		bq.pop();
		if (current->right - current->left > 1) {
			int mid = (current->left + current->right) / 2;
			BIT *ltemp = new BIT(), *rtemp = new BIT();
			ltemp->parent = current;
			ltemp->left = current->left;
			ltemp->right = mid;
			ltemp->counter = 0;
			current->lchild = ltemp;
			bq.push(ltemp);

			rtemp->parent = current;
			rtemp->left = mid;
			rtemp->right = current->right;
			rtemp->counter = 0;
			current->rchild = rtemp;
			bq.push(rtemp);
		}
		else {
			current->lchild = NULL;
			current->rchild = NULL;
		}
	}

	BIT* X = root;
	int t = 1;
	while (t < m) {
		int Xmid = (X->left + X->right) / 2;
		if (X->left == 0 && X->right == n + 1) {
			double pij = Pis[Xmid];
			if (KLCompare(pij, 1 - q, complexity))
				X = X->rchild;
			else
				X = X->lchild;
		}
		else if (X->right - X->left == 1) { //i.e., leaf node
			if (KLCompare(Pis[X->left], 1 - q2, complexity) && !KLCompare(Pis[X->right], 1 - q2, complexity)) {
				X->counter++;
				//if (X->counter > log(((double)n))) {
				//S.insert(S.begin() + X->left, i);
				//break;
				//}
			}
			else {
				if (X->counter > 0)
					X->counter--;
				else
					X = X->parent;
			}
		}
		else {
			if (!KLCompare(Pis[X->left], 1 - q3, complexity) || KLCompare(Pis[X->right], 1 - q3, complexity))
				X = X->parent;
			else if (KLCompare(Pis[Xmid], 1 - q3, complexity))
				X = X->rchild;
			else
				X = X->lchild;
		}
		t++;
	}

	bool inserted = false;
	while (!bq.empty()) bq.pop();
	bq.push(root);
	while (!bq.empty()) {
		BIT *current = bq.front();
		if (current->right - current->left > 1) {
			bq.push(current->lchild);
			bq.push(current->rchild);
		}
		else if (current->right - current->left == 1 && current->counter > 0) {
			S.insert(S.begin() + X->left, i);
			return true;
		}
		delete current;
		bq.pop();
	}
	return false;
}

vector<int>KL_LUCB_BinarySearch(const vector<vector<double>> &P, const double delta, double *complexity) {
	const double gamma = 2.0;
	int n = P.size();
	vector<int> S;
	S.push_back(1);
	if (n <= 1) {
		return S;
	}
	double delta1 = delta / (double)n;

	for (int i = 2; i <= n; i++) {
		bool inserted = KLLUCBInsertion(i, S, P, delta1, complexity);
		if (!inserted) {
			vector<int> S2(n, 0);
			return S2;
		}
	}
	return S;
}

