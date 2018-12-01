// PhysarumSolver.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include "pch.h"
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdbool>
#include <chrono>
#include <vector>

#define N (10) /* 地点の数 */
#define I (1) /* 出発地点から流れる量 */
#define LOOP_NUMBER (100) /* フィゾルム・ソルバーを繰り返す回数 */
#define DELTA (0.1) // 全体の重み付け def:0.1
#define ALPHA (1.0) // 流量に対する重み付け def:1.0
#define GAMMA (1.0) // 伝導率に対する重み付け def:1.0
#define TIME_INTERVAL (1.0) // 時間間隔 def:1.0

using namespace std;
using namespace std::chrono;

/* 地点間距離 */
/* ※1 地点aと地点bが同一の地点の場合は0が格納されている */
/* ※2 地点aと地点bが隣接地点でない場合は-1が格納されている */
static int DISTANCE_MAP[N][N] = {
	/* a->b 0   1   2   3   4   5   6   7   8   9 */
	/*0*/ {+0, +8, -1, -1, +4, +6, -1, -1, -1, -1},
	/*1*/ {+8, +0, +9, -1, +3, -1, -1, -1, -1, -1},
	/*2*/ {-1, +9, +0, +4, +3, -1, +4, -1, -1, +7},
	/*3*/ {-1, -1, +4, +0, -1, -1, -1, -1, -1, +3},
	/*4*/ {+4, +3, +3, -1, +0, +3, +4, -1, -1, -1},
	/*5*/ {+6, -1, -1, -1, +3, +0, +5, +2, -1, -1},
	/*6*/ {-1, -1, +4, -1, +4, +5, +0, +9, +3, +8},
	/*7*/ {-1, -1, -1, -1, -1, +2, +9, +0, +4, -1},
	/*8*/ {-1, -1, -1, -1, -1, -1, +3, +4, +0, +2},
	/*9*/ {-1, -1, +7, +3, -1, -1, +8, -1, +2, +0},
};

/* 伝導率 */
static double CONDUCTIVITY_MAP[N][N] = {
	/* a->b 0   1   2   3   4   5   6   7   8   9 */
	/*0*/ {+0, +1, +0, +0, +1, +1, +0, +0, +0, +0},
	/*1*/ {+1, +0, +1, +0, +1, +0, +0, +0, +0, +0},
	/*2*/ {+0, +1, +0, +1, +1, +0, +1, +0, +0, +1},
	/*3*/ {+0, +0, +1, +0, +0, +0, +0, +0, +0, +1},
	/*4*/ {+1, +1, +1, +0, +0, +1, +1, +0, +0, +0},
	/*5*/ {+1, +0, +0, +0, +1, +0, +1, +1, +0, +0},
	/*6*/ {+0, +0, +1, +0, +1, +1, +0, +1, +1, +1},
	/*7*/ {+0, +0, +0, +0, +0, +1, +1, +0, +1, +0},
	/*8*/ {+0, +0, +0, +0, +0, +0, +1, +1, +0, +1},
	/*9*/ {+0, +0, +1, +1, +0, +0, +1, +0, +1, +0},
};

typedef const char* String; // String型を定義
typedef int Point;
typedef int Distance;
typedef double Quantity;

struct Route { //通信経路
	Point start, destination; // start:通信元 , destination:通信先
	Distance distance = 0; // 通信距離
	vector<Point> point;
	Quantity quantity[N][N] = {}; // ネットワーク全体の流量の値
};

/* ガウスの削除法計算ルーチン */
void gaus(double a[N][N], double x[N], double b[N]) {

	for (int k = 0; k < N; k++) {
		for (int i = k + 1; i < N; i++) {
			double alpha = a[i][k] / a[k][k];
			for (int j = k + 1; j < N; j++) {
				a[i][j] -= alpha * a[k][j];
			}
			b[i] -= alpha * b[k];
		}
	}

	x[N - 1] = b[N - 1] / a[N - 1][N - 1];

	for (int k = N - 2; k >= 0; k--) {
		for (int j = k + 1; j < N; j++) {
			b[k] -= a[k][j] * x[j];
		}
		x[k] = b[k] / a[k][k];
	}
}

void gauss(double a[N][N], double x[N], double b[N]) {
	double copy, akk, aik;
	int k, max;

	/* ピポット選択 */
	for (k = 0; k <= N - 2; k++) {
		max = k;

		for (int i = k + 1; i <= N - 1; i++) {
			if (fabs(a[i][k]) > fabs(a[max][k])) max = i;
		}

		/* 例外処理（対角成分が小さな値になったとき） */
		if (fabs(a[max][k]) < 1.0e-7) exit(0);

		/* 行の入れ替え */
		if (max != k) {
			for (int j = 0; j <= N - 1; j++) {
				copy = a[k][j];
				a[k][j] = a[max][j];
				a[max][j] = copy;
			}
			copy = b[k];
			b[k] = b[max];
			b[max] = copy;
		}

		/* 前進消去 */
		akk = a[k][k];

		for (int j = k; j <= N - 1; j++) {
			a[k][j] /= akk;
		}
		b[k] /= akk;

		for (int i = k + 1; i <= N - 1; i++) {
			aik = a[i][k];

			for (int j = k; j <= N - 1; j++) {
				a[i][j] -= aik * a[k][j];
			}
			b[i] -= aik * b[k];
		}
	}

	/* 後退代入 */
	x[k] = b[k] / a[k][k];

	for (k = N - 2; k >= 0; k--) {
		for (int j = 0; j <= N - 1; j++) {
			x[k] -= a[k][j] * x[j];
		}
		x[k] += b[k];
	}

}

/* 2次元配列の表示 */
void print_matrix(double mat[N][N], int count) {
	cout << "count : " << setw(3) << right << count + 1 << endl;

	cout << "    ";
	for (int i = 0; i < N; i++) {
		cout << setw(14) << right << i << "\t";
	}
	cout << endl;

	for (int i = 0; i < N; i++) {
		cout << setw(3) << right << i << "|";
		for (int j = 0; j < N; j++) {
			cout << setw(14) << right << mat[i][j] << "\t";
		}
		cout << endl;
	}
	cout << endl;
}

/* 配列の表示 */
void print_list(double list[N]) {
	cout << endl;
	for (int i = 0; i < N; i++) {
		cout << setw(14) << right << list[i] << "\t";
	}
	cout << endl;
}

/* ノード間の流量を表示 */
void print_quantity(Route *route) {
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			if (DISTANCE_MAP[i][j] <= 0) continue;
			cout << "Q[" << i << "][" << j << "] ";
			cout << fixed << setprecision(15) << route->quantity[i][j] << endl;
		}
	}
}

void find_shortest(Route *route);

/* 最短経路を粘菌アルゴリズムによって計算 */
void physarum_solver(Route *route) {

	/// <image url="$(SolutionDir)\Images\pi_evaluate.png" scale = "0.6"/>

	double a[N][N] = {}; //ノードiにおけるp[j]の係数

	for (int n = 0; n < LOOP_NUMBER; n++) {

		/* ノードiにおけるp[j]の係数を定義にしたがって計算 */
		/// <image url="$(SolutionDir)\Images\pj_coefficient.png" scale = "0.4"/>
		for (int i = 0; i < N; i++) {
			a[i][i] = 0;
			for (int j = 0; j < N; j++) {
				double distance = DISTANCE_MAP[i][j];
				if (distance <= 0) continue;
				double ql = CONDUCTIVITY_MAP[i][j] / distance;
				a[i][i] += ql;
				a[i][j] = -ql;
			}
		}

		/// <image url="$(SolutionDir)\Images\pi_init.png" scale = "0.4"/>
		double p[N] = {}; // 各ノードの圧力
		/* 出発地点のノードの圧力を1.0に初期化 */
		p[route->start] = 1.0;

		/// <image url="$(SolutionDir)\Images\p_matrix.png" scale = "0.4"/>
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				double distance = DISTANCE_MAP[i][j];
				if (distance <= 0) continue;
				a[i][j] += CONDUCTIVITY_MAP[i][j] / distance * p[i];
			}
		}

		/* 圧力を圧力を求める行列の商を定義 */
		// 最初に流れている流量の値
		Quantity beginig_quantity[N] = {};
		beginig_quantity[route->start] = -I;
		beginig_quantity[route->destination] = I;

		/* 行列の行と列の入れ替え */
		double trans[N][N] = {};
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				trans[j][i] = a[i][j];
			}
		}

		/*ガウスの消去法によって、p（圧力）を求める*/
		gaus(trans, p, beginig_quantity);

		/*ノード間の流量を算出し、それを元に伝導率を更新*/
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				double distance = DISTANCE_MAP[i][j];
				if (distance <= 0) continue;
				double *quantity_address = &route->quantity[i][j]; // quantity[i][j]のアドレスを代入
				double *conductivity = &CONDUCTIVITY_MAP[i][j]; // conductivityとCONDUCTIVITY_MAP[i][j]のアドレスを同じにする *conductivityとはアドレスに存在する数値
				*quantity_address = *conductivity / distance * (p[i] - p[j]); /// <image url="$(SolutionDir)\Images\quantity_evaluate.png" scale = "0.4"/>
				*conductivity += DELTA * TIME_INTERVAL * (ALPHA * abs(*quantity_address) - GAMMA * *conductivity);/// <image url="$(SolutionDir)\Images\conductivity_renew.png" scale = "0.4"/>
			}
		}

	}

	/*現在のルートを表示*/
	find_shortest(route);
}

//流量を元に最短経路をroute.pointに格納
void find_shortest(Route *route) {

	Point point = route->start; //次に通過するノード
	Point destination = route->destination; //ゴール地点
	Point now_point = point; // 現在創作している地点
	Point before_point = point; // 探索する一つ前の地点

	// スタート地点を経路に格納
	route->point.emplace_back(point);
	// 流量の大きさが大きい方を経路に格納
	while (point != destination) { // ゴール地点に達したら探索終了
		double max = 0.0; // 流量の最大値
		for (int i = 0; i < N; i++) { // リンクしているノードから流量の大きいノードを選択し次のノードを探索
			if (DISTANCE_MAP[now_point][i] <= 0) continue; // リンクしていないものは除外
			if (i == before_point) continue; // 1つ前にの地点を除外
			// 最大値よりも流量が大きいノード間があったときpointを更新する
			if (max < abs(route->quantity[now_point][i])) {
				max = abs(route->quantity[now_point][i]);
				point = i;
			}
		}

		route->distance += DISTANCE_MAP[now_point][point];

		route->point.emplace_back(point);

		before_point = now_point;
		now_point = point;
	}
	cout << endl;
}

void print_route(Route *route) {
	cout << route->start;
	vector<Point>::iterator itr = route->point.begin();
	itr++;
	while (itr != route->point.end()) {
		cout << "," << *itr;
		itr++;
	}
}

void print_shortest(Route *route) {
	cout << "出発地点から目的地までの最短経路" << endl;
	print_route(route);
	cout << "\n出発地点から目的地までの最短経路 : " << route->distance << endl;
}

Point input_point(String prompt) {
	Point point;

	cout << prompt;
	cin >> point;
	return point;
}

int main() {
	Route route;
	route.start = input_point("出発地の地点番号を入力してください ===> ");
	route.destination = input_point("目的地の地点番号を入力してください ===> ");

	auto start = system_clock::now();
	physarum_solver(&route);
	auto end = system_clock::now();
	auto dur = end - start;
	auto msec = duration_cast<chrono::microseconds>(dur).count();
	//print_matrix(route.quantity,0);

	print_shortest(&route);
	cout << "処理速度 : " << msec << " micro sec" << endl;
}
