
#ifndef PCH_H
#define PCH_H

#define N 100 /* 地点の数 */
#define TTL 25 //TTL init:25
#define SIMULATION_TIME 10000 //シミュレーションする時間
#define COMMUNICATION_DISTANCE 1.5	//各ノードの通信可能距離
#define INTERVAL_PACKETS_SEND 2	//パケットを送る間隔 ※基本的に変えない
#define STUDY_TIME 500 //データファイルを送り始めるまでの学習時間
#define UNIT_TIME 50 // データを記録する単位時間
#define BEFOR_MOBILE_POSITION 8 // コンテンツの初期位置
#define AFTER_MOBILE_POSITION 36 // コンテンツの移動後の初期位置
#define MOVING_TIME 1000 // コンテンツ端末が移動を開始する時刻
#define SEEDS 10　//乱数の個数

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdbool>
#include <vector>
#include <queue>
#include <ctime>
#include <numeric>
#include <fstream>
#include <iostream>

using namespace std;

typedef const char* String; // String型を定義
typedef int Point;
typedef int Distance;
typedef float Quantity;
typedef float Pheromon;

struct Route { //通信経路
	Point start, destination,goal; // start:通信元 , destination:通信先 ,　goal:移動後の通信先
	Distance distance = 0; // 通信距離
	vector<Point> point;
	vector<Point> destination2goal_points;
};

// ノード
struct Node {
	static bool link[N][N]; // 各ノードの接続状態
	float x, y; // 位置座標
	bool have_mobile = false; // コンテンツを持ったモバイル端末がいるかどうか
	bool did_reach = false; // 到達済かどうか
	vector<Point> points;
	Point destination_node;
};

struct DataBace {
	int finish_packets = 0;
	int total_packets = 0;
	float probability_packets = 0;
	float *probability_packets_list = (float *)malloc(sizeof(float) * (int)(SIMULATION_TIME / UNIT_TIME));
	int *packets_hop_list = (int *)malloc(sizeof(int) * (int)(SIMULATION_TIME / INTERVAL_PACKETS_SEND));
};

#endif //PCH_H
