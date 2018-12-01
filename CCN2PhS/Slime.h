#pragma once

// 粘菌アルゴリズムの定数
//　経路情報だけでなく、流量の値が、前回のものを初期値したときと比べてどうなるか１回目の流量、１０回目の流量、１００回目の流量
#define I 1 /* 出発地点から流れる量 */
#define LOOP_NUMBER 2 /* フィゾルム・ソルバーを繰り返す回数 */
#define DELTA 0.1 // 全体の重み付け def:0.1
#define ALPHA 1.0 // 流量に対する重み付け def:1.0
#define GAMMA 1.0 // 伝導率に対する重み付け def:1.0
#define TIME_INTERVAL 1.0 // 時間間隔 def:1.0

class Slime{
private:
	static void gauss(float a[N][N], float x[N], Quantity b[N]);
	static void find_shortest();
	void print_quantity();
	static float conductivity_map[N][N];
	//static float **conductivity_map;
public:
	Slime();
	static Quantity quantity[N][N];// ネットワーク全体の流量の値
	static void physarum_solver();
	static Route *route;
	static void init_conductivity_map();
	static void decay_conductivity_map();
	static void leave_conductivity(Point point);
};

