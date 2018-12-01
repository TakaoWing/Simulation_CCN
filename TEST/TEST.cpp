#include "pch.h"
#include <iostream>

#define N 10000
#define COMMUNICATION_DISTANCE 1.5	//各ノードの通信可能距離

// ノード
struct Node {
	static int link[N][N]; // 各ノードの接続状態
	double x, y; // 位置座標
	bool have_mobile = false; // コンテンツを持ったモバイル端末がいるかどうか
};
int Node::link[N][N] = {};

void init_nodes(Node *nodes) {
	double root_N = sqrt(N);
	int ceil_root_N = int(ceil(root_N));
	int floor_root_N = int(floor(root_N));
	Node *n = nodes;
	for (int i = 0; i < floor_root_N; i++) {
		int position_y = floor_root_N * i;
		for (int j = 0; j < floor_root_N; j++) {
			n->x = j;
			n->y = i;
			*n++;
		}
	}
	if (floor_root_N != ceil_root_N) {
		int square_floor_root_N = floor_root_N * floor_root_N;
		for (int i = 0; i < N - square_floor_root_N;i++) {
			n->x = i;
			n->y = ceil_root_N - 1;
			*n++;
		}
	}

	Node *node_i = nodes, *node_j;
	for (int i = 0; i < N; i++, *node_i++) {
		node_j = nodes;
		for (int j = 0; j < N; j++, *node_j++) {
			double difference_x = node_i->x - node_j->x;
			double difference_y = node_i->y - node_j->y;
			double distance = sqrt(difference_x * difference_x + difference_y * difference_y);
			if (distance <= 0)continue;
			if (distance > COMMUNICATION_DISTANCE)continue;
			Node::link[i][j] = 1;
		}
	}
}

int main()
{
	Node nodes[N];
	init_nodes(nodes);
    std::cout << "Hello World!\n"; 
}
