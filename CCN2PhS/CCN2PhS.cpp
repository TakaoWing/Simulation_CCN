// CCN2PhS.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include "pch.h"
#include "Ant.h"
#include "Slime.h"
#include "Packet.h"
#include <limits>

bool Node::link[N][N] = {};
Quantity Slime::quantity[N][N] = {};
float Slime::conductivity_map[N][N] = {};

Route *Ant::route;
Route *Packet::route;
Route *Slime::route;
Pheromon Ant::nodes_pheromon[N][N] = {};

vector<int> Packet::trace;
int Packet::total_finish_packets = 0;
int Packet::total_delay = 0;
int Packet::total_hop = 0;

/* 2次元配列の表示 */
void print_matrix(float mat[N][N]) {
	
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
void print_matrix(int mat[N][N]) {

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

void print_route(vector<int> points) {
	cout << points[0];
	auto itr = points.begin();
	++itr;
	while (itr != points.end()) {
		cout << "," << *itr++;
	}
	cout << endl;
}

void print_shortest(Route *route) {
	cout << "出発地点から目的地までの最短経路" << endl;
	print_route(route->point);
	cout << "\n出発地点から目的地までの最短経路 : " << route->distance << endl;
}

Point input_point(String prompt) {
	Point point;

	cout << prompt;
	cin >> point;
	return point;
}

void init_nodes(Node *nodes) {
	float root_N = (float)sqrt(N);
	int ceil_root_N = int(ceil(root_N));
	int floor_root_N = int(floor(root_N));
	Node *n = nodes;
	for (float i = 0; i < floor_root_N; i++) {
		for (float j = 0; j < floor_root_N; j++) {
			n->x = j;
			n->y = i;
			*n++;
		}
	}
	if (floor_root_N != ceil_root_N) {
		int square_floor_root_N = floor_root_N * floor_root_N;
		for (float i = 0; i < N - square_floor_root_N;i++) {
			n->x = i;
			n->y = (float)ceil_root_N - 1;
			*n++;
		}
	}

	Node *node_i = nodes, *node_j;
	for (int i = 0; i < N; i++, *node_i++) {
		node_j = nodes;
		for (int j = 0; j < N; j++, *node_j++) {
			float difference_x = node_i->x - node_j->x;
			float difference_y = node_i->y - node_j->y;
			float distance = (float)sqrt(difference_x * difference_x + difference_y * difference_y);
			if (distance <= 0)continue;
			if (distance > COMMUNICATION_DISTANCE)continue;
			Node::link[i][j] = true;
		}
	}
}

void flooding(Route *route, Node *nodes) {
	queue<int> q; // 動作するノードの順番
	Point goal;
	Point start = route->start;
	Node *node = &nodes[start];
	q.push(start);
	node->did_reach = 1;
	node->points.push_back(start);
	do {
		int i = q.front();
		q.pop();
		for (int j = 0; j < N; j++) {
			if (!Node::link[i][j]) continue;
			node = &nodes[j];
			if (node->did_reach) continue;
			node->did_reach = 1;
			node->points = nodes[i].points;
			node->points.push_back(j);
			if (!node->have_mobile) {
				q.push(j);
			}
			else {
				node->destination_node = j;
				goal = j;
				goto LAVEL10;
			}
		}
	} while (!q.empty());
LAVEL10:
	auto itr = nodes[goal].points.end();
	--itr;
	do {
		--itr;
		nodes[*itr].destination_node = goal;
	} while (itr != nodes[goal].points.begin());
	route->destination = nodes[*itr].destination_node;
	for (int i = 0; i < N; i++) {
		node = &nodes[i];
		if (!node->did_reach) continue;
		node->did_reach = false;
		node->points.clear();
	}
}

void route_serch(Route *route, Node *nodes) {
	queue<int> q; // 動作するノードの順番
	Point start = route->destination;
	Node *node = &nodes[start];
	q.push(start);
	node->did_reach = 1;
	node->points.push_back(start);
	do {
		int i = q.front();
		q.pop();
		for (int j = 0; j < N; j++) {
			if (!Node::link[i][j]) continue;

			node = &nodes[j];
			if (node->did_reach) continue;
			node->did_reach = true;
			node->points = nodes[i].points;
			node->points.push_back(j);

			if (j != route->goal) q.push(j);
			else goto LAVEL10;
		}
	} while (!q.empty());
LAVEL10:
	route->destination2goal_points = nodes[route->goal].points;
	for (int i = 0; i < N; i++) {
		node = &nodes[i];
		if (!node->did_reach) continue;
		node->did_reach = false;
		node->points.clear();
	}
}

void ant_simulation(Route *route,DataBace *databace) {
	auto next_node = next(route->destination2goal_points.begin(), 1);
	Ant::route = route;
	Ant::init_nodes_pheromon();
	Ant ants[ANTS_AT_A_TIME];
	srand((unsigned)time(NULL));
	srand(1);
	for (int t = 1; t <= SIMULATION_TIME; t++) {

		// 一定時間ごとにフェロモンを衰退
		if (t % INTERVAL_DECAY_PHEROMON == 0) Ant::decay_pheromone();

		// アリ学習
		Ant::ants_send(ants);

		// モバイル端末の移動
		if (t == MOVING_TIME) {
			while (true) {
				if (next_node == route->destination2goal_points.end())break;
				Ant::put_pheromone(*next_node);
				route->destination = *next_node++;
			}
		}

		// パケットの送信
		if (t%INTERVAL_PACKETS_SEND == 0) {
			databace->total_packets++;
			Packet::packets_send_aco();
			databace->packets_hop_list[(int)(t / INTERVAL_PACKETS_SEND) - 1] = Packet::total_hop;
			Packet::total_hop = 0;
		}
		route->point.clear();

		// 単位時間ごとの情報を出力
		if (t % UNIT_TIME != 0) continue;
		cout << "time = " << t << endl;
		Packet::unit_data(databace);
		databace->probability_packets_list[(int)(t / UNIT_TIME) - 1] = databace->probability_packets;
		cout << databace->probability_packets * 100 << endl;
	}
}

void phs_simulation(Route *route, Node *nodes,DataBace *databace) {
	auto next_node = next(route->destination2goal_points.begin(), 1);
	Slime::route = route;
	Slime::init_conductivity_map();
	for (int t = 1; t <= SIMULATION_TIME; t++) {

		cout << "time = " << t << endl;

		// 粘菌学習
		//Slime::init_conductivity_map();
		flooding(route, nodes);
		Slime::physarum_solver();

		// モバイル端末の移動
		if (t == MOVING_TIME) {
			while (true) {
				if (next_node == route->destination2goal_points.end())break; // 目的地がゴールになったら移動完了
				Slime::put_conductivity_flowRate(*next_node); // 流量と伝導率の値を移動前から移動後に付加する．
				// 目的地の更新
				Point *destination = &route->destination;
				nodes[*destination].have_mobile = false;
				*destination = *next_node++;
				nodes[*destination].have_mobile = true;
			}
		}

		// パケットの送信
		if (t%INTERVAL_PACKETS_SEND == 0) {
			databace->total_packets++;
			Packet::packets_send_phs();
			databace->packets_hop_list[(int)(t / INTERVAL_PACKETS_SEND) - 1] = Packet::total_hop;
			cout << "Packet hop : " << Packet::total_hop << endl;
			print_route(Packet::trace);
			Packet::total_hop = 0;
		}
		route->point.clear();

		// 単位時間ごとの情報を出力
		if (t % UNIT_TIME != 0) continue;
		Packet::unit_data(databace);
		databace->probability_packets_list[(int)(t/UNIT_TIME)-1] = databace->probability_packets;
		cout << "Packet reception rate : " << databace->probability_packets * 100 << endl;

	}
}

int main() {
	DataBace databace;

	Node nodes[N];
	init_nodes(nodes);

	Route route;
	route.start = 0;// input_point("出発地の地点番号を入力してください ===> ");
	route.destination = 33;//nput_point("目的地の地点番号を入力してください ===> ");
	route.goal = 9;//input_point("目的地の到着地点を入力してください ===> ");
	route_serch(&route, nodes);
	nodes[route.destination].have_mobile = true;

	Packet::route = &route;
	Packet packet;

	// 蟻コロニー　シミュレーション
	//ant_simulation(&route,&databace);

	// 粘菌アルゴリズム　シミュレーション
	phs_simulation(&route, nodes,&databace);

	/*ofstream ofs("PhS-probablity_packets.csv");
	if (!ofs) {
		cerr << "ファイルオープンに失敗" << endl;
	}

	for (int i = 1; i <= (int)(SIMULATION_TIME / UNIT_TIME); i++) {
		ofs << i*50 << ",";
	}
	ofs << endl;
	for (int i = 0; i < (int)(SIMULATION_TIME / UNIT_TIME); i++) {
		ofs << databace.probability_packets_list[i] << ",";
	}

	ofs << endl;*/

	/*ofstream ofs("PhS-packets_hop.csv");
	if (!ofs) {
		cerr << "ファイルオープンに失敗" << endl;
	}
	for (int i = 1; i <= (int)(SIMULATION_TIME / INTERVAL_PACKETS_SEND);i++) {
		ofs << i* INTERVAL_PACKETS_SEND << ",";
	}
	ofs << endl;
	for (int i = 0; i < (int)(SIMULATION_TIME / INTERVAL_PACKETS_SEND); i++) {
		ofs << databace.packets_hop_list[i] << ",";
	}
	ofs << endl;*/

	free(databace.probability_packets_list);
	free(databace.packets_hop_list);
}

