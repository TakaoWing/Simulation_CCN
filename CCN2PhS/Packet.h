#pragma once
class Packet{
	bool did_reach[N];
	vector<int> link_nodes; // 接続ノードを格納
	int position_node;
	int living_time;
	int next_node;
	int hop;
public:
	Packet();
	void setup(); // 初期化
	void time2live(); // 寿命を迎えたかどうか
	void search_node(); // 接続できるノードを探す
	void choice_node_aco(); // 次に接続するノードを選択
	void choice_node_phs();
	void communication_node(); // 選択したノードと通信を行う
	void update(); // 最短経路などを更新
	static Route *route;
	static int total_finish_packets;
	static int total_delay;
	static int total_hop;
	static void unit_data(DataBace *databace);
	static void packets_send_aco();
	static void packets_send_phs();
	static vector<int> trace;;
	bool alive;
};