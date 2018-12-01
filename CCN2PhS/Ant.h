#pragma once

// アントコロニーの定数
#define RHO 0.6 //フェロモン減衰率ρ
#define INIT_P 0.1 //フェロモン初期値
#define S 30.0 //移動フェロモン上昇値定数s
#define CONTENT_PHEROMONE 30.0 //コンテンツが移動時に残すフェロモン
#define INTERVAL_ANTS_SEND 1 // アリを送り出す間隔
#define ANTS_AT_A_TIME 10 // 間隔ごとに巣から出てくるアリの数
#define CONCENTRATION_OF_PHEROMONE 2 // ノードに残すフェロモンの濃度
#define MINIMUM_PHEROMONE 0.1 // ノードのフェロモンの最小値
#define INTERVAL_DECAY_PHEROMON 10 // フェロモンを衰退させる単位時間

class Ant{
	bool did_reach[N]; //次接続するノードが通過したノードがどうか　true:通過済み false:未通過
	vector<int> link_nodes; // 接続ノードを格納
	vector<int> trace; //辿った経路を格納
	int living_time; // 生成されてからの生存時間
	int next_node; // 次に接続するノード
	int position_node; // 現在位置するノード
	double cost; // 総コスト
public:
	Ant();
	void setup(); // 初期化
	void time2live(); // 寿命を迎えたかどうか
	void search_node(); // 接続できるノードを探す
	void choice_node(); // 次に接続するノードを選択
	void communication_node(); // 選択したノードと通信を行う
	void update(); // 最短経路などを更新
	void backforward(); //  Backforward蟻行動
	static Route *route;
	static Pheromon nodes_pheromon[N][N];
	bool alive_backword_ant; // Backword蟻が生きているか
	bool alive;
	static void init_nodes_pheromon();
	static void decay_pheromone();
	static void ants_send(Ant *ants);
	static void leave_pheromone(Point point);
};