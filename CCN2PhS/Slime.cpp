#include "pch.h"
#include "Slime.h"

//float **Slime::conductivity_map;
Slime::Slime(){
}

/* 最短経路を粘菌アルゴリズムによって計算 */
void Slime::physarum_solver() {

	/* 伝導率の定義と初期化 
	float conductivity_map[N][N] = {};
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			if (!Node::link[i][j]) continue; // ノードが接続されていないとき以下の処理を行わない
			conductivity_map[i][j] = conductivity_map[j][i] = 1.0; //　接続されたノード間の伝導率を１にする.
		}
	}*/

	/// <image url="$(SolutionDir)\Images\pi_evaluate.png" scale = "0.6"/>

	float a[N][N] = {}; //ノードiにおけるp[j]の係数

	for (int n = 0; n < LOOP_NUMBER; n++) {

		/* ノードiにおけるp[j]の係数を定義にしたがって計算 */
		/// <image url="$(SolutionDir)\Images\pj_coefficient.png" scale = "0.4"/>
		for (int i = 0; i < N; i++) {
			a[i][i] = 0;
			for (int j = 0; j < N; j++) {
				if (!Node::link[i][j]) continue; // ノードが接続されていないとき以下の処理を行わない
				float ql = conductivity_map[i][j]; // ノードの距離は通信において関係ないので常に1のため
				a[i][i] += ql;
				a[i][j] = -ql;
			}
		}

		/// <image url="$(SolutionDir)\Images\pi_init.png" scale = "0.4"/>
		float p[N] = {}; // 各ノードの圧力
		/* 出発地点のノードの圧力を1.0に初期化 */
		p[route->start] = 1.0;

		/// <image url="$(SolutionDir)\Images\p_matrix.png" scale = "0.4"/>
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (!Node::link[i][j]) continue;
				a[i][j] += (float)conductivity_map[i][j] * p[i];
			}
		}

		/* 圧力を圧力を求める行列の商を定義 */
		// 最初に流れている流量の値
		Quantity beginig_quantity[N] = {};
		beginig_quantity[route->start] = -I;
		beginig_quantity[route->destination] = I;

		/* 行列の行と列の入れ替え */
		float trans[N][N] = {};
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				trans[j][i] = a[i][j];
			}
		}

		/*ガウスの消去法によって、p（圧力）を求める*/
		gauss(trans, p, beginig_quantity);

		/*ノード間の流量を算出し、それを元に伝導率を更新*/
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (!Node::link[i][j]) continue;
				Quantity *quantity_address = &Slime::quantity[i][j]; // quantity[i][j]のアドレスを代入
				float *conductivity = &conductivity_map[i][j]; // conductivityとCONDUCTIVITY_MAP[i][j]のアドレスを同じにする *conductivityとはアドレスに存在する数値
				*quantity_address = (float)*conductivity * (p[i] - p[j]); /// <image url="$(SolutionDir)\Images\quantity_evaluate.png" scale = "0.4"/>
				*conductivity += float( DELTA * TIME_INTERVAL * (ALPHA * abs(*quantity_address) - GAMMA * *conductivity));/// <image url="$(SolutionDir)\Images\conductivity_renew.png" scale = "0.4"/>
			}
		}

	}

	/*現在のルートを表示*/
	//Slime::find_shortest();
}

void Slime::init_conductivity_map(){
	/*Slime::conductivity_map = (float**)malloc(N * sizeof(float*));
	Slime::conductivity_map[0] = (float*)malloc(N*N * sizeof(float));
	for (int i = 0; i < N; i++) {
		Slime::conductivity_map[i] = Slime::conductivity_map[0] + i * N;
	}*/

	/* 伝導率の定義と初期化*/
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			if (!Node::link[i][j]) continue; // ノードが接続されていないとき以下の処理を行わない
			Slime::conductivity_map[i][j] = Slime::conductivity_map[j][i] = 1.0; //　接続されたノード間の伝導率を１にする.
		}
	}
}

void Slime::leave_conductivity(Point next_node){
	quantity[route->destination][next_node] += 10.0;
	quantity[next_node][route->destination] += 10.0;
	conductivity_map[route->destination][next_node] = 1.0e-25f;
	conductivity_map[next_node][route->destination] = 1.0e-25f;
}

void Slime::decay_conductivity_map() {
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			if (!Node::link[i][j]) continue; // ノードが接続されていないとき以下の処理を行わない
			if (fabs(Slime::conductivity_map[i][j]) > 1.0e-4)continue; // 限りなくゼロに等しい場合以下の処理を行う
			Slime::conductivity_map[i][j] = Slime::conductivity_map[j][i] = 0.0; //　接続されたノード間の伝導率を0にする.
		}
	}
}

//流量を元に最短経路をroute.pointに格納
void Slime::find_shortest() {

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
			if (!Node::link[now_point][i]) continue; // リンクしていないものは除外
			if (i == before_point) continue; // 1つ前の地点を除外
			// 最大値よりも流量が大きいノード間があったときpointを更新する
			Quantity _quantity = abs(Slime::quantity[now_point][i]);
			if (max >= _quantity) continue;
			max = _quantity;
			point = i;
		}

		route->distance += Node::link[now_point][point] ? 1 : 0;
		route->point.emplace_back(point);

		before_point = now_point;
		now_point = point;
	}
}

/* ガウスの削除法計算ルーチン */
void Slime::gauss(float a[N][N], float x[N], float b[N]) {
	float copy, akk, aik;
	int k, max;

	/* ピポット選択 */
	for (k = 0; k <= N - 2; k++) {
		max = k;

		for (int i = k + 1; i <= N - 1; i++) {
			if (fabs(a[i][k]) > fabs(a[max][k])) max = i;
		}

		/* 例外処理（対角成分が小さな値になったとき） */
		//if (fabs(a[max][k]) < 1.0e-7) return;

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


/* ノード間の流量を表示 */
void Slime::print_quantity() {
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			if (!Node::link[i][j]) continue;
			cout << "Q[" << i << "][" << j << "] ";
			cout << fixed << setprecision(15) << Slime::quantity[i][j] << endl;
		}
	}
}