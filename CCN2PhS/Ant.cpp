#include "pch.h"
#include "Ant.h"

Ant::Ant() {

}

void Ant::setup(){
	trace.clear();
	link_nodes.clear();
	trace.push_back(route->start);
	position_node = route->start;
	for (int i = 0; i < N; i++) {
		did_reach[i] = false;
	}
	did_reach[position_node] = true;
	living_time = 0;
	cost = 0;
	alive = true;
	alive_backword_ant = false;
}

void Ant::time2live() {
	if (!alive) return;
	if (living_time < TTL) return;
	alive = false;
	alive_backword_ant = false;
	trace.clear();
}

void Ant::search_node() {
	if (!alive) return;
	for (int i = 0; i < N; i++) {
		if (!Node::link[position_node][i])continue;
		if (did_reach[i]) continue;
		link_nodes.push_back(i);
	}
	if (!link_nodes.empty())return;
	alive = false;
	trace.clear();
	alive_backword_ant = false;
}

void Ant::choice_node(){
	if (!alive) return;
	vector<float> probability; // 接続されたノードに対してそれぞれに進む確率
	float sum_pheromon = 0; // 接続されたノードのフェロモンの総数
	int total_node = (int)link_nodes.size();

	if (total_node == 1) {
		next_node = link_nodes[0];
		cost++;
		link_nodes.clear();
		return;
	}

	for (int i = 0; i < total_node; i++) {
		sum_pheromon += nodes_pheromon[position_node][link_nodes[i]];
	}
	for (int i = 0; i < total_node; i++) {
		probability.push_back(nodes_pheromon[position_node][link_nodes[i]] / sum_pheromon);
	}

	vector<float> each_probabllity(total_node,0);
	each_probabllity.push_back(1.0);

	auto itr = next(each_probabllity.begin(),1);
	auto itr_finish = prev(each_probabllity.end(), 1);
	auto next_itr = probability.begin();
	while (itr != itr_finish) {
		*itr++ = (float)accumulate(probability.begin(), ++next_itr,0.0);
	}

	float randam_choice = (float)(rand() % 1001) / 1000.0f;
	for (int i = 0; i < total_node; i++) {
		if (randam_choice <= each_probabllity[i] || each_probabllity[i + 1] < randam_choice) continue;
		next_node = link_nodes[i];
		cost++;
		link_nodes.clear();
		return;
	}
}

void Ant::communication_node(){
	if (!alive) return;
	living_time++;
	position_node = next_node;
	trace.push_back(position_node);
	did_reach[position_node] = true;
}

void Ant::update(){
	if (!alive) return;
	if (position_node != route->destination) return;
	alive = false;
	alive_backword_ant = true;
	reverse(trace.begin(), trace.end());
	for (int i = 0; i < N; i++) {
		did_reach[N] = false;
	}
	did_reach[position_node] = true;
}

void Ant::backforward(){
	if (!alive_backword_ant) return;
	if (position_node != route->destination) return;
	float rise_value = (float)S / (float)pow((trace.size() - 1), CONCENTRATION_OF_PHEROMONE);
	while (true) {
		auto next_position_node = next(trace.begin(), 1);
		nodes_pheromon[position_node][*next_position_node] += rise_value;
		nodes_pheromon[*next_position_node][position_node] += rise_value;

		position_node = *next_position_node;
		trace.erase(trace.begin());
		
		if (position_node != route->start) continue;
		alive_backword_ant = false;
		trace.clear();
		break;
	}
}

void Ant::init_nodes_pheromon(){
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			if (!Node::link[i][j]) continue;
			nodes_pheromon[i][j] = nodes_pheromon[j][i] = (Pheromon)INIT_P;
		}
	}
}

void Ant::decay_pheromone(){
	Pheromon *pheromon;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (!Node::link[i][j])continue;
			pheromon = &nodes_pheromon[i][j];
			*pheromon *= (Pheromon)RHO;	//減衰率ρをかけてフェロモン減衰
			if (MINIMUM_PHEROMONE <= *pheromon) continue;
			*pheromon = (Pheromon)MINIMUM_PHEROMONE;
		}
	}
}

void Ant::ants_send(Ant * ants){
	Ant *ant;
	for (int i = 0; i < ANTS_AT_A_TIME; i++) {
		ant = &ants[i];
		ant->setup();
		int count = 0;
		while (true) {
			ant->time2live();
			ant->search_node();
			ant->choice_node();
			ant->communication_node();
			ant->update();
			if (!ant->alive) break;
		}
	}
	for (int i = 0; i < ANTS_AT_A_TIME; i++) {
		ants[i].backforward();
	}
}

void Ant::leave_pheromone(Point point){
	nodes_pheromon[point][route->distance] += CONTENT_PHEROMONE;
	nodes_pheromon[route->distance][point] += CONTENT_PHEROMONE;
}
