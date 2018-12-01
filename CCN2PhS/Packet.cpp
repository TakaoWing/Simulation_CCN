#include "pch.h"
#include "Packet.h"
#include "Ant.h"
#include "Slime.h"

Packet::Packet(){
}

void Packet::setup(){
	trace.clear();
	link_nodes.clear();
	for (int i = 0; i < N; i++) {
		did_reach[i] = false;
	}
	trace.push_back(route->start);
	position_node = route->start;
	did_reach[position_node] = true;
	living_time = 0;
	hop = 0;
	alive = true;
}

void Packet::time2live(){
	if (!alive) return;
	if (living_time < TTL) return;
	alive = false;
	trace.clear();
}

void Packet::search_node(){
	if (!alive) return;
	for (int i = 0; i < N; i++) {
		if (!Node::link[position_node][i])continue;
		if (did_reach[i]) continue;
		link_nodes.push_back(i);
	}
	if (!link_nodes.empty())return;
	alive = false;
	trace.clear();
}

void Packet::choice_node_aco(){
	if (!alive) return;
	vector<float> probability; // 接続されたノードに対してそれぞれに進む確率
	float sum_pheromon = 0; // 接続されたノードのフェロモンの総数
	int total_node = (int)link_nodes.size();

	if (total_node == 1) {
		next_node = link_nodes[0];
		link_nodes.clear();
		return;
	}
	
	for (int i = 0; i < total_node; i++) {
		sum_pheromon += Ant::nodes_pheromon[position_node][link_nodes[i]];
	}
	for (int i = 0; i < total_node; i++) {
		probability.push_back(Ant::nodes_pheromon[position_node][link_nodes[i]] / sum_pheromon);
	}

	vector<float> each_probabllity(total_node, 0);
	each_probabllity.push_back(1.0);

	auto itr = next(each_probabllity.begin(), 1);
	auto itr_finish = prev(each_probabllity.end(), 1);
	auto next_itr = probability.begin();
	while (itr != itr_finish) {
		*itr++ = (float)accumulate(probability.begin(), ++next_itr, 0.0);
	}

	float randam_choice = (float)(rand() % 1001) / 1000.0f;
	for (int i = 0; i < total_node; i++) {
		if (randam_choice <= each_probabllity[i] || each_probabllity[i + 1] < randam_choice) continue;
		next_node = link_nodes[i];
		link_nodes.clear();
		return;
	}
}

void Packet::choice_node_phs() {
	if (!alive)return;
	int total_node = (int)link_nodes.size();
	if (total_node == 1) {
		next_node = link_nodes[0];
		link_nodes.clear();
		return;
	}
	Quantity max_quantity = abs(Slime::quantity[position_node][link_nodes[0]]);
	Quantity _quantity;
	int max_index = 0;
	for (int i = 1; i < total_node; i++) {
		_quantity = abs(Slime::quantity[position_node][link_nodes[i]]);
		if (max_quantity >= _quantity) continue;
		max_quantity = _quantity;
		max_index = i;
	}
	next_node = link_nodes[max_index];
	link_nodes.clear();
}

void Packet::communication_node(){
	if (!alive) return;
	living_time++;
	hop++;
	position_node = next_node;
	trace.push_back(position_node);
	did_reach[position_node] = true;
}

void Packet::update(){
	if (!alive) return;
	if (position_node != route->destination) return;
	total_finish_packets++;
	total_delay += living_time;
	total_hop += hop;;
	alive = false;
}

void Packet::packets_send_aco() {
	Packet packet;
	packet.setup();
	while (true) {
		packet.time2live();
		packet.search_node();
		packet.choice_node_aco();
		packet.communication_node();
		packet.update();
		if (!packet.alive) break;
	}
}

void Packet::packets_send_phs() {
	Packet packet;
	packet.setup();
	while (true) {
		packet.time2live();
		packet.search_node();
		packet.choice_node_phs();
		packet.communication_node();
		packet.update();
		if (!packet.alive) break;
	}
}

void Packet::unit_data(DataBace *databace) {
	databace->finish_packets = total_finish_packets;
	databace->probability_packets = (float)databace->finish_packets / (UNIT_TIME*0.5f);
	total_finish_packets = 0;
}