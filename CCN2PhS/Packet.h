#pragma once
class Packet{
	bool did_reach[N];
	vector<int> link_nodes; // �ڑ��m�[�h���i�[
	int position_node;
	int living_time;
	int next_node;
	int hop;
public:
	Packet();
	void setup(); // ������
	void time2live(); // �������}�������ǂ���
	void search_node(); // �ڑ��ł���m�[�h��T��
	void choice_node_aco(); // ���ɐڑ�����m�[�h��I��
	void choice_node_phs();
	void communication_node(); // �I�������m�[�h�ƒʐM���s��
	void update(); // �ŒZ�o�H�Ȃǂ��X�V
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