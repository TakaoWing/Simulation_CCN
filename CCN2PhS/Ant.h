#pragma once

// �A���g�R���j�[�̒萔
#define RHO 0.6 //�t�F��������������
#define INIT_P 0.1 //�t�F�����������l
#define S 30.0 //�ړ��t�F�������㏸�l�萔s
#define CONTENT_PHEROMONE 30.0 //�R���e���c���ړ����Ɏc���t�F������
#define INTERVAL_ANTS_SEND 1 // �A���𑗂�o���Ԋu
#define ANTS_AT_A_TIME 10 // �Ԋu���Ƃɑ�����o�Ă���A���̐�
#define CONCENTRATION_OF_PHEROMONE 2 // �m�[�h�Ɏc���t�F�������̔Z�x
#define MINIMUM_PHEROMONE 0.1 // �m�[�h�̃t�F�������̍ŏ��l
#define INTERVAL_DECAY_PHEROMON 10 // �t�F�������𐊑ނ�����P�ʎ���

class Ant{
	bool did_reach[N]; //���ڑ�����m�[�h���ʉ߂����m�[�h���ǂ����@true:�ʉߍς� false:���ʉ�
	vector<int> link_nodes; // �ڑ��m�[�h���i�[
	vector<int> trace; //�H�����o�H���i�[
	int living_time; // ��������Ă���̐�������
	int next_node; // ���ɐڑ�����m�[�h
	int position_node; // ���݈ʒu����m�[�h
	double cost; // ���R�X�g
public:
	Ant();
	void setup(); // ������
	void time2live(); // �������}�������ǂ���
	void search_node(); // �ڑ��ł���m�[�h��T��
	void choice_node(); // ���ɐڑ�����m�[�h��I��
	void communication_node(); // �I�������m�[�h�ƒʐM���s��
	void update(); // �ŒZ�o�H�Ȃǂ��X�V
	void backforward(); //  Backforward�a�s��
	static Route *route;
	static Pheromon nodes_pheromon[N][N];
	bool alive_backword_ant; // Backword�a�������Ă��邩
	bool alive;
	static void init_nodes_pheromon();
	static void decay_pheromone();
	static void ants_send(Ant *ants);
	static void leave_pheromone(Point point);
};