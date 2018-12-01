#pragma once

// �S�ۃA���S���Y���̒萔
//�@�o�H��񂾂��łȂ��A���ʂ̒l���A�O��̂��̂������l�����Ƃ��Ɣ�ׂĂǂ��Ȃ邩�P��ڂ̗��ʁA�P�O��ڂ̗��ʁA�P�O�O��ڂ̗���
#define I 1 /* �o���n�_���痬���� */
#define LOOP_NUMBER 2 /* �t�B�]�����E�\���o�[���J��Ԃ��� */
#define DELTA 0.1 // �S�̂̏d�ݕt�� def:0.1
#define ALPHA 1.0 // ���ʂɑ΂���d�ݕt�� def:1.0
#define GAMMA 1.0 // �`�����ɑ΂���d�ݕt�� def:1.0
#define TIME_INTERVAL 1.0 // ���ԊԊu def:1.0

class Slime{
private:
	static void gauss(float a[N][N], float x[N], Quantity b[N]);
	static void find_shortest();
	void print_quantity();
	static float conductivity_map[N][N];
	//static float **conductivity_map;
public:
	Slime();
	static Quantity quantity[N][N];// �l�b�g���[�N�S�̗̂��ʂ̒l
	static void physarum_solver();
	static Route *route;
	static void init_conductivity_map();
	static void decay_conductivity_map();
	static void leave_conductivity(Point point);
};

