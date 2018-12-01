// PhysarumSolver.cpp : ���̃t�@�C���ɂ� 'main' �֐����܂܂�Ă��܂��B�v���O�������s�̊J�n�ƏI���������ōs���܂��B
//

#include "pch.h"
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdbool>
#include <chrono>
#include <vector>

#define N (10) /* �n�_�̐� */
#define I (1) /* �o���n�_���痬���� */
#define LOOP_NUMBER (100) /* �t�B�]�����E�\���o�[���J��Ԃ��� */
#define DELTA (0.1) // �S�̂̏d�ݕt�� def:0.1
#define ALPHA (1.0) // ���ʂɑ΂���d�ݕt�� def:1.0
#define GAMMA (1.0) // �`�����ɑ΂���d�ݕt�� def:1.0
#define TIME_INTERVAL (1.0) // ���ԊԊu def:1.0

using namespace std;
using namespace std::chrono;

/* �n�_�ԋ��� */
/* ��1 �n�_a�ƒn�_b������̒n�_�̏ꍇ��0���i�[����Ă��� */
/* ��2 �n�_a�ƒn�_b���אڒn�_�łȂ��ꍇ��-1���i�[����Ă��� */
static int DISTANCE_MAP[N][N] = {
	/* a->b 0   1   2   3   4   5   6   7   8   9 */
	/*0*/ {+0, +8, -1, -1, +4, +6, -1, -1, -1, -1},
	/*1*/ {+8, +0, +9, -1, +3, -1, -1, -1, -1, -1},
	/*2*/ {-1, +9, +0, +4, +3, -1, +4, -1, -1, +7},
	/*3*/ {-1, -1, +4, +0, -1, -1, -1, -1, -1, +3},
	/*4*/ {+4, +3, +3, -1, +0, +3, +4, -1, -1, -1},
	/*5*/ {+6, -1, -1, -1, +3, +0, +5, +2, -1, -1},
	/*6*/ {-1, -1, +4, -1, +4, +5, +0, +9, +3, +8},
	/*7*/ {-1, -1, -1, -1, -1, +2, +9, +0, +4, -1},
	/*8*/ {-1, -1, -1, -1, -1, -1, +3, +4, +0, +2},
	/*9*/ {-1, -1, +7, +3, -1, -1, +8, -1, +2, +0},
};

/* �`���� */
static double CONDUCTIVITY_MAP[N][N] = {
	/* a->b 0   1   2   3   4   5   6   7   8   9 */
	/*0*/ {+0, +1, +0, +0, +1, +1, +0, +0, +0, +0},
	/*1*/ {+1, +0, +1, +0, +1, +0, +0, +0, +0, +0},
	/*2*/ {+0, +1, +0, +1, +1, +0, +1, +0, +0, +1},
	/*3*/ {+0, +0, +1, +0, +0, +0, +0, +0, +0, +1},
	/*4*/ {+1, +1, +1, +0, +0, +1, +1, +0, +0, +0},
	/*5*/ {+1, +0, +0, +0, +1, +0, +1, +1, +0, +0},
	/*6*/ {+0, +0, +1, +0, +1, +1, +0, +1, +1, +1},
	/*7*/ {+0, +0, +0, +0, +0, +1, +1, +0, +1, +0},
	/*8*/ {+0, +0, +0, +0, +0, +0, +1, +1, +0, +1},
	/*9*/ {+0, +0, +1, +1, +0, +0, +1, +0, +1, +0},
};

typedef const char* String; // String�^���`
typedef int Point;
typedef int Distance;
typedef double Quantity;

struct Route { //�ʐM�o�H
	Point start, destination; // start:�ʐM�� , destination:�ʐM��
	Distance distance = 0; // �ʐM����
	vector<Point> point;
	Quantity quantity[N][N] = {}; // �l�b�g���[�N�S�̗̂��ʂ̒l
};

/* �K�E�X�̍폜�@�v�Z���[�`�� */
void gaus(double a[N][N], double x[N], double b[N]) {

	for (int k = 0; k < N; k++) {
		for (int i = k + 1; i < N; i++) {
			double alpha = a[i][k] / a[k][k];
			for (int j = k + 1; j < N; j++) {
				a[i][j] -= alpha * a[k][j];
			}
			b[i] -= alpha * b[k];
		}
	}

	x[N - 1] = b[N - 1] / a[N - 1][N - 1];

	for (int k = N - 2; k >= 0; k--) {
		for (int j = k + 1; j < N; j++) {
			b[k] -= a[k][j] * x[j];
		}
		x[k] = b[k] / a[k][k];
	}
}

void gauss(double a[N][N], double x[N], double b[N]) {
	double copy, akk, aik;
	int k, max;

	/* �s�|�b�g�I�� */
	for (k = 0; k <= N - 2; k++) {
		max = k;

		for (int i = k + 1; i <= N - 1; i++) {
			if (fabs(a[i][k]) > fabs(a[max][k])) max = i;
		}

		/* ��O�����i�Ίp�����������Ȓl�ɂȂ����Ƃ��j */
		if (fabs(a[max][k]) < 1.0e-7) exit(0);

		/* �s�̓���ւ� */
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

		/* �O�i���� */
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

	/* ��ޑ�� */
	x[k] = b[k] / a[k][k];

	for (k = N - 2; k >= 0; k--) {
		for (int j = 0; j <= N - 1; j++) {
			x[k] -= a[k][j] * x[j];
		}
		x[k] += b[k];
	}

}

/* 2�����z��̕\�� */
void print_matrix(double mat[N][N], int count) {
	cout << "count : " << setw(3) << right << count + 1 << endl;

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

/* �z��̕\�� */
void print_list(double list[N]) {
	cout << endl;
	for (int i = 0; i < N; i++) {
		cout << setw(14) << right << list[i] << "\t";
	}
	cout << endl;
}

/* �m�[�h�Ԃ̗��ʂ�\�� */
void print_quantity(Route *route) {
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			if (DISTANCE_MAP[i][j] <= 0) continue;
			cout << "Q[" << i << "][" << j << "] ";
			cout << fixed << setprecision(15) << route->quantity[i][j] << endl;
		}
	}
}

void find_shortest(Route *route);

/* �ŒZ�o�H��S�ۃA���S���Y���ɂ���Čv�Z */
void physarum_solver(Route *route) {

	/// <image url="$(SolutionDir)\Images\pi_evaluate.png" scale = "0.6"/>

	double a[N][N] = {}; //�m�[�hi�ɂ�����p[j]�̌W��

	for (int n = 0; n < LOOP_NUMBER; n++) {

		/* �m�[�hi�ɂ�����p[j]�̌W�����`�ɂ��������Čv�Z */
		/// <image url="$(SolutionDir)\Images\pj_coefficient.png" scale = "0.4"/>
		for (int i = 0; i < N; i++) {
			a[i][i] = 0;
			for (int j = 0; j < N; j++) {
				double distance = DISTANCE_MAP[i][j];
				if (distance <= 0) continue;
				double ql = CONDUCTIVITY_MAP[i][j] / distance;
				a[i][i] += ql;
				a[i][j] = -ql;
			}
		}

		/// <image url="$(SolutionDir)\Images\pi_init.png" scale = "0.4"/>
		double p[N] = {}; // �e�m�[�h�̈���
		/* �o���n�_�̃m�[�h�̈��͂�1.0�ɏ����� */
		p[route->start] = 1.0;

		/// <image url="$(SolutionDir)\Images\p_matrix.png" scale = "0.4"/>
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				double distance = DISTANCE_MAP[i][j];
				if (distance <= 0) continue;
				a[i][j] += CONDUCTIVITY_MAP[i][j] / distance * p[i];
			}
		}

		/* ���͂����͂����߂�s��̏����` */
		// �ŏ��ɗ���Ă��闬�ʂ̒l
		Quantity beginig_quantity[N] = {};
		beginig_quantity[route->start] = -I;
		beginig_quantity[route->destination] = I;

		/* �s��̍s�Ɨ�̓���ւ� */
		double trans[N][N] = {};
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				trans[j][i] = a[i][j];
			}
		}

		/*�K�E�X�̏����@�ɂ���āAp�i���́j�����߂�*/
		gaus(trans, p, beginig_quantity);

		/*�m�[�h�Ԃ̗��ʂ��Z�o���A��������ɓ`�������X�V*/
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				double distance = DISTANCE_MAP[i][j];
				if (distance <= 0) continue;
				double *quantity_address = &route->quantity[i][j]; // quantity[i][j]�̃A�h���X����
				double *conductivity = &CONDUCTIVITY_MAP[i][j]; // conductivity��CONDUCTIVITY_MAP[i][j]�̃A�h���X�𓯂��ɂ��� *conductivity�Ƃ̓A�h���X�ɑ��݂��鐔�l
				*quantity_address = *conductivity / distance * (p[i] - p[j]); /// <image url="$(SolutionDir)\Images\quantity_evaluate.png" scale = "0.4"/>
				*conductivity += DELTA * TIME_INTERVAL * (ALPHA * abs(*quantity_address) - GAMMA * *conductivity);/// <image url="$(SolutionDir)\Images\conductivity_renew.png" scale = "0.4"/>
			}
		}

	}

	/*���݂̃��[�g��\��*/
	find_shortest(route);
}

//���ʂ����ɍŒZ�o�H��route.point�Ɋi�[
void find_shortest(Route *route) {

	Point point = route->start; //���ɒʉ߂���m�[�h
	Point destination = route->destination; //�S�[���n�_
	Point now_point = point; // ���ݑn�삵�Ă���n�_
	Point before_point = point; // �T�������O�̒n�_

	// �X�^�[�g�n�_���o�H�Ɋi�[
	route->point.emplace_back(point);
	// ���ʂ̑傫�����傫�������o�H�Ɋi�[
	while (point != destination) { // �S�[���n�_�ɒB������T���I��
		double max = 0.0; // ���ʂ̍ő�l
		for (int i = 0; i < N; i++) { // �����N���Ă���m�[�h���痬�ʂ̑傫���m�[�h��I�������̃m�[�h��T��
			if (DISTANCE_MAP[now_point][i] <= 0) continue; // �����N���Ă��Ȃ����̂͏��O
			if (i == before_point) continue; // 1�O�ɂ̒n�_�����O
			// �ő�l�������ʂ��傫���m�[�h�Ԃ��������Ƃ�point���X�V����
			if (max < abs(route->quantity[now_point][i])) {
				max = abs(route->quantity[now_point][i]);
				point = i;
			}
		}

		route->distance += DISTANCE_MAP[now_point][point];

		route->point.emplace_back(point);

		before_point = now_point;
		now_point = point;
	}
	cout << endl;
}

void print_route(Route *route) {
	cout << route->start;
	vector<Point>::iterator itr = route->point.begin();
	itr++;
	while (itr != route->point.end()) {
		cout << "," << *itr;
		itr++;
	}
}

void print_shortest(Route *route) {
	cout << "�o���n�_����ړI�n�܂ł̍ŒZ�o�H" << endl;
	print_route(route);
	cout << "\n�o���n�_����ړI�n�܂ł̍ŒZ�o�H : " << route->distance << endl;
}

Point input_point(String prompt) {
	Point point;

	cout << prompt;
	cin >> point;
	return point;
}

int main() {
	Route route;
	route.start = input_point("�o���n�̒n�_�ԍ�����͂��Ă������� ===> ");
	route.destination = input_point("�ړI�n�̒n�_�ԍ�����͂��Ă������� ===> ");

	auto start = system_clock::now();
	physarum_solver(&route);
	auto end = system_clock::now();
	auto dur = end - start;
	auto msec = duration_cast<chrono::microseconds>(dur).count();
	//print_matrix(route.quantity,0);

	print_shortest(&route);
	cout << "�������x : " << msec << " micro sec" << endl;
}
