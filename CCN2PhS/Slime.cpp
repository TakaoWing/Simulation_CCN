#include "pch.h"
#include "Slime.h"

//float **Slime::conductivity_map;
Slime::Slime(){
}

/* �ŒZ�o�H��S�ۃA���S���Y���ɂ���Čv�Z */
void Slime::physarum_solver() {

	/* �`�����̒�`�Ə����� 
	float conductivity_map[N][N] = {};
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			if (!Node::link[i][j]) continue; // �m�[�h���ڑ�����Ă��Ȃ��Ƃ��ȉ��̏������s��Ȃ�
			conductivity_map[i][j] = conductivity_map[j][i] = 1.0; //�@�ڑ����ꂽ�m�[�h�Ԃ̓`�������P�ɂ���.
		}
	}*/

	/// <image url="$(SolutionDir)\Images\pi_evaluate.png" scale = "0.6"/>

	float a[N][N] = {}; //�m�[�hi�ɂ�����p[j]�̌W��

	for (int n = 0; n < LOOP_NUMBER; n++) {

		/* �m�[�hi�ɂ�����p[j]�̌W�����`�ɂ��������Čv�Z */
		/// <image url="$(SolutionDir)\Images\pj_coefficient.png" scale = "0.4"/>
		for (int i = 0; i < N; i++) {
			a[i][i] = 0;
			for (int j = 0; j < N; j++) {
				if (!Node::link[i][j]) continue; // �m�[�h���ڑ�����Ă��Ȃ��Ƃ��ȉ��̏������s��Ȃ�
				float ql = conductivity_map[i][j]; // �m�[�h�̋����͒ʐM�ɂ����Ċ֌W�Ȃ��̂ŏ��1�̂���
				a[i][i] += ql;
				a[i][j] = -ql;
			}
		}

		/// <image url="$(SolutionDir)\Images\pi_init.png" scale = "0.4"/>
		float p[N] = {}; // �e�m�[�h�̈���
		/* �o���n�_�̃m�[�h�̈��͂�1.0�ɏ����� */
		p[route->start] = 1.0;

		/// <image url="$(SolutionDir)\Images\p_matrix.png" scale = "0.4"/>
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (!Node::link[i][j]) continue;
				a[i][j] += (float)conductivity_map[i][j] * p[i];
			}
		}

		/* ���͂����͂����߂�s��̏����` */
		// �ŏ��ɗ���Ă��闬�ʂ̒l
		Quantity beginig_quantity[N] = {};
		beginig_quantity[route->start] = -I;
		beginig_quantity[route->destination] = I;

		/* �s��̍s�Ɨ�̓���ւ� */
		float trans[N][N] = {};
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				trans[j][i] = a[i][j];
			}
		}

		/*�K�E�X�̏����@�ɂ���āAp�i���́j�����߂�*/
		gauss(trans, p, beginig_quantity);

		/*�m�[�h�Ԃ̗��ʂ��Z�o���A��������ɓ`�������X�V*/
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (!Node::link[i][j]) continue;
				Quantity *quantity_address = &Slime::quantity[i][j]; // quantity[i][j]�̃A�h���X����
				float *conductivity = &conductivity_map[i][j]; // conductivity��CONDUCTIVITY_MAP[i][j]�̃A�h���X�𓯂��ɂ��� *conductivity�Ƃ̓A�h���X�ɑ��݂��鐔�l
				*quantity_address = (float)*conductivity * (p[i] - p[j]); /// <image url="$(SolutionDir)\Images\quantity_evaluate.png" scale = "0.4"/>
				*conductivity += float( DELTA * TIME_INTERVAL * (ALPHA * abs(*quantity_address) - GAMMA * *conductivity));/// <image url="$(SolutionDir)\Images\conductivity_renew.png" scale = "0.4"/>
			}
		}

	}

	/*���݂̃��[�g��\��*/
	//Slime::find_shortest();
}

void Slime::init_conductivity_map(){
	/*Slime::conductivity_map = (float**)malloc(N * sizeof(float*));
	Slime::conductivity_map[0] = (float*)malloc(N*N * sizeof(float));
	for (int i = 0; i < N; i++) {
		Slime::conductivity_map[i] = Slime::conductivity_map[0] + i * N;
	}*/

	/* �`�����̒�`�Ə�����*/
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			if (!Node::link[i][j]) continue; // �m�[�h���ڑ�����Ă��Ȃ��Ƃ��ȉ��̏������s��Ȃ�
			Slime::conductivity_map[i][j] = Slime::conductivity_map[j][i] = 1.0; //�@�ڑ����ꂽ�m�[�h�Ԃ̓`�������P�ɂ���.
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
			if (!Node::link[i][j]) continue; // �m�[�h���ڑ�����Ă��Ȃ��Ƃ��ȉ��̏������s��Ȃ�
			if (fabs(Slime::conductivity_map[i][j]) > 1.0e-4)continue; // ����Ȃ��[���ɓ������ꍇ�ȉ��̏������s��
			Slime::conductivity_map[i][j] = Slime::conductivity_map[j][i] = 0.0; //�@�ڑ����ꂽ�m�[�h�Ԃ̓`������0�ɂ���.
		}
	}
}

//���ʂ����ɍŒZ�o�H��route.point�Ɋi�[
void Slime::find_shortest() {

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
			if (!Node::link[now_point][i]) continue; // �����N���Ă��Ȃ����̂͏��O
			if (i == before_point) continue; // 1�O�̒n�_�����O
			// �ő�l�������ʂ��傫���m�[�h�Ԃ��������Ƃ�point���X�V����
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

/* �K�E�X�̍폜�@�v�Z���[�`�� */
void Slime::gauss(float a[N][N], float x[N], float b[N]) {
	float copy, akk, aik;
	int k, max;

	/* �s�|�b�g�I�� */
	for (k = 0; k <= N - 2; k++) {
		max = k;

		for (int i = k + 1; i <= N - 1; i++) {
			if (fabs(a[i][k]) > fabs(a[max][k])) max = i;
		}

		/* ��O�����i�Ίp�����������Ȓl�ɂȂ����Ƃ��j */
		//if (fabs(a[max][k]) < 1.0e-7) return;

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


/* �m�[�h�Ԃ̗��ʂ�\�� */
void Slime::print_quantity() {
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			if (!Node::link[i][j]) continue;
			cout << "Q[" << i << "][" << j << "] ";
			cout << fixed << setprecision(15) << Slime::quantity[i][j] << endl;
		}
	}
}