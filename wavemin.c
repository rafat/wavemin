/*
Copyright (c) 2017, Rafat Hussain
*/

#include "wavemin.h"

wave_object wave_init(char* wname) {
	wave_object obj = NULL;
	int retval;
	retval = 0;

	if (wname != NULL) {
		retval = filtlength(wname);
		//obj->filtlength = retval;
		//strcopy(obj->wname, wname);
	}

	obj = (wave_object)malloc(sizeof(struct wave_set) + sizeof(double)* 4 * retval);

	obj->filtlength = retval;
	obj->lpd_len = obj->hpd_len = obj->lpr_len = obj->hpr_len = obj->filtlength;
	strcpy(obj->wname, wname);
	if (wname != NULL) {
		filtcoef(wname, obj->params, obj->params + retval, obj->params + 2 * retval, obj->params + 3 * retval);
	}
	obj->lpd = &obj->params[0];
	obj->hpd = &obj->params[retval];
	obj->lpr = &obj->params[2 * retval];
	obj->hpr = &obj->params[3 * retval];
	return obj;
}

wt_object wt_init(wave_object wave, char* method, int siglength, int J) {
	int size, i, MaxIter;
	wt_object obj = NULL;

	size = wave->filtlength;

	if (J > 100) {
		printf("\n The Decomposition Iterations Cannot Exceed 100. Exiting \n");
		exit(-1);
	}

	MaxIter = wmaxiter(siglength, size);

	if (J > MaxIter) {
		printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", MaxIter);
		exit(-1);
	}

	if (method == NULL) {
		obj = (wt_object)malloc(sizeof(struct wt_set) + sizeof(double)* (siglength + 2 * J * (size + 1)));
		obj->outlength = siglength + 2 * J * (size + 1); // Default
		strcpy(obj->ext, "sym"); // Default
	}
	else if (!strcmp(method, "dwt") || !strcmp(method, "DWT")) {
		obj = (wt_object)malloc(sizeof(struct wt_set) + sizeof(double)* (siglength + 2 * J * (size + 1)));
		obj->outlength = siglength + 2 * J * (size + 1); // Default
		strcpy(obj->ext, "sym"); // Default
	}

	obj->wave = wave;
	obj->siglength = siglength;
	obj->J = J;
	obj->MaxIter = MaxIter;
	strcpy(obj->method, method);

	if (siglength % 2 == 0) {
		obj->even = 1;
	}
	else {
		obj->even = 0;
	}


	strcpy(obj->cmethod, "direct"); // Default
	obj->cfftset = 0;
	obj->lenlength = J + 2;
	obj->output = &obj->params[0];
	if (!strcmp(method, "dwt") || !strcmp(method, "DWT")) {
		for (i = 0; i < siglength + 2 * J * (size + 1); ++i) {
			obj->params[i] = 0.0;
		}
	}
	//wave_summary(obj->wave);

	return obj;
}


static void wconv(wt_object wt, double *sig, int N, double *filt, int L, double *oup) {
	if (!strcmp(wt->cmethod, "direct")) {
		conv_direct(sig, N, filt, L, oup);
	}	
	else {
		printf("Convolution Only accepts direct convolution");
		exit(-1);
	}
}


static void dwt_per(wt_object wt, double *inp, int N, double *cA, int len_cA, double *cD, int len_cD) {
	int l, l2, isodd, i, t, len_avg;

	len_avg = wt->wave->lpd_len;
	l2 = len_avg / 2;
	isodd = N % 2;

	for (i = 0; i < len_cA; ++i) {
		t = 2 * i + l2;
		cA[i] = 0.0;
		cD[i] = 0.0;
		for (l = 0; l < len_avg; ++l) {
			if ((t - l) >= l2 && (t - l) < N) {
				cA[i] += wt->wave->lpd[l] * inp[t - l];
				cD[i] += wt->wave->hpd[l] * inp[t - l];
			}
			else if ((t - l) < l2 && (t - l) >= 0) {
				cA[i] += wt->wave->lpd[l] * inp[t - l];
				cD[i] += wt->wave->hpd[l] * inp[t - l];
			}
			else if ((t - l) < 0 && isodd == 0) {
				cA[i] += wt->wave->lpd[l] * inp[t - l + N];
				cD[i] += wt->wave->hpd[l] * inp[t - l + N];
			}
			else if ((t - l) < 0 && isodd == 1) {
				if ((t - l) != -1) {
					cA[i] += wt->wave->lpd[l] * inp[t - l + N + 1];
					cD[i] += wt->wave->hpd[l] * inp[t - l + N + 1];
				}
				else {
					cA[i] += wt->wave->lpd[l] * inp[N - 1];
					cD[i] += wt->wave->hpd[l] * inp[N - 1];
				}
			}
			else if ((t - l) >= N && isodd == 0) {
				cA[i] += wt->wave->lpd[l] * inp[t - l - N];
				cD[i] += wt->wave->hpd[l] * inp[t - l - N];
			}
			else if ((t - l) >= N && isodd == 1) {
				if (t - l != N) {
					cA[i] += wt->wave->lpd[l] * inp[t - l - (N + 1)];
					cD[i] += wt->wave->hpd[l] * inp[t - l - (N + 1)];
				}
				else {
					cA[i] += wt->wave->lpd[l] * inp[N - 1];
					cD[i] += wt->wave->hpd[l] * inp[N - 1];
				}
			}

		}
	}



}


static void dwt_sym(wt_object wt, double *inp, int N, double *cA, int len_cA, double *cD, int len_cD) {
	int i, l, t, len_avg;

	len_avg = wt->wave->lpd_len;

	for (i = 0; i < len_cA; ++i) {
		t = 2 * i + 1;
		cA[i] = 0.0;
		cD[i] = 0.0;
		for (l = 0; l < len_avg; ++l) {
			if ((t - l) >= 0 && (t - l) < N) {
				cA[i] += wt->wave->lpd[l] * inp[t - l];
				cD[i] += wt->wave->hpd[l] * inp[t - l];
			}
			else if ((t - l) < 0) {
				cA[i] += wt->wave->lpd[l] * inp[-t + l - 1];
				cD[i] += wt->wave->hpd[l] * inp[-t + l - 1];
			}
			else if ((t - l) >= N) {
				cA[i] += wt->wave->lpd[l] * inp[2 * N - t + l - 1];
				cD[i] += wt->wave->hpd[l] * inp[2 * N - t + l - 1];
			}

		}
	}


}


void dwt(wt_object wt, double *inp) {
	int i, J, temp_len, iter, N, lp;
	int len_cA;
	double *orig, *orig2;

	temp_len = wt->siglength;
	J = wt->J;
	wt->length[J + 1] = temp_len;
	wt->outlength = 0;
	wt->zpad = 0;
	orig = (double*)malloc(sizeof(double)* temp_len);
	orig2 = (double*)malloc(sizeof(double)* temp_len);
	/*
	if ((temp_len % 2) == 0) {
	wt->zpad = 0;
	orig = (double*)malloc(sizeof(double)* temp_len);
	orig2 = (double*)malloc(sizeof(double)* temp_len);
	}
	else {
	wt->zpad = 1;
	temp_len++;
	orig = (double*)malloc(sizeof(double)* temp_len);
	orig2 = (double*)malloc(sizeof(double)* temp_len);
	}
	*/

	for (i = 0; i < wt->siglength; ++i) {
		orig[i] = inp[i];
	}

	if (wt->zpad == 1) {
		orig[temp_len - 1] = orig[temp_len - 2];
	}

	N = temp_len;
	lp = wt->wave->lpd_len;

	if (!strcmp(wt->ext, "per")) {
		i = J;
		while (i > 0) {
			N = (int)ceil((double)N / 2.0);
			wt->length[i] = N;
			wt->outlength += wt->length[i];
			i--;
		}
		wt->length[0] = wt->length[1];
		wt->outlength += wt->length[0];
		N = wt->outlength;

		for (iter = 0; iter < J; ++iter) {
			len_cA = wt->length[J - iter];
			N -= len_cA;
			
			dwt_per(wt, orig, temp_len, orig2, len_cA, wt->params + N, len_cA);
			
			temp_len = wt->length[J - iter];
			if (iter == J - 1) {
				for (i = 0; i < len_cA; ++i) {
					wt->params[i] = orig2[i];
				}
			}
			else {
				for (i = 0; i < len_cA; ++i) {
					orig[i] = orig2[i];
				}
			}
		}
	}
	else if (!strcmp(wt->ext, "sym")) {
		//printf("\n YES %s \n", wt->ext);
		i = J;
		while (i > 0) {
			N = N + lp - 2;
			N = (int)ceil((double)N / 2.0);
			wt->length[i] = N;
			wt->outlength += wt->length[i];
			i--;
		}
		wt->length[0] = wt->length[1];
		wt->outlength += wt->length[0];
		N = wt->outlength;

		for (iter = 0; iter < J; ++iter) {
			len_cA = wt->length[J - iter];
			N -= len_cA;
			dwt_sym(wt, orig, temp_len, orig2, len_cA, wt->params + N, len_cA);
			temp_len = wt->length[J - iter];

			if (iter == J - 1) {
				for (i = 0; i < len_cA; ++i) {
					wt->params[i] = orig2[i];
				}
			}
			else {
				for (i = 0; i < len_cA; ++i) {
					orig[i] = orig2[i];
				}
			}
		}
	}
	else {
		printf("Signal extension can be either per or sym");
		exit(-1);
	}

	free(orig);
	free(orig2);
}


static int ipow2(int n) {
	int p, i;
	p = 1;

	for (i = 0; i < n; ++i) {
		p *= 2;
	}

	return p;
}


static void idwt_per(wt_object wt, double *cA, int len_cA, double *cD, int len_cD, double *X) {
	int len_avg, i, l, m, n, t, l2;

	len_avg = (wt->wave->lpr_len + wt->wave->hpr_len) / 2;
	l2 = len_avg / 2;
	m = -2;
	n = -1;

	for (i = 0; i < len_cA + l2 - 1; ++i) {
		m += 2;
		n += 2;
		X[m] = 0.0;
		X[n] = 0.0;
		for (l = 0; l < l2; ++l) {
			t = 2 * l;
			if ((i - l) >= 0 && (i - l) < len_cA) {
				X[m] += wt->wave->lpr[t] * cA[i - l] + wt->wave->hpr[t] * cD[i - l];
				X[n] += wt->wave->lpr[t + 1] * cA[i - l] + wt->wave->hpr[t + 1] * cD[i - l];
			}
			else if ((i - l) >= len_cA && (i - l) < len_cA + len_avg - 1) {
				X[m] += wt->wave->lpr[t] * cA[i - l - len_cA] + wt->wave->hpr[t] * cD[i - l - len_cA];
				X[n] += wt->wave->lpr[t + 1] * cA[i - l - len_cA] + wt->wave->hpr[t + 1] * cD[i - l - len_cA];
			}
			else if ((i - l) < 0 && (i - l) > -l2) {
				X[m] += wt->wave->lpr[t] * cA[len_cA + i - l] + wt->wave->hpr[t] * cD[len_cA + i - l];
				X[n] += wt->wave->lpr[t + 1] * cA[len_cA + i - l] + wt->wave->hpr[t + 1] * cD[len_cA + i - l];
			}
		}
	}
}

static void idwt_sym(wt_object wt, double *cA, int len_cA, double *cD, int len_cD, double *X) {
	int len_avg, i, l, m, n, t, v;

	len_avg = (wt->wave->lpr_len + wt->wave->hpr_len) / 2;
	m = -2;
	n = -1;

	for (v = 0; v < len_cA; ++v) {
		i = v;
		m += 2;
		n += 2;
		X[m] = 0.0;
		X[n] = 0.0;
		for (l = 0; l < len_avg / 2; ++l) {
			t = 2 * l;
			if ((i - l) >= 0 && (i - l) < len_cA) {
				X[m] += wt->wave->lpr[t] * cA[i - l] + wt->wave->hpr[t] * cD[i - l];
				X[n] += wt->wave->lpr[t + 1] * cA[i - l] + wt->wave->hpr[t + 1] * cD[i - l];
			}
		}
	}
}


void idwt(wt_object wt, double *dwtop) {
	int J, U, i, lf, N, iter, k;
	int app_len, det_len;
	double *X_lp, *out;

	J = wt->J;
	U = 2;
	app_len = wt->length[0];
	out = (double*)malloc(sizeof(double)* (wt->siglength + 1));
	if (!strcmp(wt->ext, "per") && !strcmp(wt->cmethod, "direct")) {
		app_len = wt->length[0];
		det_len = wt->length[1];
		N = 2 * wt->length[J];
		lf = (wt->wave->lpr_len + wt->wave->hpr_len) / 2;

		X_lp = (double*)malloc(sizeof(double)* (N + 2 * lf - 1));
		iter = app_len;

		for (i = 0; i < app_len; ++i) {
			out[i] = wt->output[i];
		}

		for (i = 0; i < J; ++i) {

			//idwt1(wt, temp, cA_up, out, det_len, wt->output + iter, det_len, X_lp, X_hp, out);

			idwt_per(wt, out, det_len, wt->output + iter, det_len, X_lp);
			for (k = lf / 2 - 1; k < 2 * det_len + lf / 2 - 1; ++k) {
				out[k - lf / 2 + 1] = X_lp[k];
			}

			iter += det_len;
			det_len = wt->length[i + 2];
		}

		free(X_lp);

	}
	else if (!strcmp(wt->ext, "sym") && !strcmp(wt->cmethod, "direct")) {
		app_len = wt->length[0];
		det_len = wt->length[1];
		N = 2 * wt->length[J] - 1;
		lf = (wt->wave->lpr_len + wt->wave->hpr_len) / 2;

		X_lp = (double*)malloc(sizeof(double)* (N + 2 * lf - 1));
		iter = app_len;

		for (i = 0; i < app_len; ++i) {
			out[i] = wt->output[i];
		}

		for (i = 0; i < J; ++i) {

			//idwt1(wt, temp, cA_up, out, det_len, wt->output + iter, det_len, X_lp, X_hp, out);

			idwt_sym(wt, out, det_len, wt->output + iter, det_len, X_lp);
			for (k = lf - 2; k < 2 * det_len; ++k) {
				out[k - lf + 2] = X_lp[k];
			}

			iter += det_len;
			det_len = wt->length[i + 2];
		}

		free(X_lp);

	}
	else {
		printf("Signal extension can be either per or sym");
		exit(-1);
	}

	for (i = 0; i < wt->siglength; ++i) {
		dwtop[i] = out[i];
	}


	free(out);

}



void setDWTExtension(wt_object wt, char *extension) {
	if (!strcmp(extension, "sym")) {
		strcpy(wt->ext, "sym");
	}
	else if (!strcmp(extension, "per")) {
		strcpy(wt->ext, "per");
	}
	else {
		printf("Signal extension can be either per or sym");
		exit(-1);
	}
}


void setWTConv(wt_object wt, char *cmethod) {
	if (!strcmp(cmethod, "direct")) {
		strcpy(wt->cmethod, "direct");
	}
	else {
		printf("Convolution Only accepts direct convolution");
		exit(-1);
	}
}

void wave_summary(wave_object obj) {
	int i, N;
	N = obj->filtlength;
	printf("\n");
	printf("Wavelet Name : %s \n", obj->wname);
	printf("\n");
	printf("Wavelet Filters \n\n");
	printf("lpd : [");
	for (i = 0; i < N - 1; ++i) {
		printf("%g,", obj->lpd[i]);
	}
	printf("%g", obj->lpd[N - 1]);
	printf("] \n\n");
	printf("hpd : [");
	for (i = 0; i < N - 1; ++i) {
		printf("%g,", obj->hpd[i]);
	}
	printf("%g", obj->hpd[N - 1]);
	printf("] \n\n");
	printf("lpr : [");
	for (i = 0; i < N - 1; ++i) {
		printf("%g,", obj->lpr[i]);
	}
	printf("%g", obj->lpr[N - 1]);
	printf("] \n\n");
	printf("hpr : [");
	for (i = 0; i < N - 1; ++i) {
		printf("%g,", obj->hpr[i]);
	}
	printf("%g", obj->hpr[N - 1]);
	printf("] \n\n");
}

void wt_summary(wt_object wt) {
	int i;
	int J, t;
	J = wt->J;
	wave_summary(wt->wave);
	printf("\n");
	printf("Wavelet Transform : %s \n", wt->method);
	printf("\n");
	printf("Signal Extension : %s \n", wt->ext);
	printf("\n");
	printf("Convolutional Method : %s \n", wt->cmethod);
	printf("\n");
	printf("Number of Decomposition Levels %d \n", wt->J);
	printf("\n");
	printf("Length of Input Signal %d \n", wt->siglength);
	printf("\n");
	printf("Length of WT Output Vector %d \n", wt->outlength);
	printf("\n");
	printf("Wavelet Coefficients are contained in vector : %s \n", "output");
	printf("\n");
	printf("Approximation Coefficients \n");
	printf("Level %d Access : output[%d] Length : %d \n", 1, 0, wt->length[0]);
	printf("\n");
	printf("Detail Coefficients \n");
	t = wt->length[0];
	for (i = 0; i < J; ++i) {
		printf("Level %d Access : output[%d] Length : %d \n", i + 1, t, wt->length[i + 1]);
		t += wt->length[i + 1];
	}
	printf("\n");

}


void wave_free(wave_object object) {
	free(object);
}

void wt_free(wt_object object) {
	free(object);
}

