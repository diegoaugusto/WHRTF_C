/*
 *  whrtf.c
 *  WHRTF
 *
 *  Created by Diego Gomes on 24/05/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ReadHrtf.h"
#include "MathUtil.h"

#include <sys/time.h>
#include "runtime.h"

#define DB8_CONST "db8"
// Constants
const char* DB8 = "db8";
const int MIN_SAMPLE_NUMBER = 24;
const int MAX_SAMPLE_NUMBER = 174;
const int MAX_HRTF_SIZE = 512;
const int NUM_COEF_WITHOUT_BEGINNING = 488;	// 512-24 = 488
const int NUM_FILTROS = 4;
const char* WAVELET[] = {"db8", "db8", "db8", "db8"};


const double db8[2][8] = {
	{-0.0105, 0.0328, 0.0308, -0.1870, -0.0279, 0.6308, 0.7148, 0.2303},
	{-0.2303, 0.7148, -0.6308, -0.0279, 0.1870, 0.0308, -0.0328, -0.0105}};


/* Inclusão do respectivo módulo de definição */
#define WHRTF_SER
#include "whrtf.h"
#undef WHRTF_SER

// ########################################
// prototypes
double** getCoefSpars(int elev, int azim, char ear, int* G_size);
float* getRespImp(int numFiltros, double** G, int* G_size, int* resultLength);
float* whrtf (int elev, int azim, char ear, int* whrtfLength);


// Memory management
void cleanHostMemory(void* h);


// Find delay
double sumSquaresOfArrayElements(float* hrtf, int vecLength);
double sumArrayElements(float* vec0, int vecLength);
float** getSquaresOfArrayElements(float* vec0, float* vec1, int vecLength);
float** sumColumns(float* vec0, float* vec1, int vecLength);
short* findDelay(float** hrtf, int length);


// Sparse coefficients
double** leFiltros(char* filtro, int* filtroLength);
float*** computeIpAux(float*** Hp, int hlength, float*** Fp, int flength);
float** calculateGp(float** Rp, int rpLength, float*** Fp, int fpLength, float*** Ip, int ipLength);
void dec_poly(int hlength, float* sist, int sistLength, float** G, int* G_size, float*** Hp, float*** Fp);
double* getPartialVector(double* vec, int numOfElements);
float* convHost(float* Gl_old, int filtroLength, float* sparsArray1, int resultLength);
float* convHost1(float* Gl_old, int filtroLength, float* sparsArray1, int sparsLength);

float* shiftInParallel(float* vec, int vecLength, short delay, int maxLength);
void coef_spars(char* filtro[], int filtroLength, float* ho1d, int ho1dLength, float** G_aux, int* G_size);
float* resp_imp(char* filtros[], int numFiltros, double** G, int* G_size, int* resultLength);


void printDelay(short* delay) {
	for (int i = 0; i < 2; i++) {
		printf("delay[%d] = %d\n", i, delay[i]);
	}
}

void printDelayedHrtf(float* h, int vecLength) {
	for (int i = 0; i < vecLength; i++) {
		printf("h[%d] = %1.15f\n", i, h[i]);
	}
}


// Find Delay
/**
 Retorna o atraso (ITD) para cada ouvido usando algoritmo de busca diretamente 
 sobre a HRIR da direcao.
 D[0] = left;
 D[1] = right;
 */
short* findDelay(float** hrtf, int vecLength) {
	float** squaredElements = getSquaresOfArrayElements(hrtf[0], hrtf[1], vecLength);;
	
	double* et = (double*) malloc(2 * sizeof(double));
	et[0] = sumArrayElements(squaredElements[0], vecLength);
	et[1] = sumArrayElements(squaredElements[1], vecLength);
	
	float** ek = sumColumns(squaredElements[0], squaredElements[1], vecLength);;
	
	short* d = (short*) malloc(2*sizeof(short));
	for (int ear = 0; ear < 2; ear++) {	// ear
		for (int j = 0; j < vecLength; j++) {
			ek[ear][j] = ek[ear][j]/et[ear];
		}
		short* indexes = findIndexesGreaterThan(ek[ear], vecLength, 0.002f);
		d[ear] = indexes[0]-1;	
	}
	
	cleanHostMemory(ek);
	cleanHostMemory(et);
	cleanHostMemory(squaredElements);
	
	return d;
}

/**
 *	Soma paralelamente as colunas do vetor @vec de tamanho @vecLength.
 */
float** sumColumns(float* vec0, float* vec1, int vecLength) {
	// Variables
	float** h_result;
	
	int size = vecLength*sizeof(float);
	h_result = (float**)malloc(2 * sizeof(float*));
	h_result[0] = (float*)malloc(size);
	h_result[1] = (float*)malloc(size);
	
	for (int i = 0; i < vecLength; i++) {
		float sum0 = 0.0;
		float sum1 = 0.0;
		for (int j = 0; j < i; j++) {
			sum0 += vec0[j];
			sum1 += vec1[j];
		}		
		h_result[0][i] = (i == 0) ? vec0[i] : sum0;
		h_result[1][i] = (i == 0) ? vec1[i] : sum1;
	}
	
	return h_result;
}


/**
 *	Soma todos os elementos de um array.
 */
double sumArrayElements(float* vec, int vecLength) {
	// Sum result
	int i;
	double result = 0.0;
    for (i = 0; i < vecLength; ++i) {
		result += vec[i];
    }	
	
	return result;
}

/**
 *	Recupera os quadrados de cada elemento do array.
 */
float** getSquaresOfArrayElements(float* vec0, float* vec1, int vecLength) {
	// Variables
	float** h_result;
	int size = vecLength*sizeof(float);
	
	h_result = (float**)malloc(2 * sizeof(float*));
	h_result[0] = (float*) malloc(size);
	h_result[1] = (float*) malloc(size);
	
	for (int i = 0; i < vecLength; i++) {
		h_result[0][i] = (vec0[i] * vec0[i]);
		h_result[1][i] = (vec1[i] * vec1[i]);
	}
	
	return h_result;
}




// Shift in parallel
/**
 *	Desloca os elementos de um array em paralelo.
 */
float* shiftInParallel(float* vec, int vecLength, short delay, int maxLength) {
	// Variables
	float* h_result = NULL;
	
	h_result = (float*) calloc(maxLength, sizeof(float));
	
	for (int i = 0; i < (vecLength-delay); i++) {
		h_result[i] = vec[i+delay];
	}
	
	return h_result;
}



// Sparse coefficients

/**
 *	Função que lê os coeficientes do filtro daubechies 8
 */
double** leFiltros(char* filtro, int* filtroLength) {
	double** h = NULL;
	
	if (strcmp(DB8_CONST, filtro) == 0) {
		double db8[2][8] = {
			{-0.0105, 0.0328, 0.0308, -0.1870, -0.0279, 0.6308, 0.7148, 0.2303},
			{-0.2303, 0.7148, -0.6308, -0.0279, 0.1870, 0.0308, -0.0328, -0.0105}};
		
		int channels = 2;	// número de canais nesse filtro
		int filterLength = 8;	// tamanho do maior filtro
		*filtroLength = filterLength;
		
		h = (double**) calloc(channels, sizeof(double*));
		h[0] = (double*) calloc(filterLength, sizeof(double));
		h[1] = (double*) calloc(filterLength, sizeof(double));
		
		for (int row = 0; row < channels; row++) {
			for (int col = 0; col < filterLength; col++) {
				h[row][col] = db8[row][col];
			}
		}
	}
	return h;
}


/**
 *	Função que retorna o vetor com elementos esparsos (separados por 0),
 *	de acordo com um coeficiente de esparsidade.
 */
float* spars(double* vec, int vecLength, int sparsity, int* resultLength) {
	*resultLength = (((vecLength - 1) * sparsity) + 1);
	float* y = (float*) calloc((((vecLength - 1) * sparsity) + 1) , sizeof(float));
	
	for (int i = 0; i < vecLength; i++) {
		y[(i * sparsity)] = vec[i];
	}
	
	return y;
}


/**
 *	Operação de convolução executada no HOST.
 */
float* convHost(float* Gl_old, int filtroLength, float* sparsArray1, int resultLength) {
	float* convResult = (float*) calloc(resultLength, sizeof(float));
	for (int k = 0; k < resultLength; k++ ) {
		for (int j = 0; j < filtroLength; j++ ) {
			if ((k - j) >= 0) {
				convResult[k] += (Gl_old[j] * sparsArray1[k-j]);    // convolve: multiply and accumulate
			}
		}
	}
	return convResult;
}

float* convHost1(float* Gl_old, int filtroLength, float* sparsArray1, int sparsLength) {	
	int resultLength = filtroLength + sparsLength - 1;	
	int maxLength = (filtroLength >= sparsLength ? filtroLength : sparsLength);
	int minLength = (filtroLength <= sparsLength ? filtroLength : sparsLength);
	
	float* convResult = (float*) calloc(resultLength, sizeof(float));
	for (int k = 0; k < resultLength; k++ ) {
		for (int j = 0; j < maxLength; j++ ) {
			if ((k - j) >= 0 && (k-j) < minLength) {
				if (maxLength == filtroLength) {
					convResult[k] += (Gl_old[j] * sparsArray1[k-j]);    // convolve: multiply and accumulate
				} else {
					//printf("[%d] sparsArray1[%d] = %1.15f, Gl_old[%d] = %1.15f\n", k, j, sparsArray1[j], (k-j), Gl_old[k-j]);
					convResult[k] += (sparsArray1[j] * Gl_old[k-j]);    // convolve: multiply and accumulate
				}				
			}
		}
	}
	return convResult;
}



/**
 *	Operação de cascateamento dos filtros utilizando a convolução
 *	com FFT implementada em CUDA.
 */
void cascata(char* filtros[], int numFiltros, float** filterBank, int* filterBankLength) {
	float *convResult;
	
	int J = numFiltros;
	int filtroLength;
	double** h;
	int hLength;
	
	h = leFiltros(filtros[0], &hLength);
	
	float** filterBankAux = (float**) malloc((numFiltros + 1) * sizeof(float*));
	int* filterBankLengthAux = (int*) malloc((numFiltros + 1) * sizeof(int));
	filterBankAux[0] = (float*) malloc(hLength * sizeof(float));	// passa-altas
	filterBankLengthAux[0] = hLength;
	float* Gl = (float*) malloc(hLength * sizeof(float));	// passa-baixas
	
	for (int i = 0; i < hLength; i++) {
		filterBankAux[0][i] = h[1][i];
		Gl[i] = h[0][i];
	}
	
	float* Gl_old = NULL;
	filtroLength = hLength;
	
	for (int i = 1; i < numFiltros; i++) {
		Gl_old = (float*) calloc(filtroLength, sizeof(float));
		for (int j = 0; j < filtroLength; j++) {
			Gl_old[j] = Gl[j];
		}
		
		int sparsArray0Length, sparsArray1Length;
		float* sparsArray0 = spars(h[0], hLength, pow(2.0, i), &sparsArray0Length);
		float* sparsArray1 = spars(h[1], hLength, pow(2.0, i), &sparsArray1Length);
		
		int resultLength = (filtroLength + sparsArray0Length - 1);
		int resultSize = resultLength * sizeof(float);
		
		convResult = (float*) calloc(resultLength, sizeof(float));
		convResult = convHost1(Gl_old, filtroLength, sparsArray0, sparsArray0Length);
		
		filterBankAux[i] = (float*) calloc(resultLength, sizeof(float));
		filterBankAux[i] = convHost1(Gl_old, filtroLength, sparsArray1, sparsArray1Length);
		
		filterBankLengthAux[i] = resultLength;
		
		if ((i+1) != numFiltros) {
			free(Gl_old);	// Gl_old é usado fora do for após a última iteração
			free(Gl);
			
			filtroLength = resultLength;
			Gl = (float*) malloc(resultSize);
			for (int i = 0; i < resultLength; i++) {
				Gl[i] = convResult[i];
			}
			
			free(convResult);
		}
	}
	
	int sparsArray0Length;
	float* sparsArray0 = spars(h[0], hLength, pow(2.0, J-1), &sparsArray0Length);
	
	int resultLength = (filtroLength + sparsArray0Length - 1);
	
	filterBankAux[J] = (float*) calloc(resultLength, sizeof(float));
	filterBankAux[J] = convHost1(Gl_old, filtroLength, sparsArray0, sparsArray0Length);
	
	filterBankLengthAux[J] = resultLength;
	
	free(Gl_old);
	free(Gl);
		
	int maxLength = max(filterBankLengthAux, (numFiltros + 1));
	
	for (int i = numFiltros; i >= 0; i--) {
		filterBank[i] = (float*) calloc(maxLength, sizeof(float));
		filterBankLength[i] = filterBankLengthAux[numFiltros-i];
		for (int j = 0; j < filterBankLengthAux[numFiltros - i]; j++) {
			filterBank[i][j] = filterBankAux[numFiltros - i][j];
		}
	}
}

float* sumPolynomialCoefficients(float* dev_aux1, float* dev_aux2, int resultSize) {
	float* dev_aux_sum = (float*) calloc(resultSize, sizeof(float));
	for (int i = 0; i < resultSize; i++) {
		dev_aux_sum[i] = dev_aux1[i] + dev_aux2[i];
	}
	return dev_aux_sum;
}


float* multiplyPolynomialCoefficients(float* vecA, int vecASize, float* vecB, int vecBSize, int resultLength) {
	float* vecC = (float*) calloc(resultLength, sizeof(float));
	for (int i = 0; i < vecASize; i++) {
		for (int j = 0; j < vecBSize; j++) {
			vecC[i+j] = vecC[i+j] + (vecA[i] * vecB[j]);
		}
	}
	return vecC;
}


/**
 * Essa função computa a multiplicação das matrizes Fp e Hp. Ip_aux = Fp*Hp.
 * Resultado: [2][2][7]
 *
 *	Fp =	a	b
 *			c	d
 *
 *	Hp =	A	B
 *			C	D
 *
 *	Resultado mais rápido que método anterior com streams.
 */
float*** computeIpAux(float*** Hp, int hlength, float*** Fp, int flength) {
	float *dev_aux1, *dev_aux2;
	int resultLength = (hlength + flength - 1);
	
	float*** Ip_aux = (float***) malloc(2 * sizeof(float**));
	Ip_aux[0] = (float**) malloc(2 * sizeof(float*));
	Ip_aux[1] = (float**) malloc(2 * sizeof(float*));
	
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			Ip_aux[i][j] = (float*) malloc(resultLength*sizeof(float*));
		}
	}
	
	/*
	 Ip_aux = Fp*Hp;
	 O resultado é uma matriz 2x2x7.
	 
	 aA+bC	aB+bD
	 cA+dC	cB+dD
	 */
	
	/* 00 */	
	dev_aux1 = multiplyPolynomialCoefficients(Fp[0][0], hlength, Hp[0][0], hlength, resultLength);
	dev_aux2 = multiplyPolynomialCoefficients(Fp[0][1], hlength, Hp[1][0], hlength, resultLength);
	Ip_aux[0][0] = sumPolynomialCoefficients(dev_aux1, dev_aux2, resultLength);
	
	
	/* 01 */
	// clean auxiliary variables
	free(dev_aux1); dev_aux1 = NULL;
	free(dev_aux2); dev_aux2 = NULL;
	
	dev_aux1 = multiplyPolynomialCoefficients(Fp[0][0], hlength, Hp[0][1], hlength, resultLength);
	dev_aux2 = multiplyPolynomialCoefficients(Fp[0][1], hlength, Hp[1][1], hlength, resultLength);
	Ip_aux[0][1] = sumPolynomialCoefficients(dev_aux1, dev_aux2, resultLength);
	
	
	/* 10 */
	// clean auxiliary variables
	free(dev_aux1); dev_aux1 = NULL;
	free(dev_aux2); dev_aux2 = NULL;
	
	dev_aux1 = multiplyPolynomialCoefficients(Fp[1][0], hlength, Hp[0][0], hlength, resultLength);
	dev_aux2 = multiplyPolynomialCoefficients(Fp[1][1], hlength, Hp[1][0], hlength, resultLength);
	Ip_aux[1][0] = sumPolynomialCoefficients(dev_aux1, dev_aux2, resultLength);
	
	
	/* 11 */
	// clean auxiliary variables
	free(dev_aux1); dev_aux1 = NULL;
	free(dev_aux2); dev_aux2 = NULL;
	
	dev_aux1 = multiplyPolynomialCoefficients(Fp[1][0], hlength, Hp[0][1], hlength, resultLength);
	dev_aux2 = multiplyPolynomialCoefficients(Fp[1][1], hlength, Hp[1][1], hlength, resultLength);
	Ip_aux[1][1] = sumPolynomialCoefficients(dev_aux1, dev_aux2, resultLength);
	
	free(dev_aux1); dev_aux1 = NULL;
	free(dev_aux2); dev_aux2 = NULL;
	
	return Ip_aux;
}


/**
 * Essa função computa a multiplicação das matrizes Rp, Fp e Ip_aux. Ip_aux = Fp*Hp.
 * Resultado: [2][237]
 *
 *	Rp =	X
 *			Y
 *
 *	Fp =	a	b
 *			c	d
 *
 *	Ip =	A	B
 *			C	D
 *
 */
float** calculateGp(float** Rp, int rpLength, float*** Fp, int fpLength, float*** Ip, int ipLength) {
	float *dev_aux1, *dev_aux2, *dev_aux3, *dev_aux4, *dev_aux_sum1, *dev_aux_sum2 ;
	float **finalResult;
	
	int lengthX = (rpLength % 2 == 0 ? rpLength/2 : (rpLength/2+1));
	int lengthY = rpLength/2;
	int partialResultLengthX = lengthX + fpLength - 1;
	int partialResultLengthY = lengthY + fpLength - 1;
	int finalResultLength = (partialResultLengthX + ipLength -1);
	int finalResultSize = finalResultLength * sizeof(float);
	
	/*
	 Gp = Rp*Fp*Ip;
	 O resultado é uma matriz 2x237.
	 
	 ----- Parte 1
	 Rp	  Fp
	 (x	y)	(a	c)
	 (b	d)
	 
	 Rp*Fp
	 xa+yb	xc+yd
	 
	 */
	// Invoke kernel
	dev_aux1 = multiplyPolynomialCoefficients(Rp[0], lengthX, Fp[0][0], fpLength, partialResultLengthX);
	dev_aux2 = multiplyPolynomialCoefficients(Rp[1], lengthY, Fp[0][1], fpLength, partialResultLengthX);
	dev_aux3 = multiplyPolynomialCoefficients(Rp[0], lengthX, Fp[1][0], fpLength, partialResultLengthX);
	dev_aux4 = multiplyPolynomialCoefficients(Rp[1], lengthY, Fp[1][1], fpLength, partialResultLengthX);
	
	dev_aux_sum1 = sumPolynomialCoefficients(dev_aux1, dev_aux2, partialResultLengthX);
	dev_aux_sum2 = sumPolynomialCoefficients(dev_aux3, dev_aux4, partialResultLengthX);
	
	// Re-Inicialização de array no dispositivo para evitar erros em leituras
	free(dev_aux1); dev_aux1 = NULL;
	free(dev_aux2); dev_aux2 = NULL;
	free(dev_aux3); dev_aux3 = NULL;
	free(dev_aux4); dev_aux4 = NULL;
	
	/*
	 Gp = Rp*Fp*Ip;
	 O resultado é uma matriz 2x237.
	 
	 ----- Parte 1
	 Rp	  Fp
	 (x	y)	(a	c)
	 (b	d)
	 
	 Rp*Fp
	 xa+yb	xc+yd
	 u		  v
	 
	 ----- Parte 2	  	  
	 Rp*Fp	  Hp
	 (u	v)	(A	C)
	 (B	D)
	 
	 Rp*Fp*Ip
	 uA+vB	uC+vD
	 
	 */	
	finalResult = (float**) malloc(2 * sizeof(float*));
	finalResult[0] = (float*) malloc(finalResultSize);
	finalResult[1] = (float*) malloc(finalResultSize);
	
	dev_aux1 = multiplyPolynomialCoefficients(dev_aux_sum1, partialResultLengthX, Ip[0][0], ipLength, finalResultLength);
	dev_aux2 = multiplyPolynomialCoefficients(dev_aux_sum2, partialResultLengthX, Ip[0][1], ipLength, finalResultLength);
	dev_aux3 = multiplyPolynomialCoefficients(dev_aux_sum1, partialResultLengthX, Ip[1][0], ipLength, finalResultLength);
	dev_aux4 = multiplyPolynomialCoefficients(dev_aux_sum2, partialResultLengthX, Ip[1][1], ipLength, finalResultLength);
	
	finalResult[0] = sumPolynomialCoefficients(dev_aux1, dev_aux2, finalResultLength);
	finalResult[1] = sumPolynomialCoefficients(dev_aux3, dev_aux4, finalResultLength);
	
	// Cleaning data of auxiliary vectors
	free(dev_aux_sum1); dev_aux_sum1 = NULL;
	free(dev_aux_sum2); dev_aux_sum2 = NULL;
	free(dev_aux1); dev_aux1 = NULL;
	free(dev_aux2); dev_aux2 = NULL;
	free(dev_aux3); dev_aux3 = NULL;
	free(dev_aux4); dev_aux4 = NULL;
	
	return finalResult;
}

/**
 *	Função que realiza a operação de decimação polinomial.
 *
 *  Retorna os coeficientes G (matriz 2xK) que correspondem ao 
 *	sistema R decomposto pelo banco de 2 canais H.
 */
void dec_poly(int hlength, float* sist, int sistLength, float** G, int* G_size, float*** Hp, float*** Fp) {	
	float** Rp = (float**) calloc(2, sizeof(float*));
	Rp[0] = (float*) calloc((sistLength % 2 == 0 ? (sistLength/2) : (sistLength/2 + 1)), sizeof(float));
	Rp[1] = (float*) calloc((sistLength/2), sizeof(float));
	
	for (int i = 0; i < sistLength; i++) {
		if (i % 2 == 0) {
			Rp[0][i/2] = sist[i];
		} else {
			Rp[1][i/2] = sist[i];
		}
	}
	
	/*
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 % decomposicao polifasica do sistema R - matriz Rp
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 */
	int m = 2;
	int n = hlength;
	
	/*
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 % decomposicao polifasica dos bancos de sintese e analise
	 % matrizes Fp e Hp
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 */
	
	/*
	 Ip_aux = Fp*Hp;
	 O resultado é uma matriz 2x2x7.
	 
	 H[0][1]*F[0][0] + H[0][0]*F[1][0]	H[0][1]*F[0][1] + H[0][0]*F[1][1]
	 H[1][1]*F[0][0] + H[1][0]*F[1][0]	H[1][1]*F[0][1] + H[1][0]*F[1][1]
	 */
	// Hp e Fp são matrizes (cubos) de dimensão: [2][2][7]	
	float*** ip_aux = computeIpAux(Hp, hlength/2, Fp, hlength/2);
	
	// Gp = Rp*Fp*Ip
	float** Gp = calculateGp(Rp, sistLength, Fp, (hlength/2), ip_aux, (hlength-1));
	
	int atraso = round((n-m)/2.0);
	double esparsidade = 2.0;
	
	*G_size = ceil( ((sistLength + atraso +1)/esparsidade) +1);
	
	G[0] = (float*) malloc(*G_size * sizeof(float));
	G[1] = (float*) malloc(*G_size * sizeof(float));
	
	for (int i = 0; i < *G_size; i++) {
		G[0][i] = Gp[0][i+atraso];
		G[1][i] = Gp[1][i+atraso];
	}
}

/**
 *	Obtém os coeficientes esparsos que equivalem o sistema ho1d.
 */
void coef_spars(char* filtro[], int numFiltros, float* ho1d, int ho1dLength, float** G_aux, int* G_size) {
	int mesmofiltro = 0;
	float* sist = ho1d;
	double** h;
	float** G;
	float** coefSpars;
	int filtroLength;
	int sistLength = ho1dLength;
	
	coefSpars = (float**) malloc((numFiltros+1) * sizeof(float*));
	
	G = (float**) malloc(2 * sizeof(float*));	// sempre tem tamanho 2
	
	if (!mesmofiltro) {
		h = leFiltros(filtro[0], &filtroLength);
		mesmofiltro = (0 < numFiltros-1 && strcmp(filtro[0], filtro[1]) == 0);
	}
	
	// Banco de análise
	double** H = (double**) malloc(2 * sizeof(double*));
	H[0] = (double*) malloc(filtroLength * sizeof(double));
	H[1] = (double*) malloc(filtroLength * sizeof(double));
	
	// Banco de síntese
	double** F = (double**) malloc(2 * sizeof(double*));
	F[0] = (double*) malloc(filtroLength * sizeof(double));
	F[1] = (double*) malloc(filtroLength * sizeof(double));
	
	/*
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 % Obtencao do banco de filtros
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 */
	int power = 1;
	for (int i = filtroLength-1; i >= 0; i--) {
		int reverseIndex = filtroLength-i-1;
		double powerValue = pow((double)-1, (double)power++);
		if (h[0][reverseIndex] + h[1][i]*powerValue != 0) {
			break;	// Não é Daubechies
		} else {
			H[0][reverseIndex] = h[0][reverseIndex];
			H[1][reverseIndex] = h[0][i]*powerValue;
			F[0][reverseIndex] = h[1][reverseIndex]*powerValue;
			F[1][reverseIndex] = h[1][i];
		}
	}
	
	float*** Hp = (float***) malloc(2 * sizeof(float**));
	Hp[0] = (float**) malloc(2 * sizeof(float*));
	Hp[1] = (float**) malloc(2 * sizeof(float*));
	
	float*** Fp = (float***) malloc(2 * sizeof(float**));
	Fp[0] = (float**) malloc(2 * sizeof(float*));
	Fp[1] = (float**) malloc(2 * sizeof(float*));
	
	for (int i = 0; i < 2; i++) {
		Hp[0][i] = (float*) malloc((filtroLength/2) * sizeof(float));
		Hp[1][i] = (float*) malloc((filtroLength/2) * sizeof(float));
		Fp[0][i] = (float*) malloc((filtroLength/2) * sizeof(float));
		Fp[1][i] = (float*) malloc((filtroLength/2) * sizeof(float));
		for (int j = 0; j < (filtroLength/2); j++) {
			Hp[0][i][j] = H[i][j*2];
			Hp[1][i][j] = H[i][j*2+1];
			
			if (i == 0) {
				Fp[0][i][j] = F[0][j*2+1];
				Fp[1][i][j] = F[1][j*2+1];
			} else {
				Fp[0][i][j] = (-1) * F[1][filtroLength*i - (j*2+1)];
				Fp[1][i][j] = F[0][filtroLength*i - (j*2+1)];
			}
		}
	}
	
	for (int j = 0; j < numFiltros; j++) {		
		dec_poly(filtroLength, sist, sistLength, G, &G_size[numFiltros-j], Hp, Fp);
		
		sistLength = G_size[numFiltros-j];
		coefSpars[numFiltros-j] = G[1];
		
		sist = NULL;
		sist = G[0];
	}
	
	G_size[0] = G_size[1];
	coefSpars[0] = (float*) malloc(G_size[0] * sizeof(float));
	for (int k = 0; k < G_size[0]; k++) {
		coefSpars[0][k] = sist[k];
	}
	
	int maxG_size = max(G_size, (numFiltros+1));
	for (int k = 0; k < (numFiltros+1); k++) {
		G_aux[k] = (float*) malloc(maxG_size * sizeof(float));
		int cont = 0;
		while (cont < G_size[k]) {
			G_aux[k][cont] = coefSpars[k][cont];
			cont++;
		}
		while (cont < maxG_size) {
			G_aux[k][cont] = 0.0;
			cont++;
		}
	}	
}

/**
 *	Função que retorna os primeiros X elementos de um vetor.
 */
double* getPartialVector(double* vec, int numOfElements) {
    double* partialVec = (double*) malloc(numOfElements * sizeof(double));
	for (int i = 0; i < numOfElements; i++) {
		partialVec[i] = vec[i];
	}
	return partialVec;
}


float* sumColumnsKernel(float* dev_aux1, int aux1Length, float* dev_aux_sum, int resultSize) {
	for (int i = 0; i < resultSize && i < aux1Length; i++) {
		dev_aux_sum[i] += dev_aux1[i];
	}
    return dev_aux_sum;
}


/**
 *	Função que retorna a resposta impulsiva a partir dos coeficientes
 *	esparsos.
 */
float* resp_imp(char* filtros[], int numFiltros, double** G, int* G_size, int* resultLength) {
	float** filterBank = (float**) calloc((numFiltros + 1), sizeof(float));
	int* filterBankLength = (int*) calloc((numFiltros + 1), sizeof(int));
	
	cascata(filtros, numFiltros, filterBank, filterBankLength);
	
	// TODO implementar calc_delta
	int atrasos[5] = {1, 1, 8, 22, 50};
	
	int* L = (int*) calloc(numFiltros + 1, sizeof(int));
	L[0] = (int) pow(2.0, numFiltros);
	for (int i = 1; i < (numFiltros + 1); i++) {
		L[i] = (int) pow(2.0, (numFiltros+1)-i);
	}
	
	float** gaux = (float**) calloc((numFiltros+1), sizeof(float*));
	float** r = (float**) calloc((numFiltros+1), sizeof(float*));
	int* r_sizes = (int*) calloc((numFiltros+1), sizeof(int));
	
	for (int i = 0; i < numFiltros+1; i++) {
		int partialVecLength = (G_size[i]+atrasos[i]);
		int resultLength;
		gaux[i] = spars(getPartialVector(G[i], partialVecLength ), partialVecLength, L[i], &resultLength);
		
		// conv(filterBank[i], gaux[i])	
		int convLength = (filterBankLength[i] + resultLength - 1);
		
		r[i] = convHost1(filterBank[i], filterBankLength[i], gaux[i], resultLength);
		
		r_sizes[i] = convLength;	
	}
	

	int maxR = max(r_sizes, numFiltros+1);
	float* res = (float*) calloc(maxR, sizeof(float));
	float** raux = (float**) malloc((numFiltros + 1) * sizeof(float*));
	
	for (int i = 0; i < (numFiltros + 1); i++) {
		raux[i] = (float*) malloc(maxR * sizeof(float));
		for (int j = 0; j < maxR; j++) {
			if (j < r_sizes[i]) {
				raux[i][j] = r[i][j];
			} else {
				raux[i][j] = 0.0;
			}
		}
	}
		
	for (int i = 0; i < (numFiltros + 1); i++) {
		res = sumColumnsKernel(r[i], r_sizes[i], res, maxR);
	}
	
	int maxHlength = max(filterBankLength, numFiltros+1);
	float* resFinal = (float*) calloc((maxR - maxHlength), sizeof(float));
	
	for (int i = 0; i < (maxR-maxHlength); i++) {
		resFinal[i] = res[maxHlength + i];
	}
	
	*resultLength = (maxR-maxHlength);
	
	return resFinal;
}




// Memory management
/**
 *	Libera memória do Host.
 */
void cleanHostMemory(void* h) {
    // Free host memory
    if (h) {
        free(h);
	}
}











// ##################################################################
/**
 *	Esta função retorna os coeficientes esparsos de uma hrtf com elevação e azimute conhecidos.
 */
double** getCoefSpars(int elev, int azim, char ear, int* G_size) {
	float** hrtf = readHrtf(elev, azim, 'L', MIN_SAMPLE_NUMBER, MAX_SAMPLE_NUMBER);	
	short* delay = findDelay(hrtf, MAX_SAMPLE_NUMBER - MIN_SAMPLE_NUMBER);
	free(hrtf);
	
	//printDelay(delay);
	
	// hrtf sem os primeiros 24 coeficientes
	float** ho = readHrtf(elev, azim, 'R', MIN_SAMPLE_NUMBER, MAX_HRTF_SIZE);
	float* ho1d = NULL;
	
	int length = NUM_COEF_WITHOUT_BEGINNING - ((delay[0] > delay[1]) ? delay[0] : delay[1]);
	
	if (ear == 'L') {
		ho1d = shiftInParallel(ho[0], NUM_COEF_WITHOUT_BEGINNING, delay[0], length);
	} else if (ear == 'R') {
		ho1d = shiftInParallel(ho[1], NUM_COEF_WITHOUT_BEGINNING, delay[1], length);
	}
	
	free(ho);
	
	//printDelayedHrtf(ho1d, length);
	
	float** G_aux = NULL;
	G_aux = (float**) malloc((NUM_FILTROS+1) * sizeof(float*));
	
	coef_spars(WAVELET, NUM_FILTROS, ho1d, length, G_aux, G_size);
	
	// TODO: implementar atraso = calc_delta(Wavelet);
	int atraso[5] = {1, 1, 8, 22, 50};
	
	int novoG_size[NUM_FILTROS+1];
	for (int i = 0; i < NUM_FILTROS+1; i++) {
		novoG_size[i] = atraso[i] + G_size[i];
	}
	
	int maxValueOfGSize = max(novoG_size, NUM_FILTROS+1);
	double** G = (double**) malloc((NUM_FILTROS + 1) * sizeof(double*));
	
	for (int i = 0; i < (NUM_FILTROS + 1); i++) {
		G[i] = (double*) calloc(maxValueOfGSize, sizeof(double));
		
		int j = 0;
		for (j; j < atraso[i]; j++) {
			G[i][j] = 0.0;
		}
		for (j; j < (G_size[i]+atraso[i]); j++) {
			G[i][j] = G_aux[i][j-atraso[i]];
		}
		for (j; j < (maxValueOfGSize - (G_size[i] + atraso[i])); j++) {
			G[i][j] = 0.0;
		}
	}
	
	return G;
}

/**
 *	Esta função retorna a resposta impulsiva tendo como entrada os coeficientes esparsos
 */
float* getRespImp(int numFiltros, double** G, int* G_size, int* resultLength) {
	int auxWhrtfLength;
	float* respImp = resp_imp(WAVELET, NUM_FILTROS, G, G_size, &auxWhrtfLength);
	*resultLength = auxWhrtfLength;
	return respImp;
}

// Host code
float* whrtf (int elev, int azim, char ear, int* whrtfLength) {	
	int* G_size = (int*) calloc((NUM_FILTROS+1), sizeof(int));
	double** G = getCoefSpars(elev, azim, ear, G_size);
	
	float* whrtf = getRespImp(NUM_FILTROS, G, G_size, &*whrtfLength);
	return whrtf;
}