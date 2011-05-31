/*
*  main.cpp
*  WHRTF
*
*  Created by Diego Gomes on 24/05/11.
*  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include <fftw3.h>

#include "runtime.h"

/* Inclusão do respectivo módulo de definição */
#define WHRTF_SER
#include "whrtf.h"
#undef WHRTF_SER

#define SPARSECOEF_SER
#include "sparseCoefficients.h"
#undef SPARSECOEF_SER


#define NUM_FILTROS 4

float* convHost2(float* Gl_old, int filtroLength, float* sparsArray1, int sparsLength) {
	int resultLength = filtroLength + sparsLength - 1;
	int maxLength = (filtroLength > sparsLength ? filtroLength : sparsLength);
	float* convResult = (float*) calloc(resultLength, sizeof(float));
	for (int k = 0; k < resultLength; k++ ) {
		for (int j = 0; j < maxLength; j++ ) {
			if ((k - j) >= 0) {
				if (maxLength == filtroLength) {
					convResult[k] += (Gl_old[j] * sparsArray1[k-j]);    // convolve: multiply and accumulate
				} else {
					convResult[k] += (sparsArray1[j] * Gl_old[k-j]);    // convolve: multiply and accumulate
				}				
			}
		}
	}
	return convResult;
}

int main0 (int argc, char * const argv[]) {
	float b[5] = {1.0,2.0,3.0,4.0,5.0};
	float a[7] = {1.0, 0.0, 2.0,0.0, 3.0, 0.0, 4.0};
	
	fftw_complex *in, *out; 
	fftw_plan p;
	
	float* c = NULL;
	
	c = convHost1(a, 7, b, 5);
	
	for (int i = 0; i < 11; i++) {
		printf("c[%d] = %1.15f\n", i, c[i]);
	}
}


int main (int argc, char * const argv[]) {
    int elev = 0;
	int azim = 90;
	char ear = 'L';
	
	INIT_VARIABLES;
	INIT_RUNTIME;
	int whrtfLength;
	float* whr = whrtf(elev, azim, ear, &whrtfLength);
	END_RUNTIME;
	printf("\n[whrtf]: ");
	PRINT_RUNTIME;
	
	/*for (int i = 0; i < whrtfLength; i=i+3) {
		printf("\nColunas %d a %d \n", i, i+2);
		printf("%1.15f, %1.15f, %1.15f\n", whr[i], whr[i+1], whr[i+2]);
	}*/
	// OK - REVISADO ATÉ AQUI - 24/05/2011
}


int main2 (int argc, char * const argv[]) {
	int elev = 0;
	int azim = 0;
	char ear = 'L';
	
	double** Gl = NULL;
	int* Gl_size = NULL;
	double** Gr = NULL;
	int* Gr_size = NULL;
	
	INIT_VARIABLES;
	INIT_RUNTIME;
	
	for (azim = 0; azim < 360; azim++) {
		int flipAzim = 360 - azim;
		if (flipAzim == 360) {
			flipAzim = 0;
		}
		
		//printf("\n\n\n\nAzim = %d\n", azim);
		
		Gl_size = (int*) calloc((NUM_FILTROS+1), sizeof(int));
		Gl = getSparseCoefficients(elev, azim, ear, Gl_size);
		
		Gr_size = (int*) calloc((NUM_FILTROS+1), sizeof(int));
		Gr = getSparseCoefficients(elev, flipAzim, ear, Gr_size);
		
		/*printf("Gl_size\n");
		for (int i = 0; i < NUM_FILTROS+1; i++) {
			printf("%d\t\t", Gl_size[i]);
		}
		printf("\n\nGr_size\n");
		for (int i = 0; i < NUM_FILTROS+1; i++) {
			printf("%d\t\t", Gr_size[i]);
		}*/
		
		int whrtfLengthL, whrtfLengthR;
		float* whrtfL = getRespImp(NUM_FILTROS, Gl, Gl_size, &whrtfLengthL);
		float* whrtfR = getRespImp(NUM_FILTROS, Gr, Gr_size, &whrtfLengthR);
		
		free(Gl);
		free(Gr);
		free(Gl_size);
		free(Gr_size);
		free(whrtfL);
		free(whrtfR);
		
		Gl = NULL;
		Gr = NULL;
		Gl_size = NULL;
		Gr_size = NULL;
		whrtfL = NULL;
		whrtfR = NULL;
 	}
	
	END_RUNTIME;
	printf("\n[whrtf]: ");
	PRINT_RUNTIME;
	
	// TODO implementar calc_delta
	int atrasos[5] = {1, 1, 8, 22, 50};

		
}