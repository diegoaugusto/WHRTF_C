/*
 *  sparseCoefficients.c
 *  WHRTF
 *
 *  Created by Diego Gomes on 27/05/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "whrtf.h"

/* Inclusão do respectivo módulo de definição */
#define SPARSECOEF_SER
#include "sparseCoefficients.h"
#undef SPARSECOEF_SER

#define NUM_FILTROS 4

// Constants
const int NUM_KNOWN_ELEVATIONS = 14;
const int elevations[] = {-40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
const int atrasos[5] = {1, 1, 8, 22, 50};


// Typedef
enum boolean {
    true = 1, false = 0
};
typedef enum boolean bool;


// ########################################
// prototypes
bool isKnownElev(int elev);
void getElevationBoundaries(int elev, int* minElev, int* maxElev);
int* getAzimuths(int elev, int* numAzimuths);
void getMinMax(int elev, int azim, int* azimuths, int numAzimuths, int* min, int* max);
void getAzimuthBoundaries(int elev, int azim, int* minBound, int* maxBound);
double** getGp(double** Ga, int* Ga_size, double** Gb, int* Gb_size, int minBound, int maxBound, int azim, int* Gp_size);
double** getWhrtfWithKnownElev(int elev, int azim, char ear, int* Gp_size);
double** getSparseCoefficients(int elev, int azim, int ear, int* Gp_size);


/**
 *	Verifica se é uma elevação conhecida.
 */
bool isKnownElev(int elev) {
	bool knowElev = false;
	
	for (int i = 0; i < NUM_KNOWN_ELEVATIONS; i++) {
		if (elevations[i] == elev) {
			knowElev = true;
			break;
		}
	}
	
	return knowElev;
}

/**
 *	Recupera o intervalo de elevações.
 */
void getElevationBoundaries(int elev, int* minElev, int* maxElev) {
	for (int i = 0; i < NUM_KNOWN_ELEVATIONS - 1; i++) {
		if (elev == elevations[i]) {
			*minElev = elevations[i];
			*maxElev = elevations[i];
			break;
		} else if ((elevations[i] < elev) && (elevations[i+1] > elev)) {
			*minElev = elevations[i];
			*maxElev = elevations[i+1];
			break;
		}
	}
}

int* getAzimuths(int elev, int* numAzimuths) {
	int intervalo;
	int* azimuths = NULL;
	int j = 0;
	
	switch (elev) {
		case -20:
		case -10:
		case 0:
		case 10:
		case 20: {
			intervalo = 5;
			*numAzimuths = (int) round(360/intervalo) + 1;
			azimuths = (int*) calloc(*numAzimuths, sizeof(int));
			for (int i = 0; i <= 360; i = i + intervalo) {
				azimuths[j++] = i;
			}
			break;
		}
		case -30:
		case 30: {
			intervalo = 6;
			*numAzimuths = (int) round(360/intervalo) + 1;
			azimuths = (int*) calloc(*numAzimuths, sizeof(int));
			for (int i = 0; i <= 360; i = i + intervalo) {
				azimuths[j++] = i;
			}
			break;
		}
		case -40:
		case 40: {
			float interval = 6.43f;
			*numAzimuths = (int) round(360/interval) + 1;
			azimuths = (int*) calloc(*numAzimuths, sizeof(int));
			
			for (int i = 0; i <= 360; i = (int) round(j*interval - 0.01)) {	// problema quando j = 50. i == 321.5 e arredondava para 322.
				azimuths[j++] = i;
			}
			break;
		}
		case 50: {
			intervalo = 8;
			*numAzimuths = (int) round(360/intervalo) + 1;
			azimuths = (int*) calloc(*numAzimuths, sizeof(int));
			for (int i = 0; i <= 360; i = i + intervalo) {
				azimuths[j++] = i;
			}
			break;
		}
		case 60: {
			intervalo = 10;
			*numAzimuths = (int) round(360/intervalo) + 1;
			azimuths = (int*) calloc(*numAzimuths, sizeof(int));
			for (int i = 0; i <= 360; i = i + intervalo) {
				azimuths[j++] = i;
			}
			break;
		}
		case 70: {
			intervalo = 15;
			*numAzimuths = (int) round(360/intervalo) + 1;
			azimuths = (int*) calloc(*numAzimuths, sizeof(int));
			for (int i = 0; i <= 360; i = i + intervalo) {
				azimuths[j++] = i;
			}
			break;
		}
		case 80: {
			intervalo = 30;
			*numAzimuths = (int) round(360/intervalo) + 1;
			azimuths = (int*) calloc(*numAzimuths, sizeof(int));
			for (int i = 0; i <= 360; i = i + intervalo) {
				azimuths[j++] = i;
			}
			break;
		}
		case 90: {
			*numAzimuths = 1;
			azimuths = (int*) calloc(*numAzimuths, sizeof(int));
			azimuths[0] = 0;
			break;
		}
		default:
			break;
	}
	
	return azimuths;
}

void getMinMax(int elev, int azim, int* azimuths, int numAzimuths, int* min, int* max) {
	for (int i = 0; i < numAzimuths-1; i++) {
		if (azim == azimuths[i]) {
			*min = azimuths[i];
			*max = azimuths[i];
			break;
		} else if ((azim > azimuths[i]) && (azim < azimuths[i+1])) {
			*min = azimuths[i];
			*max = azimuths[i+1];
			break;
		}
	}
	if (*max == 360) {
		*max = 0;
	}
}

void getAzimuthBoundaries(int elev, int azim, int* minBound, int* maxBound) {
	int min, max;
	if (isKnownElev(elev)) {
		int numAzimuths = 0;
		int* azimuths = getAzimuths(elev, &numAzimuths);
		if (numAzimuths > 0) {
			getMinMax(elev, azim, azimuths, numAzimuths, &min, &max);
			*minBound = min;
			*maxBound = max;
		}
	} else {
		*minBound = -1;
		*maxBound = -1;
	}
}

/**
 *	Função que gera coeficientes esparsos ponderados a partir de dois outros conjuntos 
 *	de coeficientes esparsos Ga e Gb.
 */
double** getGp(double** Ga, int* Ga_size, double** Gb, int* Gb_size, int minBound, int maxBound, int azim, int* Gp_size) {
	if (maxBound == 0) {
		maxBound = 360;
	}
	double pesoA = ((double)(maxBound - azim) / (maxBound - minBound));
	double pesoB = 1 - pesoA;
	
	for (int i = 0; i < 5; i++) {
		//printf("Ga_size[%d] = %d, Gb_size[%d] = %d\n", i, Ga_size[i], i, Gb_size[i]);
	}
	
	int minLength = (Ga_size[NUM_FILTROS] < Gb_size[NUM_FILTROS]) ? Ga_size[NUM_FILTROS] : Gb_size[NUM_FILTROS];
	int minLengthWoAtraso = minLength;
	
	minLength += atrasos[NUM_FILTROS];
	
	double** Gp = (double**) calloc(NUM_FILTROS+1, sizeof(double*));
	
	// TODO implementar em CUDA
	for (int i = 0; i < (NUM_FILTROS + 1); i++) {
		Gp[i] = (double*) calloc(minLength, sizeof(double));
		
		if (Ga_size[i] > minLengthWoAtraso && Gb_size[i] > minLengthWoAtraso) {
			Gp_size[i] = minLengthWoAtraso;
		} else if (Ga_size[i] > Gb_size[i]) {
			Gp_size[i] = (Ga_size[i] > minLengthWoAtraso) ? minLengthWoAtraso : Ga_size[i];
		} else if (Gb_size[i] > Ga_size[i]) {
			Gp_size[i] = (Gb_size[i] > minLengthWoAtraso) ? minLengthWoAtraso : Gb_size[i];
		} else if (Ga_size[i] == Gb_size[i]) {
			Gp_size[i] = Ga_size[i];
		}


		for (int j = 0; j < minLength; j++) {
			Gp[i][j] = (pesoA * Ga[i][j]) + (pesoB * Gb[i][j]);
		}
	}
	
	return Gp;
}

/**
 *	Recupera os coeficientes esparsos para qualquer azimute em uma elevação conhecida.
 */
double** getWhrtfWithKnownElev(int elev, int azim, char ear, int* Gp_size) {
	double** Gp = NULL;
	int minBound, maxBound;
	getAzimuthBoundaries(elev, azim, &minBound, &maxBound);
	
	if (minBound != maxBound) {
		int* Ga_size = (int*) calloc((NUM_FILTROS+1), sizeof(int));
		double** Ga = getCoefSpars(elev, minBound, ear, Ga_size);
		
		int* Gb_size = (int*) calloc((NUM_FILTROS+1), sizeof(int));
		double** Gb = getCoefSpars(elev, maxBound, ear, Gb_size);
		
		Gp = getGp(Ga, Ga_size, Gb, Gb_size, minBound, maxBound, azim, Gp_size);
		
	} else if (minBound == maxBound) {
		Gp = getCoefSpars(elev, minBound, ear, Gp_size);
	}
	
	return Gp;
}


// Returning Gp e Gp_size
double** getSparseCoefficients(int elev, int azim, int ear, int* Gp_size) {
	double** Gp = NULL;
	bool knownElev = isKnownElev(elev);
	if (knownElev) {
		//printf("------> Elevação [%d] é conhecida.\n", elev);
		Gp = getWhrtfWithKnownElev(elev, azim, ear, Gp_size);
	} else {
		//printf("------> Elevação [%d] NÃO é conhecida.\n", elev);
		// Recuperar max e min das elevations
		int minElev, maxElev;
        getElevationBoundaries(elev, &minElev, &maxElev);
		
		// Recuperar Hrtfs neste azimute para estas elevations
		int* Ga_size = (int*) calloc((NUM_FILTROS+1), sizeof(int));
		double** Ga = getWhrtfWithKnownElev(minElev, azim, ear, Ga_size);
		
		int* Gb_size = (int*) calloc((NUM_FILTROS+1), sizeof(int));
		double** Gb = getWhrtfWithKnownElev(maxElev, azim, ear, Gb_size);
		
		Gp = getGp(Ga, Ga_size, Gb, Gb_size, minElev, maxElev, elev, Gp_size);
	}
	return Gp;
}