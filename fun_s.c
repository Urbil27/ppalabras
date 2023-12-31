/***************************************************
 AC - OpenMP -- SERIE
 fun_s.c
 rutinas que se utilizan en el modulo grupopal_s.c
****************************************************/
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include "defineg.h" // definiciones


/*******************************************************************
 1 - Funcion para calcular la distancia euclidea entre dos vectores
 Entrada: 2 elementos con NDIM caracteristicas (por referencia)
 Salida:  distancia (double)
********************************************************************/
double gendist (float *vec1, float *vec2){

	double distancia = 0.0;
	for(int i =0;i<NDIM;i++){
		distancia += pow((vec1[i]-vec2[i]),2);
	}
	return sqrt(distancia);
}

/***********************************************************************************
 2 - Funcion para calcular el grupo (cluster) mas cercano (centroide mas cercano)
 Entrada: nvec  numero de vectores, int
          mvec  vectores, una matriz de tamanno MAXV x NDIM, por referencia
          cent  centroides, una matriz de tamanno ngrupos x NDIM, por referencia
 Salida:  popul grupo mas cercano a cada elemento, vector de tamanno MAXV, por ref.
************************************************************************************/
void grupo_cercano (int nvec, float mvec[][NDIM], float cent[][NDIM],
		int *popul){
	for(int i=0;i<nvec;i++){
		double distanciaMinima=DBL_MAX;
		int grupoCercano;
		for(int j=0;j<ngrupos;j++){
			double distanciaActual = gendist(mvec[i],cent[j]);
			if(distanciaMinima>distanciaActual){
				distanciaMinima = distanciaActual;
				grupoCercano=j;
			}
		}
		popul[i]=grupoCercano;
	}
}


/***************************************************************************************
 3 - Funcion para calcular la calidad de la particion de clusteres.
     Ratio entre a y b.
     El termino a corresponde a la distancia intra-cluster.
     El termino b corresponde a la distancia inter-cluster.
 Entrada: mvec    vectores, una matriz de tamanno MAXV x NDIM, por referencia
          listag  vector de ngrupos structs (informacion de grupos generados), por ref.
          cent    centroides, una matriz de tamanno ngrupos x NDIM, por referencia
 Salida:  valor del CVI (double): calidad/bondad de la particion de clusters
****************************************************************************************/

double silhouette_simple(float mvec[][NDIM], struct lista_grupos *listag, float cent[][NDIM],
		float a[]){
		double b[ngrupos];
		double s[ngrupos];
		double CVI=0.0;
		
		//Distancia intra-cluster
		for(int i=0;i<ngrupos;i++){
			double distancia=0.0;
			int tamanio=listag[i].nvecg;
			if(tamanio>1){
				for(int j=0;j<tamanio;j++){
					for(int k=j+1;k<tamanio;k++){
						distancia += gendist(mvec[listag[i].vecg[j]],mvec[listag[i].vecg[k]]);
					}
				}
				a[i]=distancia/(tamanio*(tamanio-1)/2);
			}
			else{
				a[i] = 0.0;
			}	
		}
		//Distancia inter-cluster
		for(int i=0;i<ngrupos;i++){
			double distancia=0.0;
			int cont=0;
			for(int j=0;j<ngrupos;j++){
				if(i != j){
					distancia += gendist(cent[i],cent[j]);
					cont++;
				}
			}
			b[i] = distancia/cont;
		}

		for(int i=0;i<ngrupos;i++){
			s[i]=(b[i]-a[i])/(fmax(a[i],b[i]));
		}

		for(int i=0;i<ngrupos;i++){
			CVI += s[i];
		}
		return CVI/ngrupos;
}


//Algoritmo que ordena un vector por inserción
void ordenar_vector(double *vector,int tamanio){
	int pos;
	double aux;
	for(int i=1;i<tamanio;i++){
		pos=i;
		aux=vector[i];

		while(pos>0 && aux<vector[pos-1]){
			  vector[pos]=vector[pos-1];
			  pos--;
		}
		vector[pos]=aux;
	}
}


/********************************************************************************************
 4 - Funcion para relizar el analisis de campos UNESCO
 Entrada:  listag   vector de ngrupos structs (informacion de grupos generados), por ref.
           mcam     campos, una matriz de tamaño MAXV x NCAM, por referencia
 Salida:   info_cam vector de NCAM structs (informacion del analisis realizado), por ref.
*****************************************************************************************/
void analisis_campos(struct lista_grupos *listag, float mcam[][NCAM],
		struct analisis *info_cam){
	double* datos;
	for(int i=0;i<NCAM;i++){
		info_cam[i].mmax=-1;
		info_cam[i].mmin=2;
		info_cam[i].gmax=0;
		info_cam[i].gmin=0;

		for(int j = 0;j<ngrupos;j++){
			if(listag[j].nvecg>0){
				int tamanio=listag[j].nvecg;
				datos = malloc(tamanio * sizeof(double));
				for(int k=0; k<tamanio;k++){
					datos[k]=mcam[k][i];
				} 
				ordenar_vector(datos,tamanio);
				int index = tamanio/2;
				double mediana = datos[index];
				if(mediana>info_cam[i].mmax){
					info_cam[i].mmax = mediana;
					info_cam[i].gmax = j;
				}
				else if(mediana<info_cam[i].mmin){
					info_cam[i].mmin= mediana;
					info_cam[i].gmin= j;
				}
				free(datos);
			}
		}
	}
}


/*************************************
   OTRAS FUNCIONES DE LA APLICACION
**************************************/
void inicializar_centroides(float cent[][NDIM]){
	int i, j;
	float rand_val;
	srand (147);
	for (i=0; i<ngrupos; i++)
		for (j=0; j<NDIM/2; j++){
			rand_val = ((rand() % 10000) / 10000.0)*2-1;
			cent[i][j] = rand_val;
			cent[i][j+(NDIM/2)] = cent[i][j];
		}
}

int nuevos_centroides(float mvec[][NDIM], float cent[][NDIM], int popul[], int nvec){
	int i, j, fin;
	double discent;
	double additions[ngrupos][NDIM+1];
	float newcent[ngrupos][NDIM];

	for (i=0; i<ngrupos; i++)
		for (j=0; j<NDIM+1; j++)
			additions[i][j] = 0.0;

	// acumular los valores de cada caracteristica; numero de elementos al final
	for (i=0; i<nvec; i++){
		for (j=0; j<NDIM; j++) additions[popul[i]][j] += mvec[i][j];
		additions[popul[i]][NDIM]++;
	}

	// calcular los nuevos centroides y decidir si el proceso ha finalizado o no (en funcion de DELTA)
	fin = 1;
	for (i=0; i<ngrupos; i++){
		if (additions[i][NDIM] > 0) { // ese grupo (cluster) no esta vacio
			// media de cada caracteristica
			for (j=0; j<NDIM; j++)
				newcent[i][j] = (float)(additions[i][j] / additions[i][NDIM]);

			// decidir si el proceso ha finalizado
			discent = gendist (&newcent[i][0], &cent[i][0]);
			if (discent > DELTA1){
				fin = 0;  // en alguna centroide hay cambios; continuar
			}

			// copiar los nuevos centroides
			for (j=0; j<NDIM; j++)
				cent[i][j] = newcent[i][j];
		}

	}
	return fin;
}
