/******
	CENTRO DE INVESTIGACION EN MATEMATICAS
	MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

	FERNANDO CERVANTES SANCHEZ
	JUN - 2015
*****/

#include "pars_fercer.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>




/*	FUNCION:	revisar_pars

	DESCRIPCION:
				Busca entre los parametros de entrada dados en la ejecucion del programa, la lista de parametros que se requierern.
*/
void revisar_pars(PARS_ENTRADA* pars, const int n_pars, int *argc, char **argv){
	// Para quitar el nombre con el que se ejecuto el programa:
	*argc = *argc - 1;

	int retirados[*argc];
	memset(retirados, 0, *argc * sizeof(int));
	
	for(int i = 0; i < n_pars; i++){

		int n = 1;
		int encontrado = 0;
		do{
			if(!retirados[n-1]){
				if(!strcmp(pars[i].tag, argv[n])){	// El tag de la pregunta si se encontro:
					retirados[n-1] = 1;
					encontrado = 1;
					switch(pars[i].mi_tipo){
						case CHAR:
								sprintf(pars[i].mi_valor.par_s, "%s", argv[n+1]);
								break;

						case INT:
								pars[i].mi_valor.par_i = atoi(argv[n+1]);
								break;

						case DOUBLE:
								pars[i].mi_valor.par_d = atof(argv[n+1]);
								break;
					}
				}
			}
			n++;
		}while(!encontrado && n <= *argc);
		// No se encontro el tag de la pregunta:
		if(!encontrado){
			// se verifica si la pregunta es opcional:
			if(pars[i].opcional){
				// Si lo es, se copia el valor por default directamente:
				switch(pars[i].mi_tipo){
					case CHAR:
						sprintf(pars[i].mi_valor.par_s, "%s", pars[i].mi_default.par_s);
						break;

					case INT:
						pars[i].mi_valor.par_i = pars[i].mi_default.par_i;
						break;

					case DOUBLE:
						pars[i].mi_valor.par_d = pars[i].mi_default.par_d;
						break;
				}
			}else{
				// De lo contrario se pregunta al usuario por el parametro:
				printf("\n %s " ANSI_COLOR_BLUE "(default: ", pars[i].pregunta);
				switch(pars[i].mi_tipo){
					case CHAR:
							printf("\'%s\')" ANSI_COLOR_RESET "\n", pars[i].mi_default.par_s);
							printf("\n Para dejarlo en su valor por default presione " ANSI_COLOR_YELLOW "[Enter]" ANSI_COLOR_RESET ", de lo contrario, defina el valor y presione " ANSI_COLOR_YELLOW "[Enter]: " ANSI_COLOR_RESET "\n");	
							fgets(pars[i].mi_valor.par_s, 128, stdin);
							if((int)pars[i].mi_valor.par_s[0]==10){ //Default seleccionado:
								sprintf(pars[i].mi_valor.par_s, "%s", pars[i].mi_default.par_s);
							}else{ // Otro valor definido:
								pars[i].mi_valor.par_s[strlen(pars[i].mi_valor.par_s)-1] = '\0';
							}

							break;

					case INT:
							printf("%i )" ANSI_COLOR_RESET "\n", pars[i].mi_default.par_i);
							printf("\n Para dejarlo en su valor por default presione " ANSI_COLOR_YELLOW "[Enter]" ANSI_COLOR_RESET ", de lo contrario, defina el valor y presione " ANSI_COLOR_YELLOW "[Enter]: " ANSI_COLOR_RESET "\n");
							char tmp_i[128];
							fgets(tmp_i, 128, stdin);
							if((int)tmp_i[0]==10){ //Default seleccionado:
								pars[i].mi_valor.par_i = pars[i].mi_default.par_i;
							}else{
								pars[i].mi_valor.par_i = atoi(tmp_i);
							}
			
							break;

					case DOUBLE:
							printf("%f )" ANSI_COLOR_RESET "\n", pars[i].mi_default.par_d);
							printf("\n Para dejarlo en su valor por default presione " ANSI_COLOR_YELLOW "[Enter]" ANSI_COLOR_RESET ", de lo contrario, defina el valor y presione " ANSI_COLOR_YELLOW "[Enter]: " ANSI_COLOR_RESET "\n");
							char tmp_d[128];
							fgets(tmp_d, 128, stdin);
							if((int)tmp_d[0]==10){ //Default seleccionado:
								pars[i].mi_valor.par_d = pars[i].mi_default.par_d;
							}else{
								pars[i].mi_valor.par_d = atof(tmp_d);
							}

							break;
				}
			}
		}
	}
}






