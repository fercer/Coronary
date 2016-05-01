/******
	CENTRO DE INVESTIGACION EN MATEMATICAS
	MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

	FERNANDO CERVANTES SANCHEZ
	JUN - 2015
*****/

#ifndef PARS_FERCER_H_INCLUDED
#define PARS_FERCER_H_INCLUDED


#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA //"\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"




typedef enum{ CHAR, INT, DOUBLE } TIPO_ENTRADA;
typedef union{
	char par_s[128];
	int par_i;
	double par_d;
} ENTRADA;

typedef struct{
	char pregunta[128];
    char short_tag[7], long_tag[28];
	TIPO_ENTRADA mi_tipo;
	ENTRADA mi_default;
	ENTRADA mi_valor;
    char opcional;
}PARS_ENTRADA;

void revisar_pars(PARS_ENTRADA* pars, const int n_pars, int *argc, char **argv);
void mostrar_ayuda(PARS_ENTRADA* pars, const int n_pars, const char *nombre_programa);

#endif //PARS_FERCER_H_INCLUDED
