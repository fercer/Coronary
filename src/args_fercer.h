/******
	CENTRO DE INVESTIGACION EN MATEMATICAS
	MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

	FERNANDO CERVANTES SANCHEZ
	JUN - 2015
*****/

#ifndef ARGS_FERCER_H_INCLUDED
#define ARGS_FERCER_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/***********************************************************************************************************/
#if defined(_WIN32) || defined(_WIN64)
	#define ANSI_COLOR_RED
	#define ANSI_COLOR_GREEN
	#define ANSI_COLOR_YELLOW
	#define ANSI_COLOR_BLUE
	#define ANSI_COLOR_MAGENTA
	#define ANSI_COLOR_CYAN
	#define ANSI_COLOR_RESET
#else
	#define ANSI_COLOR_RED     "\x1b[31m"
	#define ANSI_COLOR_GREEN   "\x1b[32m"
	#define ANSI_COLOR_YELLOW  "\x1b[33m"
	#define ANSI_COLOR_BLUE    "\x1b[34m"
	#define ANSI_COLOR_MAGENTA "\x1b[35m"
	#define ANSI_COLOR_CYAN    "\x1b[36m"
	#define ANSI_COLOR_RESET   "\x1b[0m"
#endif

/***********************************************************************************************************/
typedef enum { MY_CHAR, MY_INT, MY_DOUBLE, MY_HELP } INPUT_TYPE;

typedef union {
	char my_value_s[128];
	int my_value_i;
	double my_value_d;
} INPUT_VAL;


typedef struct {
    char question[512];
    char short_tag[7], long_tag[28];
	INPUT_TYPE my_type;
	INPUT_VAL my_default_value;
	INPUT_VAL my_value;
    char is_optional;
} INPUT_ARGS;


/***********************************************************************************************************/
int checkArgs(INPUT_ARGS* my_args, int *argc, char **argv);
void showHelp(INPUT_ARGS* my_args, const char *program_name);

#endif //ARGS_FERCER_H_INCLUDED
