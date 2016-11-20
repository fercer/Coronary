/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* FILE NAME: args_fercer.h                                                                                  *
*                                                                                                           *
* DEVELOPMENT HISTORY:                                                                                      *
* Date        Author        Change Id    Release    Description Of Change                                   *
* Jun/2015    Fernando C.   0            1.0        Creation                                                *
* Nov/2016    Fernando C.   1            2.0        Modification of coding standard and languaje            *
*                                                   revisar_pars  ---> checkArgs                            *
*                                                   mostrar_ayuda ---> showHelp                             *
* Nov/2016    Fernando C.   2            2.1        Change the methods for a class ARGUMENTS                *
*                                                                                                           *
************************************************************************************************************/

#ifndef ARGS_FERCER_H_INCLUDED
#define ARGS_FERCER_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <vector>

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

/************************************************************************************************************
*                                                                                                           *
* CLASS NAME: ARGUMENTS                                                                                     *
*                                                                                                           *
************************************************************************************************************/
class ARGUMENTS {

	private:
		/******* STRUCTURES, ENUMS AND UNIONS ****************************************************************/
		typedef enum INPUT_TYPE  { MY_CHAR, MY_INT, MY_DOUBLE, MY_HELP } INPUT_TYPE;

		typedef union INPUT_VAL {
			char my_value_s[513];
			int my_value_i;
			double my_value_d;
		} INPUT_VAL;

		typedef struct BIN_TREE_ARGS {
			BIN_TREE_ARGS *my_upper_long;
			BIN_TREE_ARGS *my_lower_long;

			BIN_TREE_ARGS *my_upper_short;
			BIN_TREE_ARGS *my_lower_short;

			char question[513];
			char short_tag[9], long_tag[33];
			INPUT_TYPE my_type;
			INPUT_VAL my_default_value;
			INPUT_VAL my_value;
			bool is_optional;
			bool was_defined;

		} BIN_TREE_ARGS;


		/******* PRIVATE MEMBERS *****************************************************************************/
		char **my_argv;
		int my_argc;

		int n_arguments;
		BIN_TREE_ARGS *my_args_root;


		/****** PRIVATE METHODS*******************************************************************************/
		void flushLine();
		void cleanArgumentsTree(BIN_TREE_ARGS *my_root);
		BIN_TREE_ARGS *binarySearchLongTag(const char * searched_tag, BIN_TREE_ARGS *my_root);
		BIN_TREE_ARGS *binarySearchLongTagExact(const char * searched_tag, BIN_TREE_ARGS * my_root);

		BIN_TREE_ARGS *binarySearchShortTag(const char * searched_tag, BIN_TREE_ARGS *my_root);
		BIN_TREE_ARGS *binarySearchShortTagExact(const char * searched_tag, BIN_TREE_ARGS * my_root);

		void askForValue(BIN_TREE_ARGS *searched_input);
		void checkIfPassed(BIN_TREE_ARGS *searched_input);
		BIN_TREE_ARGS * getArgumentTree(const char *searched_tag);

		void displayArgument(BIN_TREE_ARGS * my_root);


	public:
		/******* CONSTRUCTORS AND DESTRUCTORS ****************************************************************/
		ARGUMENTS();
		ARGUMENTS(int input_argc, char ** input_argv);
		~ARGUMENTS();

		/******* PUBLIC METHODS ******************************************************************************/
		void setArgc(const int input_argc);
		void setArgv(char ** input_argv);

		void newArgument(
			const char *input_question,
			const char *input_short_tag,
			const char *input_tag,
			const int input_default_value,
			bool input_is_optional
		);

		void newArgument(
			const char * input_question,
			const char * input_short_tag,
			const char * input_tag,
			const double input_default_value,
			bool input_is_optional
		);

		void newArgument(
			const char * input_question,
			const char * input_short_tag,
			const char * input_tag,
			const char * input_default_value,
			bool input_is_optional
		);

		int getArgumentINT(const char *searched_tag);
		char* getArgumentCHAR(const char * searched_tag);
		double getArgumentDOUBLE(const char * searched_tag);

		void setArgumentINT(const char *searched_tag, const int input_value);
		void setArgumentCHAR(const char * searched_tag, const char *input_value);
		void setArgumentDOUBLE(const char * searched_tag, const double input_value);



		void showHelp();
};

#endif //ARGS_FERCER_H_INCLUDED
