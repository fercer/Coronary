/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* FILE NAME: nondom_fercer.c                                                                                *
*                                                                                                           *
* PURPOSE: Request arguments for the program execution and usage display.                                   *
*                                                                                                           *
* FILE REFERENCES:                                                                                          *
* Name        I/O        Description                                                                        *
* None        ----       ----------                                                                         *
*                                                                                                           *
* ABNORMAL TERMINATION CONDITIONS, ERROR AND WARNING MESSAGES:                                              *
* None                                                                                                      *
*                                                                                                           *
* NOTES:                                                                                                    *
*                                                                                                           *
* DEVELOPMENT HISTORY:                                                                                      *
* Date        Author        Change Id    Release    Description Of Change                                   *
* Jun/2015    Fernando C.   0            1.0        Creation                                                *
* Nov/2016    Fernando C.   1            2.0        Modification of coding standard and languaje            *
*                                                   revisar_pars  ---> checkArgs                            *
*                                                   mostrar_ayuda ---> showHelp                             *
*                                                                                                           *
************************************************************************************************************/

#include "args_fercer.h"


/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: checkArgs                                                                                  *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* my_args         INPUT_ARGS*         input      The definition of each requested argument                  *
* n_args          const int           input      The number of requested arguments                          *
* argc            int*                input      The number of arguments passed                             *
* argv            char**              input      The arguments passed                                       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* 1 if the display of the usage was requested                                                               *
* 0 otherwise                                                                                               *
*                                                                                                           *
************************************************************************************************************/
int checkArgs(INPUT_ARGS* my_args, int *argc, char **argv){
	const int n_args = my_args->my_default_value.my_value_i;

    int was_help_requested = 0;

	// Para quitar el nombre con el que se ejecuto el programa:
	*argc = *argc - 1;
    
	int *retired_args = NULL;
    
    if( *argc > 0 ){
        retired_args = (int *) calloc(*argc, sizeof(int));
	}
	else {
		was_help_requested = 1;
	}

	for (int i = 0; i < n_args; i++){
		int n = 1;
		int was_found = 0;

        if ( retired_args ){
            do{
                if(!retired_args[n-1]){
                    if(!strcmp(my_args[i].short_tag, argv[n]) || !strcmp(my_args[i].long_tag, argv[n])){ 
						/* The question tag was found: */
                        retired_args[n-1] = 1;
                        was_found = 1;
                        switch (my_args[i].my_type){
                            case MY_CHAR:
                                    sprintf(my_args[i].my_value.my_value_s, "%s", argv[n+1]);
                                    break;

                            case MY_INT:
                                    my_args[i].my_value.my_value_i = atoi(argv[n+1]);
                                    break;

                            case MY_DOUBLE:
                                    my_args[i].my_value.my_value_d = atof(argv[n+1]);
                                    break;

                            case MY_HELP:
                                    was_help_requested = n_args;
                                    break;
                        }
                    }
                }
                n++;
            }while (!was_found && n <= *argc);
		}

		/* The question tag was not found */
		if (!was_found){
			/* Check if the argument is optional: */
			if (my_args[i].is_optional){
				/* If it is optional, the default value is assigned directly: */
				switch (my_args[i].my_type){
					case MY_CHAR:
						sprintf(my_args[i].my_value.my_value_s, "%s", my_args[i].my_default_value.my_value_s);
						break;

					case MY_INT:
						my_args[i].my_value.my_value_i = my_args[i].my_default_value.my_value_i;
						break;

					case MY_DOUBLE:
						my_args[i].my_value.my_value_d = my_args[i].my_default_value.my_value_d;
						break;
				}
			}else if( was_help_requested == 1 ){
				/* Otherwise, the parameter value is requested to the user if the help was not requested: */
				printf("\n %s " ANSI_COLOR_BLUE "(default: ", my_args[i].question);
                
                char *resp;
                
				switch (my_args[i].my_type){
					case MY_CHAR:
						printf("\'%s\')" ANSI_COLOR_RESET "\n", my_args[i].my_default_value.my_value_s);
						printf("\n To left the default value, press " ANSI_COLOR_YELLOW "[Enter]" ANSI_COLOR_RESET ",\
							otherwise, please define the argument value and press " ANSI_COLOR_YELLOW "[Enter]: "\
							ANSI_COLOR_RESET "\n");

						resp = fgets(my_args[i].my_value.my_value_s, 128, stdin);

						if ((int)my_args[i].my_value.my_value_s[0]==10){ /* Default selected: */
							sprintf(my_args[i].my_value.my_value_s, "%s", my_args[i].my_default_value.my_value_s);
						}else{ /* Other value defined: */
							my_args[i].my_value.my_value_s[strlen(my_args[i].my_value.my_value_s)-1] = '\0';
						}
						was_help_requested--;
						break;

					case MY_INT:
						printf("%i)" ANSI_COLOR_RESET "\n", my_args[i].my_default_value.my_value_i);
						printf("\n To left the default value, press " ANSI_COLOR_YELLOW "[Enter]" ANSI_COLOR_RESET "\
							otherwise, please define the argument value and press " ANSI_COLOR_YELLOW "[Enter]: "\
							ANSI_COLOR_RESET "\n");
						char tmp_i[128];
						resp = fgets(tmp_i, 128, stdin);
						if ((int)tmp_i[0]==10){ /* Default selected: */
							my_args[i].my_value.my_value_i = my_args[i].my_default_value.my_value_i;
						}else{ /* Other value defined: */
							my_args[i].my_value.my_value_i = atoi(tmp_i);
						}
						was_help_requested--;
						break;

					case MY_DOUBLE:
						printf("%f)" ANSI_COLOR_RESET "\n", my_args[i].my_default_value.my_value_d);
						printf("\n To left the default value, press " ANSI_COLOR_YELLOW "[Enter]" ANSI_COLOR_RESET ",\
							otherwise, please define the argument value and press " ANSI_COLOR_YELLOW "[Enter]: " \
							ANSI_COLOR_RESET "\n");
						char tmp_d[128];
						resp = fgets(tmp_d, 128, stdin);
						if ((int)tmp_d[0]==10){ /* Default selected: */
							my_args[i].my_value.my_value_d = my_args[i].my_default_value.my_value_d;
						}else{ /* Other value defined: */
							my_args[i].my_value.my_value_d = atof(tmp_d);
						}
						was_help_requested--;
						break;
				}
			}
		}
	}

	if( retired_args ){
        free(retired_args);
    }
    
    if( was_help_requested > 0 ){
        showHelp(my_args, n_args, argv[0]);
    }

    return was_help_requested;
}






/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: showHelp                                                                                   *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* my_args         INPUT_ARGS*         input      The definition of each requested argument                  *
* n_args          const int           input      The number of requested arguments                          *
* program_name    const char*         input      The name of the program at the time it was executed        *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Displays the requested arguments                                                                          *
*                                                                                                           *
************************************************************************************************************/
void showHelp(INPUT_ARGS* my_args, const char *program_name){
	const int n_args = my_args->my_default_value.my_value_i;
    /* Show usage of the program: */

#if defined(_WIN32) || defined(_WIN64)
	printf("Usage: %s [OPTION [obligatory]/{optional}] Value\n", program_name);
#else
    printf( "Usage: " ANSI_COLOR_GREEN "%s" ANSI_COLOR_RESET " [OPTION " ANSI_COLOR_MAGENTA "obligatory"\
		ANSI_COLOR_RESET "/" ANSI_COLOR_YELLOW "optional" ANSI_COLOR_RESET "] Value\n", program_name);
#endif
    for( int i = 0; i < n_args; i++){

#if defined(_WIN32) || defined(_WIN64)
		if (my_args[i].is_optional) {
			printf("\t{%s %s}\t\t\t", my_args[i].short_tag, my_args[i].long_tag);
		}
		else {
			printf("\t[%s %s]\t\t\t", my_args[i].short_tag, my_args[i].long_tag);
		}
		printf("%s (", my_args[i].question);
#else
        if( my_args[i].is_optional ){
            printf( ANSI_COLOR_YELLOW "\t");
        }else{
            printf( ANSI_COLOR_MAGENTA "\t");
        }
        printf("%s %s\t\t\t" ANSI_COLOR_GREEN "%s " ANSI_COLOR_RESET "(", my_args[i].short_tag, my_args[i].long_tag, my_args[i].question);
#endif

        switch( my_args[i].my_type ){
            case MY_CHAR:
                printf(ANSI_COLOR_CYAN "%s", my_args[i].my_default_value.my_value_s );
                break;
            case MY_INT:
                printf(ANSI_COLOR_CYAN "%i", my_args[i].my_default_value.my_value_i );
                break;
            case MY_DOUBLE:
                printf(ANSI_COLOR_CYAN "%f", my_args[i].my_default_value.my_value_d );
                break;
        }
        printf(ANSI_COLOR_RESET ")\n");
    }
    printf(ANSI_COLOR_RESET "\n\n");
}
