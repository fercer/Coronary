/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* FILE NAME: args_fercer.cpp                                                                                *
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
* Nov/2016    Fernando C.   2            2.1        Change the methods for a class ARGUMENTS                *
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
void ARGUMENTS::flushLine() {
	char c;
	while ((c = getchar()) != '\n' && c != EOF)
		/* discard */;
}





/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: cleanArgumentsTree                                                                         *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* my_root         BIN_TREE_ARGS       input      A root of a binary tree part to be freed recursively       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing, frees the arguments tree                                                                         *
*                                                                                                           *
************************************************************************************************************/
void ARGUMENTS::cleanArgumentsTree( BIN_TREE_ARGS *my_root )
{
	if (!my_root) {
		return;
	}

	/* Recursively free the memory allocated to store the arguments information */
	cleanArgumentsTree(my_root->my_lower_short);
	cleanArgumentsTree(my_root->my_upper_short);

	delete my_root;
}





/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: binarySearchLongTag                                                                        *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* searched_tag    const char*         input      The searched Long tag                                      *
* my_root         BIN_TREE_ARGS       input      A root of a binary tree part to search for the searched tag*
*                                                                                                           *
* RETURNS:                                                                                                  *
* The node where the argument is located.                                                                   *
*                                                                                                           *
************************************************************************************************************/
ARGUMENTS::BIN_TREE_ARGS* ARGUMENTS::binarySearchLongTag(const char *searched_tag, BIN_TREE_ARGS *my_root)
{
	const int comp_resp = strcmp(searched_tag, my_root->long_tag);

	if (comp_resp > 0) {/* my_root < searched tag: Continue searching in the upper side */
		if (my_root->my_upper_long) {
			return binarySearchLongTag(searched_tag, my_root->my_upper_long);
		}
		else {
			return my_root;
		}
	}
	else if (comp_resp < 0) {/* my_root > searched tag: Continue searching in the lower side */
		if (my_root->my_lower_long) {
			return binarySearchLongTag(searched_tag, my_root->my_lower_long);
		}
		else {
			return my_root;
		}
	}

	/* my_root == searched tag: Returns this position */
	return my_root;
}





/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: binarySearchLongTagExact                                                                   *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* searched_tag    const char*         input      The searched Long tag                                      *
* my_root         BIN_TREE_ARGS       input      A root of a binary tree part to search for the searched tag*
*                                                                                                           *
* RETURNS:                                                                                                  *
* The node where the argument is located.                                                                   *
*                                                                                                           *
************************************************************************************************************/
ARGUMENTS::BIN_TREE_ARGS* ARGUMENTS::binarySearchLongTagExact(const char *searched_tag, ARGUMENTS::BIN_TREE_ARGS *my_root)
{
	const int comp_resp = strcmp(searched_tag, my_root->long_tag);
	
	if (comp_resp > 0) {/* my_root < searched tag: Continue searching in the upper side */
		if (my_root->my_upper_long) {
			return binarySearchLongTagExact(searched_tag, my_root->my_upper_long);
		}
		else {
			/* Tag was not found */
			return NULL;
		}
	}
	else if (comp_resp < 0) {/* my_root > searched tag: Continue searching in the lower side */
		if (my_root->my_lower_long) {
			return binarySearchLongTagExact(searched_tag, my_root->my_lower_long);
		}
		else {
			/* Tag was not found */
			return NULL;
		}
	}

	/* my_root == searched tag: Returns this position */
	return my_root;
}





/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: binarySearchShortTag                                                                       *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* searched_tag    const char*         input      The searched Short tag                                     *
* my_root         BIN_TREE_ARGS       input      A root of a binary tree part to search for the searched tag*
*                                                                                                           *
* RETURNS:                                                                                                  *
* The node where the argument is located.                                                                   *
*                                                                                                           *
************************************************************************************************************/
ARGUMENTS::BIN_TREE_ARGS* ARGUMENTS::binarySearchShortTag(const char * searched_tag, ARGUMENTS::BIN_TREE_ARGS * my_root)
{
	const int comp_resp = strcmp(searched_tag, my_root->short_tag);

	if (comp_resp > 0) {/* my_root < searched tag: Continue searching in the upper side */
		if (my_root->my_upper_short) {
			return binarySearchShortTag(searched_tag, my_root->my_upper_short);
		}
		else {
			return my_root;
		}
	}
	else if (comp_resp < 0) {/* my_root > searched tag: Continue searching in the lower side */
		if (my_root->my_lower_short) {
			return binarySearchShortTag(searched_tag, my_root->my_lower_short);
		}
		else {
			return my_root;
		}
	}

	/* my_root == searched tag: Returns this position */
	return my_root;
}





/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: binarySearchShortTagExact                                                                  *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* searched_tag    const char*         input      The searched Short tag                                     *
* my_root         BIN_TREE_ARGS       input      A root of a binary tree part to search for the searched tag*
*                                                                                                           *
* RETURNS:                                                                                                  *
* The node where the argument is located.                                                                   *
*                                                                                                           *
************************************************************************************************************/
ARGUMENTS::BIN_TREE_ARGS* ARGUMENTS::binarySearchShortTagExact(const char * searched_tag, ARGUMENTS::BIN_TREE_ARGS * my_root)
{
	const int comp_resp = strcmp(searched_tag, my_root->short_tag);

	if (comp_resp > 0) {/* my_root < searched tag: Continue searching in the upper side */
		if (my_root->my_upper_short) {
			return binarySearchShortTagExact(searched_tag, my_root->my_upper_short);
		}
		else {
			/* Tag was not found */
			return NULL;
		}
	}
	else if (comp_resp < 0) {/* my_root > searched tag: Continue searching in the lower side */
		if (my_root->my_lower_short) {
			return binarySearchShortTagExact(searched_tag, my_root->my_lower_short);
		}
		else {
			/* Tag was not found */
			return NULL;
		}
	}

	/* my_root == searched tag: Returns this position */
	return my_root;
}





/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: ARGUMENTS                                                                                  *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing, creates and initializes a new ARGUMENTS object.                                                  *
*                                                                                                           *
************************************************************************************************************/
ARGUMENTS::ARGUMENTS()
{
	n_arguments = 0;
	my_args_root = new BIN_TREE_ARGS;
	my_args_root->my_lower_long = NULL;
	my_args_root->my_upper_long = NULL;
	my_args_root->my_lower_short = NULL;
	my_args_root->my_upper_short = NULL;

	my_args_root->is_optional = true;
	my_args_root->my_type = MY_HELP;
	
#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(my_args_root->question, 512 * sizeof(char), "Show input arguments necesary to execute this program");
	sprintf_s(my_args_root->long_tag, 32 * sizeof(char), "help");
	sprintf_s(my_args_root->short_tag, 8 * sizeof(char), "h");
	
	my_args_root->my_default_value.my_value_i = 0;
#else
	sprintf(my_args_root->question, "Show input arguments necesary to execute this program");
	sprintf(my_args_root->long_tag, "help");
	sprintf(my_args_root->short_tag, "h");
	
	my_args_root->my_default_value.my_value_i = 0;
#endif

}





/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: ARGUMENTS                                                                                  *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* input_argc      int                 input      Number of arguments passed to the program at execution     *
* input_argv      char**              input      Arguments passed tothe program at execution                *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing, creates and initializes a new ARGUMENTS object.                                                  *
*                                                                                                           *
************************************************************************************************************/
ARGUMENTS::ARGUMENTS(int input_argc, char **input_argv)
{
	n_arguments = 0;
	my_args_root = new BIN_TREE_ARGS;
	my_args_root->my_lower_long = NULL;
	my_args_root->my_upper_long = NULL;
	my_args_root->my_lower_short = NULL;
	my_args_root->my_upper_short = NULL;

	my_args_root->is_optional = true;
	my_args_root->my_type = MY_HELP;

#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(my_args_root->question, 512 * sizeof(char), "Show input arguments necesary to execute this program");
	sprintf_s(my_args_root->long_tag, 32 * sizeof(char), "help");
	sprintf_s(my_args_root->short_tag, 8 * sizeof(char), "h");

	my_args_root->my_default_value.my_value_i = 0;
#else
	sprintf(my_args_root->question, "Show input arguments necesary to execute this program");
	sprintf(my_args_root->long_tag, "help");
	sprintf(my_args_root->short_tag, "h");

	my_args_root->my_default_value.my_value_i = 0;
#endif

	my_argc = input_argc;
	my_argv = input_argv;
}





/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: ARGUMENTS                                                                                  *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing, frees the memory used to store this object members                                               *
*                                                                                                           *
************************************************************************************************************/
ARGUMENTS::~ARGUMENTS()
{
	cleanArgumentsTree(my_args_root);
}






/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: setArgc                                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* input_argc      int                 input      Number of arguments passed to the program at execution     *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing                                                                                                   *
*                                                                                                           *
************************************************************************************************************/
void ARGUMENTS::setArgc(const int input_argc) 
{
	my_argc = input_argc;
}






/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: setArgv                                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* input_argv      char**              input      Arguments passed tothe program at execution                *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing                                                                                                   *
*                                                                                                           *
************************************************************************************************************/
void ARGUMENTS::setArgv(char **input_argv)
{
	my_argv = input_argv;
}






/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: newArgument                                                                                *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* input_question         char*        input      Question to display if help is requested                   *
* input_tag              char*        input      Tag identifier for current argument                        *
* input_short_tag        char*        input      Short tag identifier for current argument                  *
* input_default_value    int          input      Default value for the current Integer argument             *
* input_is_optional      bool         input      If the argument is strictly necesary                       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing                                                                                                   *
*                                                                                                           *
************************************************************************************************************/
void ARGUMENTS::newArgument( const char *input_question, const char *input_short_tag, const char *input_tag,
	const int input_default_value, bool input_is_optional)
{
	BIN_TREE_ARGS *my_parent_short = binarySearchShortTag(input_short_tag, my_args_root);

	/* Check if the new argument goes to the lower or the upper direction of its parent: */
	const int my_comp_resp_short = strcmp(input_short_tag, my_parent_short->short_tag);

	BIN_TREE_ARGS *curr_arg;

	if (my_comp_resp_short > 0) {
		my_parent_short->my_upper_short = new BIN_TREE_ARGS;
		curr_arg = my_parent_short->my_upper_short;
	}
	else if (my_comp_resp_short < 0) {
		my_parent_short->my_lower_short = new BIN_TREE_ARGS;
		curr_arg = my_parent_short->my_lower_short;
	}
	else {
		/* If the argument already exists an error is messaged */
		std::cout << std::endl <<
			ANSI_COLOR_RED "<< Error: The new input already exist >>" ANSI_COLOR_RESET << std::endl;
		return;
	}

	curr_arg->my_lower_short = NULL;
	curr_arg->my_upper_short = NULL;

	n_arguments++;
	
#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(curr_arg->question, 512 * sizeof(char), "%s", input_question);
	sprintf_s(curr_arg->long_tag, 32 * sizeof(char), "%s", input_tag);
	sprintf_s(curr_arg->short_tag, 8 * sizeof(char), "%s", input_short_tag);
#else
	sprintf(curr_arg->question, "%s", input_question);
	sprintf(curr_arg->long_tag, "%s", input_tag);
	sprintf(curr_arg->short_tag, "%s", input_short_tag);
#endif

	curr_arg->my_type = MY_INT;
	curr_arg->my_default_value.my_value_i = input_default_value;
	curr_arg->is_optional = input_is_optional;
	curr_arg->was_defined = false;

	BIN_TREE_ARGS *my_parent_long = binarySearchLongTag(input_tag, my_args_root);
	const int my_comp_resp_long = strcmp(input_tag, my_parent_long->long_tag);

	if (my_comp_resp_long > 0) {
		my_parent_long->my_upper_long = curr_arg;
	}
	else if (my_comp_resp_long < 0) {
		my_parent_long->my_lower_long = curr_arg;
	}

	curr_arg->my_lower_long = NULL;
	curr_arg->my_upper_long = NULL;
}






/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: newArgument                                                                                *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* input_question         char*        input      Question to display if help is requested                   *
* input_tag              char*        input      Tag identifier for current argument                        *
* input_short_tag        char*        input      Short tag identifier for current argument                  *
* input_default_value    double       input      Default value for the current Double argument              *
* input_is_optional      bool         input      If the argument is strictly necesary                       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing                                                                                                   *
*                                                                                                           *
************************************************************************************************************/
void ARGUMENTS::newArgument(const char *input_question, const char *input_short_tag, const char *input_tag,
	const double input_default_value, bool input_is_optional)
{
	BIN_TREE_ARGS *my_parent_short = binarySearchShortTag(input_short_tag, my_args_root);

	/* Check if the new argument goes to the lower or the upper direction of its parent: */
	const int my_comp_resp_short = strcmp(input_short_tag, my_parent_short->short_tag);

	BIN_TREE_ARGS *curr_arg;

	if (my_comp_resp_short > 0) {
		my_parent_short->my_upper_short = new BIN_TREE_ARGS;
		curr_arg = my_parent_short->my_upper_short;
	}
	else if (my_comp_resp_short < 0) {
		my_parent_short->my_lower_short = new BIN_TREE_ARGS;
		curr_arg = my_parent_short->my_lower_short;
	}
	else {
		/* If the argument already exists an error is messaged */
		std::cout << std::endl <<
			ANSI_COLOR_RED "<< Error: The new input already exist >>" ANSI_COLOR_RESET << std::endl;
		return;
	}

	curr_arg->my_lower_short = NULL;
	curr_arg->my_upper_short = NULL;

	n_arguments++;

#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(curr_arg->question, 512 * sizeof(char), "%s", input_question);
	sprintf_s(curr_arg->long_tag, 32 * sizeof(char), "%s", input_tag);
	sprintf_s(curr_arg->short_tag, 8 * sizeof(char), "%s", input_short_tag);
#else
	sprintf(curr_arg->question, "%s", input_question);
	sprintf(curr_arg->long_tag, "%s", input_tag);
	sprintf(curr_arg->short_tag, "%s", input_short_tag);
#endif

	curr_arg->my_type = MY_DOUBLE;
	curr_arg->my_default_value.my_value_d = input_default_value;
	curr_arg->is_optional = input_is_optional;
	curr_arg->was_defined = false;

	BIN_TREE_ARGS *my_parent_long = binarySearchLongTag(input_tag, my_args_root);
	const int my_comp_resp_long = strcmp(input_tag, my_parent_long->long_tag);

	if (my_comp_resp_long > 0) {
		my_parent_long->my_upper_long = curr_arg;
	}
	else if (my_comp_resp_long < 0) {
		my_parent_long->my_lower_long = curr_arg;
	}

	curr_arg->my_lower_long = NULL;
	curr_arg->my_upper_long = NULL;

}






/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: newArgument                                                                                *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* input_question         char*        input      Question to display if help is requested                   *
* input_tag              char*        input      Tag identifier for current argument                        *
* input_short_tag        char*        input      Short tag identifier for current argument                  *
* input_default_value    char*        input      Default value for the current String argument              *
* input_is_optional      bool         input      If the argument is strictly necesary                       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing                                                                                                   *
*                                                                                                           *
************************************************************************************************************/
void ARGUMENTS::newArgument(const char *input_question, const char *input_short_tag, const char *input_tag,
	const char *input_default_value, bool input_is_optional)
{
	BIN_TREE_ARGS *my_parent_short = binarySearchShortTag(input_short_tag, my_args_root);

	/* Check if the new argument goes to the lower or the upper direction of its parent: */
	const int my_comp_resp_short = strcmp(input_short_tag, my_parent_short->short_tag);

	BIN_TREE_ARGS *curr_arg;

	if (my_comp_resp_short > 0) {
		my_parent_short->my_upper_short = new BIN_TREE_ARGS;
		curr_arg = my_parent_short->my_upper_short;
	}
	else if (my_comp_resp_short < 0) {
		my_parent_short->my_lower_short = new BIN_TREE_ARGS;
		curr_arg = my_parent_short->my_lower_short;
	}
	else {
		/* If the argument already exists an error is messaged */
		std::cout << std::endl <<
			ANSI_COLOR_RED "<< Error: The new input already exist >>" ANSI_COLOR_RESET << std::endl;
		return;
	}

	curr_arg->my_lower_short = NULL;
	curr_arg->my_upper_short = NULL;

	n_arguments++;

#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(curr_arg->question, 512 * sizeof(char), "%s", input_question);
	sprintf_s(curr_arg->long_tag, 32 * sizeof(char), "%s", input_tag);
	sprintf_s(curr_arg->short_tag, 8 * sizeof(char), "%s", input_short_tag);
	sprintf_s(curr_arg->my_default_value.my_value_s, 512 * sizeof(char), "%s", input_default_value);
#else
	sprintf(curr_arg->question, "%s", input_question);
	sprintf(curr_arg->long_tag, "%s", input_tag);
	sprintf(curr_arg->short_tag, "%s", input_short_tag);
	sprintf(curr_arg->my_default_value.my_value_s, "%s", input_default_value);
#endif

	curr_arg->my_type = MY_CHAR;
	curr_arg->is_optional = input_is_optional;
	curr_arg->was_defined = false;

	BIN_TREE_ARGS *my_parent_long = binarySearchLongTag(input_tag, my_args_root);
	const int my_comp_resp_long = strcmp(input_tag, my_parent_long->long_tag);

	if (my_comp_resp_long > 0) {
		my_parent_long->my_upper_long = curr_arg;
	}
	else if (my_comp_resp_long < 0) {
		my_parent_long->my_lower_long = curr_arg;
	}

	curr_arg->my_lower_long = NULL;
	curr_arg->my_upper_long = NULL;
}






/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: askForValue                                                                                *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* searched_input  INPUT_ARGS*         input      The argument checked                                       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Stores the passed argument in the searched argument structure                                             *
*                                                                                                           *
************************************************************************************************************/
void ARGUMENTS::askForValue(BIN_TREE_ARGS *searched_input)
{
	std::cout << std::endl << searched_input->question << ANSI_COLOR_BLUE " (default: ";
	switch (searched_input->my_type) {
		case MY_CHAR:
			std::cout << searched_input->my_default_value.my_value_s;
			break;
		case MY_INT:
			std::cout << searched_input->my_default_value.my_value_i;
			break;
		case MY_DOUBLE:
			std::cout << searched_input->my_default_value.my_value_d;
			break;
	}

	std::cout << ")" << ANSI_COLOR_RESET std::endl;
	std::cout << "To left the default value, press " ANSI_COLOR_YELLOW "[Enter]" ANSI_COLOR_RESET
		", otherwise, please define the argument value and press " ANSI_COLOR_YELLOW "[Enter]: " ANSI_COLOR_RESET << std::endl;

	char input_value[513] = "";
	char *resp;
	resp = gets_s(input_value, 512*sizeof(char));
	flushLine();
	
	if ( (int)(*input_value) == 0 ){ /* Default selected: */
		switch (searched_input->my_type) {
		case MY_CHAR:
#if defined(_WIN32) || defined(_WIN64)
			sprintf_s(searched_input->my_value.my_value_s, 512 * sizeof(char), "%s",
				searched_input->my_default_value.my_value_s);
#else
			sprintf(searched_input->my_value.my_value_s, "%s", searched_input->my_default_value.my_value_s);
#endif
			break;
		case MY_INT:
			searched_input->my_value.my_value_i = searched_input->my_default_value.my_value_i;
			break;
		case MY_DOUBLE:
			searched_input->my_value.my_value_d = searched_input->my_default_value.my_value_d;
			break;
		}
	}
	else {
		switch (searched_input->my_type) {
		case MY_CHAR:
#if defined(_WIN32) || defined(_WIN64)
			sprintf_s(searched_input->my_value.my_value_s, 512 * sizeof(char), "%s", input_value);
#else
			sprintf(searched_input->my_value.my_value_s, "%s", input_value);
#endif
			break;
		case MY_INT:
			searched_input->my_value.my_value_i = atoi(input_value);
			break;
		case MY_DOUBLE:
			searched_input->my_value.my_value_d = atof(input_value);
			break;
		}
	}
}







/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: checkIfPassed                                                                              *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* searched_input  INPUT_ARGS*         input      The argument checked                                       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Stores the passed argument in the searched argument structure                                             *
*                                                                                                           *
************************************************************************************************************/
void ARGUMENTS::checkIfPassed(BIN_TREE_ARGS *searched_input) {
	/* If the parameter was not defined yet, search wether it was passed by the user: */
	int arg_idx = 1;

	while (arg_idx < my_argc) {
		if ((strcmp(*(my_argv + arg_idx) + 1, searched_input->short_tag) == 0) ||
			(strcmp(*(my_argv + arg_idx) + 2, searched_input->long_tag)  == 0)) {
			break;
		}
		arg_idx++;
	}
	
	if (arg_idx == my_argc) {
		/* The argument was not passed by the user: */
		if (searched_input->is_optional) {
			/* If the argument is optional, copy the default value: */
			switch (searched_input->my_type) {
				case MY_CHAR:
	#if defined(_WIN32) || defined(_WIN64)
					sprintf_s(searched_input->my_value.my_value_s, 512 * sizeof(char), "%s",
						searched_input->my_default_value.my_value_s);
	#else
					sprintf(searched_input->my_value.my_value_s, "%s", searched_input->my_default_value.my_value_s);
	#endif
					break;
				case MY_HELP:
				case MY_INT:
					searched_input->my_value.my_value_i = searched_input->my_default_value.my_value_i;
					break;
				case MY_DOUBLE:
					searched_input->my_value.my_value_d = searched_input->my_default_value.my_value_d;
					break;
			}
		}
		else {
			/* If the argument is not optional, ask for its value: */
			askForValue(searched_input);
		}
	}
	else {
		/* The argument was passed by the user: */
		switch (searched_input->my_type) {
			case MY_CHAR:
#if defined(_WIN32) || defined(_WIN64)
				sprintf_s(searched_input->my_value.my_value_s, 512 * sizeof(char),
					"%s", *(my_argv + arg_idx + 1));
#else
				sprintf(searched_input->my_value.my_value_s, "%s", *(my_argv + arg_idx + 1));
#endif
				break;
			case MY_HELP:
				searched_input->my_value.my_value_i = 1;
				break;
			case MY_INT:
				searched_input->my_value.my_value_i = atoi(*(my_argv + arg_idx + 1));
				break;
			case MY_DOUBLE:
				searched_input->my_value.my_value_d = atof(*(my_argv + arg_idx + 1));
				break;
		}
	}

	searched_input->was_defined = true;
}




/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: getArgumentTree                                                                            *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* search_tag      const char*         input      The tag (short or long) that is looked for                 *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Return the pointer to the searched tree                                                                   *
*                                                                                                           *
************************************************************************************************************/
ARGUMENTS::BIN_TREE_ARGS* ARGUMENTS::getArgumentTree(const char *searched_tag) {
	
	if (*(searched_tag + 1) == '-') {
		/* The tag is long */
		return binarySearchLongTagExact(searched_tag + 2, my_args_root);
	}
	else {
		/* The tag is short */
		return binarySearchShortTagExact(searched_tag + 1, my_args_root);
	}
}




/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: getArgumentINT                                                                             *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* search_tag      const char*         input      The tag (short or long) that is looked for                 *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Returns the Integer argument passed by the user at the execution time                                     *
*                                                                                                           *
************************************************************************************************************/
int ARGUMENTS::getArgumentINT(const char *searched_tag) 
{
	BIN_TREE_ARGS *found_arg = getArgumentTree(searched_tag);

	if (!found_arg) {
		std::cout << ANSI_COLOR_RED "<< Error: The tag " << searched_tag << " was not defined >>" ANSI_COLOR_RESET << std::endl;
		return -1;
	}

	if ( found_arg->was_defined == false ) {
		checkIfPassed(found_arg);
	}

	return found_arg->my_value.my_value_i;
}




/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: getArgumentCHAR                                                                            *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* search_tag      const char*         input      The tag (short or long) that is looked for                 *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Returns the String argument passed by the user at the execution time                                      *
*                                                                                                           *
************************************************************************************************************/
char* ARGUMENTS::getArgumentCHAR(const char *searched_tag)
{
	BIN_TREE_ARGS *found_arg = getArgumentTree(searched_tag);
	
	if (!found_arg) {
		std::cout << ANSI_COLOR_RED "<< Error: The tag " << searched_tag << " was not defined >>" ANSI_COLOR_RESET << std::endl;
		return NULL;
	}

	if (!found_arg->was_defined) {
		checkIfPassed(found_arg);
	}
	
	return found_arg->my_value.my_value_s;
}





/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: getArgumentDOUBLE                                                                          *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* search_tag      const char*         input      The tag (short or long) that is looked for                 *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Returns the Double argument passed by the user at the execution time                                      *
*                                                                                                           *
************************************************************************************************************/
double ARGUMENTS::getArgumentDOUBLE(const char *searched_tag)
{
	BIN_TREE_ARGS *found_arg = getArgumentTree(searched_tag);

	if (!found_arg) {
		std::cout << ANSI_COLOR_RED "<< Error: The tag " << searched_tag << " was not defined >>" ANSI_COLOR_RESET << std::endl;
		return NULL;
	}

	if (!found_arg->was_defined) {
		checkIfPassed(found_arg);
	}

	return found_arg->my_value.my_value_d;
}





/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: setArgumentINT                                                                             *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* search_tag      const char*         input      The tag (short or long) that is looked for                 *
* input_value     const int           input      The value to assign to the selected tag                    *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing                                                                                                   *
*                                                                                                           *
************************************************************************************************************/
void ARGUMENTS::setArgumentINT(const char * searched_tag, const int input_value)
{
	BIN_TREE_ARGS *found_arg = getArgumentTree(searched_tag);

	if (!found_arg) {
		std::cout << ANSI_COLOR_RED "<< Error: The tag " << searched_tag << " was not defined >>" ANSI_COLOR_RESET << std::endl;
		return;
	}

	found_arg->was_defined = true;

	found_arg->my_value.my_value_i = input_value;
}





/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: setArgumentINT                                                                             *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* search_tag      const char*         input      The tag (short or long) that is looked for                 *
* input_value     const double        input      The value to assign to the selected tag                    *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing                                                                                                   *
*                                                                                                           *
************************************************************************************************************/
void ARGUMENTS::setArgumentCHAR(const char * searched_tag, const char *input_value)
{
	BIN_TREE_ARGS *found_arg = getArgumentTree(searched_tag);

	if (!found_arg) {
		std::cout << ANSI_COLOR_RED "<< Error: The tag " << searched_tag << " was not defined >>" ANSI_COLOR_RESET << std::endl;
		return;
	}

	found_arg->was_defined = true;

#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(found_arg->my_value.my_value_s, 512 * sizeof(char), "%s", input_value);
#else
	sprintf(found_arg->my_value.my_value_s, "%s", input_value);
#endif
}





/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: setArgumentINT                                                                             *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* search_tag      const char*         input      The tag (short or long) that is looked for                 *
* input_value     const double        input      The value to assign to the selected tag                    *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing                                                                                                   *
*                                                                                                           *
************************************************************************************************************/
void ARGUMENTS::setArgumentDOUBLE(const char * searched_tag, const double input_value)
{
	BIN_TREE_ARGS *found_arg = getArgumentTree(searched_tag);

	if (!found_arg) {
		std::cout << ANSI_COLOR_RED "<< Error: The tag " << searched_tag << " was not defined >>" ANSI_COLOR_RESET << std::endl;
		return;
	}

	found_arg->was_defined = true;

	found_arg->my_value.my_value_d = input_value;
}





/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: displayArguments                                                                           *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* my_root         BIN_TREE_ARGS*      input      The current root of an arguments tree part                 *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Displays the current argument information                                                                 *
*                                                                                                           *
************************************************************************************************************/
void ARGUMENTS::displayArgument(ARGUMENTS::BIN_TREE_ARGS* my_root) {

#if defined(_WIN32) || defined(_WIN64)
	if (my_root->is_optional) {
		std::cout << "\t{-" << my_root->short_tag << "\t--"
			<< my_root->long_tag << "}\t\t\t";
	}
	else {
		std::cout << "\t[-" << my_root->short_tag << "\t--"
			<< my_root->long_tag << "]\t\t\t";
	}
	std::cout << my_root->question;
#else
	if (my_root->is_optional) {
		std::cout << ANSI_COLOR_YELLOW "\t";
	}
	else {
		std::cout << ANSI_COLOR_MAGENTA "\t";
	}

	std::cout << "-" << my_root->short_tag << "\t--"
		<< my_root->long_tag << "\t\t\t"
		<< ANSI_COLOR_GREEN my_root->my_question ANSI_COLOR_RESET;
#endif

	std::cout << " (default: ";

	switch (my_root->my_type) {
	case MY_CHAR:
		std::cout << ANSI_COLOR_CYAN my_root->my_default_value.my_value_s;
		break;
	case MY_HELP:
		std::cout << ANSI_COLOR_CYAN "no";
		break;
	case MY_INT:
		std::cout << ANSI_COLOR_CYAN my_root->my_default_value.my_value_i;
		break;
	case MY_DOUBLE:
		std::cout << ANSI_COLOR_CYAN my_root->my_default_value.my_value_d;
		break;
	}
	std::cout << ANSI_COLOR_RESET ")" << std::endl;

	if (my_root->my_lower_short) {
		displayArgument(my_root->my_lower_short);
	}

	if (my_root->my_upper_short) {
		displayArgument(my_root->my_upper_short);
	}
}






/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: showHelp                                                                                   *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Displays the requested arguments                                                                          *
*                                                                                                           *
************************************************************************************************************/
void ARGUMENTS::showHelp() {

	/* Check if the help was requested by the user: */
	checkIfPassed(my_args_root);

	if (my_args_root->my_value.my_value_i == 0) {
		return;
	}

	/* Show usage of the program: */
#if defined(_WIN32) || defined(_WIN64)
	std::cout << "Usage: " << *(my_argv) << " [OPTION [obligatory]/{optional}] Value" << std::endl;
#else
	std::cout << "Usage: " ANSI_COLOR_GREEN << *(my_argv) << ANSI_COLOR_RESET " [OPTION " ANSI_COLOR_MAGENTA "obligatory" ANSI_COLOR_RESET "/" ANSI_COLOR_YELLOW "optional" ANSI_COLOR_RESET "] Value" << std::endl;
#endif

	displayArgument(my_args_root);

	std::cout << ANSI_COLOR_RESET std::endl << std::endl;
}