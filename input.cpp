#include <iostream>
#include <string>
#include <cstring>
#include <getopt.h>

#include "defs.h"

static struct option long_options[] =
{
	{ "solid-string",			required_argument,	NULL,	'r' },
	{ "weighted-string",		required_argument,	NULL,	'l' },
	{ "output-file",			required_argument,	NULL,	'o' },
	{ "cumulative-threshold",	required_argument,	NULL,	'z'	},
	{ "help",					0,					NULL,	'h' },
};

int decode_switches ( int argc, char * argv[], struct TSwitch * sw )
{
	int opt;
	char *ep;
	double val;
	int args_counter;

	/* initialisation */
	sw -> z						=	1;

	args_counter = 0;

	while ( ( opt = getopt_long ( argc, argv, "r:l:o:z:h", long_options, NULL ) ) != -1 )
	{
		switch ( opt )
		{
			case 'r':
				sw -> right_filename= optarg;
				args_counter ++;
				break;
			case 'l':
				sw -> left_filename = optarg;
				args_counter ++;
				break;
			case 'o':
				sw -> output_filename = optarg;
				args_counter ++;
				break;
			case 'z':
				val = strtod ( optarg, &ep );
				if ( optarg == ep )
				{
					return 0;
				}
				sw -> z = val;
				args_counter ++;
				break;
			case 'h':
				return 0;
		}
	}

	if ( args_counter != 4 )
	{
		usage();
		exit ( 1 );
	}
	else
		return ( optind );
}

void usage ( void )
{
	cout << "Please provide appropriate arguments:"  << endl;
	cout << "	-l, --weighted-string\t<str>\tFilename for Left FASTQ String." << endl;
	cout << "	-r, --solid-string\t<str>\tFilename for Right FASTQ String." << endl;
	cout << "	-o, --output-file\t<str>\tFilename for result output." << endl;
	cout << "	-z, --cumulative-threshold\t<dbl>\tCumulative weight threshold."<<endl;
	cout << "	-h, --help\t<dbl>\tHelp!"<<endl;
}






