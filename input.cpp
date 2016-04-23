/*
 *  Author: Liu, C., (2015), WPT (Version 2.2)
 *  Available from: https://github.com/YagaoLiu/WPT
 *  Adapted by: Yordan Petrov Yordanov
 */

#include <iostream>
#include <string>
#include <cstring>
#include <getopt.h>

#include "defs.h"

static struct option long_options[] =
{
	{ "right-string",				required_argument,	NULL,	'r' },
	{ "left-string",				required_argument,	NULL,	'l' },
	{ "output-file",				required_argument,	NULL,	'o' },
	{ "cumulative-threshold",		required_argument,	NULL,	'z'	},
	{ "exclude-read-threshold",		required_argument,	NULL,	'x'	},
	{ "bigram-window",				required_argument,	NULL,	'w'	},
	{ "quality-score-coefficient",	required_argument,	NULL,	'q'	},
	{ "print",						required_argument,	NULL,	'p'	},
	{ "help",						0,					NULL,	'h' },
};

int decode_switches ( int argc, char * argv[], struct TSwitch * sw )
{
	int opt;
	char *ep;
	double val;
	int args_counter;

	/* initialisation */
	sw -> z							=	10;
	sw -> x							=	10;
	sw -> bigramWindow				=	1000;
	sw -> qualityScoreCoefficient	=	100;
	sw -> print						=	0;

	args_counter = 0;

	while ( ( opt = getopt_long ( argc, argv, "r:l:o:z:x:w:q:p:h", long_options, NULL ) ) != -1 )
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
			case 'x':
				val = strtod ( optarg, &ep );
				if ( optarg == ep )
				{
					return 0;
				}
				sw -> x = val;
				args_counter ++;
				break;
			case 'w':
				val = strtod ( optarg, &ep );
				if ( optarg == ep )
				{
					return 0;
				}
				sw -> bigramWindow = val;
				args_counter ++;
				break;
			case 'q':
				val = strtod ( optarg, &ep );
				if ( optarg == ep )
				{
					return 0;
				}
				sw -> qualityScoreCoefficient = val;
				args_counter ++;
				break;
			case 'p':
				val = strtod ( optarg, &ep );
				if ( optarg == ep )
				{
					return 0;
				}
				sw -> print = val;
				args_counter ++;
				break;
			case 'h':
				return 0;
		}
	}

	if ( args_counter < 3 )
	{
		usage();
		exit ( 1 );
	}
	else
		return ( optind );
}

void usage ( void )
{
	cout << "\nPlease provide compulsory  arguments:"  << endl;
	cout << "	-l, --left-string\t\t<str>\tFilename for Left FASTQ String." << endl;
	cout << "	-r, --right-string\t\t<str>\tFilename for Right FASTQ String." << endl;
	cout << "	-o, --output-file\t\t<str>\tFilename for result output." << endl;
	cout << "\nOptional arguments:"  << endl;
	cout << "	-z, --cumulative-threshold\t<dbl>\tBy default set to 10. Cumulative weight threshold. The bigger the threshold, the higher the prefix table construction tolerance."<<endl;
	cout << "	-x, --exclude-read-threshold\t<dbl>\tBy default set to 10. Exclude reads with insignificant overlap. Higher exculsion percentage generates fewer reads."<<endl;
	cout << "	-q, --quality-score-coefficient\t<int>\tBy default set to 100. Coefficent of quality score : bigram probability distribution. To disregard quality score probability distribution set to 0."<<endl;
	cout << "	-w, --bigram-window\t\t<int>\tBy default set to 1000. Length of bigram window. To disregard bigram probability distribution set to 0."<<endl;
	cout << "	-p, --print\t\t\t<bool>\tBy default set to False. Printout all probability matrixes, prefix tables and border tables. The print functionality will dramatically decrease the performance therefore, it is used only for very small data sets."<<endl;
	cout << "	-h, --help\t\t\tHelp!"<<endl;
}






