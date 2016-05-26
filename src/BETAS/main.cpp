// Copyright 2011 Hamed Shateri Najafabadi

/********************************************************************

This file is part of COSMOS.

COSMOS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

COSMOS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with COSMOS.  If not, see <http://www.gnu.org/licenses/>.

********************************************************************/


// main.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "declarations.h"

extern const char *__metadata_file;
extern const char *__counts_file;
extern const int *__distribution;
extern const int *__learn;
extern const int *__burn;
extern const int *__samp;
extern const char *__output_file;

void welcome_message()
{
	cout << endl
		<< "**************************** DEStableE version 1.0 ****************************" << endl
		<< "*                                                                             *" << endl
		<< "* Copyright 2016 Hamed S. Najafabadi                                          *" << endl
		<< "*                                                                             *" << endl
		<< "*******************************************************************************" << endl
		<< endl;
}


int main( int argc, char* argv[] )
{	
	//******************* receive and analyze the arguments

	welcome_message();

	if( argc <= 1 )
	// no argument is profided
	{
		print_commandline_format();
		return 0;
	}
	
	if( !read_arguments( argc, argv ) )
		return 1;
		
	print_arguments();

	//******************* open output

	ofstream ofs;
	if( !open_output( ofs, __output_file, ".out.txt" ) )
		return 1;

	//******************* open tsv file

	s_array *arrays[ MAX_ARRAYS ];
	int num_arrays = 0;
	s_gene *genes[ MAX_GENES ];
	int num_genes = 0;

	cout << "Opening read count file..." << endl;
	if( !open_tsv( __counts_file, genes, &num_genes, arrays, &num_arrays ) )
		return 1;
	
	//******************* open metadata file

	s_group *groups[ MAX_GROUPS ];
	int num_groups = 0;
	
	cout << "Opening metadata file..." << endl;
	if( !open_metadata( __metadata_file, arrays, num_arrays, groups, &num_groups ) )
		return 1;
		

	//******************* run the algorithm

	cout << "Initializing MCMC..." << endl;
	Metropolis_Hastings( genes, num_genes, arrays, num_arrays, groups, num_groups,
		*__distribution, *__learn, *__burn, *__samp );


	//******************* write the output
	cout << "Writing the outputs...";
	write_estimates(
		ofs, genes, num_genes, groups, num_groups );
	
	cout << endl << "Job finished successfully." << endl; 

	return 0;
}
