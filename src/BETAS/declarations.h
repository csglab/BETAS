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

#ifndef _H_DECLARATIONS // this is to make sure that this file will not be included twice
//*******************************************************************
#define _H_DECLARATIONS // this macro indicates that this file is now included


#include "structures.h" // this header file is required for definition of structures


bool read_arguments( int argc, char* argv[] );
void print_arguments();
void print_commandline_format();

void print_progress_range();
void update_progress_report( int *prev_percent, double total_count, double *progress );
double generate_normal_random( double mean, double stdev );
char* extract_phrase( const char *line, int phrase_index, char *phrase, const char *delimiters,	int *guide_index = NULL, int *guide_pos = NULL );

bool open_output( ofstream &ofs, const char *path, const char *extension );
bool open_tsv( const char *file, s_gene *genes[], int *num_genes, s_array *arrays[], int *num_arrays );
bool open_metadata( const char *file, s_array *arrays[], int num_arrays, s_group *groups[], int *num_groups );

bool read_tsv( ifstream &ifs, s_gene *genes[], int *num_genes, s_array *arrays[], int *num_arrays, char this_EOL );
bool read_metadata( ifstream &ifs, s_array *arrays[], int num_arrays, s_group *groups[], int *num_groups, char this_EOL );

int optimize_cor( s_gene *genes[], int num_genes, int Ni1, int Ne1, int Ni2, int Ne2, double *ri_average, double *ri_stdev, double *re_average, double *re_stdev, double *rei_cor, double *rei_var );
double sample_f( double old_fa1, double fa2, double fb1, double fb2, double ra_average, double ra_stdev, double rb_average, double rb_stdev, double rab_cor, double rab_var, int na1, int Na1, int distribution );
void Metropolis_Hastings( s_gene *genes[], int num_genes, s_array *arrays[], int num_arrays, s_group *groups[], int num_groups, int distribution, int learn, int burn_in, int sampling );

void write_estimates( ofstream &ofs, s_gene *genes[], int num_genes, s_group *groups[], int num_groups );


//////////////////////// Functions that need to be complied on-demand /////////////////
template <typename s_type>
int find_object(
	const char *name,
	s_type *list[],
	int list_size )
// returns num_arrays if the requested array is not found
{
	int i;
	for( i = 0; i < list_size; i ++ )
		if( strcmp( name, list[ i ] ->name ) == 0 )
			break;

	return i;
}

//*******************************************************************
#endif // this is to make sure that this file will not be included twice
