// Copyright 2010 Hamed Shateri Najafabadi

/********************************************************************

This file is part of HyperMotif.

HyperMotif is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

HyperMotif is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with HyperMotif.  If not, see <http://www.gnu.org/licenses/>.

********************************************************************/


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

using namespace std;

#include "structures.h"
#include "declarations.h"

#define _EOFCHECK(x)	if( !x || x.eof() ){ cout << "ERROR: Unexpected end of file." << endl << endl; return false; }
#define _PTRCHECK(x)	if( !x ){ cout << "ERROR: Unexpected end of line." << endl; return false; }

static char _delimiters[] = { char(9), 0 };

//////////////////////////////////////////////////////////////////////////////////////////
int _count_chr( char *string, char ch )
{
	int count = 0;
	int len = strlen( string );
	
	int i;
	for( i = 0; i < len; i ++ )
		if( string[ i ] == ch )
			count ++;
		
	return count;
}

//////////////////////////////////////////////////////////////////////////////////////////
bool read_tsv(
	ifstream &ifs,
	s_gene *genes[],
	int *num_genes,
	s_array *arrays[],
	int *num_arrays,
	char this_EOL )
// returns true if successful, false otherwise
{
	*num_genes = 0;
	*num_arrays = 0;
	
	// read the header, containing array names
	char string[ MAX_LINE_LENGTH + 1 ];
	ifs.getline( string, MAX_LINE_LENGTH, this_EOL );
	_EOFCHECK( ifs );
	
	// determine the number of arrays based on the number of delimiting characters in the header
	*num_arrays = _count_chr( string, _delimiters[ 0 ] );
	if( *num_arrays > MAX_ARRAYS )
	{
		cout << "ERROR: Too many columns (maximum allowed: " << MAX_ARRAYS << ")" << endl;
		return false;
	}
	cout << *num_arrays << " columns are being read..." << endl;
	
	// extract the array names
	int i;
	for( i = 0; i < *num_arrays; i ++ )
	{
		char array_name[ MAX_STRING_LENGTH + 1 ];
		_PTRCHECK( extract_phrase( string, i + 1, array_name, _delimiters ) );
		
		arrays[ i ] = new s_array;
		arrays[ i ] ->index = i;
		
		_COPY_STR( arrays[ i ] ->name, array_name );
	}	
	
	// read the matrix
	while( true )
	{
		ifs.getline( string, MAX_LINE_LENGTH, this_EOL ); // read the line
		if( !ifs || ifs.eof() )
			break; // this is the end of file
			
		if( *num_genes >= MAX_GENES )
		{
			cout << "ERROR: Too many genes (maximum allowed: " << MAX_GENES << ")" << endl;
			return false;
		}
			
		i = *num_genes; // the index of the gene to be read
		
		genes[ i ] = new s_gene;
		genes[ i ] ->index = i;
		
		// extract and copy the gene name
		char gene_name[ MAX_STRING_LENGTH + 1 ];
		_PTRCHECK( extract_phrase( string, 0, gene_name, _delimiters ) );		
		_COPY_STR( genes[ i ] ->name, gene_name );
		
		// allocate memory to the profile
		genes[ i ] ->profile = new double[ *num_arrays ];
				
		// extract and store the profile
		int j;
		int guide_index = -1;
		int guide_pos = -1;
		for( j = 0; j < *num_arrays; j ++ )
		{
			char value[ MAX_STRING_LENGTH + 1 ];
			_PTRCHECK( extract_phrase( string, j + 1, value, _delimiters,
				&guide_index, &guide_pos ) );
				
			if( strcmp( value, "na" ) == 0 ||
				strcmp( value, "NA" ) == 0 ||
				strcmp( value, "nan" ) == 0 ||
				strcmp( value, "NAN" ) == 0 ) // this is not a number
				genes[ i ] ->profile[ j ] = 0;
			else
				genes[ i ] ->profile[ j ] = atof( value ); // this is a number
		}
		
		// update the number of genes
		(*num_genes) ++;
	}
	
	return true;
}
