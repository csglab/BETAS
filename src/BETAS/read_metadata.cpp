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
bool read_metadata(
	ifstream &ifs,
	s_array *arrays[],
	int num_arrays,
	s_group *groups[],
	int *num_groups,
	char this_EOL )
// returns true if successful, false otherwise
{
	*num_groups = 0;

	while( true )
	{
		char string[ MAX_LINE_LENGTH + 1 ];
		
		ifs.getline( string, MAX_LINE_LENGTH, this_EOL );
		if( !ifs || ifs.eof() )
			break;
			
		if( string[ 0 ] == '#' )
			continue; // this is a comment line
			
		//////////// extract the first column, corresponding to the sample ID
		char sample_id[ MAX_STRING_LENGTH + 1 ];
		if( !extract_phrase( string, 0, sample_id, _delimiters ) )
			continue; // the line is empty
	
		// find the array index in the list with this array name
		int array_index = find_object( sample_id, arrays, num_arrays );
		if( array_index >= num_arrays ) // this array was not found
		{
			cout << "WARNING: Sample " << sample_id << " was not found in the read count table." << endl;
			continue; // ignore this line
		}
		if( arrays[ array_index ] ->group_index >= 0 ) // this array was already assigned to a group
		{
			cout << "ERROR: Duplicate entry for sample " << sample_id << "." << endl;
			return false;
		}
			
		//////////// extract the second column, corresponding to the group ID
		char group_id[ MAX_STRING_LENGTH + 1 ];
		_PTRCHECK( extract_phrase( string, 1, group_id, _delimiters ) );
		
		// find the group index in the list with this group name
		int group_index = find_object( group_id, groups, *num_groups );
		if( group_index >= *num_groups ) // this is a new group
		{
			if( group_index >= MAX_GROUPS )
			{
				cout << "ERROR: Too many groups (maximum allowed: " << MAX_GROUPS << ")" << endl;
				return false;
			}
			
			groups[ group_index ] = new s_group; // create the new group
			_COPY_STR( groups[ group_index ] ->name, group_id ); // copy the name
			groups[ group_index ] ->index = group_index;
			*num_groups = group_index + 1; // update the group list size
		}
		arrays[ array_index ] ->group_index = group_index;

		//////////// extract the third column, corresponding to the read type
		char read_type[ MAX_STRING_LENGTH + 1 ];
		_PTRCHECK( extract_phrase( string, 2, read_type, _delimiters ) );
		if( strcmp( read_type, "intronic" ) == 0 )
		{
			arrays[ array_index ] ->ei = INTRONIC;
			groups[ group_index ] ->i = 1; // intronic reads are found for this group
		}
		else if( strcmp( read_type, "exonic" ) == 0 )
		{
			arrays[ array_index ] ->ei = EXONIC;
			groups[ group_index ] ->e = 1; // exonic reads are found for this group
		}
		else
		{
			cout << "ERROR: Unrecognized read type for sample " << sample_id << "." << endl;
			return false;
		}
	}
	
	// check if all the information is complete
	
	int i;
	for( i = 0; i < num_arrays; i ++ )
		if( arrays[ i ] ->group_index < 0 )
		{
			cout << "ERROR: Sample "<< arrays[ i ] ->name << " does not have metadata." << endl;
			return false;
		}

	for( i = 0; i < *num_groups; i ++ )
		if( groups[ i ] ->i <= 0 )
		{
			cout << "ERROR: Group "<< groups[ i ] ->name << " does not have intronic reads." << endl;
			return false;
		}
		else if( groups[ i ] ->e <= 0 )
		{
			cout << "ERROR: Group "<< groups[ i ] ->name << " does not have exonic reads." << endl;
			return false;
		}

	
	return true;
}
