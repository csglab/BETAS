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


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "declarations.h"

///////////////////////////////////////////////////////////////////////////////////////////
void print_progress_range()
{
	int i;
	for( i = 0; i < 100; i ++ )
		cout << "_";
	cout << endl;
}

void update_progress_report(
	int *prev_percent,
	double total_count,
	double *progress )
{
	// update the progress report
	(*progress) ++;
	int percent = (*progress) / total_count * 100;
	for( ; (*prev_percent) < percent; (*prev_percent) ++ )
	{
		cout << "^"; cout.flush();
	}
}

////////////////////////////////////////////////////////////////////////////////////
double generate_normal_random(
	double mean,
	double stdev )
// This function returns a random number generated based on a normal distribution
{
	static double spare = 0; // the spare number genreated previously, but not used
	static bool spare_ready = false; // indicates whether a spare number exists
	
	if( spare_ready ) // a spare number exists
	{
		spare_ready = false; // unflag the spare_ready indicator
		// return the random number, scaled to the requested distribution
		return spare * stdev + mean;
	}
	
	double u, v, s;
	
	// choose two random variables u and v from uniform distributino [-1,1], so that
	// 0 < u*u + v*v < 1
	do
	{
		// generate two random variables from uniform distributino [-1,1]
		u = drand48() * 2 - 1;
		v = drand48() * 2 - 1;
		
		s = u*u + v*v;

	} while( s >= 1 || s <= 0 );

	// store the spare random number for the next call	
	spare = v * sqrt( -2.0 * log(s) / s );
	spare_ready = true;
	
	// return the random number
	return stdev * u * sqrt( -2.0 * log(s) / s ) + mean;
}

///////////////////////////////////////////////////////////////////////////////////////////
bool _is_delimiter(
	char ch,
	const char *delimiters,
	int num_delimiters )
// determines whether a character is a delimiter or not
// returns true if ch is a delimiter, false otherwise
//
// ch: the query character
// delimiters: the list of delimiters
// num_delimiters: the size of the delimiter list
{
	// examine the delimiters one by one to find a match
	for( int i = 0; i < num_delimiters; i ++ )
		if( ch == delimiters[ i ] )
			return true; // a matching delimiter is found

	return false;
}


///////////////////////////////////////////////////////////////////////////////////////////
char* extract_phrase(
	const char *line,
	int phrase_index,
	char *phrase,
	const char *delimiters,
	int *guide_index,
	int *guide_pos )
// extracts a phrase from a delimitted text
// returns NULL if unsuccessful, 'phrase' otherwise
//
// line: the line from which the phrase will be extracted
// phrase_index: the index of the phrase that will be extracted; 0 means that first phrase, and so on
// phrase: the string in which the extracted phrase will be stored
// delimiters: the list of valid delimiters, separating different phrases within the text
// guide_index: if this parameter is set, then curr_index will start from this, and
//  start_pos will start from guide_pos
//  In other words, these guide variables mean that the 'guide_index'th phrase is found when start_pos=guide_pos
{
	int len = strlen( line );
	int num_delimiters = strlen( delimiters );

	// find the first position after 'phrase_index'th occurrance of a delimiter
	int curr_index = 0;
	int start_pos = 0;
	if( guide_index && // guide_index is set
		guide_pos && // guide_pos is set
		*guide_index >= 0 && // guide_index is valid
		*guide_pos >= 0 && // guide_pos is valid
		*guide_index <= phrase_index ) // guide_index is smaller than the query index
	// therefore, it is not necessary to start from the beginning of the line
	{
		curr_index = *guide_index;
		start_pos = *guide_pos;
	}
	for( ; start_pos < len && curr_index < phrase_index; start_pos ++ )
		if( _is_delimiter( line[ start_pos ], delimiters, num_delimiters ) )
			curr_index ++;
			
	if( guide_index && guide_pos ) // guide_index and guide_pos are requested
	{
		*guide_index = curr_index;
		*guide_pos = start_pos;
	}

	// return NULL if the requested phrase is not found
	if( start_pos >= len )
		return NULL;

	// extract the phrase
	int pos;
	for( pos = start_pos; pos < len; pos ++ )
		if( _is_delimiter( line[ pos ], delimiters, num_delimiters ) )
			break;
		else
			phrase[ pos - start_pos ] = line[ pos ];

	phrase[ pos - start_pos ] = 0;

	// return the extracted phrase
	return phrase;
}
