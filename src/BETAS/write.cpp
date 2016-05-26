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

void write_estimates(
	ofstream &ofs,
	s_gene *genes[],
	int num_genes,
	s_group *groups[],
	int num_groups )
{
	// write header line 1
	ofs << "Group";
	int i, j;
	for( i = 0; i < num_groups; i ++ )
		for( j = 0; j < 10; j ++ )
			ofs << char(9) << groups[ i ] ->name;
	ofs << endl;
	
	// write header line 2
	ofs << "Parameter";
	for( i = 0; i < num_groups; i ++ )
		ofs << char(9) << "fi_average" << char(9) << "fi_stdev"
			<< char(9) << "fe_average" << char(9) << "fe_stdev"
			<< char(9) << "dfi_average" << char(9) << "dfi_stdev"
			<< char(9) << "dfe_average" << char(9) << "dfe_stdev"
			<< char(9) << "ddf_average" << char(9) << "ddf_stdev";
	ofs << endl;
	
	// print the data
	
	for( i = 0; i < num_genes; i ++ )
	{
		ofs << genes[ i ] ->name;
		
		for( j = 0; j < num_groups; j ++ )
			ofs << char(9) << genes[ i ] ->log_fi_average[ j ] << char(9) << genes[ i ] ->log_fi_stdev[ j ]
				<< char(9) << genes[ i ] ->log_fe_average[ j ] << char(9) << genes[ i ] ->log_fe_stdev[ j ]
				<< char(9) << genes[ i ] ->log_dfi_average[ j ] << char(9) << genes[ i ] ->log_dfi_stdev[ j ]
				<< char(9) << genes[ i ] ->log_dfe_average[ j ] << char(9) << genes[ i ] ->log_dfe_stdev[ j ]
				<< char(9) << genes[ i ] ->log_ddf_average[ j ] << char(9) << genes[ i ] ->log_ddf_stdev[ j ];
		ofs << endl;
	}
}
