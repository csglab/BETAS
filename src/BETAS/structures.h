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

#ifndef _H_STRUCTURES // this is to make sure that this file will not be included twice
//*******************************************************************
#define _H_STRUCTURES // this macro indicates that this file is now included

#ifndef _NAMESPACE
using namespace std;
#define _NAMESPACE
#endif


#include <string.h>
#include <stdlib.h>

// ************************************************ macro definitions

// constants
#define MAX_GENES				100000
#define MAX_ARRAYS				10000
#define MAX_GROUPS				256
#define MAX_LINE_LENGTH			100000
#define MAX_STRING_LENGTH		10000

#define INTRONIC				0
#define EXONIC					1

#define NORMAL					0
#define LAPLACE					1
#define LINEARGRADIENT			2
#define UNIFORM					3

// macros
#define _RELEASE(x)		do{ if(x) delete (x); }while(false) // this do-while loop is called only once
#define _CALL(x)		do{ if(!(x)) return false; }while(false)
#define _COPY_STR(x,y)	do{ x=new char[strlen(y)+1]; strcpy((x),(y)); }while(false)
#define MAX(x,y)		(((x)>(y))?(x):(y))
#define MIN(x,y)		(((x)<(y))?(x):(y))

// ********************************************** typedef definitions

typedef unsigned char BYTE;

// ******************************************** structure definitions

// this structure will hold the information for the array groups
struct s_group
{
	s_group()
	{
		name = NULL;
		index = -1;
		e = i = 0;
		
		ri_average = -0.33496165;
		ri_stdev = 1.13643655;
		re_average = -0.38766325;
		re_stdev = 1.1894084;
		rei_cor = 0.85447095;
		rei_var = 0.4202925;
		
	}
	~s_group()
	{
		_RELEASE( name );
	}
	
	char *name; // the name associated with this group (group ID)
	int index;
	
	int e; // whether exonic read counts are present for this group
	int i; // whether intronic read counts are present for this group
	
	double ri_average, ri_stdev, re_average, re_stdev, rei_cor, rei_var; // the parameters of the prior distribution
	// These variables will only be used if they are not learned from data. If present in the metadata file, they will be read from that file,
	// otherwise the default values that are set in s_group() will be used
};

// this structure will hold the information for the input arrays
struct s_array
{
	s_array()
	{
		name = NULL;
		index = -1;
		ei = -1;
		group_index = -1;
	}
	~s_array()
	{
		_RELEASE( name );
	}
	
	char *name; // the name associated with this array (sample ID)
	int index; // the index of this array in the s_array object
	int group_index; // the group this array belongs to
	int ei; // whether this array represents intronic read counts (0) or exonic read counts (1)	
};

// this structure will hold the information for the input genes
struct s_gene
{
	s_gene()
	{
		// initialize all the variables, setting them to zero
		
		name = NULL; // no memory allocated yet
		index = -1;
		
		profile = NULL;
		
		ni1 = ne1 = ni2 = ne2 = 0;

		log_fi_average = NULL;
		log_fi_stdev = NULL;
		log_fe_average = NULL;
		log_fe_stdev = NULL;

		log_dfi_average = NULL;
		log_dfi_stdev = NULL;
		log_dfe_average = NULL;
		log_dfe_stdev = NULL;

		log_ddf_average = NULL;
		log_ddf_stdev = NULL;

	}
	~s_gene()
	{
		// release all the allocated memory
	
		_RELEASE( name );
		_RELEASE( profile );

		_RELEASE( log_fi_average );
		_RELEASE( log_fi_stdev );
		_RELEASE( log_fe_average );
		_RELEASE( log_fe_stdev );

		_RELEASE( log_dfi_average );
		_RELEASE( log_dfi_stdev );
		_RELEASE( log_dfe_average );
		_RELEASE( log_dfe_stdev );

		_RELEASE( log_ddf_average );
		_RELEASE( log_ddf_stdev );

	}

	char *name; // the name of this gene, as read from the input FASTA file
	int index; // the index of this gene in the s_gene object
	
	int *profile; // the read count profile, as read from the input file; the number and order match that of s_array object


	// read counts for intronic (i) and exonic (e) reads, for the active group (2) and other groups (1)
	int ni1;
	int ne1;
	int ni2;
	int ne2;
	int min_n;
	double min_f;
	double var_log_f;
	
	double ri_prior;
	double re_prior;

	// MCMC results
	// the number and order of values in these arrays match that of s_group object
	
	double *log_fi_average; // the average of distribution of log(f_i) (fraction of intronic reads)
	double *log_fi_stdev; // the standard deviation of distribution of log(f_i)
	double *log_fe_average; // the average of distribution of log(f_e) (fraction of exonic reads)
	double *log_fe_stdev; // the standard deviation of distribution of log(f_e)

	double *log_dfi_average; // the average of distribution of log(f_i)-log(f_i_other_groups)
	double *log_dfi_stdev; // the standard deviation of distribution of log(f_i)-log(f_i_other_groups)
	double *log_dfe_average; // the average of distribution of log(f_e)-log(f_e_other_groups)
	double *log_dfe_stdev; // the standard deviation of distribution of log(f_e)-log(f_e_other_groups)

	double *log_ddf_average; // the average of distribution of log(delta_f_e)-log(delta_f_i)
	double *log_ddf_stdev; // the standard deviation of distribution of log(delta_f_e)-log(delta_f_i)
};

//*******************************************************************
#endif // this is to make sure that this file will not be included twice
