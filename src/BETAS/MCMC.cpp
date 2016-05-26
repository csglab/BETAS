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

extern int __MCMC_accepted;

		
////////////////////////////////////////////////////////////////////////////////////
void _Metropolis_Hastings_per_gene(
	double ri_average,
	double ri_stdev,
	double re_average,
	double re_stdev,
	double rei_cor,
	double rei_var,
	int ni1,
	int ne1,
	int ni2,
	int ne2,
	int Ni1,
	int Ne1,
	int Ni2,
	int Ne2,
	int distribution,
	int burn_in,
	int sampling,
	double *log_fi1_average,
	double *log_fi1_stdev,
	double *log_fi2_average,
	double *log_fi2_stdev,
	double *log_fe1_average,
	double *log_fe1_stdev,
	double *log_fe2_average,
	double *log_fe2_stdev,
	double *log_fi21_average,
	double *log_fi21_stdev,
	double *log_fe21_average,
	double *log_fe21_stdev,
	double *log_fei21_average,
	double *log_fei21_stdev )
// approximates the unknown parameters for distribution of fi1, fi2, fe1, and fe2
// See the document MCMC.docx for details
// returns the probability that the f_t <= f_0
//
// ri_average: the prior mean of log(f_i2)-log(f_i1)
// ri_stdev: the prior standard deviation of log(f_i2)-log(f_i1)
// re_average: the prior mean of log(f_e2)-log(f_e1)
// re_stdev: the prior standard deviation of log(f_e2)-log(f_e1)
// rei_cor: the prior correlation of log(f_e2)-log(f_e1) vs. log(f_i2)-log(f_i1)
// rei_var: the prior variance of (log(f_e2)-log(f_e1))-(log(f_i2)-log(f_i1))
// ni1: the number of intronic reads in group 1 for this gene
// ne1: the number of exonic reads in group 1 for this gene
// ni2: the number of intronic reads in group 2 for this gene
// ne2: the number of exonic reads in group 2 for this gene
// Ni1: the number of intronic reads in group 1 for all genes
// Ne1: the number of exonic reads in group 1 for all genes
// Ni2: the number of intronic reads in group 2 for all genes
// Ne2: the number of exonic reads in group 2 for all genes
// distribution: NORMAL or LAPLACE
//
// burn_in: the number of iterations for converging of the distribution
// num_iterations: the number of iterations for sampling the distribution
//
// log_fi1_average: the average of sampled distribution for log(f_i1)
// log_fi1_stdev: the standard deviation of sampled distribution for log(f_i1)
// log_fi2_average: the average of sampled distribution for log(f_i2)
// log_fi2_stdev: the standard deviation of sampled distribution for log(f_i2)
// log_fe1_average: the average of sampled distribution for log(f_e1)
// log_fe1_stdev: the standard deviation of sampled distribution for log(f_e1)
// log_fe2_average: the average of sampled distribution for log(f_e2)
// log_fe2_stdev: the standard deviation of sampled distribution for log(f_e2)
// log_fi21_average: the average of sampled distribution for log(f_i2)-log(f_i1)
// log_fi21_stdev: the standard deviation of sampled distribution for log(f_i2)-log(f_i1)
// log_fe21_average: the average of sampled distribution for log(f_e2)-log(f_e1)
// log_fe21_stdev: the standard deviation of sampled distribution for log(f_e2)-log(f_e1)
// log_fei21_average: the average of sampled distribution for [log(f_e2)-log(f_e1)] - [log(f_i2)-log(f_i1)]
// log_fei21_stdev: the standard deviation of sampled distribution for [log(f_e2)-log(f_e1)] - [log(f_i2)-log(f_i1)]

{
	*log_fi1_average = *log_fi1_stdev = *log_fi2_average = *log_fi2_stdev =
	*log_fe1_average = *log_fe1_stdev = *log_fe2_average = *log_fe2_stdev =
	*log_fi21_average = *log_fi21_stdev = *log_fe21_average = *log_fe21_stdev =
	*log_fei21_average = *log_fei21_stdev = 0;
	
	// initialize the unknown parameters
	double fi1 = double(ni1+1)/double(Ni1+2);
	double fi2 = double(ni2+1)/double(Ni2+2);
	double fe1 = double(ne1+1)/double(Ne1+2);
	double fe2 = double(ne2+1)/double(Ne2+2);
	
	sampling += burn_in; // the total number of sampling
	
	int n = 0;
	// perform the MCMC
	int i;
	for( i = 0; i < sampling; i ++ )
	{
		fi1 = sample_f( fi1, fi2, fe1, fe2, ri_average, ri_stdev, re_average, re_stdev, rei_cor, rei_var, ni1, Ni1, distribution );
		fi2 = sample_f( fi2, fi1, fe2, fe1, -ri_average, ri_stdev, -re_average, re_stdev, rei_cor, rei_var, ni2, Ni2, distribution );
		fe1 = sample_f( fe1, fe2, fi1, fi2, re_average, re_stdev, ri_average, ri_stdev, rei_cor, rei_var, ne1, Ne1, distribution );
		fe2 = sample_f( fe2, fe1, fi2, fi1, -re_average, re_stdev, -ri_average, ri_stdev, rei_cor, rei_var, ne2, Ne2, distribution );
		
		if( i >= burn_in ) // the burn_in period has passed
		{
			// for debugging purposes
			//cout << ".";

			// in this version, the results are reported based on natural logarithm
			// we can change this letter to log2 by dividing the next four lines by log(2)
			double log_fi1 = log(fi1);
			double log_fi2 = log(fi2);
			double log_fe1 = log(fe1);
			double log_fe2 = log(fe2);
			
			(*log_fi1_average) += log_fi1;
			(*log_fi1_stdev) += ( log_fi1*log_fi1 );
			(*log_fi2_average) += log_fi2;
			(*log_fi2_stdev) += ( log_fi2*log_fi2 );
			(*log_fe1_average) += log_fe1;
			(*log_fe1_stdev) += ( log_fe1*log_fe1 );
			(*log_fe2_average) += log_fe2;
			(*log_fe2_stdev) += ( log_fe2*log_fe2 );
			(*log_fi21_average) += ( log_fi2-log_fi1 );
			(*log_fi21_stdev) += ( ( log_fi2-log_fi1 ) * ( log_fi2-log_fi1 ) );
			(*log_fe21_average) += ( log_fe2-log_fe1 );
			(*log_fe21_stdev) += ( ( log_fe2-log_fe1 ) * ( log_fe2-log_fe1 ) );
			(*log_fei21_average) += ( ( log_fe2-log_fe1 ) - ( log_fi2-log_fi1 ) );
			(*log_fei21_stdev) += ( ( ( log_fe2-log_fe1 ) - ( log_fi2-log_fi1 ) ) * ( ( log_fe2-log_fe1 ) - ( log_fi2-log_fi1 ) ) );
		
			n ++;
		}

		// for debugging purposes
		//cout << __MCMC_accepted << endl;
	}

	// for debugging purposes
	cerr << __MCMC_accepted << endl;
	
	(*log_fi1_average) /= n;
	(*log_fi1_stdev) /= n;
	*log_fi1_stdev = sqrt( (*log_fi1_stdev) - (*log_fi1_average) * (*log_fi1_average) );

	(*log_fi2_average) /= n;
	(*log_fi2_stdev) /= n;
	*log_fi2_stdev = sqrt( (*log_fi2_stdev) - (*log_fi2_average) * (*log_fi2_average) );

	(*log_fe1_average) /= n;
	(*log_fe1_stdev) /= n;
	*log_fe1_stdev = sqrt( (*log_fe1_stdev) - (*log_fe1_average) * (*log_fe1_average) );


	(*log_fe2_average) /= n;
	(*log_fe2_stdev) /= n;
	*log_fe2_stdev = sqrt( (*log_fe2_stdev) - (*log_fe2_average) * (*log_fe2_average) );

	(*log_fi21_average) /= n;
	(*log_fi21_stdev) /= n;
	*log_fi21_stdev = sqrt( (*log_fi21_stdev) - (*log_fi21_average) * (*log_fi21_average) );

	(*log_fe21_average) /= n;
	(*log_fe21_stdev) /= n;
	*log_fe21_stdev = sqrt( (*log_fe21_stdev) - (*log_fe21_average) * (*log_fe21_average) );

	(*log_fei21_average) /= n;
	(*log_fei21_stdev) /= n;
	*log_fei21_stdev = sqrt( (*log_fei21_stdev) - (*log_fei21_average) * (*log_fei21_average) );

}

////////////////////////////////////////////////////////////////////////////////////
void _Metropolis_Hastings_per_group(
	s_gene *genes[],
	int num_genes,
	s_array *arrays[],
	int num_arrays,
	s_group *groups[],
	int group_index,
	int distribution,
	int learn,
	int burn_in,
	int sampling )
// performs the MH algorithm for group_index
{
	int Ni1 = 0;
	int Ne1 = 0;
	int Ni2 = 0;
	int Ne2 = 0;

	// initialize ni1, ne1, ni2 and ne2 parameters for each gene	
	int i, j;
	for( i = 0; i < num_genes; i ++ )
	{
		genes[ i ] ->ni1 =
		genes[ i ] ->ne1 =
		genes[ i ] ->ni2 =
		genes[ i ] ->ne2 = 0;
		
		for( j = 0; j < num_arrays; j ++ )
			if( arrays[ j ] ->group_index == group_index ) // update ni2 and ne2
			{
				if( arrays[ j ] ->ei == INTRONIC )
					genes[ i ] ->ni2 += genes[ i ] ->profile[ j ];
				else
					genes[ i ] ->ne2 += genes[ i ] ->profile[ j ];
			}
			else
			{
				if( arrays[ j ] ->ei == INTRONIC ) // update ni1 and ne1
					genes[ i ] ->ni1 += genes[ i ] ->profile[ j ];
				else
					genes[ i ] ->ne1 += genes[ i ] ->profile[ j ];
			}
			
		Ni1 += genes[ i ] ->ni1;
		Ne1 += genes[ i ] ->ne1;
		Ni2 += genes[ i ] ->ni2;
		Ne2 += genes[ i ] ->ne2;
	}
	
	// optimize the correlation of dfi and dfe, and obtain the parameters of their joint distribution
	double ri_average, ri_stdev, re_average, re_stdev, rei_cor, rei_var;
	if( learn ) // the parameters of the prior distribution will be learned from data
	{
		int max_pos = optimize_cor( genes, num_genes, Ni1, Ne1, Ni2, Ne2, &ri_average, &ri_stdev, &re_average, &re_stdev, &rei_cor, &rei_var );
		cout << "Distribution statistics, learned from the top " << max_pos << " most abundant genes:" << endl
			<< "Mean  log(fi2/fi)             = " << ri_average << endl
			<< "Stdev log(fi2/fi1)            = " << ri_stdev << endl
			<< "Mean  log(fe2/fe1)            = " << re_average << endl
			<< "Stdev log(fe2/fe1)            = " << re_stdev << endl
			<< "Pearson correlation           = " << rei_cor << endl
			<< "Var log(fe2/fe1)-log(fi2/fi1) = " << rei_var << endl;
	}
	else // the parameters of the prior distribution have been given in the metadata file
	{
		ri_average = groups[ group_index ] ->ri_average;
		ri_stdev = groups[ group_index ] ->ri_stdev;
		re_average = groups[ group_index ] ->re_average;
		re_stdev = groups[ group_index ] ->re_stdev;
		rei_cor = groups[ group_index ] ->rei_cor;
		rei_var = groups[ group_index ] ->rei_var;

		cout << "Distribution statistics, from the metadata file or default values:" << endl
			<< "Mean  log(fi2/fi)             = " << ri_average << endl
			<< "Stdev log(fi2/fi1)            = " << ri_stdev << endl
			<< "Mean  log(fe2/fe1)            = " << re_average << endl
			<< "Stdev log(fe2/fe1)            = " << re_stdev << endl
			<< "Pearson correlation           = " << rei_cor << endl
			<< "Var log(fe2/fe1)-log(fi2/fi1) = " << rei_var << endl;
	}
	
	// now, for each gene, sample the distribution parameters of f-related values
	for( i = 0; i < num_genes; i ++ )
	{	
		double log_fi1_average, log_fi1_stdev, log_fi2_average, log_fi2_stdev,
			log_fe1_average, log_fe1_stdev, log_fe2_average, log_fe2_stdev,
			log_fi21_average, log_fi21_stdev, log_fe21_average, log_fe21_stdev,
			log_fei21_average, log_fei21_stdev;

		_Metropolis_Hastings_per_gene( ri_average, ri_stdev, re_average, re_stdev, rei_cor, rei_var,
			genes[ i ] ->ni1, genes[ i ] ->ne1, genes[ i ] ->ni2, genes[ i ] ->ne2,
			Ni1, Ne1, Ni2, Ne2,
			distribution, burn_in, sampling,
			&log_fi1_average, &log_fi1_stdev, &log_fi2_average, &log_fi2_stdev,
			&log_fe1_average, &log_fe1_stdev, &log_fe2_average, &log_fe2_stdev,
			&log_fi21_average, &log_fi21_stdev, &log_fe21_average, &log_fe21_stdev,
			&log_fei21_average, &log_fei21_stdev );
			
		genes[ i ] ->log_fi_average[ group_index ] = log_fi2_average;
		genes[ i ] ->log_fi_stdev[ group_index ] = log_fi2_stdev;
		genes[ i ] ->log_fe_average[ group_index ] = log_fe2_average;
		genes[ i ] ->log_fe_stdev[ group_index ] = log_fe2_stdev;

		genes[ i ] ->log_dfi_average[ group_index ] = log_fi21_average;
		genes[ i ] ->log_dfi_stdev[ group_index ] = log_fi21_stdev;
		genes[ i ] ->log_dfe_average[ group_index ] = log_fe21_average;
		genes[ i ] ->log_dfe_stdev[ group_index ] = log_fe21_stdev;

		genes[ i ] ->log_ddf_average[ group_index ] = log_fei21_average;
		genes[ i ] ->log_ddf_stdev[ group_index ] = log_fei21_stdev;
	}
}

////////////////////////////////////////////////////////////////////////////////////
void Metropolis_Hastings(
	s_gene *genes[],
	int num_genes,
	s_array *arrays[],
	int num_arrays,
	s_group *groups[],
	int num_groups,
	int distribution,
	int learn,
	int burn_in,
	int sampling )
{
	// initialize the memory
	int i;
	for( i = 0; i < num_genes; i ++ )
	{
		genes[ i ] ->log_fi_average = new double[ num_groups];
		genes[ i ] ->log_fi_stdev = new double[ num_groups];
		genes[ i ] ->log_fe_average = new double[ num_groups];
		genes[ i ] ->log_fe_stdev = new double[ num_groups];

		genes[ i ] ->log_dfi_average = new double[ num_groups];
		genes[ i ] ->log_dfi_stdev = new double[ num_groups];
		genes[ i ] ->log_dfe_average = new double[ num_groups];
		genes[ i ] ->log_dfe_stdev = new double[ num_groups];

		genes[ i ] ->log_ddf_average = new double[ num_groups];
		genes[ i ] ->log_ddf_stdev = new double[ num_groups];
	}

	
	// run the MCMC
	for( i = 0; i < num_groups; i ++ )
	{
		cout << endl << "## Running MCMC for group " << groups[ i ] ->name << " ..." << endl;
		_Metropolis_Hastings_per_group( genes, num_genes, arrays, num_arrays, groups, i,
			distribution, learn, burn_in, sampling );
	}
}
