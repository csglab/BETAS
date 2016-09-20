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

static int __MCMC_samples[ 100 ];
static int __MCMC_pos = 0;
int __MCMC_accepted = 0;

////////////////////////////////////////////////////////////////////////////////////
double _get_f_log_likelihood(
	double new_fa1,
	double old_fa1,
	double fa2,
	double fb1,
	double fb2,
	double ra_average,
	double ra_stdev,
	double rb_average,
	double rb_stdev,
	double rab_cor,
	double rab_var,
	double na1,
	double Na1,
	int distribution )
// new_fa1: the newly proposed value for f
// old_fa1: the previous value for f
// fa2: the value of f for the "same" read type in "other" groups
// fb1: the value of f for the "other" read type in "same" group
// fb2: the value of f for the "other" read type in "other" groups
// ra_average: the prior mean of log(f_a2)-log(f_a1)
// ra_stdev: the prior standard deviation of log(f_a2)-log(f_a1)
// rb_average: the prior mean of log(f_b2)-log(f_b1)
// rb_stdev: the prior standard deviation of log(f_b2)-log(f_b1)
// rab_cor: the prior correlation of log(f_a2)-log(f_a1) vs. log(f_b2)-log(f_b1)
// rab_var: the prior variance of (log(f_a2)-log(f_a1))-(log(f_b2)-log(f_b1))
// na1: the number of reads from the "same" read type in the "same" group for this gene
// Na1: the number of reads from the "same" read type in the "same" group for all genes
// distribution: NORMAL or LAPLACE
{
	double log_likelihood =
		na1 * ( log(new_fa1) - log(old_fa1) ) +
		(Na1-na1) * ( log(1-new_fa1) - log(1-old_fa1) );
		
	double new_ra= log(fa2)-log(new_fa1);
	double old_ra= log(fa2)-log(old_fa1);
	double rb=log(fb2)-log(fb1);
		
	if( distribution == NORMAL )
	{
		log_likelihood += ( -0.5 / (1-rab_cor*rab_cor) ) *
			( ( (new_ra-ra_average)*(new_ra-ra_average) - (old_ra-ra_average)*(old_ra-ra_average) ) / (ra_stdev*ra_stdev)
			- 2 * rab_cor * (new_ra-old_ra) * (rb-rb_average) / (ra_stdev*rb_stdev) );
	}
	else if( distribution == LAPLACE )
	{
		double det= ra_stdev*ra_stdev * rb_stdev*rb_stdev * (1-rab_cor*rab_cor); // the determinant of sigma
		
		double ta = old_ra - ra_average; // the first number in the t vector
		double tb = rb - rb_average; // the second number in the t vector
		double px = ta*(rb_stdev*rb_stdev) - tb*ra_stdev*rb_stdev*rab_cor; // the first number in t^T * sigma^(-1)
		double py = tb*(ra_stdev*ra_stdev) - ta*ra_stdev*rb_stdev*rab_cor; // the second number in t^T * sigma^(-1)
		double old_z = (px*ta + py*tb) / det; // t^T * sigma^(-1) * t

		ta = new_ra - ra_average;
		px = ta*(rb_stdev*rb_stdev) - tb*ra_stdev*rb_stdev*rab_cor;
		py = tb*(ra_stdev*ra_stdev) - ta*ra_stdev*rb_stdev*rab_cor;
		double new_z = (px*ta + py*tb) / det;
		
		log_likelihood += sqrt(2*old_z) - sqrt(2*new_z);
	}
	else if( distribution == LINEARGRADIENT )
	{
		log_likelihood += (
			( (old_ra-rb)*(old_ra-rb) - (new_ra-rb)*(new_ra-rb) ) / ( 2 * rab_var )
			);
	}
	// if none of the above, uniform distribution is assumed, which does not affect the likelihood
		
	return exp( log_likelihood );
}

////////////////////////////////////////////////////////////////////////////////////
double sample_f(
	double old_fa1,
	double fa2,
	double fb1,
	double fb2,
	double ra_average,
	double ra_stdev,
	double rb_average,
	double rb_stdev,
	double rab_cor,
	double rab_var,
	double na1,
	double Na1,
	int distribution )
// returns a newly sampled fa1
//
// See the document MCMC.docx for details
//
// old_fa1: the previous value for f
// fa2: the value of f for the "same" read type in "other" groups
// fb1: the value of f for the "other" read type in "same" group
// fb2: the value of f for the "other" read type in "other" groups
// ra_average: the prior mean of log(f_a2)-log(f_a1)
// ra_stdev: the prior standard deviation of log(f_a2)-log(f_a1)
// rb_average: the prior mean of log(f_b2)-log(f_b1)
// rb_stdev: the prior standard deviation of log(f_b2)-log(f_b1)
// rab_cor: the prior correlation of log(f_a2)-log(f_a1) vs. log(f_b2)-log(f_b1)
// rab_var: the prior variance of (log(f_a2)-log(f_a1))-(log(f_b2)-log(f_b1))
// na1: the number of reads from the "same" read type in the "same" group for this gene
// Na1: the number of reads from the "same" read type in the "same" group for all genes
// distribution: NORMAL or LAPLACE
{
	/////// update the statistics for acceptance rate
	// increase position by 1, reset to zero if exceeds 99
	// then update the status of the new position
	__MCMC_pos ++;
	if( __MCMC_pos > 99 )
		__MCMC_pos = 0;
	__MCMC_accepted -= __MCMC_samples[ __MCMC_pos ];
	__MCMC_samples[ __MCMC_pos ] = 0;
	/////// by default, assumes that this sample will be rejected
		
	double E_x; // the mean of the proposal distribution
	double stdev_x; // the standard deviation of the proposal distribution
	double new_fa1;
	
	/*if( distribution != UNIFORM && drand48() < 0.5 ) // use the first proposal distribution half of the times, if the ri vs. re distribution is not set to uniform
	{
		E_x = log(old_fa1); // the proposal distribution must be symmetric
		stdev_x = ra_stdev * sqrt( 1 - rab_cor*rab_cor );

		// choose a new fa1
		new_fa1 = exp( generate_normal_random( E_x, 2.38 * stdev_x ) );
	}*/
	//else // use the second proposal distribution
	//{
		E_x = old_fa1; // the proposal distribution must be symmetric
		stdev_x = sqrt( double(na1+1)/double(Na1+2) * ( double(na1+2)/double(Na1+3) - double(na1+1)/double(Na1+2) ) );

		// choose a new fa1
		new_fa1 = generate_normal_random( E_x, 2.38 * stdev_x );
	//}
	
	if( new_fa1 <= 0 || new_fa1 >= 1 ) // if the new fa1 is out of boundaries, it is automatically rejected (its prior is zero); zero is not allowed, in order to prevent taking log of zero
		return old_fa1;

	// calculate the ratio of the probabilities of the new vs. old value	
	double likelihood =
		_get_f_log_likelihood( new_fa1, old_fa1,
			fa2, fb1, fb2, ra_average, ra_stdev, rb_average, rb_stdev, rab_cor, rab_var, na1, Na1,
			distribution );
		
	if( drand48() < likelihood )
	{
		// update the statistics for acceptance rate
		__MCMC_samples[ __MCMC_pos ] = 1;
		__MCMC_accepted ++;
		
		return new_fa1; // the new value is accepted
		
	} // since drand48 always returns a value <=1, likelihoods >1 are always accepted
		
	return old_fa1; // the new value is rejected
}
