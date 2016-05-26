#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#include "structures.h"
#include "declarations.h"

#define SLIDING_WINDOW	1000

////////////////////////////////////////////////////////////////////////////////////
int _compare_gene( const void *a, const void *b )
{
	s_gene *gene1 = *( (s_gene **) a );
	s_gene *gene2 = *( (s_gene **) b );
	
	if( gene2 ->min_f > gene1 ->min_f )
		return 1;
	else if( gene2 ->min_f < gene1 ->min_f )
		return -1;
	
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////
int optimize_cor(
	s_gene *genes[],
	int num_genes,
	int Ni1,
	int Ne1,
	int Ni2,
	int Ne2,
	double *ri_average,
	double *ri_stdev,
	double *re_average,
	double *re_stdev,
	double *rei_cor,
	double *rei_var )
// optimizes the correlation
// returns the number of genes included in the optimal gene set
{
	*ri_average = *ri_stdev = *re_average = *re_stdev = *rei_cor = 0;
	
	//////////////// First, store an initial assessment of statistics to be used for gene filtering and estimating the prior
	int i;
	for( i = 0; i < num_genes; i ++ )
	{
		// find and store min_n
		genes[ i ] ->min_n =
			MIN( genes[ i ] ->ni1, MIN ( genes[ i ] ->ni2, MIN( genes[ i ] ->ne1, genes[ i ] ->ne2 ) ) );

		if( genes[ i ] ->min_n > 0 )
		{
			genes[ i ] ->ri_prior =
				( log( genes[ i ] ->ni2 ) - log( Ni2 ) ) -
				( log( genes[ i ] ->ni1 ) - log( Ni1 ) );
			genes[ i ] ->re_prior =
				( log( genes[ i ] ->ne2 ) - log( Ne2 ) ) -
				( log( genes[ i ] ->ne1 ) - log( Ne1 ) );
		}
		else
			genes[ i ] ->ri_prior = genes[ i ] ->re_prior = 0;

		// find and store min_f
		double fi1 = double(genes[ i ] ->ni1+1)/double(Ni1+2);
		double fe1 = double(genes[ i ] ->ne1+1)/double(Ne1+2);
		double fi2 = double(genes[ i ] ->ni2+1)/double(Ni2+2);
		double fe2 = double(genes[ i ] ->ne2+1)/double(Ne2+2);
		genes[ i ] ->min_f =
			MIN( fi1, MIN ( fe1, MIN ( fi2, fe2 ) ) );
			
		// find and store the expected variance of initial estimate of log(fe2/fe1)-log(fi2/fi1), assuming no correlation
		double var_log_fi1 = ( double(genes[i]->ni1+2)/double(Ni1+3) - fi1 ) / fi1;
		double var_log_fi2 = ( double(genes[i]->ni2+2)/double(Ni2+3) - fi2 ) / fi2;
		double var_log_fe1 = ( double(genes[i]->ne1+2)/double(Ne1+3) - fe1 ) / fe1;
		double var_log_fe2 = ( double(genes[i]->ne2+2)/double(Ne2+3) - fe2 ) / fe2;		
		genes[ i ] ->var_log_f = var_log_fi1 + var_log_fi2 + var_log_fe1 + var_log_fe2;
			
	}
	
	// sort the genes by decreasing order of min_f
	qsort( genes, num_genes, sizeof(s_gene*), _compare_gene );

	//////////////// Next, calculate the overall correlation
	double sumX = 0;
	double sumX2 = 0;
	double sumY = 0;
	double sumY2 = 0;
	double sumXY = 0;
	double sumD = 0;
	double sumD2 = 0;
	int n = 0;
	for( i = 0; i < num_genes; i ++ )
		if( genes[ i ] ->min_n > 0 )
		{
			n ++;
			double x = genes[ i ] ->ri_prior;
			double y = genes[ i ] ->re_prior;
		
			sumX += x;
			sumX2 += (x*x);
			sumY += y;
			sumY2 += (y*y);
			sumXY += (x*y);
			sumD += (x-y);
			sumD2 += ( (x-y)*(x-y) );
		}
		
	double EX = sumX/n;
	double EX2 = sumX2/n;
	double EY = sumY/n;
	double EY2 = sumY2/n;
	double EXY = sumXY/n;
	double cor_all = ( EXY - EX*EY ) / sqrt( EX2 - EX*EX ) / sqrt( EY2 - EY*EY );
	cout << "Overall correlation is " << cor_all << " (n=" << n << ")" << endl;
	cout << "Overall variance of Ri-Re is " << sumD2/n - sumD/n * sumD/n << endl;
	
	// by default, take the overall statistics
	*ri_average = EX;
	*ri_stdev = sqrt( EX2 - EX*EX);
	*re_average = EY;
	*re_stdev = sqrt( EY2 - EY*EY);
	*rei_cor = cor_all;
	
	if( num_genes < 500 ) // if there are less than 500 genes, the correlation cannot be optimized
		return num_genes;


	cerr << "Calculating cutoff statistics ..." << endl
		<< "n" << char(9) << "i" << char(9) << "MinF" << char(9) << "Correlation" << char(9)
		<< "Mean(Ri)" << char(9) << "Stdev(Ri)" << char(9)
		<< "Mean(Re)" << char(9) << "Stdev(Re)" << char(9)
		<< "Var(Ri-Re)" << endl;

	////////////////  now, calculate the running correlation for sliding windows of size SLIDING_WINDOW, and find the maximum
	sumX = sumX2 = sumY = sumY2 = sumXY = sumD = sumD2 = 0;
	n = 0;
	double max_cor = -10;
	int max_pos = -1;
	for( i = 0; i < num_genes; i ++ )
		if( genes[ i ] ->min_n > 0 )
		{
			n ++;
			double x = genes[ i ] ->ri_prior;
			double y = genes[ i ] ->re_prior;
		
			sumX += x;
			sumX2 += (x*x);
			sumY += y;
			sumY2 += (y*y);
			sumXY += (x*y);
			sumD += (x-y);
			sumD2 += ( (x-y)*(x-y) );
		
			if( n >= 500 ) // at least 500 genes are processed, and its time to calculate sliding window statistics
			{
				EX = sumX/n;
				EX2 = sumX2/n;
				EY = sumY/n;
				EY2 = sumY2/n;
				EXY = sumXY/n;
			
				if( EX2 > EX*EX && EY2 > EY*EY ) // the standard deviations are not zero
				{
					double cor = ( EXY - EX*EY ) / sqrt( EX2 - EX*EX ) / sqrt( EY2 - EY*EY );
				
					cerr << n << char(9) << i << char(9) << genes[ i ] ->min_f << char(9) << cor << char(9)
						<< EX << char(9) << sqrt( EX2 - EX*EX) << char(9)
						<< EY << char(9) << sqrt( EY2 - EY*EY) << char(9)
						<< sumD2/n - sumD/n * sumD/n << endl;

					if( max_cor < cor ) // this is the maximum correlation so far
						max_cor = cor;
					
					//if( cor >= cor_all + (max_cor-cor_all) * 0.99 ) // this is within 99% of the maximum correlation improvement
					if( cor >= cor_all + (max_cor-cor_all) * 0.90 ) // this is within 90% of the maximum correlation improvement
					{
						*ri_average = EX;
						*ri_stdev = sqrt( EX2 - EX*EX);
						*re_average = EY;
						*re_stdev = sqrt( EY2 - EY*EY);
						*rei_cor = cor;
						*rei_var = sumD2/n - sumD/n * sumD/n;
						max_pos = i+1;
					}
				}
			}
		}
	
	return max_pos;
}
