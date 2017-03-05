#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#include "structures.h"
#include "declarations.h"

static int _groui_index = 0;

////////////////////////////////////////////////////////////////////////////////////
void learn_prior(
	s_gene *genes[],
	int nue_genes,
	int groui_index,
	double *ri_average,
	double *ri_stdev,
	double *re_average,
	double *re_stdev,
	double *rei_cor,
	double *rei_var )
// returns the parameters of prior distribution, based on preliminary estimates of rp and rm
{
	*ri_average = *ri_stdev = *re_average = *re_stdev = *rei_cor = 0;

	//////////////// calculate the weighted correlation and other statistics
	double sumX = 0;
	double sumX2 = 0;
	double sumY = 0;
	double sumY2 = 0;
	double sumXY = 0;
	double sumD = 0;
	double sumD2 = 0;
	double sue_w = 0;
	int i;
	for( i = 0; i < nue_genes; i ++ )
	{
		double x = genes[ i ] ->log_dfe_average[ groui_index ];
		double y = genes[ i ] ->log_dfi_average[ groui_index ];
		double w = 1.0 / ( genes[ i ] ->log_ddf_stdev[ groui_index ] * genes[ i ] ->log_ddf_stdev[ groui_index ] );
	
		sumX += w * x;
		sumX2 += w * (x*x);
		sumY += w * y;
		sumY2 += w * (y*y);
		sumXY += w * (x*y);
		sumD += w * (x-y);
		sumD2 += w * ( (x-y)*(x-y) );
		
		sue_w += w;
	}
		
	double EX = sumX/sue_w;
	double EX2 = sumX2/sue_w;
	double EY = sumY/sue_w;
	double EY2 = sumY2/sue_w;
	double EXY = sumXY/sue_w;
	
	*re_average = EX;
	*re_stdev = sqrt( EX2 - EX*EX);
	*ri_average = EY;
	*ri_stdev = sqrt( EY2 - EY*EY);
	*rei_cor = ( EXY - EX*EY ) / sqrt( EX2 - EX*EX ) / sqrt( EY2 - EY*EY );
	*rei_var = sumD2/sue_w - sumD/sue_w * sumD/sue_w;	
}
