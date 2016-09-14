{
	if( NF==3 ) # this is the metadata file
	{
		label[ $1 ] = $2;
		type[ $1 ] = $3;
	}
	else # this is the count file
	{
		lnum ++;
		if( lnum == 1 ) # the first line of the count file
		{
			printf( "%s", $1 );			
			for( i = 2; i <= NF; i ++ )
			{
				if( label[ $i ] == sampleLabel ) # this column corresponds to the specified sample
				{
					print_as_is[ i ] = 1;
					printf( "\t%s", $i );
				}
					
				# store the type of read this column represents
				read_type[ i ] = type[ $i ];
			}
			
			printf( "\tAll.exon\tAll.intron\n" );
		}
		else
		{
			printf( "%s", $1 );
			intron_count = 0;
			exon_count = 0;
			for( i = 2; i <= NF; i ++ )
			{
				if( print_as_is[ i ] == 1 ) # this column corresponds to the specified sample
				{
					printf( "\t%s", $i );
				}
				else
				{
					if( read_type[ i ] == "intronic" )
						intron_count += $i;
					else if( read_type[ i ] == "exonic" )
						exon_count += $i;
					else
						printf( "Error: unknown read count type for column %s\n", $i );
				}
			}

			printf( "\t%i\t%i\n", exon_count, intron_count );
		}
	}	
}
