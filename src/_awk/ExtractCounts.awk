{
	if( ReadType == "exonic" && substr($1,1,5) == "exon_" )
		printf( "%s\t%s\n", substr($1,6), $2 );
	else if( ReadType == "intronic" && substr($1,1,7) == "intron_" )
		printf( "%s\t%s\n", substr($1,8), $2 );
}
