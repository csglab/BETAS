#!/usr/bin/perl

# This Perl script is courtesy of Harm van Bakel <hvbakel@gmail.com>

# GLOBALS
$ENV{TMPDIR}     ||= '/tmp';     # location for tmp file storage
$ENV{PASTE}      ||= 'paste';    # Unix paste utility
$ENV{BATCH_SIZE} ||= 300;        # Number of temp files to paste together in one batch
$SIG{'INT'}=$SIG{'HUP'}=$SIG{'ABRT'}=$SIG{'QUIT'}=$SIG{'TRAP'}=$SIG{'STOP'}=\&INTERRUPT;

# MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Temp qw(tempfile tempdir);

# GET PARAMETERS
my $sHelp           = 0;
my $nIDcol          = 1;
my $nDataCol        = 2;
my $flCaseSensitive = 0;
my $flCommonIDsOnly = 0;
GetOptions("help!"      => \$sHelp,
           "sensitive!" => \$flCaseSensitive,
           "id-col:n"   => \$nIDcol,
           "data-col:n" => \$nDataCol,
           "common!"    => \$flCommonIDsOnly);

# PRINT HELP
$sHelp = 1 unless(@ARGV>1);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [-c -s] <file1> .. <fileN>

   Joins one column in a set of files together in a large matrix
   based on identifiers shared between all files.
    
    -id-col <integer>
      Column containing the identifier
      default: $nIDcol
    -data-col <integer>
      Column containing the data values to merge
      default: $nDataCol
    -s
      Do case-sensitive matching (default is to ignore case)
    -c
      Only consider IDs shared between all files. The default
      is to replace missing values by 'NA'.
    -help
      This help message
   
HELP
}


###########
## START ##
###########

# Make column IDs zero-based
$nIDcol--;
$nDataCol--;

# Get a list of IDs common to all files
my @asJoinIDs = get_ids($nIDcol, $flCaseSensitive, $flCommonIDsOnly, @ARGV);
die "Error: could not find any join IDs\n" unless(@asJoinIDs);

# Create the first output temp file
my ($fhTmpOut, $sTmpOut) = tempfile('multi-join-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
$fhTmpOut->print("GeneID\n");
for my $sID (@asJoinIDs){
   $fhTmpOut->print("$sID\n");
}
$fhTmpOut->close();

# Now start adding the other files one by one
my @asFileBatch;
for my $sFile (@ARGV){
   
   # Read the ID and data column into a hash
   my %hFileData;
   open IN, $sFile or die "Error: can't open '$sFile': $!\n";
   while (<IN>){
      next if /^\s*$/;
      next if /^\s*#/;
      s/[\n\r]$//g;
      my (@asLine) = split /\t/, $_, -1;
      die "Error: data column out of range in file '$sFile', line $.\n" if ($nDataCol > $#asLine);
      $asLine[$nIDcol] = uc($asLine[$nIDcol]) unless ($flCaseSensitive);
      $hFileData{$asLine[$nIDcol]} = $asLine[$nDataCol];
   }
   close IN;
   
   # Print the file data to a temporary file
   my ($fhCurOut, $sCurOut) = tempfile('multi-join-batchfile-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
   push @asFileBatch, $sCurOut;
   print $fhCurOut basename($sFile), "\n";
   foreach my $sID (@asJoinIDs){
      $sID = uc($sID) unless ($flCaseSensitive);
      if (exists $hFileData{$sID} ){
         print $fhCurOut $hFileData{$sID}, "\n";
      }
      else{
         print $fhCurOut "NA\n";
      }
   }
   
   # Check if we hit the maximum batch size
   if (@asFileBatch >= $ENV{BATCH_SIZE}){
      $sTmpOut = merge_batch($sTmpOut, @asFileBatch);
      @asFileBatch = ();
   }
}

# Add the last batch of files if necessary
$sTmpOut = merge_batch($sTmpOut, @asFileBatch) if (@asFileBatch);

# Finally, print the merged matrix
open OUT, $sTmpOut or die "Error: can't open merged temporary matrix\n";
while (<OUT>){
   print;
}
close OUT;
unlink($sTmpOut);



#################
## SUBROUTINES ##
#################


# get_ids
#
# Extract the identifiers from the set of files
sub get_ids{
   my ($nColID, $flCaseSensitive, $flCommonIDsOnly, @asFiles) = @_;
   my %hIDs;
   my @asIDs;
   my $nSortOrder = 0;
   
   for my $sFile (@asFiles){
      my %hUniqueCheck;
      open IN, $sFile or die "Error: can't open '$sFile': $!\n";
      while (<IN>){
         next if /^\s*$/;
         next if /^\s*#/;
         s/[\n\r]$//g;
         my (@asLine) = split /\t/, $_, -1;
         die "Error: ID column out of range in file '$sFile', line $.\n" if ($nColID > $#asLine);
         my $sKey = $flCaseSensitive ? $asLine[$nColID] : uc($asLine[$nColID]);
         die "Error: found non-unique identifiers in file '$sFile', line $.\n" if (exists $hUniqueCheck{$sKey});
         $hIDs{$sKey}{'identifier'} = $asLine[$nColID] unless (exists $hIDs{$sKey}{'identifier'});
         $hIDs{$sKey}{'sortorder'}  = $nSortOrder++ unless (exists $hIDs{$sKey}{'sortorder'});
         $hIDs{$sKey}{'count'}++;
         $hUniqueCheck{$sKey}++;
      }
      close IN;
   }
   
   # Return common identifiers while maintaining input sort order
   my $nMissing = 0;
   for my $sKey (sort {$hIDs{$a}{'sortorder'} <=> $hIDs{$b}{'sortorder'}} keys %hIDs){
      if ($flCommonIDsOnly){
         if ($hIDs{$sKey}{'count'} == scalar(@asFiles) ){
            push @asIDs, $sKey;
         }
         else{
            $nMissing++;
         }
      }
      else{
         push @asIDs, $sKey;
      }
   }
   warn("Warning: $nMissing IDs were skipped in the combined matrix since they were not present in all files\n") if ($nMissing and $flCommonIDsOnly);
   return @asIDs;
}


# merge_batch
#
# Merges a batch of files
sub merge_batch{
   my ($sReference, @asTmpFiles) = @_;
   
   # Create new temp file to hold the merged matrix
   my ($fhNewOut, $sNewOut) = tempfile('multi-join-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
   my $sPasteArgs = join(' ', $sReference, @asTmpFiles);
   open PASTE, "$ENV{PASTE} $sPasteArgs |" or die "Error: could not run 'paste' to combine files\n";
   while (<PASTE>){
      print $fhNewOut $_;
   }
   close PASTE;
   close $fhNewOut;
         
   # Delete all files in this batch and reset the other variables
   unlink($sReference, @asFileBatch);
   return $sNewOut;
}

# INTERRUPT
#
# Interrupt routine, make sure we exit gracefully for tmp file cleanup
sub INTERRUPT{
   exit(1); # This will call END
}
