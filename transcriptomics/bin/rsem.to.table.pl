#!/usr/bin/perl -w

#########################################################################################
#                                       rsem.to.table.pl
#########################################################################################
# 
#  This program runs summarizes RSEM runs into a single table
#
#########################################################################################
# AUTHORS:
#
# Alper Kucukural, PhD 
# 
#########################################################################################

############## LIBRARIES AND PRAGMAS ################

 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use POSIX 'floor';
 
#################### VARIABLES ######################
 my $gene_iso         = "genes";
 my $quantType         = "tpm";
 my $out              = "";
 my $indir            = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $cmd=$0." ".join(" ",@ARGV); ####command line copy

GetOptions(
	'out=s'           => \$out,
	'indir=s'         => \$indir,
	'gene_iso=s'      => \$gene_iso,
	'quantType=s'      => \$quantType,
	'help'            => \$help, 
	'version'         => \$print_version,
) or die("Unrecognized optioins.\nFor help, run this script with -help option.\n");

if($help){
    pod2usage( {
		'-verbose' => 2, 
		'-exitval' => 1,
	} );
}

if($print_version){
  print "Version ".$version."\n";
  exit;
}

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($out eq "") );	

if (!($quantType  eq "tpm" || $quantType eq "fpkm" || $quantType eq "expected_count"))
{
    $quantType="tpm";
}

my %tf = (
        expected_count => 5,
        tpm => 5,
        fpkm => 6,
    );
 
################### MAIN PROGRAM ####################++
#    maps the reads to the the genome and put the files under $out/after_ribosome/tophat directory

opendir I,$indir or die "Could not open $indir\n";
my @files = grep /$gene_iso.results$/, readdir(I);
closedir I; 
print "quantType $quantType ... col: $tf{$quantType}\n";
my @a=();
my %b=();
my %c=();
my $i=0;
foreach my $f (@files){ 
  my $in = "${indir}/$f";
  my $libname=$f;
  $libname=~s/\.$gene_iso.results//;

  $i++;
  $a[$i]=$libname;
  print "${in}\n";
  open IN,"$in" or die "Could not open $in\n";
  my $lineNum = 0;
  while(<IN>)
    {
      next if ($lineNum++ == 0);
      my @v=split; 
      #print "Using $quantType $tf{$quantType}\n";
      my $quant = $v[$tf{$quantType}];
      $quant = floor ($quant + 0.5); #if ($quantType eq "expected_count");
      $b{$v[0]}{$i}= $quant;
      $c{$v[0]}=$v[1];
    }
  close IN;
}

open OUT, ">${out}" or die "Could not open $out file for writting\n";

print OUT "Id";


for(my $j=1;$j<=$i;$j++)
{
 print OUT "\t$a[$j]";
}
print OUT "\n";

foreach my $key (keys %b)
{
 if ($gene_iso ne "isoforms") {
   #print OUT "$key\t$c{$key}"; 
   print OUT "$key"; #only print out the gene name;
 }
 else
 {
    #print OUT "$c{$key}\t$key";
   print OUT "$c{$key}_${key}"; #only print out the transcript id
 }
 for(my $j=1;$j<=$i;$j++)
 {
  print OUT "\t$b{$key}{$j}";
 }
 print OUT "\n";
}

close OUT;

__END__


=head1 NAME

rsem.to.table.pl

=head1 SYNOPSIS  

rsem.to.table.pl 
            -o out <output summary table> 
            -i indir <input directory>
            -t quantType <quantification type to summarize tpm or  expected_count>
	    -g gene_iso <gene or isoform>


rsem.to.table.pl -help

rsem.to.table.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -c conversion <ucsc conversion> 

Tab delimited ucsc id gene name conversion file

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

This program runs the Cuffdiff after cufflinks

=head1 AUTHORS

 Alper Kucukural, PhD

 
=head1 LICENSE AND COPYING

 This program is free software; you can redistribute it and / or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, a copy is available at
 http://www.gnu.org/licenses/licenses.html

