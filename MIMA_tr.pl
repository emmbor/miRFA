#!/usr/bin/perl

use warnings;

#Open files
open(my $mRNA_file, '<', "mirb_IDANDname.csv");

#Open .csv file, only first line

open my $file, '<', "MIMA_t_mirnas.csv"; 
my $firstLine = <$file>; 
close $file;


my %hash;

#Loop over each MIMA-hsa pair and save in hash
while (<$mRNA_file>){
  my $line = $_;
  my @split = split(";", $line);
	$hash{$split[0]} = $split[1];

}



#Create file
if(!open(WRITE,">pat1.csv")){
die "Could not open file for writing\n";
}


#foreach key change to its corresponding value
foreach (keys %hash){
$firstLine =~ s/$_/$hash{$_}/;

}

print (WRITE $firstLine);

#Dissociate filehandles
close(WRITE);
close $mRNA_file;
