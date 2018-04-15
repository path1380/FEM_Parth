#!/usr/apps/bin/perl
#This script computes the log error vs log N of the 
#basic method with LU decomposition

$cmdFile="./templates/Input_function.f90.template";
$outFile="./src/Input_function.f90";

$n_pts = 10;

for($i = 2;$i <= 500; $i *= 2){
	open(FILE,"$cmdFile") || die "cannot open file $cmdFile!" ;
	open(OUTFILE,"> $outFile") || die "cannot open file!" ;
	while($line = <FILE>){
		$line =~ s/\bNNNN\b/$i/;
                print OUTFILE $line;
	}
	close( OUTFILE );
	close( FILE );
	system("make compile");
	system("make build");
	system("make run");
	system("make clean");
}
