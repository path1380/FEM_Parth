#!/usr/apps/bin/perl
#This script plots the convergence rate with the
#number of divisions. Run with "make clean".

$cmdFile="./templates/Input_function.f90.template";
$outFile="./src/Input_function.f90";
$cmdFile1="output.txt";
$outFile1="output_mod.txt";
$cmdFile2="output_mod.txt";
$outFile2="grid_ind_plot.txt";

$n_pts = 10;

open(OUTFILE_FINAL,"> $outFile2") || die "cannot open file!" ;

for($i = 2;$i <= 32; $i = $i + 2){
	open(FILE,"$cmdFile") || die "cannot open file $cmdFile!" ;
	open(OUTFILE,"> $outFile") || die "cannot open file!" ;
	while($line = <FILE>){
		$line =~ s/\bNNNN\b/$i/;
                print OUTFILE $line;
	}
	#print $n_lines_output,"\n";
	close( OUTFILE );
	close( FILE );
	system("make compile");
	system("make build");
	system("make run");

	$j=0;
	$n_lines_output = 0;

	open(FILE,"$cmdFile1") || die "cannot open file $cmdFile1!" ;
	open(OUTFILE,"> $outFile1") || die "cannot open file!" ;
	
	while($line = <FILE>){
		if ($j < 10) {
			$line =~ s/  $j KSP Residual norm //;
		} elsif ($j >= 10) {
			$line =~ s/ $j KSP Residual norm //;
		}
		$j = $j+1;
		print OUTFILE $line;
		$n_lines_output = $n_lines_output + 1;
	}
	close( OUTFILE );
	close( FILE );

	$j = 0;

	open(FILE,"$cmdFile2") || die "cannot open file $cmdFile2!" ;
	#open(OUTFILE,"> $outFile2") || die "cannot open file!" ;
	

	while($line = <FILE>){
		$j = $j + 1;
		if ($j == $n_lines_output - 1) {
			$err_k_minus_1 = $line;
		} elsif ($j == $n_lines_output) {
			$err_k = $line;
		}
	}

	$conv_rate = $err_k/$err_k_minus_1;
	#print $conv_rate,"\n";

	print OUTFILE_FINAL $i," ",$conv_rate,"\n";

	close( FILE );
	
	system("rm output.txt output_mod.txt");

	#system("make clean");
}
close( OUTFILE_FINAL );
