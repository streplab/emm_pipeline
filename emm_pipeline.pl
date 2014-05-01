#!/usr/bin/perl

#First argument is the Fasta database file, second argument is the results from SRST2, should start in the same directory as the scores file.
use File::Glob;
use List::Util qw(min max);
use List::MoreUtils qw{ any };
use Getopt::Long;
use Cwd 'abs_path';
use Cwd;

#Get path to script
my $script_path = abs_path($0);
my $script_dir = substr($script_path, 0, rindex($script_path, "/"));
my $working_dir = getcwd();

my $fastq_directory;
my $scoresName;
my $fasta_input;
my $definitions_file;
my $gene_fasta;
my $forward;
my $reverse;

my %opt=qw();
GetOptions(\%opt, "help|h!", "fastq_directory:s", "scoreName:s", "emm_db:s", "emm_definitions:s", "mga_genes:s", "forward:s", "reverse:s");

if (exists($opt{"help"})) {
	print "**********************************************\nImplementation of the emm pipeline:\n**********************************************\nperl emm_pipeline.pl --fastq_directory /path/to/fastq/directory --scoreName Scores_output_name --emm_db emm.fasta -- emm_definitions emm_reference.txt --mga_genes scpA_mga_mrp_enn.fasta\n\n";
	print "--fastq_directory\tPath to directory containing paired-end fastq files.\n\t\t\tMust be full path to directory, please do not use '.'\n\t\t\tor '..' to declare path\n\t\t\tIf no path is given, the current working\n\t\t\tdirectory is used\n";
	print "--scoreName\t\tName of SRST2 results file\n\t\t\t[optional: default name 'Results']\n";
	print "--emm_db\t\tMultifasta file containing the emm database\n\t\t\t(emm.fasta)\n\t\t\t[If none is provided, emm.fasta is looked\n\t\t\tfor in the directory containing the script]\n";
	print "--emm_definitions\tText file containing the definitions for the\n\t\t\temm database file\n\t\t\t(emm_reference.txt)\n\t\t\t[If none is provided, emm_reference.txt is looked\n\t\t\tfor in the directory containing the script]\n";
	print "--mga_genes\t\tMultifasta file containing mga, mrp, enn, scpA genes\n\t\t\t(mga_mrp_scpA_enn.fasta)\n\t\t\t[If none is provided, mga_mrp_scpA_enn.fasta is looked\n\t\t\tfor in the directory containing script]\n";
	print "--forward\t\tIndicator delimiting the forward reads file for\n\t\t\tpaired read fastq files\n\t\t\t[optional: default '_R1']\n";
	print "--reverse\t\tIndicator delimiting the reverse reads file for\n\t\t\tpaired read fastq files\n\t\t\t[optional: default '_R2']\n\n";
	print "Note: We recommend using paired end reads of at least 100nt in length.\nWe have not tested the efficiency of the pipeline with reads shorter than 80nt.  This pipeline requires paired-end reads.";
	exit;
}

if(exists($opt{"fastq_directory"})){
	$fastq_directory=$opt{"fastq_directory"};
}
else{
	$fastq_directory=$working_dir;
	my $any_fastqs=glob("*.fastq");

	if(!defined($any_fastqs)){	
		print "Please provide the full path to the directory containing fastq files to use. See help file [--help]";
		exit;
	}
}

if(exists($opt{"scoreName"})){
	$scoresName=$opt{"scoreName"};
}
else{
	$scoresName = "Results";
}

if(exists($opt{"emm_db"})){
	$fasta_input=$opt{"emm_db"};
}
else{
	$fasta_input=$script_dir."/emm.fasta";
	unless(-e $fasta_input){
		print "Please provide the emm database fasta file. See help file [--help]";
		exit;
	}
}

if(exists($opt{"emm_definitions"})){
	$definitions_file=$opt{"emm_definitions"};
}
else{
	$definitions_file=$script_dir."/emm_reference.txt";
	unless(-e $definitions_file){
		print "Please provide the emm definitions file. See help file [--help]";
		exit;
	}
}

if(exists($opt{"mga_genes"})){
	$gene_fasta=$opt{"mga_genes"};
}
else{
	$gene_fasta=$script_dir."/mga_mrp_scpA_enn.fasta";
	unless(-e $gene_fasta){
		print "Please provide the mga, mrp, enn, scpA multifasta file. See help file [--help]";
		exit;
	}
}
if(exists($opt{"forward"})){
	$forward = $opt{"forward"};
}
else{
	$forward = "_R1";
}
if(exists($opt{"reverse"})){
	$reverse = $opt{"reverse"};
}
else{
	$reverse = "_R2";
}

my @info;
my @sortedInfo;
my $secondAllele;
my $secondScore;
my @storeNames;

#Run SRST2
system("srst2_mosaik.py --input_pe $fastq_directory/*.fastq --forward $forward --reverse $reverse --output $scoresName --log --mlst_db $fasta_input --mlst_definitions $definitions_file --aligner mosaik --save_scores") == 0 or die "Can't use python script"; 

#Organize SRST2 output
system("mkdir pdfs");
system("mkdir pileups");
system("mkdir sam");
system("mkdir sam_mod");
system("mkdir sorted_bam");
system("mkdir unsorted_bam");
system("mkdir scores");
system("mv *.pdf ./pdfs");
system("mv *.pileup ./pileups");
system("mv *.sam ./sam");
system("mv *.mod ./sam_mod");
system("mv *.sorted.bam ./sorted_bam");
system("mv *.unsorted.bam ./unsorted_bam");
system("mv *.scores ./scores");

mkdir "Pipeline";
chdir "Pipeline";

#creates a file called reducedScores.txt with the headings Sample, First Allele, First Score, Second Allele, Second Score
open $newfile, ">", "Top3Scores.txt" or die "Can't open the output file!";
print $newfile "Sample \t First Allele \t First Score \t Second Allele \t Second Score \t Third Allele \t Third Score \n";

chdir "..";
chdir "scores";

#Gets all .score files in the directory
my @scoresFile = glob("*.scores");

#Opens each score file and organizes the scores
foreach $file(@scoresFile){
	$count = 0;
	open FILE, "<$file" or die "Can't find input file!";

	#Resets arrays to null	
	@info = ();
	@sortedInfo = ();
	#Splits each line of the score file into its elements
	foreach $line(<FILE>){
		chomp;
		($allele, $scores, $avDepth, $edge1, $edge2, $percCov, $size, $mismatch, $indels, $trun, $neighTrun, $LeastConRate, $LeastConMis, $LeastConDepth, $LeastConPVal) = split("\t", $line);

		#If the percent coverage is over 90, then the score is recorded into the info array
		if($percCov >= 90){
			$info[$count][0] = substr $file, (index($file, "__")+2), ((index($file, "."))-(index($file,"__")+2));
			$info[$count][1] = $allele;
			$info[$count][2] = $scores;
			$count = $count + 1;
		}
	}

print "$info[0][0]\n";

	#Sorts the info array by the scores
	@sortedInfo = sort{$a->[2] <=> $b->[2]}@info;

	#Records the best three scores into the output file
	if($count > 2){
		$i = 0;
		do{
			#looks for an emm type that isn't the same as the top scoring emm-type
			if((substr $sortedInfo[0][1], 7, index($sortedInfo[0][1], ".")-7) eq (substr $sortedInfo[$i][1], 7, index($sortedInfo[$i][1], ".")-7)){
				$i = $i + 1;
				$equal = 1;
			}
			else{
				$equal = 0;
			}
		}while($equal && $i < $count);
		#When a second emm-type is found, sets the second score to it, if not found "NA" is recorded
		if($equal){
			print $newfile $sortedInfo[0][0] . "\t" . $sortedInfo[0][1] . "\t" . $sortedInfo[0][2] . "\tNA\tNA\tNA\tNA\n";
		}
		else{
			#Once a second emm-type is found, continues to look for a third
			$secondAllele = $sortedInfo[$i][1];
			$secondScore = $sortedInfo[$i][2];
			$i = $i + 1;
			$equal2 = 1;
			while($equal2 && $i < $count){
				if((substr $secondAllele, 7, index($secondAllele, ".")-7) eq (substr $sortedInfo[$i][1], 7, index($sortedInfo[$i][1], ".")-7) or (substr $sortedInfo[0][1], 7, index($sortedInfo[0][1], ".")-7) eq (substr $sortedInfo[$i][1], 7, index($sortedInfo[$i][1], ".")-7)){
					$i = $i + 1;
					$equal2 = 1;
				}
				else{
					$equal2 = 0;
				}
			}
			if($equal2){
				print $newfile $sortedInfo[0][0] . "\t" . $sortedInfo[0][1] . "\t" . $sortedInfo[0][2] . "\t" . $secondAllele . "\t" . $secondScore . "\tNA\tNA\n";
			}
			else{
				print $newfile $sortedInfo[0][0] . "\t" . $sortedInfo[0][1] . "\t" . $sortedInfo[0][2] . "\t" . $secondAllele . "\t" . $secondScore . "\t" . $sortedInfo[$i][1] . "\t" . $sortedInfo[$i][2] . "\n";
			}
		}
	}
	elsif($count == 2){
		if((substr $sortedInfo[0][1], 7, index($sortedInfo[0][1], ".")-7) eq (substr $sortedInfo[1][1], 7, index($sortedInfo[1][1], ".")-7)){
			print $newfile $sortedInfo[0][0] . "\t" . $sortedInfo[0][1] . "\t" . $sortedInfo[0][2] . "\tNA\tNA\tNA\tNA\n";
		}
		else{
			print $newfile $sortedInfo[0][0] . "\t" . $sortedInfo[0][1] . "\t" . $sortedInfo[0][2] . "\t" . $sortedInfo[1][1] . "\t" . $sortedInfo[1][2] . "\tNA\tNA\n";
		}
	}
	elsif($count == 1){
		print $newfile $sortedInfo[0][0] . "\t" . $sortedInfo[0][1] . "\t" . $sortedInfo[0][2] . "\tNA\tNA\tNA\tNA\n";
	}
	else{
		print $newfile (substr $file, (index($file, "__")+2), 7) . "\tNA\tNA\tNA\tNA\tNA\tNA\n";
	}
}

chdir "..";
chdir "Pipeline";

close $newfile;

open $scoresFile, "<", "Top3Scores.txt" or die "Can't open the input file!";
open $nextNew, ">", "top3_lessThan10.txt" or die "Can't open output file!";
print $nextNew "Sample\tAllele1\tAllele2\tAllele3\n";

#reads in the alleles of those that have >1 match from the generated output file.
$count2 = 0;
foreach $line2(<$scoresFile>){
	chomp;
	($sample, $allele1, $score1, $allele2, $score2, $allele3, $score3) = split("\t", $line2);

	if($allele3 ne "NA" and $score3 < 10){
		$blasted[$count2][0] = $sample;
		$blasted[$count2][1] = $allele1;
		$blasted[$count2][2] = $allele2;
		$blasted[$count2][3] = $allele3;
		if($count2 != 0){
			print $nextNew "$blasted[$count2][0]\t$blasted[$count2][1]\t$blasted[$count2][2]\t$blasted[$count2][3]\n";
		}
		$count2 = $count2 + 1;
	}
	elsif($allele2 ne "NA" and $score2 < 10){
		$blasted[$count2][0] = $sample;
		$blasted[$count2][1] = $allele1;
		$blasted[$count2][2] = $allele2;
		$blasted[$count2][3] = "";
		print $nextNew "$blasted[$count2][0]\t$blasted[$count2][1]\t$blasted[$count2][2]\n";
		$count2 = $count2 + 1;
	}
}

open $EmmLikeFile, ">", "Emm_like_genes.txt" or die "Can't open the output file!";
print $EmmLikeFile "Sample\tGene_Type\tEmm_Label\n";

open $BlastResults, ">", "Blast_Results.txt" or die "Can't open the output file!";
print $BlastResults "Sample\tGene\tPosition\tScore\tE-Value\tScaffold\n";

#Turns the multifasta into a directory of individual fastas
mkdir "fastas";
chdir "fastas";
system("split_fasta.pl $fasta_input");
chdir "..";

mkdir "multiFastas";

open $final_emm, ">", "final_emm_Results.txt" or die "Can't open the output file!";
print $final_emm "Sample\temm\n";

open $mapping_blast, ">", "Blast_map.txt" or die "Can't open the output file!";

#Analyses each with more than one allele
for($n = 1; $n < $count2; $n++){
	##De novo assembly of those with >1 match via the A5 pipeline
	mkdir "$blasted[$n][0]";
	chdir "$blasted[$n][0]";
	system("a5_pipeline.pl $fastq_directory/$blasted[$n][0]"."*$forward"."*.fastq $fastq_directory/$blasted[$n][0]"."*$reverse"."*.fastq $blasted[$n][0] 2>&1 | tee -a log.txt");
	chdir "..";

	#Creates a blast database from the final scaffolds of de novo assembled sequences
	mkdir "BLAST_database";
	chdir "BLAST_database";
	system("formatdb -i ../$blasted[$n][0]/$blasted[$n][0].final.scaffolds.fasta -p F -o T -n $blasted[$n][0].database");
	chdir "..";

	#Creates a multifasta file with includes mga1, mga2, mrp, scpA, emmMatch1, emmMatch2, emmMatch3 (if any)
	chdir "multiFastas";
	open $newfasta, ">", "$blasted[$n][0]_emm.fasta" or die "Can't open the output file!";
	if($blasted[$n][3] eq ""){
		system("cat $gene_fasta ../fastas/$blasted[$n][1].fa ../fastas/$blasted[$n][2].fa > $blasted[$n][0]_emm.fasta");
	}
	else{
		system("cat $gene_fasta ../fastas/$blasted[$n][1].fa ../fastas/$blasted[$n][2].fa ../fastas/$blasted[$n][3].fa > $blasted[$n][0]_emm.fasta");
	}
	close $newFasta;
	chdir "..";

	#Blasts the scaffolds to the created multifasta
	mkdir "BLASTS";
	chdir "BLASTS";
	system("blastall -p blastn -d ../BLAST_database/$blasted[$n][0].database -i ../multiFastas/$blasted[$n][0]_emm.fasta -o $blasted[$n][0]_BLAST");

	##READ BLAST FILE, record blast scores and first positions
	open BLAST, "<$blasted[$n][0]_BLAST" or die "can't find input file!";

	my @bScore = ();
	my @savedlines = ();
	my @blasPos = ();
	my @scafPos = ();
	my @eVal = ();
	my $emm = ();
	my $mga = ();
	my $emmMRP = ();
	my @emmMRPwhich = ();
	my $numemm = ();
	my $mrpRecord = ();
	$count3 = 0;
	foreach $line3(<BLAST>){
		chomp;
		$savedlines[$count3] = $line3;
		$count3++;
	}
	$z = 0;
	for($r = 0; $r < $count3; $r++){
		#Check for the results, read in top score of result, look for the first position and save, continue to next result
		if($savedlines[$r] eq "Sequences producing significant alignments:                      (bits) Value\n"){
			$bScore[$z] = (split(/\s+/, $savedlines[$r+2]))[1];
			$eVal[$z] = (split(/\s+/, $savedlines[$r+2]))[2];
			$scafPos[$z] = substr($savedlines[$r+2],0,index($savedlines[$r+2],"_"));
			#check next 30 lines for the Sbjct position
			$sbj = 0;
			$num = $r;
			do{
				if(substr($savedlines[$num], 0, 6) eq "Sbjct:"){
					$blastPos[$z] = substr($savedlines[$num], 7, index($savedlines[$num], " ", 7)-7);
					$sbj = 1;
				}
				else{
					$sbj = 0;
					$num++;
				}
			}while($sbj == 0);
			$z++;
		}
	}


	#Picks mga gene for reference (either mga1 or mga2)
	if($bScore[0] > $bScore[1]){
		$mga = 0;
	}
	else{
		$mga = 1;
	}

	print $BlastResults "$blasted[$n][0]\tmga" . ($mga+1) . "\t$blastPos[$mga]\t$bScore[$mga]\t$eVal[$mga]\t$scafPos[$mga]\n";
	print $BlastResults "$blasted[$n][0]\tmrp\t$blastPos[2]\t$bScore[2]\t$eVal[2]\t$scafPos[2]\n";
	print $BlastResults "$blasted[$n][0]\tenn\t$blastPos[3]\t$bScore[3]\t$eVal[3]\t$scafPos[3]\n";
	print $BlastResults "$blasted[$n][0]\tscpA\t$blastPos[4]\t$bScore[4]\t$eVal[4]\t$scafPos[4]\n";
	print $BlastResults "$blasted[$n][0]\t$blasted[$n][1]\t$blastPos[5]\t$bScore[5]\t$eVal[5]\t$scafPos[5]\n";
	print $BlastResults "$blasted[$n][0]\t$blasted[$n][2]\t$blastPos[6]\t$bScore[6]\t$eVal[6]\t$scafPos[6]\n";
	if($z == 8){
		print $BlastResults "$blasted[$n][0]\t$blasted[$n][3]\t$blastPos[7]\t$bScore[7]\t$eVal[7]\t$scafPos[7]\n";
	}

	#check if mrp exists, if yes, flag and record any emms that are mrp
	$emmMRP = 0;
	#print $z;
	if($bScore[2] > 300){
		for($numemm = 5; $numemm < $z; $numemm++){
			if((($blastPos[$mga] < $blastPos[$numemm]) and ($blastPos[$numemm] < $blastPos[2])) or ($blastPos[2] < $blastPos[$numemm]) and ($blastPos[$numemm] < $blastPos[$mga])){
				$emmMRPwhich[$emmMRP] = $numemm;
				$emmMRP++;
				$mrpRecord = ($numemm - 4);
				print $EmmLikeFile ($blasted[$n][0] . "\tmrp\t" . $blasted[$n][$mrpRecord] . "\n");
			}
		}
	}
	
	###Get emm
	#If there are only two emms and 1 is mrp, pick the other
	if($z == 7 and $emmMRP == 1){
		if($emmMRPwhich[0] == 5){
			$emm = 6;
		}
		else{
			$emm = 5;
		}
	}
#If there are two emms and neither are mrp, checks to see which one is correct, if both map to emm region, takes the one with the best score, if they map to two separate regions, records potential enn gene
	elsif($z == 7 and $emmMRP == 0){
		if(abs($blastPos[5] - $blastPos[6]) < 1000){
			if($bScore[5] > $bScore[6]){
				$emm = 5;
			}
			else{
				$emm = 6;
			}
		}
		elsif(abs($blastPos[$mga] - $blastPos[5]) < abs($blastPos[$mga] - $blastPos[6])){
			$emm = 5;
			print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][2]\n";
		}
		else{
			$emm = 6;
			print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][1]\n";
		}
	}
	#If there are three emms and one is mrp, checks to see which one is correct, if both map to emm region, takes the one with the best score, if they map to two separate regions, records potential enn gene
	elsif($z == 8 and $emmMRP == 1){
		if($emmMRPwhich[0] == 5){
			if(abs($blastPos[6] - $blastPos[7]) < 1000){
				if($bScore[6] > $bScore[7]){
					$emm = 6;
				}
				else{
					$emm = 7;
				}
			}
			elsif(abs($blastPos[$mga] - $blastPos[6]) < abs($blastPos[$mga] - $blastPos[7])){
				$emm = 6;
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][3]\n";
			}
			elsif(abs($blastPos[$mga] - $blastPos[7]) < abs($blastPos[$mga] - $blastPos[6])){
				$emm = 7;
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][2]\n";
			}
		}
		elsif($emmMRPwhich[0] == 6){
			if(abs($blastPos[5] - $blastPos[7]) < 1000){
				if($bScore[5] > $bScore[7]){
					$emm = 5;
				}
				else{
					$emm = 7;
				}
			}
			elsif(abs($blastPos[$mga] - $blastPos[5]) < abs($blastPos[$mga] - $blastPos[7])){
				$emm = 5;
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][3]\n";
			}
			elsif(abs($blastPos[$mga] - $blastPos[7]) < abs($blastPos[$mga] - $blastPos[5])){
				$emm = 7;
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][1]\n";
			}
		}
		elsif($emmMRPwhich[0] == 7){
			if(abs($blastPos[6] - $blastPos[5]) < 1000){
				if($bScore[6] > $bScore[5]){
					$emm = 6;
				}
				else{
					$emm = 5;
				}
			}
			elsif(abs($blastPos[$mga] - $blastPos[5]) < abs($blastPos[$mga] - $blastPos[6])){
				$emm = 5;
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][2]\n";
			}
			elsif(abs($blastPos[$mga] - $blastPos[6]) < abs($blastPos[$mga] - $blastPos[5])){
				$emm = 6;
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][1]\n";
			}
		}
	}
	#If there are three emms and none are mrp, checks to see which one is correct, if all map to emm region, takes the one with the best score, if they map to two separate regions, records potential enn genes
	elsif($z == 8 and $emmMRP == 0){
		#Records gene closest to the mga region as emm
		@allEmms = (abs($blastPos[$mga] - $blastPos[5]), abs($blastPos[$mga] - $blastPos[6]), abs($blastPos[$mga] - $blastPos[7]));
		$emmMin = min @allEmms;
		if($emmMin == (abs($blastPos[$mga] - $blastPos[5]))){$emm = 5;}
		elsif($emmMin == (abs($blastPos[$mga] - $blastPos[6]))){$emm = 6;}
		elsif($emmMin == (abs($blastPos[$mga] - $blastPos[7]))){$emm = 7;}
		else{$emm = 0;}

		#double checks that two or more genes aren't mapping to emm, if so, takes the one with the best score		
		if($emm == 5){
			if(abs($blastPos[5] - $blastPos[6]) < 1000){
				if($bScore[6] > $bScore[5]){
					$emm = 6;
					if(abs($blastPos[6] - $blastPos[7]) <1000){
						if($bScore[7] > $bScore[6]){
							$emm = 7;
						}
					}
					else{
						print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][3]\n";
					}
				}
				elsif(abs($blastPos[5] - $blastPos[7]) <1000){
					if($bScore[7] > $bScore[5]){
						$emm = 7;
					}
				}
				else{
					print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][3]\n";
				}
			}
			elsif(abs($blastPos[5] - $blastPos[7]) < 1000){
				if($bScore[7] > $bScore[5]){
					$emm = 7;
				}
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][2]\n";
			}
			else{
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][2]\n";
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][3]\n";
			}
		}
		#double checks that two or more genes aren't mapping to emm, if so, takes the one with the best score	
		elsif($emm == 6){
			if(abs($blastPos[5] - $blastPos[6]) < 1000){
				if($bScore[5] > $bScore[6]){
					$emm = 5;
					if(abs($blastPos[5] - $blastPos[7]) < 1000){
						if($bScore[7] > $bScore[5]){
							$emm = 7;
						}
					}
					else{
						print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][3]\n";
					}
				}
				elsif(abs($blastPos[7] - $blastPos[6]) <1000){
					if($bScore[7] > $bScore[6]){
						$emm = 7;
					}
				}
				else{
					print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][3]\n";
				}
			}
			elsif(abs($blastPos[7] - $blastPos[6]) < 1000){
				if($bScore[7] > $bScore[6]){
					$emm = 7;
				}
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][1]\n";
			}
			else{
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][1]\n";
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][3]\n";
			}
		}
		#double checks that two or more genes aren't mapping to emm, if so, takes the one with the best score	
		elsif($emm == 7){
			if(abs($blastPos[7] - $blastPos[6]) < 1000){
				if($bScore[6] > $bScore[7]){
					$emm = 6;
					if(abs($blastPos[6] - $blastPos[5]) < 1000){
						if($bScore[5] > $bScore[6]){
							$emm = 5;
						}
					}
					else{
						print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][1]\n";
					}
				}
				elsif(abs($blastPos[5] - $blastPos[7]) <1000){
					if($bScore[7] < $bScore[5]){
						$emm = 5;
					}
				}
				else{
					print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][1]\n";
				}
			}
			elsif(abs($blastPos[5] - $blastPos[7]) < 1000){
				if($bScore[5] > $bScore[7]){
					$emm = 5;
				}
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][2]\n";
			}
			else{
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][1]\n";
				print $EmmLikeFile "$blasted[$n][0]\tenn\t$blasted[$n][2]\n";
			}
		}
	}
	#If there are three potential genes and two are mrp, records the third gene as emm
	elsif($z == 8 and $emmMRP == 2){
		if(grep(!/^5/, @emmMRPwhich)){
			$emm = 5;
		}
		elsif(grep(!/^6/,@emmMRPwhich)){
			$emm = 6;
		}
		else{
			$emm = 7;
		}
	}

	#Records results into final results sheet
	print $final_emm $blasted[$n][0] . "\t" . substr($blasted[$n][$emm-4],4) . "\n";

	$storeNames[$n] = $blasted[$n][0];
	
	chdir "..";

	#Calculate the distances between each of the alleles present and mga
	@Dist = ();
	@DistSort = ();
	$DistNum = 0;
	$Dist[$DistNum][0] = "mga";
	$Dist[$DistNum][1] = 0;
	$DistNum++;
	if($bScore[2] > 300){
		$Dist[$DistNum][0] = "mrp";
		$Dist[$DistNum][1] = abs($blastPos[$mga] - $blastPos[2]);
		$DistNum++;
	}
	if($eVal[3] < 0.0001){
		$Dist[$DistNum][0] = "enn";
		$Dist[$DistNum][1] = abs($blastPos[$mga] - $blastPos[3]);
		$DistNum++;
	}
	if($scafPos[4] == $scafPos[$mga]){
		$Dist[$DistNum][0] = "scpA";
		$Dist[$DistNum][1] = abs($blastPos[$mga] - $blastPos[4]);
		$DistNum++;
	}
	$Dist[$DistNum][0] = substr($blasted[$n][1],4);
	$Dist[$DistNum][1] = abs($blastPos[$mga] - $blastPos[5]);
	$DistNum++;
	$Dist[$DistNum][0] = substr($blasted[$n][2],4);
	$Dist[$DistNum][1] = abs($blastPos[$mga] - $blastPos[6]);
	$DistNum++;
	if($z == 8){
		$Dist[$DistNum][0] = substr($blasted[$n][3],4);
		$Dist[$DistNum][1] = abs($blastPos[$mga] - $blastPos[7]);
		$DistNum++;
	}

	#Sorts the alleles by distance from mga	
	@DistSort = sort{$a->[1] <=> $b->[1]}@Dist;

	#Prints the alleles with a dash representing 100bp in between each allele
	print $mapping_blast "$blasted[$n][0]: $DistSort[0][0]";
	$dash = int(($DistSort[1][1]/100)+0.5);
	for($d = 0; $d < $dash; $d++){
		print $mapping_blast "_";
	}
	print $mapping_blast $DistSort[1][0];
	for($b = 2; $b < $DistNum; $b++){
		$dash = int((($DistSort[$b][1]-$DistSort[($b-1)][1])/100)+0.5);
		if($dash == 0){
			print $mapping_blast "/";
		}
		for($d = 0; $d < $dash; $d++){
			print $mapping_blast "_";
		}
		print $mapping_blast $DistSort[$b][0];
	}
	print $mapping_blast "\n";
}

close $nextNew;
close $EmmLikeFile;
close $mapping_blasts;

chdir "..";

my $ResultsName = glob("$scoresName*.txt");

#Reads SRST2 results page and records any results that were not further analyzed
open SRST2, "<$ResultsName" or die "Can't open input file!";
$resultCount = 0;
foreach $results(<SRST2>){
	chomp;
	@reading = split(/\t/, $results);
	if($resultCount == 0){
	}
	elsif((any{/$reading[0]/} @storeNames) == 0){
		print $final_emm "$reading[0]\t$reading[1]\n";
	}

	$resultCount++;
}

#Deletes creates fastas and multifastas
system("rm -rf fastas/");
system("rm -rf multiFastas/");
