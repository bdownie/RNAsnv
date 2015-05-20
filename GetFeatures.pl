#!/usr/bin/perl

use strict; 
use warnings;
use Getopt::Long;



my $vcf_file;
my $gtf_file;
my $intron_bed;
my $is_repeats = "";
my $rna_edit_db;
my $reference;
my $dna_file;
my $SEQLEN = 80;
my $OPP_SEQLEN = 20;
my $dbsnp;
my $CONFIG="";

my %in_dna;
my $min_dna_depth = 10;

GetOptions('v=s'=>\$vcf_file,'g=s'=>\$gtf_file, 'a=s'=>\$is_repeats, 'e=s'=>\$rna_edit_db,'i=s'=>\$intron_bed,'r=s'=>\$reference,'d=s'=>\$dna_file,'c=s'=>\$CONFIG,'dbsnp=s'=>\$dbsnp);

unless ($gtf_file) { print STDERR "Usage: -v <vcf> -g <gtf> -i <intron bed> -a <Alu track> -e <verified edit db> -r <ref> -d <dna> -dbsnp <dbsnp>\n"; exit; }

unless (-e $vcf_file) {
	print STDERR "Couldn't find vcf: $vcf_file\n";
	print STDERR "Usage: -v <vcf> -g <gtf> -i <intron bed> -a <Alu track> -e <verified edit db> -r <ref> -d <dna>\n"; 
	exit;
}
unless ($gtf_file) {
	print STDERR "Warning: No annotation file given, annotation-based features will be empty.\n"; 
}
unless ($is_repeats) {
	print STDERR "Warning: No repeats file given, repeats features will be empty.\n"; 
}
unless ($rna_edit_db) {
	print STDERR "Warning: No RNA edit DB file given, edit feature will be empty.\n"; 
}
unless ($intron_bed) {
	print STDERR "Warning: No intron BED file given, intron-based features will be empty\n";
}
unless ($reference) {
	print STDERR "Warning: No reference file given, reference-based features will be empty\n";
}
unless ($dna_file) {
	print STDERR "Warning: No WGS-SNV file given, SNV-based features will be empty\n";
}

my %CONFIG_HASH;
if ($CONFIG) {
        open C, $CONFIG or die "Couldn't open config file $CONFIG $!";
        while (my $line = <C>) {
                $line =~ s/#.*//;
                $line =~ s/^\[.*//;
                chomp $line;
                if ($line =~ /(.*)=(.*)/) {
                        my $first = uc($1);
                        my $second = $2;
                        $first =~ s/\s+//g;
                        $second =~ s/\s+//g;
                        $CONFIG_HASH{$first} = $second;
                }
        }
}
my $bedtools_loc;
if ($CONFIG_HASH{'BEDTOOLS'}) { 
	$bedtools_loc = $CONFIG_HASH{'BEDTOOLS'};
}
else {
	$bedtools_loc = `which bedtools`;
	chomp $bedtools_loc;
}




my %closest_class;
my %repeat_dist;
my %repeat_dist_long;
my %repeat_classes;
my %site_class;
my $gene_id;
my %exons;
my %biotypes;
my $num_fields = 0;
my %gene_in_exon;
my %gene_strings;
my %vcf_lines;
my %mafs;
my %dps;
my %seqs;
my %d_values;

my %mfes;
my %meas;
my %qds;
my %ht_scores;
my %bqrs_s;
my %rprs_s;
my %var_quals;
my %med_edits;

my $gene_biotype;
my %nearest_splice_site;
my $site_type;
my @exon_ends;
my $last_site;

my %site_types;
my %gene_ids;
my %splice_dists;
my %orig_vcf_lines;
my %alt_read_count;

my %maxreadenddist;
my %medianreadenddist;

my %mutation_type;
my %SP_AVG_VAL;

my %snpclust;

my $progress = 0;
use Bio::DB::Fasta;
my $db;
if ($reference) { 
	$db = Bio::DB::Fasta->new($reference);
}
my @site_ids = `grep -v "^#" $vcf_file | awk '{print \$1 ":" \$2}'`;

my %VDB;
my %SGB;
my %MQB;
my %MQSB;
my %BQB;
my %MQ0F;
my %ICB;
my %HOB;
my %AC;
my %AN;

my %fivePrimeNucleotide;
my %threePrimeNucleotide;

open F, $vcf_file;
while (my $line = <F>) { 
	next if ($line =~ /^#/);
	chomp $line; 

	my @sp = split /\s+/, $line;
	if ($num_fields == 0) { 
		$num_fields = 1+$#sp;
	}

	my $maf = 0;
	my $dp = 0;
	my $i = 9; 
	if ($line =~ /DP4=(\d+),(\d+),(\d+),(\d+);/) {
		my $first = $1; 
		$first += $2;
		my $second = $3;
		$second += $4;
		$dp += $first + $second;
		$maf = $second;
	}
	elsif ($sp[8] =~ /:AD:/) { 
		for (; $i <= $#sp; $i++) { 
			if ($sp[$i] eq "./.") {
			}
			elsif ($sp[$i] =~ /:/) {
				my @sp2 = split ":", $sp[$i];
				my ($first,$second) = split ",", $sp2[1];
				if (defined($second)) { 
					$dp += $first + $second;
					$maf += $second;
				}
			}
			else { 
			}
		}
	}
	else { die "Couldn't calculate alternate read frequency (No DP4 Info field or AD format"; }
	
	unless ($dp) { die "$line\n"; }
	my $chr = $sp[0];
	my $site = $sp[1];
	my $site_string = "$chr:$site";

	if ($sp[8] =~ /SP$/) { 
		my $val = 0;
		for (my $j = 9; $j <= $#sp; $j++) { 
			my @tmparray = split ":", $sp[$j];
			$val += $tmparray[$#tmparray];
		}
		my $denom = $#sp - 8;
		if ($denom < 1) { print STDERR $line; exit; }

		$SP_AVG_VAL{$site_string} = int(100*$val/$denom)/100;
	}

	$mafs{$site_string}= int(100*$maf/$dp)/100;
	$dps{$site_string} = $dp;


	my $seq;
	my ($a, $b); 

	my $qd = "NA";
	my $ht_score = "NA";
	my $bqrs = "NA";
	my $rprs = "NA";
	my $var_qual = $sp[5];
	my $median_edit = -1;
	if ($line =~ /MedianEditDist=([\d\.]+?)[;\s]+/) { 
		$median_edit = $1;
	}
	if ($line =~ /MedianReadEndDist=(\d+?)[;\s]+/) { 
		$medianreadenddist{$site_string} = $1;
	}
	if ($line =~ /MaxReadEndDist=(\d+?)[;\s]+/) { 
		$maxreadenddist{$site_string} = $1;
	}
	if ($line =~ /SNVClust=(\d+?)[;\s]+/) { 
		$snpclust{$site_string} = $1;
	}
	if ($line =~ /AltReadCount=(\d+?)[;\s]+/) { 
		$alt_read_count{$site_string} = $1;
	}
	if ($line =~ /VDB=(\S+?)[;\s]+/) { 
		$VDB{$site_string} = $1;
	}
	if ($line =~ /SGB=(\S+?)[;\s]+/) { 
		$SGB{$site_string} = $1;
	}
	if ($line =~ /MQB=(\S+?)[;\s]+/) { 
		$MQB{$site_string} = $1;
	}
	if ($line =~ /MQSB=(\S+?)[;\s]+/) { 
		$MQSB{$site_string} = $1;
	}
	if ($line =~ /BQB=(\S+?)[;\s]+/) { 
		$BQB{$site_string} = $1;
	}
	if ($line =~ /MQ0F=(\S+?)[;\s]+/) { 
		$MQ0F{$site_string} = $1;
	}
	if ($line =~ /ICB=(\S+?)[;\s]+/) { 
		$ICB{$site_string} = $1;
	}
	if ($line =~ /HOB=(\S+?)[;\s]+/) { 
		$HOB{$site_string} = $1;
	}
	if ($line =~ /AC=(\S+?)[;\s]+/) { 
		$AC{$site_string} = $1;
	}
	if ($line =~ /AN=(\S+?)[;\s]+/) { 
		$AN{$site_string} = $1;
	}

	my $ref = $sp[3];
	
	if ($reference) { 
		my $up_site = $db->seq($chr,($site - 1) => ($site - 1));
		my $down_site = $db->seq($chr,($site + 1) => ($site + 1));

		# Could be reverse strand transcribed
		if (($ref eq "T") || ($ref eq "G")) { 
			$up_site =~ tr/ACGT/TGCA/;
			$down_site =~ tr/ACGT/TGCA/;
			my $tmp = $up_site;
			$up_site = $down_site;
			$down_site = $tmp;
		}
		$fivePrimeNucleotide{$site_string} = $up_site;
		$threePrimeNucleotide{$site_string} = $down_site;
	}


	$sp[4] =~ s/,.*//;

	$mutation_type{$site_string} = $sp[3] . $sp[4];

	$qds{$site_string} = $qd;
	$ht_scores{$site_string} = $ht_score;
	$bqrs_s{$site_string} = $bqrs;
	$rprs_s{$site_string} = $rprs;
	$var_quals{$site_string} = $var_qual;
	$med_edits{$site_string} = $median_edit;
	$orig_vcf_lines{$site_string} = $line;
	$orig_vcf_lines{$site_string} = $line;
}

my %closest_vars;
my %closest_vars_mutation;
my %closest_var_identical;
my %closest_var_also_edit;

my $cmd = "$bedtools_loc closest -t first -io -a $vcf_file -b $vcf_file > closest1.vcf.vcf";
my $count1 = `grep -v "^#" $vcf_file | wc -l`;
my $count2 = `wc -l closest1.vcf.vcf`;
#if ($count1 ne $count2) { 
	system($cmd);
#}
open NEAR, "closest1.vcf.vcf";

while (my $line = <NEAR>) { 
	my @sp = split /\s+/, $line;
	my $site_string = "$sp[0]:$sp[1]";

	my $site = $sp[1];
	my $site2 = $sp[$num_fields + 1];
	my $ref = $sp[$num_fields + 3];
	my $var = $sp[$num_fields + 4];
	next unless ($site && $site2);
	my $dist = log2(abs($site - $site2));
	#$dist = log2($dist);

	$closest_vars{$site_string} = $dist;
	$closest_vars_mutation{$site_string} = $ref . $var;

	my $mut_type_closest = "$ref$var";
	my $mut_type = $mutation_type{$site_string};
	if ($mut_type_closest eq $mut_type) { $closest_var_identical{$site_string} = "TRUE"; }

	# Is the closest variant also a candidate edit site?
	if (
		($mut_type eq "GA" && ($mut_type_closest eq "GA" || $mut_type_closest eq "TC")) ||
		($mut_type eq "TC" && ($mut_type_closest eq "GA" || $mut_type_closest eq "TC")) ||
		($mut_type eq "AG" && ($mut_type_closest eq "AG" || $mut_type_closest eq "CT")) ||
		($mut_type eq "CT" && ($mut_type_closest eq "AG" || $mut_type_closest eq "CT"))) { 
			$closest_var_also_edit{$site_string} = "TRUE";
	}
}
close NEAR;


if ($gtf_file) { 
	process_gtf($vcf_file,$gtf_file,\%site_types,\%biotypes,\%gene_ids,\%splice_dists,\%nearest_splice_site);
}

unless ($num_fields) { $num_fields = 9; }

if ($intron_bed) {
	process_intron_bed($vcf_file,$intron_bed,\%splice_dists,\%site_types,$num_fields);
}

my %is_dbsnp;
if ($dbsnp) { 
	process_db_snp($vcf_file, $dbsnp,\%is_dbsnp);
}


if ($is_repeats) { 
	process_repeats($vcf_file, $is_repeats,\%closest_class,\%repeat_dist,\%repeat_dist_long,\%repeat_classes);
}

if ($rna_edit_db) { 
	process_edit_db($vcf_file,$rna_edit_db,\%site_class);
}

if ($dna_file) { 
	process_SNV($vcf_file,$dna_file,\%site_class,\%in_dna);
}


my $last_site_string = "";
print "Chr\tSite\tGeneID\tBiotype\tMutationType\tClass\tMAF\tDepth\tRepeatClass\tClosestClass\tRepeatDist\tRepeatDistLong\tSpliceDist\tDist_to_Nearest_var\tClosestVarType\tClosestVarIdentical\tClosestVarAlsoEdit\tVar_qual\tMedianEdit\tMedianDistFromEnd\tMaxDistFromEnd\tSNVClusterCount\tStrandBias\tAltReadCount\tIndbSNV\tVariantDistanceBias\tSegregationBasedMetric\tMappingQualityBias\tMappingQualityStrandBias\tBaseQualityBias\tMQ0Fraction\tInbreedingCoefficientBinomial\tHOMBias\tAlleleCount\tAlleleNumber\tNuc5prime\tNuc3prime\tSiteClass\n";
foreach my $site_string (@site_ids) { 
	chomp $site_string;
	#next unless($in_dna{$site_string});
	#unless ($in_dna{$site_string}) { $sc = "Unknown"; }

	my ($chr,$site) = split ":", $site_string;
	my $class = $site_types{$site_string};
	my $biotype = $biotypes{$site_string};
	my $gene_id = $gene_ids{$site_string};
	my $maf = $mafs{$site_string};
	my $dp = $dps{$site_string};
	my $splice_dist = $splice_dists{$site_string};
	my $seq = $seqs{$site_string};
	my $mfe = $mfes{$site_string};
	my $mea = $meas{$site_string};
	my $d = $d_values{$site_string};
	my $qd = $qds{$site_string};
	my $ht_score = $ht_scores{$site_string};
	my $bqrs = $bqrs_s{$site_string};
	my $rprs = $rprs_s{$site_string};
	my $var_qual = $var_quals{$site_string};
	my $repeat_class = $repeat_classes{$site_string};
	my $mut_type = $mutation_type{$site_string};
	my $closest_class_type = $closest_class{$site_string};
	my $med_edit = $med_edits{$site_string};
	my $mred = $medianreadenddist{$site_string};
	my $maxred = $maxreadenddist{$site_string};
	my $rep = $repeat_dist{$site_string};
	my $rep_long = $repeat_dist_long{$site_string};
	my $dist = $closest_vars{$site_string};
	my $cvm = $closest_vars_mutation{$site_string};
	my $sc;
	if ($in_dna{$site_string}) { 
		$sc = $site_class{$site_string};
	}
	else { 
		$sc = "NA";
	}
	my $sp = $SP_AVG_VAL{$site_string};
	my $snp_clust = $snpclust{$site_string};
	my $arc = $alt_read_count{$site_string};

	unless ($arc) { $arc = "NA"; }
	unless ($class) { $class = "intergenic"; }
	unless ($biotype) { $biotype = "NA"; }
	unless ($gene_id) { $gene_id = "NA"; }
	unless (defined($dbsnp)) { $dbsnp = "FALSE"; }
	unless ($mut_type) { 
		$mut_type = "NA"; 
		die "$site_string has no mutation";
	}
	unless ($maf) { $maf = "NA"; }
	unless ($dp) { $dp = "NA"; }
	unless ($splice_dist) { $splice_dist = "NA"; }
	unless ($mfe) { $mfe = "NA"; }
	unless ($qd) { $qd = "NA"; }
	unless ($ht_score) { $ht_score = "NA"; }
	unless ($bqrs) { $bqrs = "NA"; }
	unless ($rprs) { $rprs = "NA"; }
	unless ($repeat_class) { $repeat_class = "NA"; }
	unless ($rep) { $rep = 0; }
	unless ($rep_long) { $rep_long = 0; }
	unless ($d) { $d = 0; }
	unless ($mea) { $mea = 0; }
	unless ($dist) { $dist = "NA"; }
	unless ($closest_class_type) { $closest_class_type = "NA"; }
	unless ($med_edit) { $med_edit = -1; }
	unless (defined($sp)) { $sp = "NA"; }

	$mut_type =~ s/,.*//; 

	unless ($cvm) { $cvm = "NA"; }
	$cvm =~ s/,.*//; 
	unless ($snp_clust) { $snp_clust = "NA"; }
	unless ($sc) { $sc = "Unknown"; }

	# Samtools calculated values
	my $vdb = $VDB{$site_string};
	my $sgb = $SGB{$site_string};
	my $mqb = $MQB{$site_string};
	my $mqsb = $MQSB{$site_string};
	my $bqb = $BQB{$site_string};
	my $mq0f = $MQ0F{$site_string};
	my $icb = $ICB{$site_string};
	my $hob = $HOB{$site_string};
	my $ac = $AC{$site_string};
	my $an = $AN{$site_string};

	my $fivepn = $fivePrimeNucleotide{$site_string};
	my $threepn = $threePrimeNucleotide{$site_string};


	unless ($fivepn) { $fivepn = "NA"; }
	unless ($threepn) { $threepn = "NA"; }

	unless ($vdb) { $vdb = "NA"; }
	unless ($sgb) { $sgb = "NA"; }
	unless ($mqb) { $mqb = "NA"; }
	unless ($mqsb) { $mqsb = "NA"; }
	unless ($bqb) { $bqb = "NA"; }
	unless ($mq0f) { $mq0f = "NA"; }
	unless ($icb) { $icb = "NA"; }
	unless ($hob) { $hob = "NA"; }
	unless ($ac) { $ac = "NA"; }
	unless ($an) { $an = "NA"; }

	unless ($maxred) { $maxred = "NA"; }
	unless ($mred) { $mred = "NA"; }

	#next unless($in_dna{$site_string});

	#if ($dbsnp eq "TRUE") { $sc = "NA"; }
	my $cv_same = $closest_var_identical{$site_string};
	unless ($cv_same) { $cv_same = "FALSE"; }

	my $cv_also_edit = $closest_var_also_edit{$site_string};
	unless($cv_also_edit) { $cv_also_edit = "FALSE"; }


	my $out_string = "$chr\t$site\t$gene_id\t$biotype\t$mut_type\t$class\t$maf\t$dp\t$repeat_class\t$closest_class_type\t$rep\t$rep_long\t$splice_dist\t$dist\t$cvm\t$cv_same\t$cv_also_edit\t$var_qual\t$med_edit\t$mred\t$maxred\t$snp_clust\t$sp\t$arc\t$dbsnp\t$vdb\t$sgb\t$mqb\t$mqsb\t$bqb\t$mq0f\t$icb\t$hob\t$ac\t$an\t$fivepn\t$threepn\t$sc\n";
	print $out_string;
	$last_site_string = $site_string;
}


sub process_edit_db {
	my ($vcf_file,$rna_edit_db,$sc_ref) = @_;

	$cmd = "$bedtools_loc intersect -wb -a $vcf_file -b $rna_edit_db | uniq";
	open VER, "-|", $cmd;

	while (my $line = <VER>) { 
		my @sp = split /\s+/, $line;

		my $dna_alt = $sp[$num_fields + 4];
		my $rna_alt = $sp[4];

		if (($sp[3] eq "T" && $sp[4] eq "C") || ($sp[3] eq "A" && $sp[4] eq "G")) { 
			my $site_string = "$sp[0]:$sp[1]";
			$$sc_ref{$site_string} = "Edit";
		}
	}
}

sub process_db_snp {
	my ($vcf_file,$dbsnp,$ref_db) = @_;

	$cmd = "$bedtools_loc intersect -wb -b $vcf_file -a $dbsnp > dbsnp.vcf";
	#system ($cmd);
	open DBSNP, "dbsnp.vcf";

	my $dbsnp_num_fields = 8;
	while (my $line = <DBSNP>) { 
		my @sp = split /\s+/, $line;
		my $site_string = "$sp[0]:$sp[1]";
		my $ref = $sp[$dbsnp_num_fields + 3];
		my $alt = $sp[$dbsnp_num_fields + 4];
	
		my $dbsnp_ref = $sp[3];
		my $dbsnp_alt = $sp[4];
	
		if ($dbsnp_ref eq $ref) { 
			my @sp2 = split ",", $dbsnp_alt;
			foreach my $dbs (@sp2) { 
				if ($dbs eq $alt) { 
					$$ref_db{$site_string} = "TRUE"
				}
			}
		}
	}
	close DBSNP;
}

sub process_SNV {
	my ($vcf_file, $dna_file, $ref_sc,$ref_idna) = @_;

	open DNA, $dna_file;
	while (my $line = <DNA>) { 
		next if ($line =~ /^#/);
		my @sp = split /\s+/, $line;
		my $site_string = "$sp[0]:$sp[1]";
		if ($line =~ /DP=(\d+)/) { 
			if ($1 > $min_dna_depth) { 
				$$ref_idna{$site_string} = 1;
			}
		}
		else { 
			$$ref_idna{$site_string} = 1;
		}
	
		if ($sp[4] ne ".") {
			$sp[4] =~ s/,.*//; 
			next unless $mutation_type{$site_string};
			my ($first,$second) = split "", $mutation_type{$site_string};
			if ($sp[4] eq $second) { 
				$$ref_sc{$site_string} = "SNV";
			}
		}
	}
	close DNA;
}


sub process_repeats {
	#process_repeats($vcf_file, $is_repeats,\%closest_class,\%repeat_dist,\%repeat_dist_long,\%repeat_classes);
	my ($vcf_file, $is_repeats,$ref_cc,$ref_rd,$ref_rdl,$ref_rc) = @_;

	$cmd = "$bedtools_loc closest -io -t first -a $vcf_file -b $is_repeats > closest2.repeat.vcf";
	#my $count = `grep -v "^#" $vcf_file | wc -l`;
	#my $count2 = `touch closest2.repeat.vcf && wc -l closest2.repeat.vcf`;
	#if ($count ne $count2) { 
		system($cmd);
	#}
	open NEAR, "closest2.repeat.vcf";
	
	while (my $line = <NEAR>) { 
		my @sp = split /\s+/, $line;
		my $site_string = "$sp[0]:$sp[1]";

		my $type = $sp[$num_fields + 3];
		$$ref_cc{$site_string} = $type;
	}
	close NEAR;


	$cmd = "$bedtools_loc intersect -wb -a $vcf_file -b $is_repeats | uniq";
	open REP, "-|", $cmd;
	
	while (my $line = <REP>) { 
		my @sp = split /\s+/, $line;

		my $site_string = "$sp[0]:$sp[1]";

		my $site = $sp[1];
		my $a = $sp[$num_fields + 1];
		my $b = $sp[$num_fields + 2];
		my $repeat_class = $sp[$num_fields+3];
		my $val = log2(abs($site - $a));
		my $val2 = log2(abs($site -  $b));
	
		if ($val > $val2) { 
			my $tmp = $val; 
			$val = $val2; 
			$val2 = $tmp;
		}
		if ($repeat_class ne "NotRepeat") { 
			$val = 0 - $val;
			$val2 = 0 - $val2;
		}
	
		$$ref_rd{$site_string} = $val;
		$$ref_rc{$site_string} = $repeat_class;
		$$ref_rdl{$site_string} = $val2;
	}
}


sub process_gtf {
	#process_gtf($vcf_file,$gtf_file,\%site_types,\%biotypes,\%gene_ids,\%splice_dists,\%nearest_splice_site);
	my ($vcf_file,$gtf_file,$ref_site_types,$ref_biotypes,$ref_gene_ids,$ref_splice_dists,$ref_nearest_splice_site) = @_;

	$cmd = "$bedtools_loc intersect -wb -a $vcf_file -b $gtf_file | uniq";
	open GTF, "-|", $cmd;
	while (my $line = <GTF>) { 
		last unless $line;
		chomp $line;
		my @sp = split /\s+/, $line;
	
		my $chr = $sp[0];
		my $site = $sp[1];
		my $site_string = "$chr:$site";

		my $subtract = 0;
		my $is_intron = 0;
		my $is_utr = 0;
		my $type = $sp[2 + $num_fields];
		my $splice_dist = "";
		if ($line =~ /gene_id "(\S+?)";/) { 
			$gene_id = $1;
		}
		$site_type = "intergenic";
	
		my $direction = $sp[$num_fields + 6];
		$splice_dist = 0;
		if ($type eq "gene") {
			if ($line =~ /gene_biotype "(\S+?)";/) { 
				$gene_biotype = $1;
			}
		}
		elsif (($type eq "exon") || ($type eq "CDS")) { 
			$subtract = $sp[$num_fields + 4];
			if ($direction eq "-") { 
				$subtract = $sp[$num_fields + 3];
			}
			$splice_dist = log2(abs($site - $subtract));
			$site_type = "exon";
		}
		elsif ($type eq "UTR") { 
			$site_type = "UTR";
		}
	
	
		if (defined($nearest_splice_site{$site_string})) {
			if ($nearest_splice_site{$site_string} > $splice_dist) {
				$$ref_nearest_splice_site{$site_string} = $splice_dist;
			}
		}
		else { 
			$$ref_nearest_splice_site{$site_string} = $splice_dist;
		}
		
		$$ref_site_types{$site_string} = $site_type;
		$$ref_biotypes{$site_string} = $gene_biotype;
		$$ref_gene_ids{$site_string} = $gene_id;
		$$ref_splice_dists{$site_string} = $splice_dist;
	}
}
sub process_intron_bed {
	my ($vcf_file_intern,$intron_bed_intern,$sd_ref,$st_ref,$num_fields_intern) = @_;

	$cmd = "$bedtools_loc intersect -wb -a $vcf_file_intern -b $intron_bed_intern | uniq";
	open INT, "-|", $cmd;

	while (my $line = <INT>) { 
		my @sp = split /\s+/, $line;

		my $site_string = "$sp[0]:$sp[1]";
		my $subtract = $sp[$num_fields_intern + 2];
		my $direction = $sp[$num_fields_intern + 5];
		if ($direction eq "-") { 
			$subtract = $sp[$num_fields_intern + 1];
		}
		my $value = log2(abs($sp[1] - $subtract));
		if (!($$sd_ref{$site_string}) || ($value < $$sd_ref{$site_string})) { 
			$$sd_ref{$site_string} = $value;
		}

		$$st_ref{$site_string} = "intron";
	}

}


sub log2 {
    my $n = shift;
	return $n;
	if ($n < 1) { return 0; }
    return log($n)/log(2);
}


