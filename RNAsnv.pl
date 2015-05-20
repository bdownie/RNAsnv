#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  RNAsnv.pl
#
#===============================================================================

use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;
use Pod::Usage;
use IO::Handle;
my  $io = IO::Handle->new();
Getopt::Long::Configure("pass_through");



my $RNAsnv_dir = abs_path($0);
$RNAsnv_dir =~ s/(.*)\/.*/$1\//;
my $HEADER;

#my $MAX_PHRED33_VALUE=74;
#my $MIN_PHRED_QUAL= 20;

# Default config values
my %CONFIG_HASH;
$CONFIG_HASH{'MINMAPQUAL'} = 255;
$CONFIG_HASH{'JAVAMEM'} = 10;
$CONFIG_HASH{'THREADS'} = 1;
$CONFIG_HASH{'MIN_DP'} = 2;
$CONFIG_HASH{'MAX_HOMONUCLEOTIDE_STRETCH'} = 5;
$CONFIG_HASH{'DIST_FROM_END'} = 6;

my $compress = 0;
my @bamfiles;
my $CONFIG = 0;
my $WGS= 0;
my $EDIT_DB=0;
my $train=0;
my $predict=0;
my $help = 0;
my $man = 0;


GetOptions('c=s'=>\$CONFIG,'train'=>\$train,'wgs=s'=>\$WGS,'predict'=>\$predict,'h|help|?' => \$help, man => \$man) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;
pod2usage("Provide config file location using -c") unless $CONFIG;
pod2usage("Choose only one of -train or -predict") if ($train && $predict);

@bamfiles = @ARGV;



if ($CONFIG) { 
	$CONFIG = abs_path($CONFIG);
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

if ($CONFIG_HASH{'BAMFILES'}) {
	@bamfiles=split /\s+/,$CONFIG_HASH{'BAMFILES'};
}
pod2usage("Provide at least one bam file") unless $bamfiles[0];

#
# File/exec locations
my $SAMTOOLS = $CONFIG_HASH{'SAMTOOLS'};
my $BCFTOOLS = $CONFIG_HASH{'BCFTOOLS'};
my $REFERENCE = $CONFIG_HASH{'REFERENCE'};
my $PICARD = $CONFIG_HASH{'PICARD'};
my $BEDTOOLS = $CONFIG_HASH{'BEDTOOLS'};

# Parameters
my $MIN_DP = $CONFIG_HASH{'MIN_DP'};
my $THREADS = $CONFIG_HASH{'THREADS'};
my $MINMAPQUAL = $CONFIG_HASH{'MINMAPQUAL'};
my $JAVAMEM = $CONFIG_HASH{'JAVAMEM'};

my $MAX_HOMONUCLEOTIDE_STRETCH = $CONFIG_HASH{'MAX_HOMONUCLEOTIDE_STRETCH'};
my $DIST_FROM_END = $CONFIG_HASH{'DIST_FROM_END'};

if ($CONFIG_HASH{'WGS'}) { 
	$WGS= $CONFIG_HASH{'WGS'};
}
if ($CONFIG_HASH{'EDIT_DB'}) { 
	$EDIT_DB= $CONFIG_HASH{'EDIT_DB'};
}
unless ($HEADER) { 
	$HEADER = $CONFIG_HASH{'HEADER'};
}

# Set up appropriate formats
#$PICARD =~ s/(.*)\/.*?$/$1\//;
#$JAVAMEM =~ s/g$//; 

# Validate data entry
unless (-e $REFERENCE) { die "Need a reference file, none given"; }
#unless (-d $SAMTOOLS ) { die "Directory $SAMTOOLS doesn't exist"; }



if ($MINMAPQUAL =~ /\D/) { die "MINMAPQUAL must be an integer"; }
if ($MIN_DP =~ /\D/) { die "MIN_DP must be an integer"; }
if ($THREADS =~ /\D/) { die "THREADS must be an integer"; }
#if ($JAVAMEM =~ /\D/) { die "JAVAMEM must be an integer"; }
#unless ($JAVAMEM =~ /g$/) { $JAVAMEM .= "g"; }

if ($MAX_HOMONUCLEOTIDE_STRETCH =~ /\D/) { die "MAX_HOMONUCLEOTIDE_STRETCH must be an integer"; }
if ($DIST_FROM_END =~ /\D/) { die "MINIMUM_DIST_FROM_5PRIME must be an integer"; }

if ($train && (!-e $WGS || -z $WGS)) {
	pod2usage("Need WGS bam or vcf file to train model");
}
if ($train && (!-e $EDIT_DB || -z $EDIT_DB)) {
	pod2usage("Need RNA edit DB (e.g. RADAR) to train model");
}



unless ($HEADER) { 
	$HEADER = getHeader(@bamfiles);
}

pod2usage("Couldn't dynamically determine header. Please provide in config file.") unless $HEADER;
my $FINAL_VCF = "$HEADER.final.vcf";


my $working_dir = "$HEADER.RNAsnv";
mkdir $working_dir;
open LOG, ">$working_dir/$HEADER.log" or die "Couldn't open $!";
select LOG;
$| = 1;
select STDOUT;
my $return_val;
my $INPUT;
my $OUTPUT = "$HEADER.raw.bam";
my $run_cmd;


my $JAVA_CMD = "ionice -c 3 java -Xmx$JAVAMEM -Djava.io.tmpdir=tmp ";
if ($THREADS) { $JAVA_CMD .= " -XX:ParallelGCThreads=$THREADS "; }
$JAVA_CMD .= "-jar ";

if ($WGS) { $WGS = abs_path($WGS); }

if ($WGS =~ /am$/i) {
	my $alt_bai = $WGS;
	$alt_bai =~ s/bam$/bai/;
        unless ((-e "$WGS.bai") || (-e $alt_bai)) {
		$run_cmd = "samtools index $WGS";
		print LOG "$run_cmd\n";
		$return_val = system ($run_cmd);
		if ($return_val) {
			print LOG "Return value: $return_val, exitting\n";
			print LOG "The bam file is probably not coordinate sorted\n";
			print STDERR "Return value: $return_val for $HEADER, exitting\n";
			print STDERR "The bam file is probably not coordinate sorted\n";
			exit;
		}
	}
}


# Merge and sort BAM files
if ($#bamfiles > 0) {
	my $MASTER_BAM_LOC = "$working_dir/$OUTPUT";
	for (my $i = 0; $i <= $#bamfiles; $i++) { 
		unless (-f $bamfiles[$i]) { print STDERR "Bad file found: $bamfiles[$i]\n"; exit; }
	}
	$run_cmd = "$JAVA_CMD $PICARD/MergeSamFiles.jar SORT_ORDER=coordinate ";
	foreach my $bamfile (@bamfiles) { 
		$run_cmd .= "INPUT=$bamfile ";
	}
	$run_cmd .= "OUTPUT=$MASTER_BAM_LOC USE_THREADING=true QUIET=TRUE VERBOSITY=ERROR > $working_dir/$HEADER.stdout 2>&1";

	unlink $MASTER_BAM_LOC;
	print LOG "$run_cmd\n";
	$return_val = system ($run_cmd);
	if ($return_val) { 
		print LOG "Return value: $return_val, exitting\n";
		print STDERR "Return value: $return_val for $HEADER, exitting\n";
		exit;
	}
	chdir "$working_dir";
}
else { 
	my $source_file = abs_path($bamfiles[0]);
	chdir "$working_dir";
	if (-e $OUTPUT) {
		if (-l $OUTPUT) { 
			unlink $OUTPUT;
		}
		else { 
			die "$OUTPUT file exists, won't overwrite";
		}
	}
	symlink $source_file, $OUTPUT;
}
	
# Set minimum mapping quality (default 255 for STAR)
# Stage 1
$INPUT = $OUTPUT;
$OUTPUT = "$HEADER.sorted.bam";
$run_cmd = "$SAMTOOLS view -hS -q $MINMAPQUAL $INPUT | awk '{OFS = \"\\t\"; if (\$1 ~ /^@/) { print} else if (\$5 > 60) { \$5 = 60;  print \$0 } else { print } }'  |  $JAVA_CMD $PICARD/SortSam.jar I=/dev/stdin O=$OUTPUT SO=coordinate VERBOSITY=ERROR QUIET=TRUE >> $HEADER.stdout  2>&1";
#$run_cmd = "$SAMTOOLS view -hS -q $MINMAPQUAL $INPUT | awk '{OFS = \"\\t\"; if (\$1 ~ /^@/) { print} else if (\$5 > 60) { \$5 = 60;  print \$0 } else { print } }' | samtools view -b - > $OUTPUT";
print LOG "$run_cmd\n";
$return_val = system ($run_cmd);
if ($return_val) { 
	print LOG "Return value: $return_val, exitting\n";
	print STDERR "Return value: $return_val for $HEADER, exitting\n";
	exit;
}

# Remove all duplicate reads using picard tools
$INPUT = $OUTPUT;
$OUTPUT = "$HEADER.bam";
$run_cmd = "$JAVA_CMD  $PICARD/MarkDuplicates.jar I=$INPUT O=$OUTPUT M=/dev/null REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE VERBOSITY=ERROR QUIET=TRUE >> $HEADER.stdout 2>&1; samtools index $OUTPUT >> $HEADER.stdout 2>&1 ";
print LOG "$run_cmd\n";
$return_val = system ($run_cmd);
if ($return_val) { 
	print LOG "Return value: $return_val, exitting\n";
	print STDERR "Return value: $return_val for $HEADER, exitting\n";
	exit;
}

$INPUT=$OUTPUT;
my $VCF = "$HEADER.vcf";

# Parallelize code done with assistance of http://www.research.janahang.com/efficient-way-to-generate-vcf-files-using-samtools/

unless (-e "$REFERENCE.fai") { 
	$run_cmd = "$SAMTOOLS faidx $REFERENCE";
	print LOG "$run_cmd\n";
	$return_val = system ($run_cmd);
	if ($return_val) { 
		print LOG "Return value: $return_val, exitting\n";
		print STDERR "Return value: $return_val for $HEADER, exitting\n";
		exit;
	}
}

$run_cmd = "$SAMTOOLS view -H $INPUT | grep \"\@SQ\" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P $THREADS  sh -c \"ionice -c 3 $SAMTOOLS mpileup -g -t DP,SP -O -r {} -f $REFERENCE $INPUT | $BCFTOOLS  filter -i 'DP >= $MIN_DP' | $BCFTOOLS view -O b > tmp.{}.bcf\"";

print LOG "$run_cmd\n";
$return_val = system ($run_cmd);
if ($return_val) { 
	print LOG "Return value: $return_val, exitting\n";
	print STDERR "Return value: $return_val for $HEADER, exitting\n";
	exit;
}

$run_cmd = "$BCFTOOLS concat -O b tmp.*.bcf | tee $HEADER.bcf | $BCFTOOLS view - |  $BCFTOOLS call -m - | $BCFTOOLS filter -i '%TYPE!=\"indel\"' - |  $BCFTOOLS filter -e '%TYPE=\"ref\"' |  $RNAsnv_dir/remove_3prime_5prime_variants -b $INPUT -d $DIST_FROM_END -m $MAX_HOMONUCLEOTIDE_STRETCH | $BEDTOOLS sort -header -i - > $VCF 2>> $HEADER.stdout";
#$run_cmd = "cat $HEADER.bcf | $BCFTOOLS view - |  $BCFTOOLS call -m - | $BCFTOOLS filter -i '%TYPE!=\"indel\"' - |  $BCFTOOLS filter -e '%TYPE=\"ref\"' |  $RNAsnv_dir/remove_3prime_5prime_variants -b $INPUT -d $DIST_FROM_END -m $MAX_HOMONUCLEOTIDE_STRETCH | $BEDTOOLS sort -header -i - > $VCF 2>> $HEADER.stdout";
print LOG "$run_cmd\n";
$return_val = system ($run_cmd);
# This may be necessary to make sure that all the files get sync'd
#$io->sync;
if ($return_val) { 
	print LOG "Return value: $return_val, exitting\n";
	print STDERR "Return value: $return_val for $HEADER, exitting\n";
	exit;
}
else { 
#exit;
	unlink glob "tmp.*.bcf";
}

if ($train || $predict) {
	my $get_features = "$RNAsnv_dir/GetFeatures.pl -c $CONFIG ";
	my $ANNOTATION = $CONFIG_HASH{'ANNOTATION'};
	my $INTRON = $CONFIG_HASH{'INTRON'};
	my $REPEAT_INTERVALS = $CONFIG_HASH{'REPEAT_INTERVALS'};


	$get_features .= "-v $VCF -r $REFERENCE ";

	if ($INTRON) { 
		$get_features .= "-i $INTRON ";
	}
	if ($REPEAT_INTERVALS) { 
		$get_features .= "-a $REPEAT_INTERVALS ";
	}
	if ($ANNOTATION) { 
		$get_features .= "-g $ANNOTATION ";
	}

	my $r_cmd = "R --file='$RNAsnv_dir/RNAsnv.R' CONFIG=$CONFIG $HEADER.features ";

	#print "($WGS)($EDIT_DB)\n"; exit;

	if ($WGS && $EDIT_DB) { 
		$get_features .= "-e $EDIT_DB -d $HEADER.WGS.vcf";
	
	
		unlink glob "tmp.*.bed";
		#open VCF, "-|", "$SAMTOOLS view -H $INPUT | grep \"\@SQ\" | sed 's/^.*SN://g' | cut -f 1";
		open VCF_FILE, $VCF;
		my $last = "0";
		while (my $line = <VCF_FILE>) { 
			next if ($line =~ /^#/);
			
			my @s= split /\s+/, $line;
			if ($s[0] ne $last) { 
				if ($last ne "0") { close OUT; }
				open OUT, ">>tmp.$s[0].bed";
			}
			my $chr = $s[0];
			my $start = $s[1] - 1;
			my $end = $s[1];
			my $one = $s[3];
			my $two = $s[4];
			print OUT "$chr\t$start\t$end\t$one/$two\n";
		}
		close OUT;
		close VCF_FILE;
		#unlink glob "tmp.*.bed";
					

		# Match cram or bam
		if ($WGS =~ /am$/i) { 
			$run_cmd = "ls tmp.*.bed | sed 's/^tmp.//' | sed 's/.bed//' | xargs -I {} -n 1 -P $THREADS  sh -c \"ionice -c 3 $SAMTOOLS mpileup -uv -t DP,SP -O -f $REFERENCE -r {} -l tmp.{}.bed $WGS | $BCFTOOLS call -m - | $BCFTOOLS filter -i 'DP > 10 && QUAL > 100' | $BCFTOOLS view -O b - > tmp.{}.WGS.bcf\"";
			#$run_cmd = "$SAMTOOLS mpileup -uv -t DP,SP -O -f $REFERENCE -l $HEADER.vcf.bed $WGS| $BCFTOOLS call -m - | $BCFTOOLS filter -i 'DP > 10 && QUAL > 100'  > $HEADER.WGS.vcf"; 
			print LOG "$run_cmd\n";
			#print "$run_cmd\n"; exit;
			$return_val = system ($run_cmd);
			if ($return_val) { 
				print LOG "Return value: $return_val, exitting\n";
				print STDERR "Return value: $return_val for $HEADER, exitting\n";
				exit;
			}

			#$run_cmd = "$SAMTOOLS mpileup -uv -t DP,SP -O -f $REFERENCE -l $HEADER.vcf.bed $WGS| $BCFTOOLS call -m - | $BCFTOOLS filter -i 'DP > 10 && QUAL > 100'  > $HEADER.WGS.vcf"; 
			#print LOG "$run_cmd\n"tmp.;
			$return_val = system ($run_cmd);
			$run_cmd = "$BCFTOOLS concat tmp.*.WGS.bcf |  bedtools sort -header -i - > $HEADER.WGS.vcf";
			#$run_cmd = "$SAMTOOLS view -H $INPUT | grep \"\@SQ\" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P $THREADS  sh -c \"ionice -c 3 $SAMTOOLS mpileup -uv -t DP,SP -O -f $REFERENCE -l $HEADER.{}.bed $WGS | $BCFTOOLS call -m - | $BCFTOOLS filter -i 'DP > 10 && QUAL > 100' > $HEADER.{}.WGS.vcf\""; 
			print LOG "$run_cmd\n";
			$return_val = system ($run_cmd);
			if ($return_val) { 
				print LOG "Return value: $return_val, exitting\n";
				print STDERR "Return value: $return_val for $HEADER, exitting\n";
				exit;
			}
			unlink glob "tmp.*.bed";
			unlink glob "tmp.*.WGS.bcf";
		}
		elsif ($WGS =~ /.vcf$/i) { 
			$run_cmd = "$BEDTOOLS intersect -a $WGS -b $HEADER.vcf | uniq > $HEADER.WGS.vcf"; 
			print LOG "$run_cmd\n";
			$return_val = system ($run_cmd);
			if ($return_val) { 
				print LOG "Return value: $return_val, exitting\n";
				print STDERR "Return value: $return_val for $HEADER, exitting\n";
				exit;
			}
		}
		else { die "Unknown WGS file"; }
	}
	if ($predict) {
		my $GBM_MODEL = $CONFIG_HASH{'GBM_MODEL'};
		$r_cmd .= "$HEADER.vcf $GBM_MODEL ";
	}
#	print "($WGS)($EDIT_DB)\n";
	$get_features .= " > $HEADER.features";
	print LOG "$get_features\n";
	$return_val = system ($get_features);
	if ($return_val) { 
		print LOG "Return value: $return_val, exitting\n";
		print STDERR "Return value: $return_val for $HEADER, exitting\n";
		exit;
	}

	print LOG "$r_cmd\n";
	$return_val = system ($r_cmd);
	if ($return_val) { 
		print LOG "Return value: $return_val, exitting\n";
		print STDERR "Return value: $return_val for $HEADER, exitting\n";
		exit;
	}
}

sub getHeader {
	my @files = @_;

	my $HEADER = $files[0];
	$HEADER =~ s/\.bam//; 
	$HEADER =~ s/.*\///; 
	my @file_names = @files;
	my @master = split "", $HEADER;
	for (my $i = 1; $i <= $#file_names; $i++) { 
		my $tmp = $file_names[$i]; 
		$tmp =~ s/\.bam//; 
		$tmp =~ s/.*\///; 
		my @seq = split "", $tmp;
		for (my $j = 0; $j <= $#seq; $j++) { 
			if ((!defined($master[$j]) || !defined($seq[$j])) || ($master[$j] ne $seq[$j])) { 
				$HEADER = substr($HEADER,0,$j);
			}
		}
	}
	$HEADER =~ s/[\.\_\-\:]$//;

	return $HEADER;
}


__END__

=head1 NAME

RNAsnv.pl - calls RNA-seq based variants from BAM files. 

            GBM model can be trained (when DNA BAM/VCF file is provided) 
            or used to predict variant site origin from a pre-existing model.

            A human-based model can be downloaded from X.

=head1 SYNOPSIS

RNAsnv.pl -c <configuration file> [options] [BAM1 BAM2 ...]

 Options:
	-help 		this help message
	-man		displays RNAsnv man page
	-c		configuration file
	-train 		train model using parameters in config file
	-predict 	predict variant site origin using parameters in config file
	-wgs		WGS BAM or VCF file

=cut
