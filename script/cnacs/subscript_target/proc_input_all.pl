#! /usr/local/bin/perl -w
use strict;

# cumulative lengths of chromosomes
my %cum_length;
my $hg_version = (defined $ENV{'HUMAN_GENOME_VERSION'}) ? $ENV{'HUMAN_GENOME_VERSION'} : 'hg38';
if ($hg_version eq 'hg38') {
	$cum_length{'chr1'}  = 0;
	$cum_length{'chr2'}  = 248956422;
	$cum_length{'chr3'}  = 491149951;
	$cum_length{'chr4'}  = 689445510;
	$cum_length{'chr5'}  = 879660065;
	$cum_length{'chr6'}  = 1061198324;
	$cum_length{'chr7'}  = 1232004303;
	$cum_length{'chr8'}  = 1391350276;
	$cum_length{'chr9'}  = 1536488912;
	$cum_length{'chr10'} = 1674883629;
	$cum_length{'chr11'} = 1808681051;
	$cum_length{'chr12'} = 1943767673;
	$cum_length{'chr13'} = 2077042982;
	$cum_length{'chr14'} = 2191407310;
	$cum_length{'chr15'} = 2298451028;
	$cum_length{'chr16'} = 2400442217;
	$cum_length{'chr17'} = 2490780562;
	$cum_length{'chr18'} = 2574038003;
	$cum_length{'chr19'} = 2654411288;
	$cum_length{'chr20'} = 2713028904;
	$cum_length{'chr21'} = 2777473071;
	$cum_length{'chr22'} = 2824183054;
	$cum_length{'chrX'}  = 2875001522;
	$cum_length{'chrY'}  = 3031042417;
} else {
	$cum_length{'chr1'}  = 0;
	$cum_length{'chr2'}  = 249250621;
	$cum_length{'chr3'}  = 492449994;
	$cum_length{'chr4'}  = 690472424;
	$cum_length{'chr5'}  = 881626700;
	$cum_length{'chr6'}  = 1062541960;
	$cum_length{'chr7'}  = 1233657027;
	$cum_length{'chr8'}  = 1392795690;
	$cum_length{'chr9'}  = 1539159712;
	$cum_length{'chr10'} = 1680373143;
	$cum_length{'chr11'} = 1815907890;
	$cum_length{'chr12'} = 1950914406;
	$cum_length{'chr13'} = 2084766301;
	$cum_length{'chr14'} = 2199936179;
	$cum_length{'chr15'} = 2307285719;
	$cum_length{'chr16'} = 2409817111;
	$cum_length{'chr17'} = 2500171864;
	$cum_length{'chr18'} = 2581367074;
	$cum_length{'chr19'} = 2659444322;
	$cum_length{'chr20'} = 2718573305;
	$cum_length{'chr21'} = 2781598825;
	$cum_length{'chr22'} = 2829728720;
	$cum_length{'chrX'}  = 2881033286;
	$cum_length{'chrY'}  = 3036303846;
}

### Load all BAFs including homozygous SNPs ###

open BAF, '<', $ARGV[2] || die "cannot open $!";

my %pos2baf;
while (<BAF>) {
	s/[\r\n]//g;
	my @curRow = split(/\t/, $_);
	my $chr = $curRow[0];
	my $pos = $curRow[1];
	my $key = $chr . "\t" . $pos;
	$pos2baf{$key} = $curRow[2];
}
close(BAF);

### End of loading all BAFs ###



### Load all the CNAs ###

open SEG, '<', $ARGV[1] || die "cannot open $!";
open OUT_SEG, '>', $ARGV[4] || die "cannot open $!";

my %start2end;

while (<SEG>) {
	s/[\r\n]//g;
	my @curRow = split(/\t/, $_);
	my $chr_num = $curRow[1];
	my $chr = 'chr' . $chr_num;
	if ( $chr_num == 23 ) {
		$chr = 'chrX';
	} elsif ( $chr_num == 24 ) {
		$chr = 'chrY';
	}
	
	my $start = $curRow[2];
	my $end = $curRow[3];
	
	my $key = $chr . "\t" . $start;
	$start2end{$key} = $end;
	
	my $ploidy = 'NA';
	if ($curRow[5] =~ /\d/ ) {
		$ploidy = $curRow[5] + 4;
	}
	my $as = $curRow[6];
	
	if ( defined $cum_length{$chr} ) {
		$start += $cum_length{$chr};
		$end += $cum_length{$chr};
	} else {
		next;
	}
	
	$start = $start / 1000000;
	$end = $end / 1000000;
	
	print OUT_SEG $start . ',' . $end . ',' . $ploidy . ',' . $as . "\n";
}
close(SEG);
close(OUT_SEG);

### End of loading all the CNAs ###



### Check whether heterozygous SNPs are present in the CNA regions ###

open SIG, '<', $ARGV[0] || die "cannot open $!";
my %seg2flag;
my %seg2as;
my $seg_id = 0;
my $processing = 0;
my $cur_end;

while (<SIG>) {
	s/[\r\n]//g;
	my @curRow = split(/\t/, $_);
	my $chr = $curRow[0];
	my $pos = $curRow[1];
	my $key = $chr . "\t" . $pos;
	
	if ( defined $start2end{$key} ) {
		$seg_id++;
		$seg2flag{$seg_id} = 1;
		$processing = 1;
		$cur_end = $start2end{$key};
	}
	
	my $as = $curRow[3];
	if ( $processing == 1 ) {
		$seg2flag{$seg_id} = 0 if ( ( $as ne 'NA' ) && ( $as > 0.12 ) );
		
		if ( $pos == $cur_end ) {
			$cur_end = "";
			$processing = 0;
		}
		
		if ( $as eq 'NA' ) {
			$as = $pos2baf{$key} if ( defined $pos2baf{$key} );
			next if ( ! defined $pos2baf{$key} );
		}
		
		if ( defined $seg2as{$seg_id} ) {
			$seg2as{$seg_id} .= ',' . $as;
		} else {
			$seg2as{$seg_id} = $as;
		}
	}
}
close(SIG);

### End of checking heterozygous SNPs ###



### Determine threshold of BAFs for drawing plots ###

my %seg2thresh;
foreach my $id ( keys %seg2as ) {
	my @nums = split(/,/, $seg2as{$id});
	my $median = 1;
	$median = &percentile(50, @nums) if ( @nums > 0 );
	$seg2thresh{$id} = $median;
}

### End of determining thresholds ###



### Output ###

open SIG, '<', $ARGV[0] || die "cannot open $!";
open OUT_SIG, '>', $ARGV[3] || die "cannot open $!";
$seg_id = 0;
while (<SIG>) {
	s/[\r\n]//g;
	my @curRow = split(/\t/, $_);
	my $chr = $curRow[0];
	my $pos = $curRow[1];
	my $key = $chr . "\t" . $pos;
	if ( defined $start2end{$key} ) {
		$seg_id++;
		if ( $seg2flag{$seg_id} == 1 ) {
			$processing = 1;
			$cur_end = $start2end{$key};
		}
	}
	
	my $ploidy = $curRow[2];
	$ploidy += 4 if ( $ploidy ne 'NA' );
	
	my $as = $curRow[3];
	my $thresh = 0.12;
	if ( $processing == 1 ) {
		$thresh = $seg2thresh{$seg_id} if ( defined $seg2thresh{$seg_id} );
		$as = $pos2baf{$key} if ( defined $pos2baf{$key} );
		if ( $pos == $cur_end ) {
			$cur_end = "";
			$processing = 0;
		}
	}
	unless ( ( $as =~ /\d/ ) && ( $as > $thresh ) ) {
		$as = 'NA';
	}
	
	if ( defined $cum_length{$chr} ) {
		$pos += $cum_length{$chr};
	} else {
		next;
	}
	
	$pos = $pos / 1000000;
	
	print OUT_SIG $pos . ',' . $ploidy . ',' . $as . "\n";
}
close(SIG);
close(OUT_SIG);

### End of output ###



sub percentile {
	my $percent = $_[0];
	my @sorted = sort { $a <=> $b } @_[ 1 .. $#_ ];
	my $idx = int( @sorted * $percent / 100 );
	my $down_dif = @sorted * $percent / 100 - $idx;
	my $up_dif = 1 - $down_dif;
	my $value = $sorted[ $idx - 1 ] * $up_dif + $sorted[ $idx ] * $down_dif;
	return $value;
}
