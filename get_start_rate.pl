##### Calculate the stop rate using mpileup file

open IN,"$ARGV[0]" or die $!;   ## The mpileup result
open STOP,">$ARGV[1]" or die $!;   ## The stop rate file output
print STOP "gene_id\tpos\tbase\tstop\ttotal\tstop_rate\n";

my %base;

while(<IN>){
	chomp;
	next if /mpileup/;
	my @sp=split /\s+/;
	$base{$sp[1]}=uc ($sp[2]);
#	if($sp[0]=~ /rRNA/){next if $sp[1]!~ /132/;}
#	next if $sp[3]<5;        ## remove the pos covered less than 10
	my $stop=0;
	while($sp[4]=~ s/\^//){
		$stop+=1;
	}
	my $stop_rate=int($stop/$sp[3]*10000)/100;
	my $tmp=$sp[1];    ## LK
	$base{$tmp}="N" if !exists $base{$tmp};
	print STOP "$sp[0]\t$tmp\t$base{$tmp}\t$stop\t$sp[3]\t$stop_rate\n";
}
