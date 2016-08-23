#!/usr/bin/perl
use strict;

my $ips_file = "Fraxinus_pennsylvanica_120313_peptides.IPSresults.tsv";

my $bed_file = "best_candidates.eclipsed_orfs_removed.bed";

my $go_file = "go-2014-10-26.obo";

my $output_file = "Fraxinus_pennsylvanica_120313";
my $species = "Fraxinus pennsylvanica";

# --------------------------------------
# build a hash to look up go term namespace
# --------------------------------------
my %go2ns;
open IN, $go_file;
my $id;
my $ns;
while(<IN>){
	if(/^id: GO:(\d+)/){
		$id = $1;
	}
	elsif(/namespace: (\S+)/){
		$ns = $1;
		$go2ns{$id} = $ns;
	}
	#print "$id => $ns\n";
}
close IN;
# --------------------------------------
# build a hash to convert aa to nt names
# --------------------------------------
my %aa2nt;
open IN, $bed_file;
<IN>;
while(<IN>){

	/^(\S+)/;
	my $nt = $1;

	/ID=(m\.\d+)/;
	my $aa = $1;

	$aa2nt{$aa} = $nt;
	#print "$aa $nt\n";
}
close IN;

# --------------------------------------
# get IPS results
# --------------------------------------
my %bp_go;
my %cc_go;
my %mf_go;
open IN, $ips_file;
while(my $line = <IN>){
    $line =~ m/^(\S+)/i;
    my $aa_full = $1;
	$aa_full =~ /(m.\d+)$/;
	my $aa_short = $1;
	my $nt = $aa2nt{$aa_short};
	#print "$aa_full => $aa_short => $nt\n";

    while($line =~ m/GO:(\d+)/gi){
        my $go = $1;
		my $ns = $go2ns{$go};

		if($ns =~ /biological/){
			$bp_go{"$nt = $go"} = 1;
		}
		elsif($ns =~ /cellular/){
			$cc_go{"$nt = $go"} = 1;
		}
		elsif($ns =~ /molecular/){
			$mf_go{"$nt = $go"} = 1;
		}
		else{
			print "ignoring namespace $ns (GO:$go)\n";
		}
		#print "\t$go\n";
    }
}
close IN;

# --------------------------------------
# build GAF file
# --------------------------------------
open OUTB, ">$output_file\_biological_process.GAF";
print OUTB "(species=$species)(type=Biological Process)(curator=GO)\n";
foreach my $go (keys %bp_go){
	print OUTB "$go\n";
}
close OUTB;
# ----
open OUTC, ">$output_file\_cellular_component.GAF";
print OUTC "(species=$species)(type=Cellular Component)(curator=GO)\n";
foreach my $go (keys %cc_go){
	print OUTC "$go\n";
}
close OUTC;
# ----
open OUTM, ">$output_file\_molecular_function.GAF";
print OUTM "(species=$species)(type=Molecular Function)(curator=GO)\n";
foreach my $go (keys %mf_go){
	print OUTM "$go\n";
}
close OUTM;
