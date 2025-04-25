#!/usr/bin/perl

# libraries
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use feature ("say");

my $blast1; my $blast2;
my $positions;
my $gbk1; my $gbk2;
my $col;
my $help;
my $input;

GetOptions(
    'b1|blast1=s' => \$blast1,
    'b2|blast2=s' => \$blast2,
    'p|positions=s' => \$positions,
    'g1|genbank1=s' => \$gbk1,
    'g2|genbank2=s' => \$gbk2,
    'c|col|colors=s' => \$col,
    'h|help' => \$help,
);


pod2usage(-verbose => 2) if $help;
pod2usage(-msg => "** ERROR: No input ** \n") unless ($blast1 && $blast2 && ($positions || ($gbk1 && $gbk2)));


unless ($positions || ($gbk1 && $gbk2)) {
    die "Provide gene positions either as 2 genbank files (the blast ids need to correspond to locus tags), or as a single tab file with structure '<gene ID>\t<start>\t<end>'\n";
}

my %pos = ();
if ($gbk1 && $gbk2) {
    %pos = &get_pos;
}
if ($positions) {
    open POSITIONS, $positions;
    while (<POSITIONS>) {
	my @split = split;
	my $gene = $split[0];
	$pos{$gene}{'start'} = $split[1];
	$pos{$gene}{'end'} = $split[2];
	$pos{$gene}{'or'} = '-' if ($pos{$gene}{'start'} > $pos{$gene}{'end'});
	$pos{$gene}{'or'} = '+' if ($pos{$gene}{'start'} < $pos{$gene}{'end'});
	

    }
    close POSITIONS;
}
my %colors = ();
if ($col) {
    %colors = &get_cols();
} else {
    %colors = ();
}

my $i = 0;
my %hits = ();
for my $file ($blast1,$blast2) {
    $i++;
    open IN, $file;
    while (<IN>) {
        next if (/^#/);
        my @split = split /\t/;
        if (! exists $hits{$i}{$split[0]}) {
            $hits{$i}{$split[0]}{'id'} = $split[1];
            $hits{$i}{$split[0]}{'evalue'} = $split[10];
        }
        else {
            if ($hits{$i}{$split[0]}{'evalue'} > $split[10]) {
#               print "$hits{$i}{$split[0]}{'evalue'} >  $split[10]\n";
                $hits{$i}{$split[0]}{'id'} = $split[1];
                $hits{$i}{$split[0]}{'evalue'} = $split[10];
            }
        }
    }
    close IN;
}


print "start1\tend1\tstart2\tend2\tcol\tgene1\tgene2\n";
for my $key (sort keys %{$hits{'1'}}) {
    
    my $best = $hits{'1'}{$key}{'id'};

    if ($hits{'2'}{$best}{'id'}) {
	if ($key eq $hits{'2'}{$best}{'id'}) {
	    my $col = "";
	    if ($colors{$key} || $colors{$best}) {
		$col = $colors{$best} if ($colors{$best});
		$col = $colors{$key} if ($colors{$key});
	    } else {
		$col = "grey"
	    }
	    print $pos{$key}{'start'}, "\t";
	    print $pos{$key}{'end'}, "\t";
	    print $pos{$best}{'start'}, "\t";
	    print $pos{$best}{'end'}, "\t";
	    print $col, "\t";
	    print $key, "\t";
	    print $best, "\n";
	}
    } else {
#       print STDERR "No best hit for $key ($hits{'1'}{$key}{'id'}, $hits{'1'}{$key}{'evalue'}) ($best, $hits{'2'}{$best}{'id'}, $hits{'2'}{$best}{'evalue'}) \n";
    }
}








#############################
sub get_pos {
    my %pos = ();
    foreach my $gb ($gbk1, $gbk2) {
	my $name = $gb;
	$name =~ s/\.gb.*//;
#	my @keys = keys %{$ids{$name}};
	
	open GEN, $gb or die "Could not open $gb\n";
	my $cds = my $start = my $end = 0;
	my $check = 0;
	my $or = my $tag = "";
	while(<GEN>) {
	    unless (($cds)||(/[<>]?\d+[<>]?\.\.[<>]?\d+/)) {
		next;
	    }
	    if (/\d+[<>]?\.\.[<>]?\d+/) {
		if (/^\s+\s\sCDS\s\s\s+[A-Za-z\(]*[<>]?(\d+)[<>]?\.\.[<>]?(\d+)/) {
		    $cds = 1;
		    $start = $1;
		    $end = $2;
		    $or = "+";
		    $or = "-" if (/complement/); 
		} else {
		    $cds = 0;
		}
		next;
	    }
	    if ($cds) {
		if (/^\s+\/locus_tag=\"(.*)\"\s*$/) {
		    $tag = $1;
		    if ($or eq "+") {
			$pos{$tag}{'start'} = $start;
			$pos{$tag}{'end'} = $end;
			$pos{$tag}{'or'} = "+";
		    } elsif ($or eq "-") {
			$pos{$tag}{'start'} = $end;
			$pos{$tag}{'end'} = $start;
			$pos{$tag}{'or'} = "-";
		    }
		}
	    }
	}
    }
    return %pos;
}

sub get_cols {
    say STDERR $col;
    my %colors = ();
    open COL, $col;
    while (my $c = <COL>) {
	chomp $c;
	my @columns = split("\t", $c); 
	$colors{$columns[0]} = $columns[2];
    }
    close COL;
    return %colors;
}


###################


=head1 NAME

    blastp2genoplotr_comp.pl - Short desc

=head1 USAGE

    blastp2genoplotr_comp.pl -blast1 -blast2    (-p | (-g1 & -g2))   [-h (show help)]

=head1 DESCRIPTION

    Longer desc

=head1 OPTIONS

=head2 REQUIRED

=over 3

=item B<-b1, --blast1>

=over 1

=item Tabular blast file for Genome 1 vs Genome 2

=back

=item B<-b2, --blast2>

=over 1

=item Tabular blast file for Genome 2 vs Genome 1

=back

=back

=head2 OPTIONAL

=over 3

=item B<-h, --help>

=over 1

=item Show help

=back 

=item B<-p, --positions>

=over 1

=item Tabular file with structure <gene ID>\t<start>\t<end>. Mandatory if no genbank files are provided

=back

=item B<-g1, --genbank1>

=over 1

=item Genbank file for Genome 1. Contig multi-genbank files are not properly processesd -- concatenation prior to using this script is necessary. Mandatory if a tabular file with positions is not provided.

=back

=item B<-g2, --genbank2>

=over 1

=item Genbank file for Genome 2. Contig multi-genbank files are not properly processesd -- concatenation prior to using this script is necessary. Mandatory if a tabular file with positions is not provided.

=back

=back

=head1 DEPENDENCIES

=head1 AUTHOR

Daniel Tamarit (L<daniel.tamarit@icm.uu.se>)

=head1 DATE

2018-05-28

=head1 COPYRIGHT AND LICENSE

Copyright 2018 Daniel Tamarit

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
