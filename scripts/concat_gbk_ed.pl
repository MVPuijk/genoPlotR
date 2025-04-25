#! /usr/bin/perl -w
use Getopt::Long;

open (FILE, $ARGV[0]) or die "Can't open $ARGV[0]: $!";

open (TMP, ">temp_xxx.gbk");


my $ignoresymbols;
GetOptions("i" => \$ignoresymbols
          );



$total_length=0;
$i=0;
while (<FILE>) {
  chomp;
  next if (/^CONTIG/);

  if(/^LOCUS\s+\S+\s+(\d+)\s+bp/) {
      $length=$1;
      $total_length+=$length;
#      print $length,"\t",$total_length,"\tLENGTH\n";
  }elsif (/^\s+source\s+(\d+)\.\.(\d+)/) {
      $start=$1+$total_length-$length;
      $end=$2+$total_length-$length;
      $i=0;
      print TMP "     contig          ",$start,"..",$end,"\n";
  }elsif (/^\s+CDS\s+(\<*)(\d+)\.\.(\>*)(\d+)/) {
      $start=$2+$total_length-$length;
      $end=$4+$total_length-$length;
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     CDS             ",$start,"..",$end,"\n";
      print TMP "                     /pseudo\n" if $_ =~ /\<|\>/;
  }elsif (/^\s+CDS\s+complement\((\<*)(\d+)\.\.(\>*)(\d+)\)/) {
      $start=$2+$total_length-$length;
      $end=$4+$total_length-$length;
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     CDS             complement(",$start,"..",$end,")\n";
      print TMP "                     /pseudo\n" if $_ =~ /\<|\>/;
  }elsif (/^\s+gene\s+(\<*)(\d+)\.\.(\>*)(\d+)/) {
      $start=$2+$total_length-$length;
      $end=$4+$total_length-$length;
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     gene            ",$start,"..",$end,"\n";
      print TMP "                     /pseudo\n" if $_ =~ /\<|\>/;
  }elsif (/^\s+gene\s+complement\((\<*)(\d+)\.\.(\>*)(\d+)\)/) {
      $start=$2+$total_length-$length;
      $end=$4+$total_length-$length;
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     gene            complement(",$start,"..",$end,")\n";
      print TMP "                     /pseudo\n" if $_ =~ /\<|\>/;
  }elsif (/^\s+tRNA\s+(\<*)(\d+)\.\.(\>*)(\d+)/) {
      $start=$2+$total_length-$length;
      $end=$4+$total_length-$length;
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     tRNA            ",$start,"..",$end,"\n";
  }elsif (/^\s+tRNA\s+complement\((\<*)(\d+)\.\.(\>*)(\d+)\)/) {
      $start=$2+$total_length-$length;
      $end=$4+$total_length-$length;
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     tRNA            complement(",$start,"..",$end,")\n";
  }elsif (/^\s+rRNA\s+(\<*)(\d+)\.\.(\>*)(\d+)/) {
      $start=$2+$total_length-$length;
      $end=$4+$total_length-$length;
#      print STDERR $length, "\n";
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     rRNA            ",$start,"..",$end,"\n";
  }elsif (/^\s+rRNA\s+complement\((\<*)(\d+)\.\.(\>*)(\d+)\)/) {
      $start=$2+$total_length-$length;
      $end=$4+$total_length-$length;
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     rRNA            complement(",$start,"..",$end,")\n";
  }elsif (/^\s+misc_RNA\s+(\<*)(\d+)\.\.(\>*)(\d+)/) {
      $start=$2+$total_length-$length;
      $end=$4+$total_length-$length;
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     misc_RNA        ",$start,"..",$end,"\n";
  }elsif (/^\s+misc_RNA\s+complement\((\<*)(\d+)\.\.(\>*)(\d+)\)/) {
      $start=$2+$total_length-$length;
      $end=$4+$total_length-$length;
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     misc_RNA        complement(",$start,"..",$end,")\n";
  }elsif (/^\s+repeat_region\s+(\<*)(\d+)\.\.(\>*)(\d+)/) {
      $start=$2+$total_length-$length;
#     print STDERR $length, "\n";
      $end=$4+$total_length-$length;
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     repeat_region   ",$start,"..",$end,"\n";
  }elsif (/^\s+repeat_region\s+complement\((\<*)(\d+)\.\.(\>*)(\d+)\)/) {
      $start=$2+$total_length-$length;
      $end=$4+$total_length-$length;
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     repeat_region   complement(",$start,"..",$end,")\n";
  }elsif (/^\s+gap\s+(\<*)(\d+)\.\.(\>*)(\d+)/) {
      $start=$2+$total_length-$length;
      $end=$4+$total_length-$length;
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     gap             ",$start,"..",$end,"\n";
  }elsif (/^\s+gap\s+complement\((\<*)(\d+)\.\.(\>*)(\d+)\)/) {
      $start=$2+$total_length-$length;
      $end=$4+$total_length-$length;
      $start=$1.$start unless ($ignoresymbols);
      $end=$3.$end unless ($ignoresymbols);
      print TMP "     gap             complement(",$start,"..",$end,")\n";
  } else {
      if (/^\/\//) {
	  $i=1;
      }elsif (/BASE COUNT/){
	  $i=1;
      }elsif (/ORIGIN/){
	  $i=2;
      }elsif ($i == 2){
	  s/\d+//g;
	  s/\s+//g;
	  $seq.=$_;
      }else{
	  print TMP $_,"\n" if $i==0;
      }
  }
}
$k=1;
print TMP "ORIGIN\n";
#print TMP $seq,"\n";
@sequences = split "", $seq;
$final_length = scalar(@sequences)-1;
for ($j=0;$j<=$final_length;$j+=60){
    if ($k < 10 ) {
	$space="       ";
    }elsif ($k < 100 ) {
	$space="      ";
    }elsif ($k < 1000 ) {
	$space="     ";
    }elsif ($k < 10000 ) {
	$space="    ";
    }elsif ($k < 100000 ) {
	$space="   ";
    }elsif ($k < 1000000 ) {
	$space="  ";
    }
    print TMP $space,$k," ",@sequences[$j..$j+9]," ",@sequences[$j+10..$j+19]," ",@sequences[$j+20..$j+29]," ",@sequences[$j+30..$j+39]," ",@sequences[$j+40..$j+49]," ",@sequences[$j+50..$j+59],"\n" if $j+58 <= $final_length-1; 

#     print TMP "  ",$k," ",@sequences[$j..$final_length-1],"\n" if $j+58 > $final_length-1;

    if ($j+58 > $final_length-1) {
	print TMP "  ",$k," ";

	if ($final_length-$k > 10) {
	    print TMP @sequences[$j..$j+9], " ";
	    if ($final_length-$k > 20) {
		print TMP @sequences[$j+10..$j+19], " ";
		if ($final_length-$k > 30) {
		    print TMP @sequences[$j+20..$j+29], " ";
		    if ($final_length-$k > 40) {
			print TMP @sequences[$j+30..$j+39], " ";
			if ($final_length-$k > 50) {
			    print TMP @sequences[$j+40..$j+49], " ";
			    print TMP @sequences[$j+50..$final_length], "\n";
			} else {
			    print TMP @sequences[$j+40..$final_length], "\n";
			}
		    } else {
			print TMP @sequences[$j+30..$final_length], "\n";
		    }
		} else {
		    print TMP @sequences[$j+20..$final_length], "\n";
		}
	    } else {
		print TMP @sequences[$j+10..$final_length], "\n";
	    }
	} else {
	    print TMP @sequences[$j..$final_length], "\n";
#	    print STDERR @sequences[$j..$final_length], "\n";
	}

    }

    $k+=60;
}
print TMP "//\n";
close FILE;
close TMP;

print "LOCUS       XX_XXXXXXXXXXXX";
print sprintf("% 13d",$total_length);
print " bp     DNA     linear   CON XX-XXX-20XX\n";

open TMP, "temp_xxx.gbk";
@in = <TMP>;
print  @in;

unlink "temp_xxx.gbk";
