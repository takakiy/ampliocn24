#!/usr/bin/perl -w

if ((@ARGV == 0) || ($ARGV[0]=~/-h/)) {
	print "OPTION -e clip -max number (-min 0/number) 
			-s single-reads fq (-pe paired-end files)  
			-max 100 -min 100 ## only 100 bp
			-max 100 -min 0  non-filter out of short read ## CUTTING CONSTANT LENGTH \n";
	print "       -e sel -lim LEN (-lgc min max ) ## Filter by READS LEN\n";
	print "       -e stat   ## STATICS OF READS LEN  \n";
	print "       -e sep    ## SEPARATE one fq (Paired-end) to two fq\n";
	print "       -e inter    ## INTERLEAVE two fq (Paired-end) to single fq\n";
	print "       -e pick -num ##   PICK UP \n";
	print "       -e get -l list ##   EXTRACT FOR LIST (PAIR) \n";
    print "       -e rm -l list ##   REMOVE FOR LIST (PAIR) \n";
    print "       -e trim -fixlen 3 3 (head tail) ##   REMOVE FOR LIST (PAIR) \n";

	exit;

}

use IO::Uncompress::Gunzip qw($GunzipError);


foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};
$|= 1;

if (exists $option->{-s}) {
    if ( $option->{-s}->[0] =~ /.gz$/ ) {
        $fh1 = IO::Uncompress::Gunzip->new( $option->{-s}->[0] ) ||
                die "Could not make uncompress object: $GunzipError";
    } else {
        $fh1="IN1";
        open ($fh1,"$option->{-s}->[0]");
    }
    $fqtype= "single";
    
} elsif ($option->{-pe}) {
     if ( $option->{-pe}->[0] =~ /.gz$/ ) {
        $fh1 = IO::Uncompress::Gunzip->new( $option->{-pe}->[0] ) ||
                die "Could not make uncompress object1: $GunzipError";
    } else {
        $fh1="IN1";
        open ($fh1,"$option->{-pe}->[0]");
    }
    if ( $option->{-pe}->[1] =~ /.gz$/ ) {
        $fh2 = IO::Uncompress::Gunzip->new( $option->{-pe}->[1] ) ||
                die "Could not make uncompress object2: $GunzipError";
    } else {
        $fh2="IN2";
        open ($fh2,"$option->{-pe}->[1]");
    }
    $fqtype= "pair";

}

# &output($fh2);

sub output {
    my ($fhh)=@_;
    $n=0;
    while(<$fhh>) {
        print;
        $n++;
        last if ($n > 10);
    }
}
#exit;

if ($option->{-e}->[0] eq "clip") {
	$lim_max= exists $option->{-max} ? $option->{-max}->[0] : 0;
	$lim_min= exists $option->{-min} ? $option->{-min}->[0] : 0;

	if ($fqtype eq "single") {
		$fqtype= "single";
		($info)= &clipping_by_len([$fh1],$lim_max,$lim_min,$fqtype);
		print "MIN-MAX: $lim_min $lim_max\n";
		print "SEQNUM: $info->{count} GET: $info->{get}\n";
		print "OUTPUT: out_cliplen_R1_OP.fq\n";
	} elsif ($fqtype eq "pair") {
		$fqtype= "pair";
		($info)= &clipping_by_len([$fh1,$fh2],$lim_max,$lim_min,$fqtype);
		print "MIN-MAX: $lim_min $lim_max\n";
		print "SEQNUM: $info->{count} GET: $info->{get}\n";
		print "OUTPUT: out_cliplen_R1_MP.fq out_cliplen_R2_MP.fq\n";
	}
    
} elsif ($option->{-e}->[0] eq "sel") {
	$lim_len= exists $option->{-lim} ? $option->{-lim}->[0] : 0;
	$lim_gc= exists $option->{-lgc} ? $option->{-lgc} : [0,1];
    
    if ($fqtype eq "single") {
        $fqtype= "single";
        ($info)= &select_by_length([$fh1],$lim_len,$fqtype,$lim_gc);
		print "LIM_LEN: $lim_len\n";
		print "SEQNUM: $info->{count} GET: $info->{get}\n";
		print "OUTPUT: out_sellen_R1_OP.fq\n";

	} elsif ($fqtype eq "pair") {
		($info)= &select_by_length([$fh1,$fh2],$lim_len,$fqtype,$lim_gc);
		print "LIM_LEN: $lim_len\n";
	#	print "SEQNUM: $info->{count} PAIR: $info->{pair} SINGLE: $info->{orph_fwd} $info->{orph_rev}\n";
		print "OUTPUT: out_sellen_R1_PE.fq out_sellen_R2_PE.fq\n";
		print "OUTPUT: out_sellen_R1_OP.fq out_sellen_R2_OP.fq\n";
        
        print "IN: $info->{count}\nAMB: $info->{amb}\nPE: $info->{pair}\n";
	}
    
} elsif ($option->{-e}->[0] eq "stat") {
	$width= 20;

##	$upper= exists $option->{-upper} ? $option->{-upper}->[0] : 400;
  $width= exists $option->{-width} ? $option->{-width}->[0] : 10;

  open (LOG,">out_manu_fastq_stat.logs");

	if ($fqtype eq "single") {
		$fqtype= "single";
		($info)= &statics_readlength([$fh1],$fqtype,$width);
		print "#READS: $info->{count1}\nTOTAL_BASE: $info->{sum1}\n";
		print "MIN: $info->{min1} MAX: $info->{max1} AVE: $info->{ave1} \n";
		print "LEN_DIST 1\n";

		print "#### READ 1\nLEN\tCOUNT\tRATE\tCUMRate\n";
		$sum_distrate= 0;
		foreach $inter (sort{$a<=>$b} keys %{$info->{len1}}) {
			$inter_count= int($info->{len1}->{$inter}/$info->{count1}*100);
			$inter_rate= sprintf("%.3f",$info->{len1}->{$inter}/$info->{count1}*100);
			$sum_distrate+= $info->{len1}->{$inter}/$info->{count1};
            $sum_distrate= sprintf("%.4f",$sum_distrate);
			print "${inter}-\t$info->{len1}->{$inter}\t".('*' x$inter_count)."  $inter_rate\t$sum_distrate\n";
		}

	} elsif (exists $option->{-pe}) {
		$fqtype= "pair";
		($info)= &statics_readlength([$fh1,$fh2],$fqtype,$width);
		print "#READS: $info->{count1} x2\nTOTAL_BASE: r1 $info->{sum1} r2 $info->{sum2}\n";
		print "MIN: $info->{min1} MAX: $info->{max1} AVE: $info->{ave1} \n";
		print "MIN2: $info->{min2} MAX2: $info->{max2} AVE2: $info->{ave2} \n";
		print "#### READ 1\nLEN\tCOUNT\tRATE\tCUMRate\n";
		$sum_distrate= 0;
		foreach $inter (sort{$a<=>$b} keys %{$info->{len1}}) {
			$inter_count= int($info->{len1}->{$inter}/$info->{count1}*100);
			$inter_rate= sprintf("%.3f",$info->{len1}->{$inter}/$info->{count1}*100);
			$sum_distrate+= $info->{len1}->{$inter}/$info->{count1};
            $sum_distrate= sprintf("%.4f",$sum_distrate);
            print "${inter}-\t$info->{len1}->{$inter}\t".('*' x$inter_count)."  $inter_rate\t$sum_distrate\n";
		}
	}

} elsif ($option->{-e}->[0] eq "sep") {
	($info)= &convert_single_separate([$fh1,$fh2]);
	print "SEQNUM: $info->{count}\n";
	print "OUTPUT: xxxx _R1_001.fastq\n";

} elsif ($option->{-e}->[0] eq "pick") {

	if ($fqtype eq "single") {
		$fqtype= "single";
		($info)= &pickup_reads([$fh1],$option->{-num}->[0],$fqtype);
		print "SEQNUM: $info->{count} GET: $info->{get}\n";
		print "OUTPUT: out_pickup_R1_MP.fq\n";
	} elsif ($fqtype eq "pair") {
		$fqtype= "pair";
		($info)= &pickup_reads([$fh1,$fh2],$option->{-num}->[0],$fqtype);
		print "SEQNUM: $info->{count} GET: $info->{get}\n";
		print "OUTPUT: out_pickup_R1_MP.fq out_pickup_R2_MP.fq\n";
	}

} elsif ($option->{-e}->[0] eq "get") {
   
	($selinfo)= &create_selectinfo($option->{-l}->[0]);
	print "SEL: $selinfo->{selnum}\n";

	if ($fqtype eq "single") {
		$fqtype= "single";
		($info)= &extract_reads([$fh1],$selinfo,$fqtype);
		print "SEQNUM: $info->{count} GET: $info->{get}\n";
		print "OUTPUT: out_extract_PE_R1.fq\n";
	} elsif ($fqtype eq "pair") {
		$fqtype= "pair";
		($info)= &extract_reads([$fh1,$fh2],$selinfo,$fqtype);
		print "SEQNUM: $info->{count} GET: $info->{get}\n";
		print "OUTPUT: out_extract_PE_R1.fq out_extract_PE_R2.fq\n";
	}
} elsif ($option->{-e}->[0] eq "rm") {
    
    ($selinfo)= &create_selectinfo($option->{-l}->[0]);
    print "SEL: $selinfo->{selnum}\n";
    
    if ($fqtype eq "single") {
        $fqtype= "single";
        ($info)= &rm_reads([$fh1],$selinfo,$fqtype);
        print "SEQNUM: $info->{count} GET: $info->{get}\n";
        print "OUTPUT: out_extract_PE_R1.fq\n";
    } elsif ($fqtype eq "pair") {
        $fqtype= "pair";
        ($info)= &rm_reads([$fh1,$fh2],$selinfo,$fqtype);
        print "SEQNUM: $info->{count} GET: $info->{get}\n";
        print "OUTPUT: out_extract_PE_R1.fq out_extract_PE_R2.fq\n";
    }

} elsif ($option->{-e}->[0] eq "inter") {
   
	($selinfo)= &convert_interleaved_fq([$fh1,$fh2]);
	print "SEL: $selinfo->{count}\n";

} elsif ($option->{-e}->[0] eq "trim") {
    @timelen= @{$option->{-fixlen}};
    
    if ($fqtype eq "single") {
        $fqtype= "single";
   #     print "SEQNUM: $info->{count} GET: $info->{get}\n";
        print "OUTPUT: out_trim_SE_R1.fq\n";
    } elsif ($fqtype eq "pair") {
        $fqtype= "pair";
        ($info)= &fixtrim_reads([$fh1,$fh2],$fqtype,[@timelen]);
   #     print "SEQNUM: $info->{count} GET: $info->{get}\n";
        print "OUTPUT: out_trim_PE_R1.fq out_trim_PE_R2.fq\n";
    }

}




#==========================================================================
sub fixtrim_reads {
    my ($filehand,$fqtype,$timelens)= @_;

    my $info= {};
    
    if ( $fqtype eq "single" ) {
        my ($fh1)= @$filehand;
        open (OUT1,">out_trim_PE_R1.fq");
    } elsif ($fqtype eq "pair") {
        my ($fh1,$fh2)= @$filehand;
        open (OUT1,">out_trim_PE_R1.fq");
        open (OUT2,">out_trim_PE_R2.fq");
    }

    ($trimhead,$trimtail)= @{$timelens};

    while (<$fh1>) {
        next if (/^(\#|\s*$)/);
        $name1= $_;
        chomp($name1);
        $base1= <$fh1>;
        chomp($base1);
        $qname1= <$fh1>;
        chomp($qname1);
        $qual1= <$fh1>;
        chomp($qual1);
        $info->{count}++;
        
        $blen= length($base1);
        $netlen= $blen-$trimhead-$trimtail;
        $subseq1= substr($base1,$trimhead,$netlen);
        $subqual1= substr($qual1,$trimhead,$netlen);

        print OUT1 "$name1\n$subseq1\n$qname1\n$subqual1\n";

        if ($fqtype eq "pair") {
            $name2= <$fh2>;
            chomp($name2);
            $base2= <$fh2>;
            chomp($base2);
            $qname2= <$fh2>;
            chomp($qname2);
            $qual2= <$fh2>;
            chomp($qual2);
            
            $subseq2= substr($base2,$trimhead,$netlen);
            $subqual2= substr($qual2,$trimhead,$netlen);
            print OUT2 "$name2\n$subseq2\n$qname2\n$subqual2\n";

        }
        &display_couter($info->{count},[100000,10,3],$info->{count});
    }
    print "\n";
    
    return ($info);


}




sub rm_reads {
    my ($filehand,$selinfo,$fqtype)= @_;
    my ($name,$base,$qname,$qual,$clipseq,$clipqual);
    
    if ( $fqtype eq "single" ) {
        my ($fh1)= @$filehand;
        open (OUT1,">out_extract_PE_R1.fq");
    } elsif ($fqtype eq "pair") {
        my ($fh1,$fh2)= @$filehand;
        open (OUT1,">out_extract_PE_R1.fq");
        open (OUT2,">out_extract_PE_R2.fq");
    }
    
    my $info= {};
    $info->{get}= 0;
    $info->{count}= 0;
    
    while (<$fh1>) {
        next if (/^(\#|\s*$)/);
        $name1= $_;
        chomp($name1);
        $base1= <$fh1>;
        chomp($base1);
        $qname1= <$fh1>;
        chomp($qname1);
        $qual1= <$fh1>;
        chomp($qual1);
        $info->{count}++;
        
        $readname= $1 if ($name1 =~ /^@(\S+)/);
        #print "$readname\n";
        
        unless (exists $selinfo->{read}->{$readname}) {
            print OUT1 "$name1\n$base1\n$qname1\n$qual1\n";
            $info->{get}++;
        }
        if ($fqtype eq "pair") {
            $name2= <$fh2>;
            chomp($name2);
            $base2= <$fh2>;
            chomp($base2);
            $qname2= <$fh2>;
            chomp($qname2);
            $qual2= <$fh2>;
            chomp($qual2);
            unless (exists $selinfo->{read}->{$readname}) {
                print OUT2 "$name2\n$base2\n$qname2\n$qual2\n";
            }
        }
        &display_couter($info->{count},[100000,10,3],$info->{count});
    }
    print "\n";
    
    return ($info);
}

sub extract_reads {
    my ($filehand,$selinfo,$fqtype)= @_;
    my ($name,$base,$qname,$qual,$clipseq,$clipqual);
    
    if ( $fqtype eq "single" ) {
        my ($fh1)= @$filehand;
        open (OUT1,">out_extract_PE_R1.fq");
    } elsif ($fqtype eq "pair") {
        my ($fh1,$fh2)= @$filehand;
        open (OUT1,">out_extract_PE_R1.fq");
        open (OUT2,">out_extract_PE_R2.fq");
    }

    my $info= {};
    $info->{get}= 0;
    $info->{count}= 0;
    
    while (<$fh1>) {
        next if (/^(\#|\s*$)/);
        $name1= $_;
        chomp($name1);
        $base1= <$fh1>;
        chomp($base1);
        $qname1= <$fh1>;
        chomp($qname1);
        $qual1= <$fh1>;
        chomp($qual1);
        $info->{count}++;
        
        $readname= $1 if ($name1 =~ /^@(\S+)/);
        #print "$readname\n";
        
        if (exists $selinfo->{read}->{$readname}) {
            print OUT1 "$name1\n$base1\n$qname1\n$qual1\n";
            $info->{get}++;
        }
        if ($fqtype eq "pair") {
            $name2= <$fh2>;
            chomp($name2);
            $base2= <$fh2>;
            chomp($base2);
            $qname2= <$fh2>;
            chomp($qname2);
            $qual2= <$fh2>;
            chomp($qual2);
            if (exists $selinfo->{read}->{$readname}) {
                print OUT2 "$name2\n$base2\n$qname2\n$qual2\n";
            }
        }
        &display_couter($info->{count},[100000,10,3],$info->{count});
    }
    print "\n";
    
    return ($info);
}



sub pickup_reads {
	my ($filehand,$num,$fqtype)= @_;
	my ($name1,$base1,$qname1,$qual1);
	my ($name2,$base2,$qname2,$qual2);
	my $info= {};
		
    if ( $fqtype eq "single" ) {
        my ($fh1)= @$filehand;
        open (OUT1,">out_pickup_R1_MP.fq");
    } elsif ($fqtype eq "pair") {
        my ($fh1,$fh2)= @$filehand;
        open (OUT1,">out_pickup_R1_MP.fq");
        open (OUT2,">out_pickup_R2_MP.fq");
    }

	while (<$fh1>) {
 		next if (/^(\#|\s*$)/);
		chomp();
		$name1= $_;
		$base1= <$fh1>;
		chomp($base1);
		$qname1= <$fh1>;
		chomp($qname1);
		$qual1= <$fh1>;
		chomp($qual1);
		$info->{count}++;

		&display_couter($info->{count},[100000,10,3],$info->{count});
	
		if ($fqtype eq "single") {
			print OUT1 "$name1\n$base1\n$qname1\n$qual1\n";
			if ($info->{count} < $num) {
				next;
			} elsif ($info->{count} == $num) {
				last;
			}
		}
		
		$name2= <$fh2>;
		chomp($name2);
		$base2= <$fh2>;
		chomp($base2);
		$qname2= <$fh2>;
		chomp($qname2);
		$qual2= <$fh2>;
		chomp($qual2);

		print OUT1 "$name1\n$base1\n$qname1\n$qual1\n";				
		print OUT2 "$name2\n$base2\n$qname2\n$qual2\n";
		$info->{get}++;

		last if ($info->{count} == $num);
	}

	return ($info);
}


#==========================================================================
sub convert_single_separate {
	my ($filehand)= @_;

    my ($fh1)= @$filehand;
    
	open (OUT_F,">out_sep_R1.fastq");
	open (OUT_R,">out_sep_R2.fastq");
	
	$info= {};

	while (<$fh1>) {
 		next if (/^(\#|\s*$)/);
		chomp();
		$name_1= $_; print OUT_F "$name_1\n";
	
		$base_1= <$fh1>;
		chomp($base_1); print OUT_F "$base_1\n";
		$qname_1= <$fh1>;
		chomp($qname_1); print OUT_F "$qname_1\n";
		$qual_1= <$fh1>;
		chomp($qual_1); print OUT_F "$qual_1\n";


		$name_2= <$fh1>;
		chomp($name_2); print OUT_R "$name_2\n";
		$base_2= <$fh1>;
		chomp($base_2); print OUT_R "$base_2\n";
		$qname_2= <$fh1>;
		chomp($qname_2); print OUT_R "$qname_2\n";
		$qual_2= <$fh1>;
		chomp($qual_2); print OUT_R "$qual_2\n";

		$info->{count}++;
	}

	print "OUTPUT: ${headname}_R1_001.fastq ${headname}_R2_001.fastq\n";

	return ($info);
}


sub select_by_length {
	my ($filehand,$lim_len,$fqtype,$lim_gc)= @_;

	my $info= {};
		
	if ($fqtype eq "single") {
        open (OUT3,">out_sellen_R1_OP.fq");
        my ($fh1)= @$filehand;
        
	} elsif ($fqtype eq "pair") {
        open (OUT1,">out_sellen_R1_PE.fq");
        open (OUT2,">out_sellen_R2_PE.fq");
        open (OUT3,">out_sellen_R1_OP.fq");
        open (OUT4,">out_sellen_R2_OP.fq");
        my ($fh1,$fh2)= @$filehand;
	}
	while (<$fh1>) {
 		next if (/^(\#|\s*$)/);
		$name1= $_;
		chomp($name1);
		$base1= <$fh1>;
		chomp($base1);
		$qname1= <$fh1>;
		chomp($qname1);
		$qual1= <$fh1>;
		chomp($qual1);
		$seqlen1= length($base1);

        $info->{count}++;
		&display_couter($info->{count},[100000,10,3],$info->{count});
		
		if ($$lim_gc[0] == 0 && $$lim_gc[1] == 1) {
			$baseinfo1->{gc_rate}= 0.5;
		} else {
			$baseinfo1= &cal_base_count($base1);
		}

		if ($fqtype eq "single") {
			
            if ($seqlen1 >= $lim_len) {
                if ($baseinfo1->{gc_rate} >= $$lim_gc[0] && $baseinfo1->{gc_rate} <= $$lim_gc[1]) {
					print OUT3 "$name1\n$base1\n$qname1\n$qual1\n";
					$info->{get}++;
				}
			}
			next;
		}

		$name2= <$fh2>;
		chomp($name2);
		$base2= <$fh2>;
		chomp($base2);
		$qname2= <$fh2>;
		chomp($qname2);
		$qual2= <$fh2>;
		chomp($qual2);

		$seqlen2= length($base2);
		if ($$lim_gc[0] == 0 && $$lim_gc[1] == 1) {
			$baseinfo2->{gc_rate}= 0.5;
		} else {
			$baseinfo2= &cal_base_count($base2);
		}

        $_=$base1;
        $num_amb_1= tr/ATGC//c;
        $_=$base2;
        $num_amb_2= tr/ATGC//c;

        
        
		if ($seqlen1 >= $lim_len && $seqlen2 >= $lim_len) {
			if ($baseinfo1->{gc_rate} >= $$lim_gc[0] && $baseinfo1->{gc_rate} <= $$lim_gc[1]) {
				if ($baseinfo2->{gc_rate} >= $$lim_gc[0] && $baseinfo2->{gc_rate} <= $$lim_gc[1]) {
                    
                    if ($num_amb_1 == 0 && $num_amb_2 == 0) {
                        print OUT1 "$name1\n$base1\n$qname1\n$qual1\n";
                        print OUT2 "$name2\n$base2\n$qname2\n$qual2\n";
                        $info->{pair}++;
                    } else {
                        $info->{amb}++;
                    }
				}
			}

		} elsif ($seqlen1 >= $lim_len && $seqlen2 < $lim_len) {
			print OUT3 "$name1\n$base1\n$qname1\n$qual1\n";
			$info->{orph_fwd}++
		} elsif ($seqlen1 < $lim_len && $seqlen2 >= $lim_len) {
			print OUT4 "$name2\n$base2\n$qname2\n$qual2\n";
			$info->{orph_rev}++
		}

	}

	return ($info)
}

sub clipping_by_len {
	my ($filehand,$lim_max,$lim_min,$fqtype)= @_;
	my ($name,$base,$qname,$qual,$clipseq,$clipqual);

 
	my $info= {};
	$info->{get}= 0;
	$info->{count}= 0;

	if ($fqtype eq "single") {
        my ($fh1)= @$filehand;
        open (OUT3,">out_cliplen_R1_OP.fq");
    } elsif ($fqtype eq "pair") {
        my ($fh1,$fh2)= @$filehand;
		open (OUT1,">out_cliplen_R1_PE.fq");
		open (OUT2,">out_cliplen_R2_PE.fq");
		open (OUT3,">out_cliplen_R1_OP.fq");
		open (OUT4,">out_cliplen_R2_OP.fq");
	}
    
	while (<$fh1>) {
 		next if (/^(\#|\s*$)/);
		$name1= $_;
		chomp($name1);
		$base1= <$fh1>;
		chomp($base1);
		$qname1= <$fh1>;
		chomp($qname1);
		$qual1= <$fh1>;
		chomp($qual1);
		$info->{count}++;
		$seqlen1= length($base1);
		$subseq1= $seqlen1 >= $lim_max ? substr($base1,0,$lim_max) : $base1;
		$subqual1= $seqlen1 >= $lim_max ? substr($qual1,0,$lim_max) : $qual1;


		if ($fqtype eq "single") {
			if ($lim_min == 0) {
				print OUT3 "$name1\n$subseq1\n$qname1\n$subqual1\n";
				$info->{get}++;
			} elsif (length($subseq1) >= $lim_min) {
				print OUT3 "$name1\n$subseq1\n$qname1\n$subqual1\n";
				$info->{get}++;
			}
			next;
		}

		$name2= <$fh2>;
		chomp($name2);
		$base2= <$fh2>;
		chomp($base2);
		$qname2= <$fh2>;
		chomp($qname2);
		$qual2= <$fh2>;
		chomp($qual2);

		$seqlen2= length($base2);
		$subseq2= $seqlen2 >= $lim_max ? substr($base2,0,$lim_max) : $base2;
		$subqual2= $seqlen2 >= $lim_max ? substr($qual2,0,$lim_max) : $qual2;

		if ($lim_min == 0) {
			print OUT1 "$name1\n$subseq1\n$qname1\n$subqual1\n";
			print OUT2 "$name2\n$subseq2\n$qname2\n$subqual2\n";
			$info->{get}++;
		} else {
			if (length($subseq1) >= $lim_min && length($subseq2) >= $lim_min) {
				print OUT1 "$name1\n$subseq1\n$qname1\n$subqual1\n";
				print OUT2 "$name2\n$subseq2\n$qname2\n$subqual2\n";
				$info->{get}++;
			} else {
				print OUT3 "$name1\n$base1\n$qname1\n$qual1\n" if (length($subseq1) >= $lim_min);
				print OUT4 "$name2\n$base2\n$qname2\n$qual2\n" if (length($subseq2) >= $lim_min);
			}
		}

	}
	return ($info);
}


sub statics_readlength {
	my ($filehand,$fqtype,$width)= @_;
	my ($info,$name1,$base1,$qname1,$qual1,$seqlen1);
	my ($name2,$base2,$qname2,$qual2,$seqlen2);
	
	$info= {};

    open (LOG,">out_manu_fastq_stat.logs");


    $info->{count}= 0;
	$info->{min1}= 100000;
	$info->{max1}= 0;
	$info->{sum1}= 0;
	$info->{min2}= 100000;
	$info->{max2}= 0;
	$info->{sum2}= 0;

    if ( $fqtype eq "single" ) {
        my ($fh1)= @$filehand;
        open (OUT1,">out_readlen_R1.txt");
    } elsif ($fqtype eq "pair") {
        my ($fh1,$fh2)= @$filehand;
        open (OUT1,">out_readlen_R1.txt");
        open (OUT2,">out_readlen_R2.txt");
    }

	while (<$fh1>) {
 		next if (/^(\#|\s*$)/);
		$name1= $_;
		chomp($name1);
		$base1= <$fh1>;
		chomp($base1);
		$qname1= <$fh1>;
		chomp($qname1);
		$qual1= <$fh1>;
		chomp($qual1);
		$info->{count1}++;
		$seqlen1= length($base1);
		$info->{min1}= $seqlen1 if ($seqlen1 < $info->{min1});
		$info->{max1}= $seqlen1 if ($seqlen1 > $info->{max1});
		$info->{sum1}+= $seqlen1;
		$interval= int($seqlen1/$width)*$width;
		$info->{len1}->{$interval}++;

		$name_ori=$1 if ($name1=~ /^@(\S+)/);
		print LOG "$name_ori\t$seqlen1\n";

		if ($fqtype eq "single") {
			next;
		}

		$name2= <$fh2>;
		chomp($name2);
		$base2= <$fh2>;
		chomp($base2);
		$qname2= <$fh2>;
		chomp($qname2);
		$qual2= <$fh2>;
		chomp($qual2);
		$info->{count2}++;
		$seqlen2= length($base2);
		$info->{min2}= $seqlen2 if ($seqlen2 < $info->{min2});
		$info->{max2}= $seqlen2 if ($seqlen2 > $info->{max2});
		$info->{sum2}+= $seqlen2;
		$interval2= int($seqlen2/$width)*$width;
		$info->{len2}->{$interval2}++;

		&display_couter($info->{count1},[100000,10,3],$info->{count1});
	}
	print "\n";

	$info->{ave1}= sprintf("%.2f",$info->{sum1}/$info->{count1});
	$info->{ave2}= sprintf("%.2f",$info->{sum2}/$info->{count2}) if ($fqtype eq "pair");

	return ($info);
}


sub convert_interleaved_fq {
	my ($filehand)= @_;
	my ($name,$base,$qname,$qual,$clipseq,$clipqual);

	my $info= {};
	$info->{count}= 0;

    if ( $fqtype eq "single" ) {
        my ($fh1)= @$filehand;
    } elsif ($fqtype eq "pair") {
        my ($fh1,$fh2)= @$filehand;
    }

    open (OUT,">out_interleaved_PE.fq");

	while (<$fh1>) {
 		next if (/^(\#|\s*$)/);
		$name1= $_;
		chomp($name1);
		$base1= <$fh1>;
		chomp($base1);
		$qname1= <$fh1>;
		chomp($qname1);
		$qual1= <$fh1>;
		chomp($qual1);
		$info->{count}++;
		
#		$readname= $1 if ($name1 =~ /^@(\S+)/);
#		$readname_F= "${readname}/1";
#		print OUT "$readname_F\n$base1\n$qname1\n$qual1\n";

		print OUT "${name1}\n$base1\n$qname1\n$qual1\n";

		$name2= <$fh2>;
		chomp($name2);
		$base2= <$fh2>;
		chomp($base2);
		$qname2= <$fh2>;
		chomp($qname2);
		$qual2= <$fh2>;
		chomp($qual2);

#		$readname_R= "${readname}/2";
#		print OUT "$readname_R\n$base2\n$qname2\n$qual2\n";

		print OUT "${name2}\n$base2\n$qname2\n$qual2\n";

		&display_couter($info->{count},[100000,10,3],$info->{count});
	}
	print "\n";

	return ($info);

}
#==========================================================================

sub display_couter {
  my ($count,$limits,$sign)=@_;
  my ($lim);

  ($subunit,$munit,$rep)= (@$limits);

  print "\-" if ($count%$subunit == 0);
  print "$sign" if ($count%($subunit*$munit) == 0);
  print "\n" if ($count%($subunit*$munit*$rep) == 0);
}


sub create_selectinfo {
	my ($infile)= @_;
	my ($selinfo,$contig,$read);

print "$infile\n";
	$selinfo= {};
	
    if ($infile =~ /.gz$/) {
        open (IN,"gzip -d $infile |") || die "cat't open pipe to $infile";
    } else {
        open (IN,"$infile") || die "cat't open pipe to $infile";
    }

	while (<IN>) {
 		chomp;
 		next if (/^(\#|\s*$)/);
		next if (/^>/);
		@item= split;
		$selinfo->{read}->{$_}++ foreach (@item);
	}

	$selinfo->{selnum}= keys %{$selinfo->{read}};

	return ($selinfo);

}


sub cal_base_count {
	my ($seq)= @_;
	my ($baseinfo);
	my ($gc_num,$gc_cont,$ag_num,$ag_cont);
	
	$baseinfo= {};
	$_=$seq;
	$baseinfo->{a}= tr/Aa//;
	$baseinfo->{t}= tr/Tt//;
	$baseinfo->{g}= tr/Gg//;
	$baseinfo->{c}= tr/Cc//;
	$baseinfo->{n}= tr/[ATGCatgc]//c;
	$baseinfo->{base}= $baseinfo->{a}+$baseinfo->{t}+$baseinfo->{g}+$baseinfo->{c};
	
   $gc_num= $baseinfo->{g}+$baseinfo->{c};
   $gc_cont= $gc_num != 0 ? $gc_num/$baseinfo->{base} : 0;
   $baseinfo->{gc_rate}= sprintf ("%.4f",$gc_cont);
   $ag_num= $baseinfo->{g}+$baseinfo->{a};
   $ag_cont= $ag_num != 0 ? $ag_num/$baseinfo->{base} : 0;
   $baseinfo->{ag_rate}= sprintf ("%.4f",$ag_cont);	

	return ($baseinfo);
}
