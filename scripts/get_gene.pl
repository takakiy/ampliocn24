#!/usr/bin/perl -w

#====================================  OPTIONS  =================================
## REV. 2015.11.12

if ($ARGV[0]=~/\-h/){
     print "-l list file -d Database -o (m)/s\n";
     print "-a seq1 seq2 .... file -d Database -o (m)/s\n";
     exit;
}

foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};

$getinfo= {};

if (exists $option->{-l}) {

} else {
    print "ERROR: SPECIFy LIST\n";
          exit;
}

$getinfo={};
$getinfo->{qrynum}=0;
$getinfo->{getnum}=0;
$getinfo->{nogetnum}=0;


if (exists $option->{-a}) {
    $getinfo->{memb}->{getseq}= $option->{-a};
    $getinfo->{tax}->{$_}++ foreach (@{$option->{-a}});

} elsif (exists $option->{-l}) {
    ($getinfo)= &create_getinfo($option->{-l}->[0],$getinfo);
}

use File::Basename;
($dbbase,$dbdirect,$dbsuf)= fileparse($option->{-d}->[0],('.nuc','.pep','.fas','.fna','.fa','.faa','.fasta'));
print "DBASE $dbbase $dbdirect $dbsuf\n";

open(LOG,">out_getgene.logs.txt");

$outfile_mode= exists $option->{-o} ? $option->{-o}->[0] : "m";

($getinfo)= &get_sequence($getinfo,$option->{-d});

($getinfo)= &output_seq($getinfo,$outfile_mode,$dbsuf);

print "  QUERY: $getinfo->{qrynum}\n  GET: $getinfo->{getnum}\n  NOGET: $getinfo->{nogetnum}\n";


 print "LOG: out_getgene.logs.txt\n";


#=== === === === === === === === === === === === === === === === ===

sub create_getinfo {
    my ($lstfile,$getinfo)= @_;
    
    open ("LI","$lstfile");
    @lines=<LI>;
    $grp= "getseq";
    
    foreach (@lines){
        chomp;
        next if (/^(\#|\s*$)/);
        if (/^\>(\S+)/) {
            $grp= $1;
            print "$grp\n";
       
        } else {
            foreach (split/\s+/) {
                $getinfo->{tax}->{$_}= $grp;
                $getinfo->{qrynum}++;
                push @{$getinfo->{memb}->{$grp}},$_;

            }
        }
    }
    return($getinfo);
}



sub get_sequence {
    my ($getinfo,$dbfiles)= @_;

    $|=1;

    # ($seqinfo)= &readseq_fasta($option->{-d},);
    $seqinfo={};
    $seqc=0;
    
    foreach $fileone (@$dbfiles) {
    #    print "$fileone\n";
        open (SEQ,"$fileone");

        while(<SEQ>){
            chomp;
            next if (/^(\#|\s*$)/);
            
            if (/^>(\S+)/){
                $seqname= $1;
                $index= $_;
                $index=~ s/^>(\S+)\s*//;
                $seqinfo->{each}->{$seqname}++;
                $seqinfo->{total}++;
                if ($seqinfo->{each}->{$seqname} == 1) {
                     $getinfo->{index}->{$seqname}= $index;
                }
                $seqc++;
                &display_counter($seqc,[100000,10,3],$seqc);
              
                if (exists $getinfo->{tax}->{$seqname}) {
                    $getinfo->{seqraw}->{$seqname}=[];
                 #   $getinfo->{getnum}++;
                }
                
                
            }else{

                if (exists $getinfo->{tax}->{$seqname}) {
                    push @{$getinfo->{seqraw}->{$seqname}},$_;
                }
            }
        }
      
    }
    
    print "\n";
    
    $duplinum= 0;
    foreach $seqname(keys %{$seqinfo->{each}}) {
        $duplinum++ if ($seqinfo->{each}->{$seqname} > 1);
        print LOG "DUP\t$seqname\t$seqinfo->{each}->{$seqname}\n" if ($seqinfo->{each}->{$seqname} > 1);

    }
    
    print "DB @$dbfiles\n";
    print "DBNUM $seqinfo->{total} DUPLI $duplinum\n";
    
    return($getinfo);

}


sub output_seq {
    
    my($getinfo,$outfile_mode,$dbsuf)=@_;

    mkdir "EXT" unless (-d "EXT");

    foreach $grp (sort keys %{$getinfo->{memb}}) {
        print "GRP $grp\n";
        
        if ($outfile_mode eq "m") {
            $outfile= "./EXT/${grp}$dbsuf";
            open (OUT,">$outfile");
        }
        
        foreach $seqname (@{$getinfo->{memb}->{$grp}}) {
     #       print "XX $seqname\n";

            if (exists $getinfo->{seqraw}->{$seqname} ) {
                $sequence= join "",@{$getinfo->{seqraw}->{$seqname}};
                $sequence=~ s/\s+//g;
                $sequence=~ tr/a-z/A-Z/;
                $mlseq= &multiline_seq($sequence);
                $seqindex= $getinfo->{index}->{$seqname};

                if ($outfile_mode eq "s") {                    
                    $outfile= "./EXT/${seqname}$dbsuf";
                    $outfile=~ tr/:;/_/;
                    open (OUT,">$outfile");
                }
                if ($seqindex eq "") {
                    print OUT ">$seqname\n$mlseq\n";
                } else {
                    print OUT ">$seqname $seqindex\n$mlseq\n";
                }
                $getinfo->{getnum}++;
                
            } else {
                $getinfo->{nogetnum}++;
                print LOG "NO\t$seqname\n";
            }
        }

    }
    
    return ($getinfo);
}









#==========================================================================
#=============================DISPLAY COUNTER =============================
sub display_counter {
  my ($count,$limits,$sign)=@_;
  my ($lim);

  ($subunit,$munit,$rep)= (@$limits);

  print "\-" if ($count%$subunit == 0);
  print "$sign" if ($count%($subunit*$munit) == 0);
  print "\n" if ($count%($subunit*$munit*$rep) == 0);
}


sub readseq_fasta {
	my ($seqfile)= @_;
	my ($seqinfo,$seqname,$fileone,@lines,$index,$dupinfo,$onefinfo);
  
	$seqinfo={};
	$dupinfo= {};
	$seqc= 0;

	foreach $fileone (@$seqfile) {
		open (SEQ,"$fileone");
		@lines= <SEQ>;
		$onefinfo= {};
		foreach(@lines){
			chomp;
			next if (/^(\#|\s*$)/);        
			if (/^>(\S+)/){
				$seqname=$1;
				$index= $_;
				$dupinfo->{$seqname}++;
				$dupinfo->{total}++;
				if ($dupinfo->{$seqname} == 1) {
					push @{$onefinfo->{tax}},$seqname;
					$index=~ s/^>(\S+)\s*//;
					$seqinfo->{index}->{$seqname}= $index;
				}
				$seqc++;
				&display_counter($seqc,[50000,10,3],$seqc);
			}else{
				push @{$seqinfo->{seqline}->{$seqname}},$_ if ($dupinfo->{$seqname} == 1);
			}
		}
		
		foreach $seqname (@{$onefinfo->{tax}}) {
			if (exists $seqinfo->{seqline}->{$seqname}) {
				$seqinfo->{taxnum}++;
				push @{$seqinfo->{tax}},$seqname;
			} else {
				delete $seqinfo->{index}->{$seqname};
				print "   NOSEQ $seqname\n";
			}							 
		}
		
		print "FILE $fileone SEQ_NUMBER  $seqinfo->{taxnum}\/$dupinfo->{total}\n";
	}
  
	foreach $seqname (sort keys %{$dupinfo}) {
		next if ($seqname=~ /total/);
		print "  DUPLI $seqname $dupinfo->{$seqname}\n" if ($dupinfo->{$seqname} > 1);
	}
	
	($seqinfo)= &join_seqline($seqinfo);
  

	return ($seqinfo);
  
}
sub join_seqline {
	my ($seqinfo)= @_;
	my ($sequence);
  
	foreach $seqname (@{$seqinfo->{tax}}) {
		$sequence= join "",@{$seqinfo->{seqline}->{$seqname}};
		$sequence=~ s/\s+//g;
		$sequence=~ tr/a-z/A-Z/;
	#	$sequence=~ tr/U/T/;
		$seqinfo->{seq}->{$seqname}= $sequence;
		$seqinfo->{len}->{$seqname}= length $sequence;
		delete $seqinfo->{seqline}->{$seqname};
		# print "$seqname $seqinfo->{len}->{$seqname}\n";
	}
	return ($seqinfo);
}

sub multiline_seq {
	my $seq= shift @_;
	my $unit = defined $_[0] ? $_[0] : 60;
	my ($len,@subs,$subseq,$x,$line_seq);
	$len=length $seq;
	@subs=(); 
	for ($x= 0; $x < $len; $x += $unit){
		$subseq=substr $seq,$x,$unit;
		push @subs,$subseq;
	}
	$line_seq=join "\n", @subs;
	return $line_seq;
}


#==========================================================================















