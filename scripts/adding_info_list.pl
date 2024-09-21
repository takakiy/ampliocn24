#!/usr/bin/perl -w

if ((@ARGV == 0) || ($ARGV[0]=~/-h/)) {
	print " -i original file -a add info file -key 0 0 (-i -a) -val 2..4 ###   If adding all, -val no specified\n";

	print " -i xx -b yyy -key 0 0 -val 1 2  # -key (original add) Simple extract of tophit\n";
	print " -i xx -bs yyy -key 0 0 -val 1 2  # Extract of tophit under selected species\n";

	exit;

}

foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};

open (LOG,">out_adding.logs");


if (exists $option->{-b}) {
	($tag_item,$key_item)= exists $option->{-key} ? @{$option->{-key}} : (0,0);

	$val_item= exists $option->{-val} ? $option->{-val} : [3,4,5,6,7,8,9,10,11];
	$val_num= @$val_item;
	($blastinfo)= &create_blastinfo_top($option->{-b}->[0],$key_item,$val_item);
	$blastinfo->{val_num}= $val_num;
	print "VAL_NUM: $val_num\n";

	&output_addinfo($option->{-i}->[0],$blastinfo,$tag_item);

} elsif (exists $option->{-bs}) {
	($tag_item,$key_item)= exists $option->{-key} ? @{$option->{-key}} : (0,0);

	($blastinfo)= &create_blastinfo_spec($option->{-bs}->[0],$key_item);

	&output_addinfo($option->{-i}->[0],$blastinfo,$tag_item);
} elsif (exists $option->{-p}) {
	($taxinfo,$continfo)= &create_taxinfo_orf($option->{-p}->[0]);
	&output_addinfo($option->{-i}->[0],$taxinfo);
	&output_addinfo($option->{-i}->[1],$continfo);
}elsif (exists $option->{-t}) {
	($taxinfo)= &create_taxinfo($option->{-t}->[0]);
	&output_addinfo($option->{-i}->[0],$taxinfo);

}elsif (exists $option->{-c}) {
	($taxinfo)= &create_tmpinfo($option->{-c}->[0]);
	&output_addinfo($option->{-i}->[0],$taxinfo);
}elsif (exists $option->{-s}) {
	($seedinfo)= &create_seedinfo($option->{-s}->[0]);
	&output_addinfo($option->{-i}->[0],$seedinfo);	
}elsif (exists $option->{-a}) {   ####   Simple adding
	($tag_item,$key_item)= exists $option->{-key} ? @{$option->{-key}} : (0,0);
    print "TEMP $tag_item ADD $key_item\n";
	$val_item= exists $option->{-val} ? $option->{-val} : [-1];
# print "XXX$option->{-val}->[3]\n";
    ($addinfo)= &create_addinfo($option->{-a}->[0],$key_item,$val_item);
	&output_addinfo($option->{-i}->[0],$addinfo,$tag_item);	
}


###================================================================================
sub create_addinfo {
	my ($infile,$key_item,$val_item)= @_;
	my ($val_num,@item,@values);

	open (IN,"$infile");

	$val_num= 0;

#@valsxx= @$val_item;
print "values_item @$val_item\n";
 

	my $addinfo= {};
	while (<IN>) {
		chomp;
 		next if (/^(\#|\s*$)/);
		@item= split/\t/;
		@values= ();

		if ($$val_item[0] >= 0) {
			foreach (@$val_item) {
				if (defined $item[$_]) {
					push @values,$item[$_];
				} else {
					push @values,' -- ';
				}
			}
			$val_num= @$val_item;
		} else {
			@values= @item;
			$val_num= @item;
		}
# $xnum= @item;
# print "$xnum\n";
		$addinfo->{each}->{$item[$key_item]}->{info}= join "\t",@values;

		$addinfo->{val_num}= $val_num;
	}

print "$val_num VAL @$val_item\n";
#exit;
	return ($addinfo);
}

sub create_seedinfo {
	my ($infile)= @_;

	open (IN,"$infile");

	my $seedinfo= {};
	while (<IN>) {
		chomp;
 		next if (/^(\#|\s*$)/);
		@item= split/\t/;
		($name)= shift @item;
		foreach (@item) {
			s/\s+//;
			$seedinfo->{each}->{$_}->{info}= $name if (/\w+/);
#print "$_\n";
		}
	}

	return ($seedinfo);
}

sub create_tmpinfo {
	my ($infile)= @_;

	open (IN,"$infile");

	my $taxinfo= {};
	while (<IN>) {
		chomp;
 		next if (/^(\#|\s*$)/);
		($name,$taxone)= ((split/\t/)[1],(join " ",(split/\t/)[1,2,3,4]));
		$taxone=~ s/^\s+//;

		$taxinfo->{each}->{$name}->{info}= $taxone;
	}

	return ($taxinfo);
}

sub create_blastinfo_top {
	my ($infile,$key_item,$val_item)= @_;

	open (IN,"$infile");
	my $blastinfo= {};
	$info= {};
	while (<IN>) {
		chomp;
 		next if (/^(\#|\s*$)/);
		($obj,$sbj)= (split/\t/)[$key_item,3];
		$info->{$obj}->{hit}++;
		$info->{$obj}->{info}= [] if ($info->{$obj}->{hit} == 1);

		push @{$info->{$obj}->{info}},$_ if ($sbj !~ /^no_/);
	}

	@val_item_s= @$val_item;

	foreach $obj (sort keys %{$info}) {

		@top_hit= ();

		if (@{$info->{$obj}->{info}} > 0) {
			@hitsort= sort{(split/\t/,$b)[6] <=> (split/\t/,$a)[6]}@{$info->{$obj}->{info}};  # sort by score
			@hitones= split/\t/,$hitsort[0];
			push @top_hit,$hitones[$_] foreach (@val_item_s);
			$blastinfo->{each}->{$obj}->{info}= join "\t",(@top_hit);
			
		} else {
			push @top_hit," --" foreach (0..$#val_item_s);
			$blastinfo->{each}->{$obj}->{info}= join "\t",(@top_hit);
		}
	}
	return ($blastinfo);
}

sub create_blastinfo_spec {
	my ($infile,$key_item)= @_;

	open (IN,"$infile");
	my $blastinfo= {};
	$info= {};
	while (<IN>) {
		chomp;
 		next if (/^(\#|\s*$)/);
		($obj,$sbj)= (split/\t/)[$key_item,3];
		$info->{$obj}->{hit}++;
		$info->{$obj}->{info}= [] if ($info->{$obj}->{hit} == 1);

		push @{$info->{$obj}->{info}},$_ if ($sbj !~ /^no_/);
	}

	foreach $obj (sort keys %{$info}) {
		($sbj,$s_len,$s_desc,$eval,$qcov,$scov,$ident,$simi,$orga)= (" --"," --"," --"," --"," --"," --"," --"," --"," --"," --");

		@top_hit= (" --"," --"," --"," --"," --"," --"," --"," --"," --"," --");
		@top_select= (" --"," --"," --"," --"," --"," --"," --"," --"," --"," --");
		$select= 0;

		$rank= 0;
print LOG "$obj\n";
		if (@{$info->{$obj}->{info}} > 0) {
			@hitsort= sort{(split/\t/,$b)[6] <=> (split/\t/,$a)[6]}@{$info->{$obj}->{info}};  # sort by score
			$hitnum= @hitsort;
print LOG "  $hitnum ==> ";

			$select= 0;
			@top_hit= ();
			@top_select= ();
			foreach $hitone (@hitsort) {

				($sbj,$s_len,$s_desc,$score,$eval,$qcov,$scov,$ident,$simi)= (split/\t/,$hitone)[3,4,5,6,7,8,9,10,11];
				$s_desc=~ s/\n+/ /g;
				$orga= " --";
				$_= $s_desc;
				@orgas= ();
				@orgas= /\[(.+?)\]/g;
				$orga= $orgas[-1] if (@orgas > 0);

				$rank++;
				@top_hit= ($sbj,$s_len,$s_desc,$eval,$qcov,$scov,$ident,$simi,$orga,$rank) if ($rank == 1);

				$select++ if ($sbj =~ /^(bse|vok|rma|cpha)\:/);
				@top_select= ($sbj,$s_len,$s_desc,$eval,$qcov,$scov,$ident,$simi,$orga,$rank) if ($select == 1);

print LOG " $select $sbj $score ";

				last if ($select == 1);
			}
			
			if ($select == 1) {
				$blastinfo->{each}->{$obj}->{info}= join "\t",(@top_hit,@top_select);
			} else {
				$blastinfo->{each}->{$obj}->{info}= join "\t",(@top_hit);
			}

		} else {
			$blastinfo->{each}->{$obj}->{info}= join "\t",(@top_hit);
		}


print LOG "\n";
	}
	return ($blastinfo);
}

sub create_taxinfo {
	my ($infile)= @_;

	open (IN,"$infile");

	my $taxinfo= {};
	while (<IN>) {
		chomp;
 		next if (/^(\#|\s*$)/);
		($name,$taxone)= split/\t/;
		$taxone=~ s/^\s+//;
#print "$taxone\n";
		$taxinfo->{each}->{$name}->{info}= $taxone;
	}

	return ($taxinfo);

}
sub create_taxinfo_orf {
	my ($infile)= @_;

	open (IN,"$infile");
	my $taxinfo= {};
	while (<IN>) {
		chomp;
 		next if (/^(\#|\s*$)/);
		($name,$taxone)= split/\t/;
		$taxone=~ s/^\s+//;
#print "$name\n";

		$taxinfo->{each}->{$name}->{info}= $taxone;
		$cont= $1 if ($name=~/^(\S+)_/);
		push @{$taxinfo->{cont}->{$cont}},$taxone;
	}
	my $continfo= {};
	foreach $cont (sort keys %{$taxinfo->{cont}}) {
		$orfnum= @{$taxinfo->{cont}->{$cont}};
		$continfo->{each}->{$cont}->{info}= $orfnum."\t".(join " \/\/ ",@{$taxinfo->{cont}->{$cont}});

	}
	return ($taxinfo,$continfo);

}


###================================================================================

sub output_addinfo {
	my ($infile,$addinfo,$tag_item)= @_;

	open (IN,"$infile");
	$fname=(split/\//,$infile)[-1];

	open (OUT,">out_${fname}");
  print "OUTFILE out_${fname}\n";
	
	$val_num= exists $addinfo->{val_num} ? $addinfo->{val_num} : 0;

    
	while (<IN>) {
		chomp;
 		if (/^(\#)/) {
			print OUT "$_\n";
		} elsif (/^(\#|\s*$)/) {
			next;
		} else {
			$line= $_;
			($name)= (split/\t/)[$tag_item];
			if (exists $addinfo->{each}->{$name}->{info}) {
				$addinfotxt= $addinfo->{each}->{$name}->{info};
			} else {
				$nulltxt= ' --  ' x$val_num;
				
				$addinfotxt= join "\t",(split/  /,$nulltxt);
			}
			print OUT "$line\t$addinfotxt\n";

		}
	}

}


