use strict;
#use lib "/xcommon/bin/_mdpa";
#use lib "/Users/dming/mdpa/_mdpa";
use lib "/opt/xcommon/bin/_mdpa";
use _FILEprocc;
use _PDBprocc;

package _gMSMSpoints;

sub generate_MSMSpoints{
    my ($pdbid0,$pdbfile,$selechain,$structuredatadir,$msmsdir,$workdir,
	$msms_density,$msms_prob,$msmslayer,$flag_msmsdense,$flag_msms_saveall,
	$msmssurff,$msmssurfpdbf,$flag_check,$flag_quite
	)=@_;
    

    my $peptidepdbf=$pdbid0.'_peptide.pdb';        #temp file: peptide chain atom crd
    my $natom= &_PDBprocc::read_pdb($structuredatadir,$pdbfile,$selechain,'getpeptidepdb',$peptidepdbf,$workdir);
    if($natom <= 0){
	print "INPUT PDB FILE CONTAINS NOTHING TO BUILD A PEPTIDE\n";
	return -9999;
    }
    system ("rm -f msms.exe pdb_to_xyzr atmtypenumbers");
    system "cp $msmsdir/msms.exe .";
    system "cp $msmsdir/pdb_to_xyzr .";
    system "cp $msmsdir/atmtypenumbers .";

    
    my $nredundance=1;
    if($flag_msmsdense eq 'double' or $flag_msmsdense eq '2') {
	$nredundance=2;
    }
    elsif($flag_msmsdense eq 'triple' or $flag_msmsdense eq '3') {
	$nredundance=3;
    }
    elsif($flag_msmsdense eq 'tetra' or $flag_msmsdense eq '4') {
	$nredundance=4;
    }
    else {
	if( $flag_msmsdense eq 'true') {
	    $nredundance=9999;
	}
	elsif ( $flag_msmsdense eq 'false') {
	}
	else {
	    print "WARNING> UNRECOGNIZABLE OPTION FOR MSMSDENSE, IGNORED\n";   
	}
    }

    my %surfcrd=(); my $nsurf=0;
    for (my $i=1; $i <= $msmslayer; $i++) {  #check multiple layers, suppose.	
	my $pdbid="_tmp_msms_".$pdbid0."_".$i;	
	my $protf=$pdbid."_peptide.pdb";
	if($i eq 1) {`cp -f $workdir/$peptidepdbf $workdir/$protf`};
	my $xyzrf=$pdbid.'.xyzr';
	my $vertf=$pdbid.'.vert';
	my $facef=$pdbid.'.face';
	my $msmsfhead=$pdbid;
	my $msmsoutf=$pdbid."_msms.out";

	print  "MSMS> ./msms.exe -if $xyzrf -de $msms_density  -prob $msms_prob -of $msmsfhead -no_header > $msmsoutf \n" if ($flag_check eq 'true');
	system "./pdb_to_xyzr $workdir/$protf  > $xyzrf";
	system "./msms.exe -if $xyzrf -de $msms_density  -prob $msms_prob -of $msmsfhead -no_header > $msmsoutf "; 
	system ("rm -f $xyzrf  $msmsoutf") if ( $flag_check ne 'true');

	my ($x,$y,$z,$nx,$ny,$nz,$index1,$index2);
	my $sid=0; my $sid0=-1;
	my $redundance='';                 #redundance check

	my %index=(); undef %index;                 #recording surfpoint index 
	undef %surfcrd;                             #possible output data
	open(vf,"<$vertf") || die "gMSMS_ERR_$pdbid0> NO VERT GENERATED!\n";
	while(<vf>){
	    $sid0++;
	    my @vl=split /\s+/,$_;
	    $x=@vl[1];$y=@vl[2];$z=@vl[3];$nx=@vl[4];$ny=@vl[5];$nz=@vl[6];
	    $index1=@vl[7];$index2=@vl[8];


	    $redundance=0;             #redundance check
	    if( exists $index{$index2}) {
		$redundance=$index{$index2}
	    }

	    next if $redundance >=$nredundance;       #filer out most of verts    
	    next if($i < $msmslayer && $redundance >=1 );

	    $index{$index2}++;
	    push (@{$surfcrd{$sid}},$x,$y,$z,$nx,$ny,$nz,$index1,$index2); 
	    $sid++; 
	}
	close(vf);
	$nsurf=$sid;
#
#----
	if($i < $msmslayer){
	    my $j=$i+1;
 	    my $nextprotf="_tmp_msms_".$pdbid0."_".$j."_peptide.pdb";
	    &append_to_pdbf(\%surfcrd,$protf,$nextprotf,$workdir);
#	    system "rm -f $xyzrf $vertf  $facef $msmsoutf $protf" if $flag_check ne 'true';
	    my $tmpsurff=$pdbid0."_".$i.'.msms';
	    &write_surfdata($tmpsurff,\%surfcrd,$workdir) if $flag_msms_saveall eq 'true';

	    if($msmssurfpdbf){
		my $tmpsurfpdbf=$pdbid0."_msms_".$i.'.pdb';
		&swritepdbf($tmpsurfpdbf,\%surfcrd,'S',$workdir) if ($flag_msms_saveall eq 'true');
	    }
	}
	system "rm -f $xyzrf $vertf $facef $workdir/$protf $msmsoutf" if $flag_check eq 'false';

    }
    system ("rm -f $workdir/$peptidepdbf") if $flag_check ne 'true';
    system ("rm -f msms.exe pdb_to_xyzr atmtypenumbers");

    &write_surfdata($msmssurff,\%surfcrd,$workdir) if $msmssurff;
    &swritepdbf($msmssurfpdbf,\%surfcrd,'S',$workdir) if $msmssurfpdbf;
    
    return $nsurf;   #Data output instead of surfcrd file output 
}



sub write_surfdata {
    
    my ($surff,$surfcrd,$output_datadir)=@_;
    open(SF,"> $output_datadir/$surff") || die "CAN NOT OPEN FOR WRITE:  $output_datadir/$surff\n";
    my $nsurf=0;
    foreach my $sid (sort {$a<=>$b} keys %$surfcrd) {
	my $x=${$$surfcrd{$sid}}[0]; my $y=${$$surfcrd{$sid}}[1]; my $z=${$$surfcrd{$sid}}[2];
	my $nx=${$$surfcrd{$sid}}[3]; my $ny=${$$surfcrd{$sid}}[4]; my $nz=${$$surfcrd{$sid}}[5]; 
	my $index1=${$$surfcrd{$sid}}[6]; my $index2=${$$surfcrd{$sid}}[7];
	printf SF " %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f %6d %6d\n",$x,$y,$z,$nx,$ny,$nz,$index1,$index2;
	$nsurf++;
    }
    close(SF);
    #print "$nsurf points were written into $output_datadir/$surff\n";
}



sub swritepdbf{
    my $wpdbfilename=@_[0];
    my $crdhash=@_[1];
    my $chainid=@_[2];
    my $output_dir=@_[3];
    open(wp,"> $output_dir/$wpdbfilename");
    my $lid=0; my $resName='SUR';
    foreach my $hid (sort {$a<=>$b} keys %$crdhash) {
	my $x=@{$$crdhash{$hid}}[0];
	my $y=@{$$crdhash{$hid}}[1];
	my $z=@{$$crdhash{$hid}}[2];
	$lid++; my $resSeq=substr($lid,0,4); 
	printf wp "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s\n", 
	'ATOM  ',$lid,' SUR',' ',$resName,$chainid,$resSeq,' ',$x,$y,$z,1.0,10.0,'','SURF','','';
    }
    close(wp);
}


sub append_to_pdbf{
    my ($surfcrd,$protf,$newprotf,$workdir)= @_;
    system "cp -f $workdir/$protf $workdir/$newprotf";

    my $pdbformat="%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%5d  %8.3f%8.3f%8.3f\n";
    my $i=&_FILEprocc::lineoffile('.',$newprotf);
    if( -e "$workdir/$newprotf") {          #wcrd -- writing crd
	open(WCRD,">> $workdir/$newprotf") || die "CAN NOT OPEN FOR WRITE:  $workdir/$newprotf\n";}
    else {
	open(WCRD,"> $workdir/$newprotf") || die "CAN NOT OPEN FOR WRITE:  $workdir/$newprotf\n";}    
    printf WCRD "TER\n" if $i > 0;
    
    foreach my $sid (sort {$a<=>$b} keys %$surfcrd) {
	my $x=${$$surfcrd{$sid}}[0]; my $y=${$$surfcrd{$sid}}[1]; my $z=${$$surfcrd{$sid}}[2];
	my $nx=${$$surfcrd{$sid}}[3]; my $ny=${$$surfcrd{$sid}}[4]; my $nz=${$$surfcrd{$sid}}[5]; 
	my $index1=${$$surfcrd{$sid}}[6]; 
	my $index2=${$$surfcrd{$sid}}[7];
	$i++;
	my $rid=substr($i,0,4); 
	printf WCRD $pdbformat,'ATOM  ',$i,' CA ',' ','ALA',' ',$rid," ",$x,$y,$z,1.0,10.0,'SURF',$index2,$nx,$ny,$nz;
	#printf $pdbformat,'ATOM  ',++$i,' CA ',' ','ALA',' ',$i," ",$x,$y,$z,1.0,10.0,'SURF',$index2,$nx,$ny,$nz;
    }
    close(WCRD);
    
#    print "---> $protf,$workdir/$newprotf \n";    exit;

}


1;
