package _PDBprocc;
use strict;

sub read_pdb{
    my ($pdbdatadir,$pdbf,$selechain,$flag_read,$recordf,$outputdir,$flag_check)=@_;
#
    $selechain='FirstChain' if ($selechain eq '' || $selechain eq 'undefine' || $selechain eq 'unselectchain');
    my $blankchain;
    if($selechain =~ /BLANK/i) {
	$blankchain=' ';
	$selechain =~ s/BLANK//gi;
    } 
    print "READ_PDB> $selechain ---- $flag_read\n" if ($flag_check eq "true");

    my ($nca,$nhet,$natom,$nenv); $nca=0;$nhet=0;$natom=0;$nenv=0;
    my (%cacrd,%hetcrd,%envcrd,@envid); %cacrd=();%hetcrd=();%envcrd=();@envid=(); 
    undef %cacrd;undef %hetcrd;undef %envcrd; undef @envid;
    my ($atom,$serial,$atomname,$altloc,$resname,$chainid,$resseq);
    my ($icode,$x,$y,$z,$occupancy,$tempfactor,$segid,$element,$charge);
    my ($resseq0,$icode0);$resseq0='null';$icode0='null'; my $atmtype='null';
    my @revisedres=('CGU','MSE','CME','CSS','KCX','TRO','SEP',' FE','OXY');
    my $hetatm_name='***';
    my $wcrdformat=' %8.3f %8.3f %8.3f %3s %1s %5d %1s %6.2f';
    #print wcrd "$x $y $z $resname $chainid $resseq $icode $tempfactor"
    open(wcrd,">  $outputdir/$recordf") if($recordf ne '' and $recordf ne 'none');
    open(pf,"<  $pdbdatadir/$pdbf");
    while (<pf>){
	$atom=substr($_,0,6);
	$serial=substr($_,6,5);
	$atomname=substr($_,12,4); 
	$atmtype=substr($_,13,1); $atmtype='X' if $atmtype eq ' ';
	$altloc=substr($_,16,1);
	$resname=substr($_,17,3);
	$chainid=substr($_,21,1);
	$resseq=substr($_,22,4);
	$icode=substr($_,26,1);
	$x=substr($_,30,8);
	$y=substr($_,38,8);
	$z=substr($_,46,8);
	$occupancy=substr($_,54,6);
	$tempfactor=substr($_,60,6);
	$segid=substr($_,72,4);
	$element=substr($_,76,2);
	$charge=substr($_,78,2);
	last if($atom eq "ENDMDL");
	next if($atom ne "ATOM  " and $atom ne 'HETATM');
	$selechain=$chainid if $selechain eq 'FirstChain';

	if(($flag_read eq 'getnca' or $flag_read eq 'getcacrd') and ($atom eq "ATOM  " or grep (/$resname/, @revisedres)) and  $atomname eq ' CA '){
	    if($resseq ne $resseq0 or $icode ne $icode0) {
		$resseq0=$resseq;
		$icode0=$icode;
		if ($selechain =~ /$chainid/i or ( $chainid eq $blankchain) or $selechain =~/ALL/i){
		    $nca++   if($flag_read eq 'getnca' or $flag_read eq 'getcacrd');
		    push (@{$cacrd{$nca}},$x,$y,$z,$resname,$chainid,$resseq,$icode,$tempfactor) if($flag_read eq 'getcacrd'); 
		    printf wcrd "$wcrdformat\n",$x,$y,$z,$resname,$chainid,$resseq,$icode,$tempfactor if($recordf);
		}
	    }
	}
	if(($flag_read eq 'getnhet' or $flag_read eq 'gethetcrd') and ($atom eq 'HETATM' and ( $hetatm_name eq '***' or  $resname eq $hetatm_name) and ! grep (/$resname/, @revisedres) and $resname ne 'HOH')){
	    if ($selechain =~ /$chainid/i or $chainid eq $blankchain or $selechain =~/ALL/i){
		$nhet++   if($flag_read eq 'getnhet' or $flag_read eq 'gethetcrd');
		push (@{$hetcrd{$nhet}},$x,$y,$z,$resname,$chainid,$resseq,$icode,$tempfactor) if($flag_read eq 'gethetcrd' or $flag_read eq 'getall');
		printf wcrd "$wcrdformat\n",$x,$y,$z,$resname,$chainid,$resseq,$icode,$tempfactor if($recordf ne '' and $recordf ne 'null');
	    }
	}
	if(($flag_read eq 'getnenv' or $flag_read eq 'getenvcrd') and ($atom eq "ATOM  " or grep (/$resname/, @revisedres)) and $atmtype ne 'H'){
	    if ($selechain =~ /$chainid/i or $chainid eq $blankchain or $selechain =~/ALL/i){
		$nenv++;
		push(@envid,$serial);
		push (@{$envcrd{$nenv}},$x,$y,$z,$resname,$chainid,$resseq,$icode,$tempfactor) if($flag_read eq 'getenvcrd');
		printf wcrd "$wcrdformat\n",$x,$y,$z,$resname,$chainid,$resseq,$icode,$tempfactor if($recordf);
	    }
	}
	if($flag_read eq 'getpeptidecrd' and ($atom eq "ATOM  " or grep (/$resname/, @revisedres))){
	    if ($selechain =~ /$chainid/i or $chainid eq $blankchain or $selechain =~/ALL/i){
		$natom++;
		printf wcrd "$wcrdformat\n",$x,$y,$z,$resname,$chainid,$resseq,$icode,$tempfactor if($recordf ne '' and $recordf ne 'null');
	    }
	}
	if($flag_read eq 'getpeptidepdb' and ($atom eq "ATOM  " or grep (/$resname/, @revisedres))){
	    if ($selechain =~ /$chainid/i or $chainid eq $blankchain or $selechain =~/ALL/i){
		print wcrd "$_";
		$natom++;
	    }
	}
    }
    close(pf);
    close(wcrd) if($recordf ne '' and $recordf ne 'none');
    return $nca  if($flag_read eq 'getnca');
    return $nca,\%cacrd  if($flag_read eq 'getcacrd');
    return $nhet  if($flag_read eq 'getnhet');
    return $nhet,\%hetcrd  if($flag_read eq 'gethetcrd');
    return $nenv  if($flag_read eq 'getnenv');
    return $nenv,\%envcrd  if($flag_read eq 'getenvcrd');
    return $natom if($flag_read eq 'getpeptidecrd' or $flag_read eq 'getpeptidepdb');
}



1;
