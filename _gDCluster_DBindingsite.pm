#!/usr/bin/perl
#
#-- Based on Backbone-enhanced ENM
#-- Generating DPA data for PDB protein structure
#     To predict Protein/Ligand interaction
#-- originally writen by DENGMING MING 
#     in Los Alamos National Lab, 2004
# 
#   ADD MULTIPLE-CHAIN PROTEINS ("-chain" option)
#   ADD MULTIPLE-LAYER OF MSMS-POINTS, LAYER 2 BUILT ON LAYER 1, AND SO ON ("-layer" option)
#   ADD REDUCED NMODES ("-nmodes" option)
#     by DENGMING MING
#     in Fudan University, 2013/11/29
#
#   ADD HIEARCHICAL CALC., I.E., GENERTING L_DPA-CLUSTERS FOR EACH (L)LEVEL,
#       AND FIND CONNECTION BETWEEN 1-, 2-,...,L-DPA-CLUSTERS
#     in Nanjing Tech University, 2021/05/21
#
#
use strict;
#use lib "/xcommon/bin/_mdpa";
#use lib "/Users/dming/mdpa/_mdpa";
use lib "/opt/xcommon/bin/_mdpa";
use Sys::Hostname;
use _FILEprocc;
use _PDBprocc;
use _gDPA;
use _anaDPA;
package _gDCluster_DBindingsite;
#
my $wformat=' %3d %3s %1s %5d %1s';
my $wformat2='%4s %3d %3s %1s %5d %1s';
#
sub protein_dpa{
    my ($pdbid,$selechain,
    $cutcc,$delta_cutcc,$delta_cutlg,$cutlg,$cutll,$intcutcc,$wcc,$wctc,$wlg,$wll,$intwcc,
    $nmodes,$modekappa,$ev_threshold,$flag_enforceDPA,$flag_udpa,
    $msms_density,$msms_prob,$msmslayer,
    $flag_msmsdense,$flag_msms_saveall,$flag_wmsms,$flag_wspdb,$flag_uMSMS,$input_msmsf,
    $epislon_dpa,$MinPts_dpa,$epislon_cluster_dpa,$cutdpaperct,$cutdpaperct2,$adjcutperct,
    $ndpavalue_section,$topdpa_num_cutoff,$dpabindingcutoff,$flag_flexible_bcut,$flag_wcpdb,
    $flag_check,$flag_task,$flag_quite,
    $home,$scrh,$workdir,$exedir,$msmsdir,$structuredatadir,$outputdatadir,$dpadir,$totproccf
    ) = @_;
    my $jobid=$pdbid."_L".$msmslayer;
    my $proccf=$jobid.'_DPAproccessing.log';
#
#-- Generate C-alpha
    my $pdbfile=$pdbid.'.pdb';  
    my $cacrdf=$pdbid.'_ca.crd';
    &processing_record($flag_quite,$outputdatadir,$proccf,$totproccf,'ini');
    &processing_record($flag_quite,$outputdatadir,$proccf,$totproccf,-9999,"DPA_ERR> \"$pdbfile\" -- not in \"$structuredatadir\", using option -p for input-file, -h for help") if(! -e "$structuredatadir/$pdbfile");
    my ($nca,$cacrd)=&_PDBprocc::read_pdb($structuredatadir,$pdbfile,$selechain,'getcacrd',$cacrdf,$workdir,$flag_check);
    my ($nchain,$chainlist)=&find_chain($cacrd);
    print "DPA_CALC_$pdbid> CHAIN $nchain :: @{$chainlist} | LAYER $msmslayer\n";
    
    my $natom=&_PDBprocc::read_pdb($structuredatadir,$pdbfile,$selechain,'getpeptidecrd',$flag_check);
    &processing_record($flag_quite,$outputdatadir,$proccf,$totproccf,-1000,"DPA_ERR>  $pdbfile -- TOO FEW RESIDUES, CHAIN MIGHT NOT EXIST")  if($nca<=5);
    &processing_record($flag_quite,$outputdatadir,$proccf,$totproccf,0,"DPA_WARNING> $jobid -- TOO FEW ATOMS: NCA_$nca NCHAIN_$nchain -- $natom -- possibly broken structure") if($natom/$nca < 2.);
#
#
#-- Generate DPA 
    my $dpafile=$jobid.'.dpa1'; 
    my $calc_info=0; my $infocalc=' using DPA'; my $cutcc_work=-1;

    if($flag_udpa eq 'true' && -e "$dpadir/$dpafile") {
	#print "DPA_$jobid> USING DPA: $dpadir/$dpafile\n" if ! $flag_quite;
        #-- write MSMS point-crds which stored in given dpa-file: $dpafile 
	&_gDPA::writemsms_w_dpa($jobid,$dpadir,$dpafile,$workdir,$flag_wmsms,$flag_wspdb);
    }
    else {
	($calc_info,$infocalc,$cutcc_work)=&_gDPA::generate_dpa(
	    $jobid,$pdbfile,$selechain,$cacrdf,$dpafile,
	    $msms_density,$msms_prob,$msmslayer,
            $flag_msmsdense,$flag_msms_saveall,$flag_wmsms,$flag_wspdb,$flag_uMSMS,$input_msmsf,
	    $cutcc,$delta_cutcc,$delta_cutlg,$cutlg,$cutll,$intcutcc,$wcc,$wctc,$wlg,$wll,$intwcc,
	    $nmodes,$modekappa,$ev_threshold,$flag_enforceDPA,
	    $flag_check,$flag_quite,
	    $structuredatadir,$msmsdir,$exedir,$workdir,$dpadir,$outputdatadir,
	    $proccf);
    }

    `rm -f $workdir/$cacrdf` if $flag_check eq 'false';    
    &processing_record($flag_quite,$outputdatadir,$proccf,$totproccf,$calc_info,"DPA__ERR> $jobid -- $infocalc") if($calc_info ne 0);
    &processing_record($flag_quite,$outputdatadir,$proccf,$totproccf,0,         "DPA_CALC> $jobid -- NCA_$nca NCHAIN_$nchain -- $dpafile -- CPU-TIME $infocalc") if($calc_info eq 0);    
    &processing_record($flag_quite,$outputdatadir,$proccf,$totproccf,0,         "DPA_CALC> $jobid -- cutcc changed $cutcc --> $cutcc_work")
	if($cutcc_work ne $cutcc && $infocalc ne ' using DPA');

    #updating processing information.....
    if ($flag_task eq 'genedpa' or $flag_task eq 'gdpa') {
	$flag_udpa='false';	
	&appendcopy($flag_quite,$outputdatadir,$proccf,$outputdatadir,$totproccf);
	`rm -f $outputdatadir/$proccf`;
	next;
    }	
#
#-- Generate DPA_CLUSTERS
#    $epislon_dpa,$MinPts_dpa,$epislon_cluster_dpa,$cutdpaperct,$cutdpaperct2,
#    $adjcutperct,$ndpavalue_section,$topdpa_num_cutoff,$dpabindingcutoff,
#    $flag_wcpdb,$file_output_dpa,
#
#    high-valued MSMS-points poured into clusters: $ndpaclu,
#
### increase bcut for higher layer DPA-cluster points
#   change bcut for different layers
#
    $dpabindingcutoff = $dpabindingcutoff + $msms_prob * ($msmslayer-1) if $flag_flexible_bcut eq 'true'; 
###
#    print ">>>>> $dpabindingcutoff <<<< $msmslayer ---- $msms_prob\n";
    my ($calc_info,$infocalc,$ndpaclu,$dpabindingsites,$dpacenters,$flag_sele_TOPDPA,
	$cutdpaperct,$cutdpaperct2,$extremefita,$extremefitb,$serr,$rcorr)=
	&_anaDPA::get_dpapredictbindsites($jobid,$dpafile,$cacrd,$dpadir,$exedir,$workdir,$outputdatadir,
					  $topdpa_num_cutoff,$dpabindingcutoff,$cutdpaperct,$cutdpaperct2,
					  $adjcutperct,$ndpavalue_section,$flag_wcpdb,
					  $epislon_dpa,$MinPts_dpa,$epislon_cluster_dpa,$flag_check);

    &processing_record($flag_quite,$outputdatadir,$proccf,$totproccf,$calc_info,"anaDPA> $jobid -- NCA_$nca NCHAIN_$nchain -- $infocalc ") if($calc_info ne 0);       
    &processing_record($flag_quite,$outputdatadir,$proccf,$totproccf,0,"anaDPA> $jobid -- NCA_$nca NCHAIN_$nchain ---- $ndpaclu DPA_cluster(s) found");
    &processing_record($flag_quite,$outputdatadir,$proccf,$totproccf,0,"anaDPA> $jobid -- EVD_FITTING cutdpaperct changed $cutdpaperct --> $cutdpaperct ") if($cutdpaperct ne $cutdpaperct &&  $flag_sele_TOPDPA eq'EVDFITTING');
    &processing_record($flag_quite,$outputdatadir,$proccf,$totproccf,0,"anaDPA> $jobid -- DIRECT_SELE cutdpaperct2 changed $cutdpaperct2 --> $cutdpaperct2 ") if($cutdpaperct ne $cutdpaperct && $flag_sele_TOPDPA eq 'DIRECT');
#############
    foreach my $rkid (sort {$a<=>$b} keys %$dpabindingsites){
        my $nn=$#{$$dpabindingsites{$rkid}} + 1;
#	print "LAYER_DPASITE> $jobid"."C$rkid : $nn --> @{$$dpabindingsites{$rkid}}\n";}
	print "DPASITES_$pdbid> L$msmslayer"."C$rkid : $nn --> @{$$dpabindingsites{$rkid}}\n" if $flag_check eq 'true';}
#     foreach my $rkid (sort {$a<=>$b} keys %{$dpacenters}){
#	foreach my $crdid (sort {$a<=>$b} keys %{$$dpacenters{$rkid}}){
#	    print "LAYER_CENTER> $jobid--$rkid--$crdid: @{$$dpacenters{$rkid}{$crdid}}\n";}}
#     exit();
#############
#    my $outputfile=$jobid.'.dpa';
#    if ($file_output_dpa){
#	$outputfile=$file_output_dpa.".dpa";
#	if ($file_output_dpa=~ /(.*)\./) {
#	    $outputfile=$1.".dpa";
#	}	
#    }
#    open(outf,"> $outputdatadir/$outputfile");
#    foreach my $rkid (sort {$a<=>$b} keys %$dpabindingsites){
#	foreach my $crdid (sort {$a<=>$b} @{$$dpabindingsites{$rkid}}){
#	    my $resName=@{$$cacrd{$crdid}}[3];
#	    my $chainID=@{$$cacrd{$crdid}}[4];  
#	    my $resSeq =@{$$cacrd{$crdid}}[5];
#	    my $iCode  =@{$$cacrd{$crdid}}[6];
#	    my $tFactor=@{$$cacrd{$crdid}}[7];
#	    printf "$wformat2\n",$jobid,$rkid,$resName,$chainID,$resSeq,$iCode if($flag_quite ne 'true');
#	    printf outf "$wformat\n",$rkid,$resName,$chainID,$resSeq,$iCode;
#	}}
#
#-----
#
    &appendcopy($flag_quite,$outputdatadir,$proccf,$outputdatadir,$totproccf);
    `rm -f $outputdatadir/$proccf`;
    return $calc_info,$infocalc,$ndpaclu,$dpabindingsites,$dpacenters,$nca,$cacrd,$dpafile;
}

sub find_chain{
    my $crd=@_[0]; my %chains=(); undef %chains;
    foreach my $i (sort {$a<=>$b} keys %$crd) {
	$chains{@{$$crd{$i}}[4]}++;
    }
    my $nchain= keys %chains;
    my @chainlist = keys %chains;
    return $nchain,\@chainlist;
}


sub appendcopy{
    my ($flag_quite,$dir1,$file1,$dir2,$file2)=@_;
#
#-- copy proccessing file to the total_proccessing:
    return if($flag_quite eq 'true');
    open(procf,">>$dir2/$file2")||die "CAN NOT OPEN WRITE $dir2/$file2" ;
    open(pp,"< $dir1/$file1")||die "CAN NOT OPEN READ $dir1/$file1" ;
    while(<pp>){
	print procf $_;
    }
    close(pp);
    close(procf);
}
#
#-- writing down running-data/informations in $procf, $totprocf
sub processing_record{
    my ($flag_quite,$outputdir,$procf,$totprocf,$info,$commend)=@_;

    my $timetitle = &timetitle();
    return if($flag_quite eq 'true');

    if($info <= -9999){   #serve error!
        print "DPA_ERR> $commend\n";}
     
    open(pp,">> $outputdir/$procf") || die "CAN NOT OPEN PRCESSING-FILE: $outputdir/$procf\n";

    if($info eq 'ini'){
	print pp "\nRecording DPA prediction at $timetitle\n";
	print pp "COMMAND: dpa.pl ";for(my $i=0;$i<=$#ARGV;$i++){print pp " $ARGV[$i]";}
	print pp "\n";}
    else{
	print pp "$commend -- recorded at $timetitle \n" if $commend ne '';}
    close(pp);
    if ($info < 0){
	&appendcopy($flag_quite,$outputdir,$procf,$outputdir,$totprocf);
	`rm -f $outputdir/$procf`;
    }
}

sub timetitle{
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime; #localtime(time);
    $year = $year - 100; # Handle century
    $mon++;
    $hour=$hour+8;

    if (length($year) < 2) { $year = "0" . $year;} # Make the years all two digits
    if (length($mon) < 2) { $mon = "0" . $mon;} # Make the months all two digits 
    if (length($mday) < 2) { $mday = "0" . $mday;} # Make the dates all two digits
    my $timetile= "$year-$mon-$mday $hour:$min";
    return $timetile;
}

1;
