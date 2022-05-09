#use lib "/xcommon/bin/_mdpa";
use lib "/Users/dming/mdpa/_mdpa";

use strict;

use _FILEprocc;
use _PDBprocc;
use _gMSMSpoints;
package _gDPA;

use Sys::Hostname;

sub generate_dpa{
    #
    #This subroutine calculate DPA-values for MSMS points around given PDB-structure
    #on return, it gives info: Success/Failed, info2: CPU time, cutcc0: the cutcc for ANM-NM calc.
    #           it also creates DPA-FILE: $dpafile in $dpadir, if info == 0
    #
    my ($jobid,$pdbfile,$selechain,$cacrdf,$dpafile,	
	$msms_density,$msms_prob,$msmslayer,
        $flag_msmsdense,$flag_msms_saveall,$flag_wmsms,$flag_wspdb,$flag_uMSMS,$input_msmsf,
	$cutcc,$delta_cutcc,$delta_cutlg,$cutlg,$cutll,$intcutcc,$wcc,$wctc,$wlg,$wll,$intwcc,
	$nmodes,$modekappa,$ev_threshold,$flag_enforceDPA,
	$flag_check,$flag_quite,
	$structuredatadir,$msmsdir,$exedir,$workdir,$dpadir,$outputdatadir,
	$proccf)=@_;
    if(-e "$dpadir/$dpafile"){ my $olddpaf='old_'.$dpafile; `mv -f $dpadir/$dpafile $dpadir/$olddpaf`};
    
#
#-- Get/Generate MSMS_points

    my ($msmssurff,$msmssurfpdbf,$nmsmspoints);
    if($flag_uMSMS eq 'true') {  #using existed MSMS-file
	$msmssurff=$input_msmsf;
	return -9995, "CAN NOT FIND INPUT MSMS-FILE: $workdir/$msmssurff"
	    if(! $msmssurff || ! -e "$workdir/$msmssurff"); 
	$nmsmspoints=&_FILEprocc::lineoffile($workdir,$msmssurff);  
    }    
    else {
	$msmssurff=$jobid.'.msms';  #if ($flag_wmsms eq 'true');
	$msmssurfpdbf=$jobid.'_msms'.'.pdb' if ($flag_wspdb eq 'true');
	
	($nmsmspoints)=&_gMSMSpoints::generate_MSMSpoints(
	    $jobid,$pdbfile,$selechain,$structuredatadir,$msmsdir,$workdir,
	    $msms_density,$msms_prob,$msmslayer,$flag_msmsdense,$flag_msms_saveall,
	    $msmssurff,$msmssurfpdbf,$flag_check,$flag_quite
	    );    # if($flag_wmsms or $flag_wspdb or $flag_check);    
	return -1,"NO MSMS-points generated with ",
	"\"-de $msms_density -prob $msms_prob -layer $msmslayer -dense $flag_msmsdense\" ",
	"try \"-msms 1. 1.5 ...\"" if($nmsmspoints <= 0); 
	return -9998,"NO MSMS_FILE CREATED in $structuredatadir" if(! -e "$workdir/$msmssurff");
    }
    return -1,"Too few MSMS_POINTS ($nmsmspoints) AVAILABLE (<10)" if ($nmsmspoints < 10);
#
#-- DPA calc.
    my $cutcc0=$cutcc;my $cutlg0=$cutlg;my $cutll0=$cutll;
    my $wcc0=$wcc;my $wctc0=$wctc;my $wlg0=$wlg; my $wll0=$wll;
    my $intcutcc0=$intcutcc; my $intwcc0=$intwcc;

    my $dpainputf=$jobid.'_dpa.inp';my $dpaoutputf=$jobid.'_dpa.out';
    my $info='null'; my $info2='null';
dpacycle:
    print "CALCDPA_jobid> $cutcc0,$cutlg0,$cutll0,$intcutcc0,$wcc0,$wctc0,$wlg0,$wll0,$intwcc0,$nmodes,$modekappa,$ev_threshold,$flag_enforceDPA\n"  if $flag_check eq 'true';
    open(dpain,">$workdir/$dpainputf");
    print dpain "\"$workdir/$cacrdf\" \"$workdir/$msmssurff\" \"$dpadir/$dpafile\"\n";
    print dpain "$wcc0 \t$wctc0 \t$wlg0 \t$wll0 \t$intwcc0\n";
    print dpain "$cutcc0 \t$cutlg0 \t$cutll0 \t$intcutcc0\n";
    print dpain "$nmodes \t$modekappa\n";
    print dpain "$ev_threshold\n";
    close(dpain);
    system "$exedir/dpa_ggsspnma.exe < $workdir/$dpainputf > $workdir/$dpaoutputf";

    ($info,$info2)=&read_calc_info($workdir,$dpaoutputf);  #info: Success/Failed, info2: CPU time
    if($flag_enforceDPA eq 'true'){  #loosen method works for even abnormal (noncompact) tructures
	if($info==-1&&$cutcc0<50){$cutcc0+=$delta_cutcc;$cutlg0=$cutcc0+$delta_cutcc; goto dpacycle;}
	if($info==-9&&$cutcc0<50){$cutcc0+=$delta_cutcc;$cutlg0=$cutcc0+$delta_cutcc; goto dpacycle;}}  
            #info=-9: singular KESSIAN, increase cutcc
    else{
	    #stringent method designed for ONLY normal (compact) structions and time-saving
	if(($info==-1||$info==-9)&&$cutcc0<16){
	    $cutcc0+=$delta_cutcc;$cutlg0=$cutcc0+$delta_cutlg;$intcutcc0=$cutcc0;	    
	    goto dpacycle; }	
	if(($info==-1||$info==-9)&&$cutcc0>=16){$info2=$info2.' try "-enforcedpa"';}
    }

    `rm -f $workdir/$dpainputf $workdir/$dpaoutputf $workdir/$cacrdf` if $flag_check eq 'false';
    `rm -f $workdir/$msmssurff` if($flag_uMSMS ne 'true' && $flag_wmsms ne 'true' && $flag_check eq 'false'); 

    return $info,$info2,$cutcc0;
}


sub dpainfo_procc{
    my ($proccf,$jobid,$info,$info2,$outputdir)=@_;
    my $timetitle = &timetitle();
    open(pp,">> $outputdir/$proccf");
    if($info==0) {
	printf pp "%4s\t%5d\t%9.4f%31s %14s\n",$jobid,$info2," (CPU seconds) -- Recording at ",$timetitle;}
    if($info==-9999){
	print pp "$jobid OPEN cacrd error   -- Recording at $timetitle\n";}
    elsif($info==-9998){
	print pp "$jobid READ cacrd error   -- Recording at $timetitle\n";}	
    elsif($info==-9997){
	print pp "$jobid OPEN surfcrd error   -- Recording at $timetitle\n";}
    elsif($info==-9996){
	print pp "$jobid READ surfcrd error   -- Recording at $timetitle\n";}
    elsif($info==-9995){
	print pp "$jobid Too few MSMS points  -- Recording at $timetitle\n";}
    elsif($info==-5000){
	print pp "$jobid protein CRD pair-distance < 0.3   -- Recording at $timetitle\n";}
    elsif($info==-5100){
	print pp "$jobid protein/ligand (surf point) CRD pair-distance < 1   -- Recording at $timetitle\n";}
    elsif($info==-2000){
	print pp "$jobid Dimension of RP (neighboring contacts) is too smal   -- Recording at $timetitle\n";}
    elsif($info==-1000){
	print pp "$jobid Dimension of HESSIAN is too smal   -- Recording at $timetitle\n";}
    elsif($info==-9){
	print pp "$jobid Kessian is singular (NO \"G^t K^-1 G\")    -- Recording at $timetitle\n";}
    elsif($info==-1){
	print pp "$jobid CUTCC for APOprotein is too small    -- Recording at $timetitle\n";}
    elsif($info eq "null"){
	print pp "$jobid\tDPA calc. segment false  -- Recording at $timetitle\n";}
    close(pp);
}

sub read_calc_info{
    my $dir=@_[0]; my $file=@_[1];
    my $info='null'; my $info2='null';

    open(tmpf,"<  $dir/$file");
    while(<tmpf>){
	my @tmp=split /\s+/,$_;
	$info=@tmp[1];  $info=@tmp[0] if @tmp[0] ne '';
	$info2=$_; chomp($info2); $info2 =~ s/$info//gi; $info2 =~ s/^\s+|\s+$//g;
	$info2='"'.$info2.'"';
	last;
    }
    close(tmpf);
    return $info,$info2;
}

sub timetitle{
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime; #localtime(time);
    $year = $year - 100; # Handle century
    $mon++;
    $hour=$hour+8;
    if (length($year) < 2) { $year = "0" . $year;} # Make the years all two digits
    if (length($mon) < 2)  { $mon  = "0" . $mon;} # Make the months all two digits 
    if (length($mday) < 2) { $mday = "0" . $mday;} # Make the dates all two digits
    my $timetile= "$year-$mon-$mday $hour:$min";
    return $timetile;
}
sub writemsms_w_dpa{
    my ($jobid,$dpadir,$dpafile,$workdir,$flag_wmsms,$flag_wspdb)=@_;

    my $surff=$jobid.".msms"; 
    my $surfpdbf=$jobid."_msms.pdb";
    
    my %surfp=(); undef %surfp; 
    my ($x,$y,$z);
    open (sdpaf,"< $dpadir/$dpafile"); 
    my $n=0;
    while (<sdpaf>){
	my @tmp=split /\s+/,$_;
	if(@tmp[0] eq ''){ $x=@tmp[1];$y=@tmp[2];$z=@tmp[3];}
	if(@tmp[0] ne ''){ $x=@tmp[0];$y=@tmp[1];$z=@tmp[2];}
	push(@{$surfp{++$n}},$x,$y,$z);
    }close(sdpaf);

    if($flag_wmsms eq "true"){	
	system "rm -f $workdir/$surff " if -e "$workdir/$surff";
	open(WCRD,"> $workdir/$surff") || die "CAN NOT OPEN FOR WRITE MSMS: $workdir/$surff\n";
	foreach my $i (sort {$a<=>$b} keys %surfp) {
	    $x=@{$surfp{$i}}[0];$y=@{$surfp{$i}}[1];$z=@{$surfp{$i}}[2];
	    printf WCRD " %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f %6d %6d\n",$x,$y,$z,
	    1.0,1.0,1.0,0,$i;
	}
	close(WCRD);
    }

    if($flag_wspdb eq "true"){	
	system "rm -f $workdir/$surfpdbf " if -e "$workdir/$surfpdbf";
	open(WPDB,"> $workdir/$surfpdbf") || die "CAN NOT OPEN FOR WRITE MSMS: $workdir/$surfpdbf\n";

	foreach my $i (sort {$a<=>$b} keys %surfp) {
	    $x=@{$surfp{$i}}[0];$y=@{$surfp{$i}}[1];$z=@{$surfp{$i}}[2];
	    my $ii=substr($i,0,4);
	    printf WPDB "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s\n", 
	    'ATOM  ',$i,' SUR',' ','SUR','A',$ii,' ',$x,$y,$z,1.0,10.0,'','SURF','','';
	}
	close(WPDB);
    }
    
    return;
}

1;
