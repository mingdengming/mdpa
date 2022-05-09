package _CLUSTER;

#
############## Clusterization of given set of points using OPTICS methods #######################
sub OPTICS{
    ($crddata,$epislon,$MinPts,$epislon_cluster,$flag_wpdb,$output_dir)=@_;               #Global PARAMETERS in Cluster
    $nbMinPts=$MinPts-1;
    $epislons=$epislon**2;$epislon_clusters=$epislon_cluster**2;
    $for_core_distanceSq='core_distanceSq';
    $if_processed='if_processed';
    $for_reach_distSq='reachable_distanceSq';
    %oset=();%infoset=();%seedlist=();$nlist=0;%listseed=();@orderlist=();      #Global LIST of database 
    undef %oset;undef %infoset;undef %seedlist; undef %listseed;undef @orderlist;
    my $nop=-1; foreach $i ( keys %$crddata){push (@{$oset{$i}},@{$$crddata{$i}}); $nop++;}
#########
#    foreach (sort {$a<=>$b} keys %oset){ print "CHECK_OSET>$_:    @{$oset{$_}}\n"; }
#########

    for(my $i=0;$i<=$nop;$i++){
	next if ($infoset{$i}{$if_processed} eq 'yes');
	&ExpandClusterOrder($i);
    }
    my $optics_cluster=&optics_cluster();

    &showCluster($flag_wpdb,$optics_cluster,$output_dir) if($flag_wpdb ne '');

#########
#    foreach (sort {$a<=>$b}keys %$optics_cluster){  print "OPTICS> $_  @{$$optics_cluster{$_}}\n";};
#    print "ORDERLIST> @orderlist\n"; # exit();
#########
    my $num_effect_cluster=0;
    foreach my $cluid (keys %$optics_cluster){ 
	$num_effect_cluster++ if $cluid > 0;  #$cluid=0 is for noise, not an eff. cluster
    }

    return $num_effect_cluster,$optics_cluster;
}
sub optics_cluster{
    my %cluster=();undef %cluster;my $clusterid=0; #0 for NOISE, 1,2,3... for nontrivial cluters
    
    foreach my $i (@orderlist){
	$rcsq=$infoset{$i}{$for_core_distanceSq};
	$rdsq=$infoset{$i}{$for_reach_distSq} ;
	if($rdsq > $epislon_clusters or $rdsq eq 'undef'){
	    if($rcsq<=$epislon_clusters and $rcsq ne 'undef' ){
		$clusterid++;
		push(@{$cluster{$clusterid}},$i);}
	    else{
		push(@{$cluster{0}},$i);}
	}
	else{
	    push(@{$cluster{$clusterid}},$i);}
    }
    return \%cluster;
}
sub ExpandClusterOrder{
    my $id=@_[0];
    my ($nb,$neighbor)=Neighbor($id);
    $infoset{$id}{$if_processed}='yes';
    $infoset{$id}{$for_core_distanceSq}=CoreDistance($neighbor);
    push (@orderlist,$id); #print "ADDING ($id) to orderlist\n";
    $infoset{$id}{$for_reach_distSq}='undef';
    return  if($infoset{$id}{$for_core_distanceSq} eq 'undef');   #NOT A core-object
    &OrderSeedUpdata($neighbor,$id);
    return if($nlist <=0);
    while ($nlist >=1) {
	($first_rds,$first_listid)=NextSeedlist(); 
	splice(@{$seedlist{$first_rds}},0,1);$nlist--;
	delete($listseed{$first_listid});
	if(scalar(@{$seedlist{$first_rds}}) eq 0){
	    delete($seedlist{$first_rds});}
	($nb,$neighbor)=Neighbor($first_listid);
	$infoset{$first_listid}{$if_processed}='yes';
	$infoset{$first_listid}{$for_core_distanceSq}=CoreDistance($neighbor);
	push (@orderlist,$first_listid); 
	#print "whileloop)ADDING ($first_listid) with CORE $infoset{$first_listid}{$for_core_distanceSq} to orderlist\n" if($id ne '');

	next if($infoset{$first_listid}{$for_core_distanceSq} eq 'undef');
	&OrderSeedUpdata($neighbor,$first_listid);
    }
}
sub CoreDistance{
    my $neighbor=@_[0];
    my $nnb=0; my $rcores='undef';
    foreach my $rs (sort {$a<=>$b} keys %$neighbor){
	$nnb+=scalar(@{$$neighbor{$rs}});
	if($nnb >= $nbMinPts){      #counting the object itself
	    $rcores=$rs;
	    last};
    }
    return $rcores;
}
sub Neighbor{     #need: %oset,$epislons
    my $id=@_[0];
    my %neighbor=(); undef %neighbor;
    my $x0,$y0,$z0,$x,$y,$z,$dx,$dy,$dz;my $nb=0;
    $x0=@{$oset{$id}}[0];$y0=@{$oset{$id}}[1];$z0=@{$oset{$id}}[2];
    foreach my $i (keys %oset){
	next if($i eq $id);
	$x=@{$oset{$i}}[0];$y=@{$oset{$i}}[1];$z=@{$oset{$i}}[2];
	$dx=$x0-$x;$dy=$y0-$y;$dz=$z0-$z;
	next if($dx>$epislon or $dy>$epislon or $dz>$epislon or $dx<-$epislon or $dy<-$epislon or $dz<-$epislon);
	$ds=$dx**2+$dy**2+$dz**2;
	next if($ds > $epislons);
	push (@{$neighbor{$ds}},$i);
	$nb++;
    }
    return $nb,\%neighbor;
}
sub OrderSeedUpdata{
    my $neighbor=@_[0];my $id=@_[1];
    my $reach_distSq,$old_rds;
    my $core_distSq=$infoset{$id}{$for_core_distanceSq};  #$core_distSq must PREDEFINED

    foreach my $ds (sort {$a<=>$b} keys %$neighbor){
	foreach my $nbid (@{$$neighbor{$ds}}){
	    $reach_distSq=$core_distSq;$reach_distSq=$ds if($ds > $core_distSq);
	    if ($infoset{$nbid}{$if_processed} eq 'yes'){
		next if ($infoset{$nbid}{$for_reach_distSq} ne 'undef' or $infoset{$nbid}{$for_core_distanceSq} ne 'undef');
		for(my $tid=0;$tid<=$#orderlist;$tid++){
		    if(@orderlist[$tid] eq $nbid){
			splice(@orderlist,$tid,1);last;}}
	    }
	    if($listseed{$nbid} eq ''){
		$listseed{$nbid}=$reach_distSq;
		$infoset{$nbid}{$for_reach_distSq}=$reach_distSq;    #renew reachable_distance
		push(@{$seedlist{$reach_distSq}},$nbid);
		$nlist++;
	    }
	    elsif($reach_distSq <$listseed{$nbid}){
		$old_rds=$listseed{$nbid};
		$listseed{$nbid}=$reach_distSq;
		$infoset{$nbid}{$for_reach_distSq}=$reach_distSq;
		my $tn=-1;
		for(my $tid=0;$tid<=$#{$seedlist{$old_rds}};$tid++){
		    if(@{$seedlist{$old_rds}}[$tid] eq $nbid){
			splice(@{$seedlist{$old_rds}},$tid,1);$nlist--;last;}}
		if(scalar(@{$seedlist{$old_rds}}) eq 0){
		    delete($seedlist{$old_rds});}
		push(@{$seedlist{$reach_distSq}},$nbid);$nlist++; #adding $nbid with new Rds
	    }
	}
    }
    return;
}
sub NextSeedlist{
    my @rds=(sort {$a<=>$b} keys %seedlist);
    my $first_rds=@rds[0];
    my $first_listid=@{$seedlist{@rds[0]}}[0];
    return $first_rds,$first_listid;
}
sub rank_of_clusters{
    my $nclu=@_[0];                   # number of clueter
    my $crdclu=@_[1];                 # ID (point) list for each cluster
    my $crd=@_[2];                    # crd for each ID (point)
    my %clurank=(); undef %clurank;
   
    push (@{$clurank{-1}},-1,-1,-1);
    if ($nclu<=0){
	push(@{$clurank{0}},-1,-1);
	return \%clurank;
    }
    if ($nclu==1){
	push(@{$clurank{1}},1,1);	
	return \%clurank;
    }
    my $cluminino=9999;               #the minimum of MSMSpoints numbers for all Clusters
    foreach my $cid  (keys %$crdclu){
	next if $cid eq 0;            #this is for noise-point-set, ignored
	my $tmpn=scalar(@{$$crdclu{$cid}});
	$cluminino=$tmpn if $tmpn < $cluminino;
    }
    my @pertavgall=();undef @pertavgall;
    foreach my $cid (sort {$a<=>$b} keys %$crdclu){    #DPADATA ordered as "pert" values
	next if $cid <=0;
	my $perttt=0.; my $pertavg=-1; my %tmpcrd=(); undef %tmpcrd;
	my $ncrd=0;                                   #take first $cluminino data
	foreach my $crdid (@{$$crdclu{$cid}}){
	    $perttt+=@{$$crd{$crdid}}[-1]; $ncrd++;last if $ncrd > $cluminino; #take first $cluminino data
	    push (@{$tmpcrd{$crdid}},@{$$crd{$crdid}});

	}
	$pertavg=$perttt/$ncrd;
	my $gr=&gyrationradio(%tmpcrd);
	push (@{$clurank{$cid}},$pertavg,$gr); push(@pertavgall,$pertavg);
    }
    my $rankid=0; my %flag_rank=();undef %flag_rank;
    foreach (sort{$b<=>$a} @pertavgall){
	$rankid++;
	foreach $cid (keys %clurank){
	    if($flag_rank{$cid} ne 1 and abs(@{$clurank{$cid}}[0]-$_) le 0.0001){
		push (@{$clurank{$cid}},$rankid);
		$flag_rank{$cid}=1;
	    }
	}
    }
########
#    foreach my $cid (sort {$a<=>$b} keys %clurank) {
#	print "ClusterRANK> $cid -- @{$clurank{$cid}}\n";}
########
    return \%clurank;   #ranking the input clusters "$crdclu"
    }

sub gyrationradio{
    my %crd=@_;
    my $gr=-1; my @center=(); my $ntt=0;
    foreach my $i (sort {$a<=>$b}keys %crd){
	$ntt++;
	for (my $k=0; $k<3; $k++){ @center[$k]+=@{$crd{$i}}[$k];}
	#print "GYARADIO> $i --  @{$crd{$i}}[0] @{$crd{$i}}[1] @{$crd{$i}}[2] \n";
    }
    for (my $k=0; $k<3; $k++){ @center[$k]=@center[$k]/$ntt}
    
    my $tmpv=0;
    foreach my $i (keys %crd){
	for (my $k=0; $k<3; $k++){ $tmpv+=(@{$crd{$i}}[$k]-@center[$k])**2;}}

    $tmpv=$tmpv/$ntt;
    $gr=sqrt($tmpv); $density=9999;
    $density=(3/4./3.14159)*$ntt/($gr**3) if $gr > 0;
    
    #print "GYARADIO> $gr, $density  -- $ntt\n";

    return $gr;
}
sub showCluster{          
    my $wpdbfilename=@_[0];
    my $cluster=@_[1];   
    my $outputdir=@_[2];

    my @chainid_list=("O".."Z");my @resName_list=('CLU','DLU','XLU','ZLU');
    open(wp,"> $outputdir/$wpdbfilename");
    my $iatom=0;
    foreach $cid (sort {$a<=>$b} keys %$cluster) {
	$chainid=X if($cid eq 0);
	$chainid=@chainid_list[($cid-1)%scalar(@chainid_list)] if($cid >=1);
	$fold_id=int(($cid-1)/scalar(@chainid_list));
	$resname=@resName_list[$fold_id];
	$resname=@resName_list[$fold_id%scalar(@resName_list)] if($fold_id >= scalar(@resName_list));;
	$resname='NOS' if($cid eq 0);
	my $ia_chain=0;
	foreach (my $i=0;$i<=$#{$$cluster{$cid}};$i++){
	    my $id=@{$$cluster{$cid}}[$i];$pid=$id;
	    my $x=@{$oset{$id}}[0]; my $y=@{$oset{$id}}[1]; my $z=@{$oset{$id}}[2];
	    my $resSeq=@{$oset{$id}}[3];
	    my $resName=@{$oset{$id}}[4];
	    my $value=10.00;
	    my $value=@{$oset{$id}}[5] if @{$oset{$id}}[5] ne '';
	    $iatom++;
	    $ia_chain++;

#	    printf wp "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s\n", 'HETATM',
#	    $iatom,' CA ',' ',$resname,$chainid,$ia_chain,' ',$x,$y,$z,1.0,10.0,'','MING','_C','LA';

	    printf wp "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%4d\n", 'HETATM',
	    $iatom,' CA ',' ',$resName,$chainid,$resSeq,' ',$x,$y,$z,1.0,10.0,'','MING',$pid;

#	    printf wp "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%4d\n", 'HETATM',
#	    $iatom,' CA ',' ',$resName,$chainid,$ia_chain,' ',$x,$y,$z,1.0,$value,'','MING',$pid;

	    }
    }
    close(wp);
}

###
############### END OF CLUSTER PROGRAM #######################

1;
