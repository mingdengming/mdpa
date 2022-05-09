use strict;
#use lib "/xcommon/bin/_mdpa";
#use lib "/Users/dming/mdpa/_mdpa";
use lib "/opt/xcommon/bin/_mdpa";
use     _CLUSTER;
package _anaDPA;

############### ANALYSIS OF DPA ######################################################################
#
#   ADD HIEARCHICAL CALC., I.E., GENERTING L_DPA-CLUSTERS FOR EACH (L)LEVEL,
#       AND FIND CONNECTION BETWEEN 1-, 2-,...,L-DPA-CLUSTERS
#       THE PREDICTION IS BASED ON THE CONNECTED-CLUSTERS 
#     in Nanjing Tech University, 2021/05/21
#
#   NOISE-DPA-MSMS-POINT Set now included by comment out "next if($targetclusterster_id le 0); "
#   at "sub getcontacts_with_ClusterPOINTS"
#   NOISE-DPA-MSMS-POINT Set ranked BY ADDING "$rank = 0 if $cid == 0;"
#   at "sub get_dpapredictbindsites"
#
######################################################################################################


sub get_layercluster_connection{
#
# the following calculates layer-cluster connetivity 
# 
#  Giving Layer1, Cluster O1, P1, Q1....    recorded in $mdpa_centers, prediction in $mdpa_bindingsites
#         Layer2, Cluster O2, P2, Q2....    recorded in $mdpa_centers, prediction in $mdpa_bindingsites
#         Layer3........................
#  Gernate Prediction Cluster C1, C2, C3....
#         
#        IF maxmsmslayer == 1: Then C1 = Q1, C2 = P1, C3 = Q1, ......, C_n = Z1
#
#    (1) IF maxmsmslayer >= 2: Then
#           consider connectivity betweey Layer 1 and 2:
#              first, set C1 = Q1, C2 = P1, C3 = Q1, ......
#       (2)    if O2 conn C_i (i=1, n), then merge O2 to C_i, 
#                 else C_(n+1)=O2, n = n+1. 
#           repeat (2) for P2, Q2.....
#           now, get new C1, C2, ..., Cn
#       repeat (1) for maxmsmslayer >= 3 (consider Layer1,2 and 3)
#       repeat (1) for maxmsmslayer >= 4 (consider Layer1,2,3 and 4)
#       ............................................................    
#       repeat (1) until maxmsmslayer == maxmsmslayer.
#
#     return n,C1,C2,....,Cn,prd1,prd2,prd3,.....
#       
#
    my ($pdbid,$maxmsmslayer,
        $connpoints,$conncutoff,$epislon_dpa,$MinPts_dpa,$epislon_cluster_dpa,
	$mdpa_nclu,$mdpa_bindingsites,$mdpa_centers,$nca,$cacrd,$flag_check)=@_;
    print "MDPA_$pdbid"."_ANA> HIEARCHICAL LAYER-CONECTIVITY ANALYSIS\n";
#############
#    foreach my $layer (sort {$a<=>$b} keys %$mdpa_bindingsites){
#        my $nclu = keys %$mdpa_bindingsites;
#        foreach my $rkid (sort {$a<=>$b} keys %{$$mdpa_bindingsites{$layer}}){
#            my $nsite = keys @{$$mdpa_bindingsites{$layer}{$rkid}};
#            print "DPA0> $pdbid -- $layer -- $nclu -- $rkid --> $nsite :: @{$$mdpa_bindingsites{$layer}{$rkid}}\n";
#        }
#    }
#       	my $dpabindingsites = $$mdpa_bindingsites{$layer};
#        my $layer =1 ;
#	my $dpacenters = $$mdpa_centers{$layer};
#                print "----> $layer ----> $$mdpa_nclu{$layer}\n";
#    	foreach my $rkid (sort {$a<=>$b} keys %$dpacenters){
#            last if $rkid >=2;
#            foreach my $crdid (sort {$a<=>$b} keys %{$$dpacenters{$rkid}}){
#                print "DPA_MSMSCTRS> $pdbid--$layer--$rkid--$crdid: @{$$dpacenters{$rkid}{$crdid}}\n";}}
#    	foreach my $rkid (sort {$a<=>$b} keys %{$$mdpa_centers{$layer}}){
#            last if $rkid >=2;
#            foreach my $crdid (sort {$a<=>$b} keys %{$$mdpa_centers{$layer}{$rkid}}){
#                print "DPA_CTRS> $pdbid--$layer--$rkid--$crdid @{$$mdpa_centers{$layer}{$rkid}{$crdid}}\n";}}
#        foreach my $rkid (sort {$a<=>$b} keys %$dpabindingsites){
#            print "GDPACTS0000> $pdbid---$layer- $rkid : @{$$dpabindingsites{$rkid}}\n";
#            foreach my $crdid (sort {$a<=>$b} @{$$dpabindingsites{$rkid}}){
#                my $resName=@{$$cacrd{$crdid}}[3]; my $chainID=@{$$cacrd{$crdid}}[4];
#                my $resSeq =@{$$cacrd{$crdid}}[5]; my $iCode  =@{$$cacrd{$crdid}}[6];
#                printf "%4s %3d %3d %3s %1s %5d %1s\n",$pdbid,$layer,$rkid,$resName,$chainID,$resSeq,$iCode;
#        }}
#    }
###############
#
    my %clusters; my %predictions; undef %clusters,%predictions; #record the converged clusters/predictions;
    my %conn; undef %conn;   #recording the cluster-connectivity of each converged cluster,
                             #    e.g.  L1-L2: a cluster at layer 1 extended to one at layer 2
    my $nclu=0;     #index for Whole-DPA-clusters, 0 for NOISE-DPA-MSMS-point set
    my %noise_centers=''; my %noise_psites=''; undef %noise_centers, %noise_psites;
    my $layer=1; 
    foreach my $rkid (sort {$a<=>$b} keys %{$$mdpa_bindingsites{$layer}}){
    #
        if($rkid <=0){
            my $nzp=0;
            foreach my $crdid (sort {$a<=>$b} keys %{$$mdpa_centers{$layer}{$rkid}}){
                push(@{$noise_centers{0}{$nzp++}}, @{$$mdpa_centers{$layer}{$rkid}{$crdid}});}
        }
        next if $rkid <=0;   #ignore noise-DPA-MSMS-point set
    #    rkid == 0 : for noise-DPA-MSMS-point set (centers), 
    #                noise-point set will be collected at ALL layers and
    #                finally created ONE total-noise-set, with rankid = 0!
    #
    #                this set is NOT merged to normal DPA-cluster/prediction
    #                after normal-DPA-cluster merging happen, 
    #                noise-set will be evalued if merged to certain WHOLE-DPA-clusters
    #   
        $nclu += 1;      
        $conn{$nclu}='L1C'.$rkid; my $np=0;
        push (@{$predictions{$nclu}},@{$$mdpa_bindingsites{$layer}{$rkid}});
        foreach my $crdid (sort {$a<=>$b} keys %{$$mdpa_centers{$layer}{$rkid}}){
            push(@{$clusters{$nclu}{$np++}}, @{$$mdpa_centers{$layer}{$rkid}{$crdid}});}
        print "CLU_CONN_$pdbid> FIRST_LAYER_CLU_ADDED $nclu: $conn{$nclu}\n" if($flag_check eq 'true');
    }

    while ($layer < $maxmsmslayer){
        $layer += 1;
        my $tmpclusters=$$mdpa_centers{$layer};           #for next layer
        my $tmppredictsites=$$mdpa_bindingsites{$layer};      #for next layer

#        foreach my $rkid (sort {$a<=>$b} keys %$tmpclusters){
#            foreach my $crdid (sort {$a<=>$b} keys %{$$tmpclusters{$rkid}}){
#                print "CHECKING_CTRS> $pdbid --> ($layer :: $rkid) --> $crdid: @{$$tmpclusters{$rkid}{$crdid}}\n";
#            }}
#        foreach my $rkid (sort {$a<=>$b} keys %$tmppredictsites){
#            my $nn =   $#{$$tmppredictsites{$rkid}};
#            print "CHECKING_SITES> $pdbid -- ($layer :: $rkid) --> $nn -->   @{$$tmppredictsites{$rkid}}[0..10]\n";}
#        next;

        foreach my $r1 (keys %$tmpclusters){              #check each cluster "r1" in new $layer
            if($r1 <=0){
                my $nzp=0;
                foreach my $crdid (sort {$a<=>$b} keys %{$$mdpa_centers{$layer}{$r1}}){
                    push(@{$noise_centers{0}{$nzp++}}, @{$$mdpa_centers{$layer}{$r1}{$crdid}});}
            }
            next if $r1 <=0;   #ignore noise-DPA-MSMS-point set
            my $nn = $#{$$tmppredictsites{$r1}} +1;
            print "CLU_CONN_$pdbid> L$layer"."C$r1 IGNORED, NO-BINDING-SITE FOUND\n" if $nn <=0 and $flag_check eq 'true';
            next if $nn <=0;   #this DPA-cluster canNOT find contact residues

            my $tmpclu = $$tmpclusters{$r1};
            my $flag_merged='false';
            foreach my $rkid ( keys %clusters){            #cluster "rkid" in WHOLE-DPA-CLUSTERS
                my $aclu = $clusters{$rkid};
                my ($isconn)=&isconn_clusterAB($tmpclu,$aclu,
                    $connpoints,$conncutoff,$epislon_dpa,$MinPts_dpa,$epislon_cluster_dpa);
                if ($isconn eq 'true')  {
                    #-- NOW cluster-r1 in current-layer merged to cluster-rkid in whole-dpa-clusters 
                    $conn{$rkid}=$conn{$rkid}."_L".$layer."C".$r1;
                    print "CLU_CONN_$pdbid> L$layer"."C$r1 -MERGEDTO $rkid: $conn{$rkid}, CONN-METHOD($connpoints,$conncutoff)\n" if $flag_check eq 'true';
                    my $mergedsites=&sum_array_clear($$tmppredictsites{$r1},$predictions{$rkid});
                    $predictions{$rkid}=$mergedsites;
#                    push(@{$predictions{$nclu}},@{$mergedsites});
                    my $np=keys %{$clusters{$rkid}}; #print "....> $np\n";
#   		    push(@{$predictions{$rkid}},@{$$tmppredictsites{$r1}});
                    foreach my $i (sort {$a<=>$b}  keys %$tmpclu){ 
			push(@{$clusters{$rkid}{$np++}}, @{$$tmpclu{$i}});
                     }
                     $flag_merged='true';
                     last;}
            } #finish checking cluster-ri in this-layer with all the clusters in Whole-DPA-clusters

            if($flag_merged eq 'true'){     #this tmpclu --merged--> to a whole-dpa-cluster,already merged
	        next;}                      # check next cluster in this new layer
            else {                          #
                $nclu++; my $np=0;          #tmpclu as a new whole-dpa-cluste, add new index
                $conn{$nclu}='L'.$layer."C".$r1;  
                push (@{$predictions{$nclu}},@{$$tmppredictsites{$r1}});
                foreach my $i (sort {$a<=>$b}  keys %$tmpclu){ 
                    push (@{$clusters{$nclu}{$np++}},@{$$tmpclu{$i}});}
                print "CLU_CONN_$pdbid> L$layer"."C$r1 -NOT-CONN-THUS-ADDED $nclu: $conn{$nclu}, CONN-METHOD ($connpoints,$conncutoff)\n" if $flag_check eq 'true';
#                foreach my $i (keys %conn){
#		    print "CONN> $i $conn{$i}\n";}

            }
        }
   }

#   foreach my $rkid (sort {$a<=>$b} keys %predictions){
#        print "DDD_CONN> $pdbid --> $rkid ---> $conn{$rkid}\n";
#        print "DDD_SITE> $pdbid --> $rkid ---> @{$predictions{$rkid}}\n";
#   }
#   print "CENTERS>>>>>>>>>>>>>>> \n";
#   foreach my $rkid (sort {$a<=>$b} keys %clusters){
#        foreach my $crdid (sort {$a<=>$b} keys %{$clusters{$rkid}}){
#            print "DDD_CTRS> $pdbid --> $rkid -- $crdid: @{$clusters{$rkid}{$crdid}}\n"; 
#        }
#   }

#    foreach my $rkid (sort {$a<=>$b} keys %$apredicts){
#        print "DPA_PRED> $pdbid $layer $rkid --- @{$$apredicts{$rkid}}\n";
#    }
#    foreach my $rkid (sort {$a<=>$b} keys %$aclusters){
#        foreach my $crdid (sort {$a<=>$b} keys %{$$aclusters{$rkid}}){
#                print "DPA_MSMSCTRS> $pdbid--$layer--$rkid--$crdid: @{$$aclusters{$rkid}{$crdid}}\n";}}

    print "MDPA_$pdbid"."_ANA> Done!\n";
    return $nclu,\%predictions,\%clusters,\%conn,\%noise_centers;
}

#$crddata,$epislon,$MinPts,$epislon_cluster,$flag_wpdb,$output_dir

sub sum_array_clear{
    my ($thisarray,$target)=@_;
    my @sum=''; undef @sum;
    foreach my $e (@{$target}) {push (@sum,$e)};
    foreach my $e (@{$thisarray}) {
        my $flag_exist='false';
        foreach my $e0 (@{$target}) {
            if ($e eq $e0){$flag_exist='true'; last}
            }
        if ($flag_exist ne 'true'){
              push (@sum,$e);}
    }
    return \@sum;
}
sub isconn_clusterAB{
    my ($clu1,$clu2,$connp,$conncut,$epislon1,$MinPts1,$epislonclu1)=@_;
    my $mergcrd; my $n = 0; 
    foreach my $i (sort {$a<=>$b}  keys %$clu1){ push (@{$$mergcrd{$n++}},@{$$clu1{$i}});}
    foreach my $i (sort {$a<=>$b}  keys %$clu2){ push (@{$$mergcrd{$n++}},@{$$clu2{$i}});}
#=head
#    print "\n\n";
#    foreach my $i (sort {$a<=>$b} keys %$clu1){
#	print "-1111--> $i --- @{$$clu1{$i}}\n";}
#    print "<><><><><><><><><><\n";
#    foreach my $i (sort {$a<=>$b} keys %$clu2){
#	print "-2222--> $i --- @{$$clu2{$i}}\n";}
#    print "\n";
#     print "<><><><><><><><><><\n";
#     foreach my $i (sort {$a<=>$b} keys %$mergcrd){
#         print "-9999--> $i --- @{$$mergcrd{$i}}\n";}
#     print "\n\n";
#    print "CONN>>>> $connp ---- $conncut \n\n"; exit;
#=cut

    if($connp <= -1) {  #using OPTICS to determine if clu1 and clu2 are connected
        my ($nclu,$topclu)=&_CLUSTER::OPTICS($mergcrd,$epislon1,$MinPts1,$epislonclu1); #,"_tmp_319.dpa1","./");
        if ($nclu == 1) {
            return 'true',$mergcrd;}
        elsif ($nclu >= 2) {
	    return 'false';}
        else {
            print "CLUANA> ERROR IN DETERMINING CLU1/2 CONN: NCLU = $nclu\n"; exit;}
    }
    else{
        my $nconn=0;
        foreach my $i (sort {$a<=>$b} keys %$clu1){
            my @p1; undef @p1; @p1 = @{$$clu1{$i}};
            foreach my $j (sort {$a<=>$b} keys %$clu2){
                my @p2; undef @p2; @p2 = @{$$clu2{$j}}; #print "COMP>> @p1[0..2] ==== @p2[0..2] \n"; 
                my $d12= &distp2p($$clu1{$i},$$clu2{$j}); 
                $nconn ++ if $d12 <= $conncut;
                last if $nconn >= $connp;
            }
            last if $nconn >= $connp;
        }       
        if ($nconn >= $connp){
	    return 'true';}
        else {
	    return 'false';}
    }
}

sub distp2p {
   my ($p1,$p2)=@_;
   my $dist = (@{$p1}[0]-@{$p2}[0])**2 + 
              (@{$p1}[1]-@{$p2}[1])**2 +
              (@{$p1}[2]-@{$p2}[2])**2;

   $dist = sqrt($dist);
#   print "$dist >> @{$p1}[0] @{$p1}[1] @{$p1}[2]     @{$p2}[0] @{$p2}[1] @{$p2}[2]\n"; exit;

   return $dist;
}

#
#--- KEY subroutine
sub get_dpapredictbindsites{ 
    my ($jobid,$dpafile,$cacrd,$dpa_dir,$exe_dir,$work_dir,$output_dir,
	$topdpa_num_cutoff,$dpabindingcutoff,$cutdpaperct0,$cutdpaperct2,
	$adjcutperct,$ndpavalue_section0,$flag_wcpdb,
	$epislon_dpa,$MinPts_dpa,$epislon_cluster_dpa,$flag_check)=@_;
    my $dpaclu=0; 
#
    print "ANADPA_PARA> $jobid,$dpafile,$cacrd,$dpa_dir,$exe_dir,$work_dir,$output_dir,	$topdpa_num_cutoff,$dpabindingcutoff,$cutdpaperct0,$cutdpaperct2,$adjcutperct,$ndpavalue_section0,$flag_wcpdb,$epislon_dpa,$MinPts_dpa,$epislon_cluster_dpa,$flag_check\n" if $flag_check eq 'true';
#
#-- List DPA values
   my %surfdpa=(); undef %surfdpa; 
    my ($pert,$x,$y,$z);
    open (SDPAF,"< $dpa_dir/$dpafile"); 
    while (<SDPAF>){
	my @tmp=split /\s+/,$_;
	if(@tmp[0] eq ''){ $pert=@tmp[4]; $x=@tmp[1];$y=@tmp[2];$z=@tmp[3];}
	if(@tmp[0] ne ''){ $pert=@tmp[3]; $x=@tmp[0];$y=@tmp[1];$z=@tmp[2];}
	push(@{$surfdpa{$pert}},$x,$y,$z);
    }close(SDPAF);
#
#--  Find TOPDPA points
    my $flag_sele_TOPDPA='EVDFITTING'; #$flag_sele_TOPDPA indicating the method of selecting TOP DPA POINTS
                                       #default method is EVD fitting methods, see Ming & Wall, JMB 2006

#
#-- select TOP_DPA_DATA: use EVD fitting
    my $cutdpaperct=$cutdpaperct0;my $ndpavalue_section=$ndpavalue_section0;
    my ($seleinfo,$topdpacrd,$extremefita,$extremefitb,$serr,$rcorr);
  decreaseperct:
    ($seleinfo,$topdpacrd,$extremefita,$extremefitb,$serr,$rcorr)=
	&sele_topdpa($jobid,\%surfdpa,$ndpavalue_section,$cutdpaperct,'notshowfit',$topdpa_num_cutoff,$exe_dir,$work_dir,$flag_check);
    return -9999,'SELE_TOPDPA> TOO FEW DPA DATA ( < 10 ) -- NO TOP_DPA' if ($seleinfo eq -9999);
    return -9001,'DPA DATA HANDLING ERR: OPEN,READ,or N>1000' if ($seleinfo <= -9000);

    print "ANADPA_EVFIT> $seleinfo,$topdpacrd,$extremefita,$extremefitb,$serr,$rcorr <<<=== $jobid,\%surfdpa,$ndpavalue_section,$cutdpaperct,'notshowfit',$topdpa_num_cutoff,$exe_dir,$work_dir,$flag_check\n" if $flag_check eq 'true';

    goto decreaseperct2 if ($seleinfo eq -1000);   #can not use EVD fitting
    if ($seleinfo eq -1 && $cutdpaperct > 0.50) {
	$cutdpaperct=$cutdpaperct-0.01;      
	goto decreaseperct;}
    goto decreaseperct2 if($rcorr < 0.80);
    goto cont_cluster if($seleinfo >=1);
#
#-- select TOP_DPA_DATA: DIRECT select a top percentage, DONOT use EVD fitting
  decreaseperct2:
    $flag_sele_TOPDPA='DIRECT';
    ($seleinfo,$topdpacrd,$extremefita,$extremefitb,$serr,$rcorr)=
	&sele_topdpa2($jobid,\%surfdpa,$cutdpaperct2,$topdpa_num_cutoff);
    if ($seleinfo eq -1 && $cutdpaperct2 > 0.90) {
	$cutdpaperct2=$cutdpaperct2-0.01;  
	goto decreaseperct2;}
    return -1,"NO enough TOP_DPA points (< $topdpa_num_cutoff) SELECTED, smaller cutdpaperct (< 0.9) is desired" if ($seleinfo < 0);
#
#-- Clusterizing TOP_DPA points
  cont_cluster:
#######
#    foreach my $id (sort {$a<=>$b} keys %$topdpacrd){print "TDPACHECK> $id -- @{$$topdpacrd{$id}}  \n";} exit();
#######
    my $topdpa_pdbf=$jobid.'_dpacluster.pdb' if($flag_wcpdb eq 'true') ;

    my ($num_effect_cluster,$topdpacluster)=&_CLUSTER::OPTICS($topdpacrd,$epislon_dpa,$MinPts_dpa,$epislon_cluster_dpa,$topdpa_pdbf,$output_dir);
#######    
#    foreach my $cid (sort {$a<=>$b} keys %$topdpacluster) {print "OPTICS> $cid -- @{$$topdpacluster{$cid}}\n";}# exit();
#######
    if($num_effect_cluster<=0) {   #can not form effect cluster with selected TOP_DPA points
	if($flag_sele_TOPDPA eq 'EVDFITTING' && $cutdpaperct > 0.70){
	    $cutdpaperct=$cutdpaperct-0.01;
	    goto decreaseperct;
	}
	elsif($flag_sele_TOPDPA eq 'DIRECT' &&  $cutdpaperct2 > 0.90){
	    $cutdpaperct2=$cutdpaperct2-0.01;
	    goto decreaseperct2;
	}
	else{
	    return -1, "NO DPA_CLUSTER is FORMED at %cutdpa=0.50, need spared_clu_option -clup 5 2 5 "
		if($flag_sele_TOPDPA eq 'EVDFITTING');
	    return -1, "NO DPA_CLUSTER is FORMED at %cutdpa=0.90, need spared_clu_option -clup 5 2 5 "
		if($flag_sele_TOPDPA eq 'DIRECT');
	}
    }
#
#-- Find the C_alphas that TOPDPA cluster contact with
    my ($ndpaclu,$dpabindingsites0,$dpacenters)=&getcontacts_with_ClusterPOINTS($topdpacluster,$topdpacrd,$cacrd,$dpabindingcutoff);
    my $rank_topclu=&_CLUSTER::rank_of_clusters($ndpaclu,$topdpacluster,$topdpacrd);
    my %dpabindingsites=(); undef %dpabindingsites;
    foreach my $cid (sort {$a<=>$b} keys %$dpabindingsites0){
	my $rank=@{$$rank_topclu{$cid}}[-1];
#
#--ADDING rank "0" for NOISE-DPA-POINT set included 
        $rank = 0 if $cid == 0;    #0 for noise cluster
#
#        print "GGGGG>> $ndpaclu -- $cid --> $rank \n";
	return -1, "NO PREDICT from DPA_CLUSTER,too small bind_cutoff $dpabindingcutoff / wired structure"
                  if (scalar(@{$$dpabindingsites0{$cid}})<=0&& $cid==1);
	push (@{${dpabindingsites{$rank}}},@{$$dpabindingsites0{$cid}});
        ######
#	print "GDPACTS0> $jobid-- $cid --RANK: $rank :: @{$$dpabindingsites0{$cid}}\n";
	######
    }
#########vvv  
#    foreach my $cid (sort {$a<=>$b} keys %$topdpacluster) {print "OPTICS> $cid -- @{$$topdpacluster{$cid}}\n";}
#    foreach $cid (sort {$a<=>$b} keys %$rank_topclu){print "TOPCLU_RANK> $cid -- @{$$rank_topclu{$cid}}\n";}
#    foreach my $cid (sort {$a<=>$b} keys %$dpabindingsites0){
#	foreach my $crdid (sort {$a<=>$b} keys %{$$dpacenters{$cid}}){
#	    print "GDPACTR> $jobid--$cid--$crdid: @{$$dpacenters{$cid}{$crdid}}\n";}}
#    foreach my $rkid (sort {$a<=>$b} keys %dpabindingsites){
#	print "GDPACTS> $jobid-- $rkid : @{$dpabindingsites{$rkid}}\n";} 
#    exit();
#########^^^
   return  0,'dpaCluster_done',$ndpaclu,\%dpabindingsites,$dpacenters,$flag_sele_TOPDPA,$cutdpaperct,$cutdpaperct2,$extremefita,$extremefitb,$serr,$rcorr; 
}

sub sele_topdpa{
    my ($jobid,$sdpa,$nsec,$cutperct,$flagshowfit,$topdpa_num_cutoff,$exe_dir,$work_dir,$flag_check)=@_;
#
    my %topdpa=(); #DATA IN topdpa :: X, Y, Z, pointID, TOP, pertValue
    my $flagseledpa=-9999;
    my $ndpadata=0;
#
#-- Extrem Value Fitting (EVF)
    my $cutdpa=-9999;my $serr=9999.; my $rcorr=0; #standard error, correlation 
    my $sdpalistf=$jobid.'_sdpadat.list_tmp0000';       
    my $extremefitinputf=$jobid.'extremefitting_input0000';
    my $efitout=$jobid.'extremefitting_out0000';
    open(SDPALIST,"> $work_dir/$sdpalistf")||die "Cannot open $work_dir/$sdpalistf";
    foreach my $pert (sort {$b<=>$a} keys %$sdpa){
	for(my $k=0; $k<=$#{$$sdpa{$pert}};$k+=3){print SDPALIST "  $pert\n";$ndpadata++;}} #MIGHT BE > 1 POINT FOR A PERT
    close(SDPALIST);
    return  -9999 if $ndpadata < 10;
    open(EXTREMEFITINPUT,"> $work_dir/$extremefitinputf");
    print EXTREMEFITINPUT " '$work_dir/$sdpalistf' $nsec $cutperct $flagshowfit\n";
    close(EXTREMEFITINPUT);
    system ("$exe_dir/extremefit.exe  < $work_dir/$extremefitinputf  > $work_dir/$efitout");
    open(TMPF,"< $work_dir/$efitout"); my $line=0;my $extremefita = my $extremefitb ='null';
    while(<TMPF>){ 
	my @tmp=split /\s+/,$_;$line++;#print 'READFIT> ',$_;
	if($line==1 && @tmp[0] eq ''){$cutdpa=@tmp[1];$serr=@tmp[2];$rcorr=@tmp[3];};
	if($line==1 && @tmp[0] ne ''){$cutdpa=@tmp[0];$serr=@tmp[1];$rcorr=@tmp[2];};
	if($line==2 && @tmp[0] eq ''){ $extremefita=@tmp[1];$extremefitb=@tmp[2];last;};
	if($line==2 && @tmp[0] ne ''){ $extremefita=@tmp[0];$extremefitb=@tmp[1];last;};
    } 
    close(TMPF);
    system ("rm -f $work_dir/$extremefitinputf $work_dir/$efitout $work_dir/$sdpalistf") if($flag_check eq 'false'); 
    
    return -9001 if $cutdpa <= -9000;  # handling SDPA data file problem
    return -1000 if $cutdpa <= 0; #$calc_info=$cutdpa in case EVD failed

    my $ntop=-1;      #DATA IN topdpa :: X, Y, Z, pointID, TOP, pertValue
    foreach my $pert (sort {$b<=>$a} keys %$sdpa){
	last if $pert < $cutdpa;
	for (my $k=0;$k<=$#{$$sdpa{$pert}};$k+=3){
	    $ntop++;
	    push(@{${topdpa{$ntop}}},@{$$sdpa{$pert}}[$k],@{$$sdpa{$pert}}[$k+1],
		 @{$$sdpa{$pert}}[$k+2],$ntop,'TOP',$pert);#print "===> $ntop >>>> @{${topdpa{$ntop}}} \n";	    
	}}
    if ($ntop >= $topdpa_num_cutoff){
	return $ntop,\%topdpa,$extremefita,$extremefitb,$serr,$rcorr;}
    else{
	return -1;
    }
}
sub sele_topdpa2{
    my ($jobid,$sdpa,$cutperct,$topdpa_num_cutoff)=@_;
#
    my %topdpa=(); my $flagseledpa=-9999;
#
    my $cutdpa=-9999;my $serr=9999.; my $rcorr=0; #standard error, correlation 
    my $ndpatotal=0;
    foreach my $pert (sort {$b<=>$a} keys %$sdpa){
	for(my $k=0; $k<=$#{$$sdpa{$pert}};$k+=3){
	    $ndpatotal++;}}
    my $ncut=$ndpatotal*(1.-$cutperct);
    return -1 if($ncut < $topdpa_num_cutoff);
    my $i=0;
    foreach my $pert (sort {$b<=>$a} keys %$sdpa){
	for(my $k=0; $k<=$#{$$sdpa{$pert}};$k+=3){
	    $i++;
	    if($i>=$ncut){$cutdpa=$pert;goto finddpa}}}
  finddpa:
    my $ntop=-1;
    foreach my $pert (sort {$b<=>$a} keys %$sdpa){
	last if $pert < $cutdpa;
	for (my $k=0;$k<=$#{$$sdpa{$pert}};$k+=3){
	    $ntop++;
	    push(@{${topdpa{$ntop}}},@{$$sdpa{$pert}}[$k],@{$$sdpa{$pert}}[$k+1],
		 @{$$sdpa{$pert}}[$k+2],$ntop,'TOP',$pert);
	}}
    if ($ntop >= $topdpa_num_cutoff){
	return $ntop,\%topdpa,-9,-9,-9,0.;}
    else{
	return -1;
    }
}
sub getcontacts_with_ClusterPOINTS{
    ######
    # get C-alpha ID which contacts target clusters, see: %targetcts
    ######
    my $targetcluster=@_[0];                     
    my $targetcrd=@_[1];
    my $cacrd=@_[2];
    my $bcut=@_[3];

    my $ndpaclu=0;                          #number of effective DPA clusters
    my %targetctrs=(); undef %targetctrs;   #ID list for C-alpha's to which each targe_cluster CONTACT
    my %targetcts=(); undef %targetcts;     #CRD infos for each MSMS/DPA_points in eqch targe_cluster
#
#    foreach my $cid (sort {$a<=>$b} keys %$targetcluster) {print "GCPCTS_CHECK_CLU> $cid -- @{$$targetcluster{$cid}}\n"; last if $cid >3; }
#    foreach my $cid (sort {$a<=>$b} keys %$targetcrd){print "GCPCTS_CHECK_TCRD>  $cid   @{$$targetcrd{$cid}}\n"; last if $cid >3;}
#    foreach my $cid (sort {$a<=>$b} keys %$cacrd) {print "GCPCTS_CHECK_CA> $cid -- @{$$cacrd{$cid}}\n"; last if $cid >3;}
#    foreach my $cid (sort {$a<=>$b} keys %$targetcluster) {print "GCPCTS_CHECK_CLU> $cid -- @{$$targetcluster{$cid}}\n";}
#

    foreach my $targetclusterster_id (sort {$a<=>$b} keys %$targetcluster){
############################################################################
#                                                                          #
# WE ADDED NOISE-PREDICTIONS                                               #
#-- NOISE-DPA-points IGNORED  for ligand-binding-site prediction           # 
#                                                                          #
#  #ID <= 0 identifies the targe_cluster formed by "noisy MSMS/DPA_points" #
#	next if($targetclusterster_id le 0);                               #
#                                                                          #
############################################################################
        $ndpaclu++;
	######
        #print "GTARGETCTS_CHECK> $targetclusterster_id -- @{$$targetcluster{$targetclusterster_id}}\n";
	######
	my %targetcenters=(); undef %targetcenters; my $np=-1;
	foreach my $crd_id (@{$$targetcluster{$targetclusterster_id}}){
	    push (@{$targetcenters{++$np}},@{$$targetcrd{$crd_id}}); #print "GTARGETCTS> $targetclusterster_id -- $crd_id -- $np -->@{$$topcrd{$crd_id}}\n";

	    push (@{$targetctrs{$targetclusterster_id}{$crd_id}},@{$$targetcrd{$crd_id}});

	}
	my $targetcontacts=&get_Ligand_contacts_Calpha($cacrd,\%targetcenters,$bcut); #print "GTARGETCTS> $targetclusterster_id --  @{$targetcontacts}   <<---\n";

	$targetcts{$targetclusterster_id}=$targetcontacts;
    }
#########vvv
#    foreach my $did (sort {$a<=>$b} keys %targetcts){
#	print "GTARGETCTS_> $did  - @{$targetcts{$did}}\n"; } 
#    foreach $cid (sort {$a<=>$b} keys %targetctrs){
#	foreach $rol (sort {$a<=>$b} keys %{$targetctrs{$cid}}){
#	    print "GTARGETCTRS_> $cid -- $rol --- @{$targetctrs{$cid}{$rol}}\n";}}
#    exit();
#########^^^

    return $ndpaclu,\%targetcts,\%targetctrs;
}

sub get_Ligand_contacts_Calpha{
    my $cacrd=@_[0];                    #Hash 1: %$cacrd   
    my $hetcrd=@_[1];                   #Hash 2: %$hetcrd
    my $bindingcutoff=@_[2];
    my $bindingcutoffsq=$bindingcutoff**2;
#########
#    foreach (sort{$a<=>$b} keys %$hetcrd){
#	print "GLGDCNT_HETCRD> $bindingcutoff || $_ --  @{$$hetcrd{$_}}\n"
#	}
#    foreach (sort{$a<=>$b} keys %$cacrd){
#	print "GLGDCNT_CA_CRD> $bindingcutoff || $_ --  @{$$cacrd{$_}}\n"
#	}
#    exit();
#########
    my @contacts=();                    #output -- Recording index in cacrd and the distance
    my @d1=();my $dtmp=''; my $dissq='';

    foreach my $caid (keys %$cacrd){              #For each CA searching ALL Ligand atoms
	foreach my $hetid (keys %$hetcrd){
	    @d1=();                              #tmp displacement between Ligand and CA
	    for (my $k=0;$k<=2;$k++){
		$dtmp=@{$$cacrd{$caid}}[$k]-@{$$hetcrd{$hetid}}[$k];
		last if($dtmp > $bindingcutoff or $dtmp < -$bindingcutoff);
		push(@d1,$dtmp);
	    }
	    next if($#d1 <2);                     # at least 1 of 3 components of displacement > cutoff
	    $dissq=(@d1[0]*@d1[0]+@d1[1]*@d1[1]+@d1[2]*@d1[2]);
	    if($dissq<$bindingcutoffsq){          #find the contact CA
		push (@contacts,$caid);   #   print "GCTS> $bindingcutoff -- $caid\n";
		last;                     #stop search when any ligand atom found within the cutoff
	    }
	}
    }

     ######
     # print "GETLGCNT_CA_ID>: $bindingcutoff -- @contacts\n"; 
     ######
    return  \@contacts;
}
###
############### END OF DPA ANALYSIS #######################

1;
