use strict;
#use lib "/xcommon/bin/_mdpa";
#use lib "/Users/dming/mdpa/_mdpa";
use lib "/opt/xcommon/bin/_mdpa";
package _cmpEXP;

############### ANALYSIS OF DPA ######################################################################
#
#   COMPARISON WITH EXPERIMENTAL DATA
#
######################################################################################################

my $vdm_gap=1.0;                 # cutoff that two atom are within VDW contact
my $covalence_cutoff=1.8;        #
my $covalence_bond_gap=0.3;
my ($atom_vdws,$atom_covradii,$pdb_positive_charge_ele,$aromatic_ring)=&constants;
my $pdbid=''; my $selechain='';  my $flag_check='true'; 
my $ligand_contact_cutoff=5.0;

sub cmp_expri{
    my ($pdbid0,$selechain0,$structuredatadir,$workdir,$clu2ligcutoff,
	$nca,$conn,$predictions,$clusters,$cacrd,$msmslayer,$flag_wsite,$flag_positive_pred,$flag_check0)
	 =@_;
    $pdbid=$pdbid0;
    $selechain=uc($selechain0); $selechain=~ s/\s//g; #CHAIN: "A", "BCD", "ALL", etc.
    $ligand_contact_cutoff = $clu2ligcutoff;  #default: 6.0A
    $flag_check=$flag_wsite;
    print "\n...........CMP_EXPERIMETAL_DATA................\n";
    my ($ac_ligand_atoms,$ac_protein_atoms,$ac_protein_resi,$ac_nonligand_protein_resi,
                $ac_list,$ac_ligand_entries,$ac_nonligand)=&get_pdbsite(
		$selechain,$structuredatadir,$workdir,$clu2ligcutoff,$flag_wsite,$flag_check0);
#
#--LBS-PDB-SITE INFO FOR THIS CHAIN:
#
    my %ac_selechain=(); undef %ac_selechain; my $tmpselechain=uc($selechain);
    foreach my $ac (keys %$ac_protein_atoms) {  #AC-ID,e.g. "RQ3 C 80"
          my $tmpchain = '';
          if( $$ac_ligand_entries{$ac} =~ /\w+\s([A-Za-z])\s.\d+/){$tmpchain=$1;}
          elsif($$ac_ligand_entries{$ac} =~ /\w+\s([A-Za-z])\d+/){$tmpchain=$1;}
	  #next if ($tmpselechain ne 'ALL' and $tmpchain ne $selechain);
	  next if ($tmpselechain ne 'ALL' and  $selechain !~ $tmpchain);
          if ($flag_check0 eq 'true') {
	      print "_GSITE_$pdbid"."_ACLIGAND_LBS-LIST_SELECT>> $ac :: @{$$ac_protein_resi{$ac}}\n";
          }
          foreach my $resi (keys %{$$ac_protein_atoms{$ac}}) {
             my $tmpchain = ''; my $ii = '';
             if($resi =~ /\w+\s([A-Za-z])\s+(\d+)/){$tmpchain=$1; $ii=$2}
             elsif($resi =~ /\w+\s([A-Za-z])(\d+)/){$tmpchain=$1; $ii=$2}
             if(! $tmpchain or ! $ii) {
                  print "_GSITE_AC_LBS> CHAIN/RESID READ-ERR>>  $ac :: $$ac_ligand_entries{$ac} ::: $resi\n";
                  next;
	     }
	     #next if ($tmpselechain ne 'ALL' and $tmpchain ne $selechain);
	     next if ($tmpselechain ne 'ALL' and  $selechain !~ $tmpchain);
             push(@{$ac_selechain{$ac}},$ii);
             foreach my $atom (@{$$ac_protein_atoms{$ac}{$resi}}){
#                 print "_GSITE_$pdbid","_ACLIG_LBS_ATOMCRD >> $ac :: $resi >>> $atom\n";
              }
          }
	  print "_GSITE_$pdbid"."_CHAIN_".$selechain."_AC_LBS> $ac :: $$ac_ligand_entries{$ac}",
                " --> @{$ac_selechain{$ac}}\n" if ($flag_check eq 'true');
    }
#
#--LBS-PREDICTION DATA
    my %mdpa=(); undef %mdpa;
    foreach my $rkid (sort {$a<=>$b} keys %$predictions){
    #    printf  "TS>#%3d :: %s\n",$rkid,$$conn{$rkid};
        foreach my $resid (sort {$a<=>$b} @{$$predictions{$rkid}}){
            my $resName=@{$$cacrd{$resid}}[3]; my $chainID=@{$$cacrd{$resid}}[4];
            my $resSeq =@{$$cacrd{$resid}}[5]; my $iCode  =@{$$cacrd{$resid}}[6];
    #        printf  "TS> %3d %3s %1s %5d %1s\n",$rkid,$resName,$chainID,$resSeq,$iCode;
	    push (@{$mdpa{$rkid}{'PREDICT'}},$resSeq);
    }}
#
#-- LBS-PREDICTION-OVERLAP
    foreach my $rkid (sort { $a cmp $b } keys %mdpa){
       print "_MDPA__$pdbid"."_CHAIN_".$selechain."_PREDICT> $rkid :: $$conn{$rkid}",
	     " >> @{$mdpa{$rkid}{'PREDICT'}}\n" if ($flag_check eq 'true');

        my $maxmcc=-999; my @tmpstat=(); undef @tmpstat;
        foreach my $ac (keys %ac_selechain){
	     my ($mcc,$precision,$recall,$tp,$np,$nexpr)=
		&overlap_two_sets(\@{$mdpa{$rkid}{'PREDICT'}},\@{$ac_selechain{$ac}},$nca);
             if ($mcc > $maxmcc) {
                 $maxmcc = $mcc;
                 undef @tmpstat; push (@tmpstat,$ac,$mcc,$precision,$recall,$tp,$np,$nexpr);
             }
        }
        $mdpa{$rkid}{'OVERLAP_AC'} = 'NULL'; 
        if ($maxmcc > -999){
            $mdpa{$rkid}{'CLUSTER'} = $$conn{$rkid};
            $mdpa{$rkid}{'OVERLAP_AC'} = @tmpstat[0]; #$ac; 
            $mdpa{$rkid}{'OVERLAP_AC_ENTRY'} = $$ac_ligand_entries{@tmpstat[0]}; 
            $mdpa{$rkid}{'OVERLAP_MCC'} = @tmpstat[1]; #$ac_ligand_entries{$ac}; 
            $mdpa{$rkid}{'OVERLAP_PRECI'} = @tmpstat[2]; #$precision;
            $mdpa{$rkid}{'OVERLAP_RECALL'} = @tmpstat[3]; #$recall;
            $mdpa{$rkid}{'OVERLAP_TP'} = @tmpstat[4]; #$tp;
            $mdpa{$rkid}{'OVERLAP_NP'} = @tmpstat[5]; #$np;
            $mdpa{$rkid}{'OVERLAP_NEXPR'} = @tmpstat[6]; #$nexpr;
        }
    }
    foreach my $rkid (sort {$a<=>$b} keys %mdpa){
        next if $mdpa{$rkid}{'OVERLAP_AC'} eq 'NULL';
	next if $flag_positive_pred eq 'true' and $mdpa{$rkid}{'OVERLAP_MCC'} <= 0;
        print "_MDPA_$pdbid"."_CHAIN_".$selechain."_LBS-OVERLAP> $rkid :: ",$mdpa{$rkid}{'CLUSTER'};
        printf " <<>> %s :: %s >> %5.2f  %5.2f %5.2f %4d %4d %4d\n",
	     $mdpa{$rkid}{'OVERLAP_AC'},$mdpa{$rkid}{'OVERLAP_AC_ENTRY'},
             $mdpa{$rkid}{'OVERLAP_MCC'},$mdpa{$rkid}{'OVERLAP_PRECI'},$mdpa{$rkid}{'OVERLAP_RECALL'},
	     $mdpa{$rkid}{'OVERLAP_TP'},$mdpa{$rkid}{'OVERLAP_NP'},$mdpa{$rkid}{'OVERLAP_NEXPR'};
	}

    foreach my $rkid (sort {$a<=>$b} keys %mdpa){
        next if $mdpa{$rkid}{'OVERLAP_AC'} eq 'NULL';
        next if $flag_positive_pred eq 'true' and $mdpa{$rkid}{'OVERLAP_MCC'} <= 0;
        my $thisac = $mdpa{$rkid}{'OVERLAP_AC'};
        my %tmpclu=(); undef %tmpclu; my $np=0;      #point-crd of predicted pocket
        foreach my $i (sort {$a<=>$b} keys %{$$clusters{$rkid}}){
            my $x = @{$$clusters{$rkid}{$i}}[0];
            my $y = @{$$clusters{$rkid}{$i}}[1];
            my $z = @{$$clusters{$rkid}{$i}}[2];
       #     print ">>> $rkid ---> $i --->>> $x, $y, $z\n"; 
            push (@{$tmpclu{++$np}} ,$x,$y,$z);
        }
        my %tmplig=(); undef %tmplig;  my $nliga=0; #ligand atom crd in predicted pocket
        foreach my $ligand (keys %{$$ac_ligand_atoms{$thisac}}) {
            foreach my $lgatm (@{$$ac_ligand_atoms{$thisac}{$ligand}}){
#                print "ATOMCRD>> $thisac :: $ligand >>> \n";
                my $lgx = substr($lgatm,30,8);        #
                my $lgy = substr($lgatm,38,8);        #ligand atom crd
                my $lgz = substr($lgatm,46,8);        # 
                push(@{$tmplig{++$nliga}},$lgx,$lgy,$lgz);
#                print ">>> $rkid --> $ligand --> $lgatm -- $lgx -- $lgy -- $lgz\n";
            }
        }
        if($nliga > 0) {
            print "0000000> ",$nliga,"<<<<<<<<<<<";
            my ($precision,$recall,$ncnt,$nfind)= &overlap_two_point_sets($nliga,$np,\%tmplig,\%tmpclu,$clu2ligcutoff);
	    print "_MDPA_$pdbid"."_CHAIN_".$selechain."_LIG-OVERLAP> $rkid :: ",$mdpa{$rkid}{'CLUSTER'};
            printf " <<>> %s :: %s >> %5.2f %5.2f  %4d %4d  %4d %4d\n",
                   $thisac,$mdpa{$rkid}{'OVERLAP_AC_ENTRY'}, $precision,$recall,$ncnt,$np,$nfind,$nliga;
        }
        else{
            print "_MDPA_$pdbid"."_CHAIN_".$selechain."_LIG-OVERLAP> NO LIGAND SITE FOUND IN PDB\n";
        }
    }

    foreach my $rkid (sort {$a<=>$b} keys %mdpa){
        foreach my $ac (keys %{$mdpa{$rkid}{'OVERLAP'}}){
	    print "_MDPA__$pdbid"."_CHAIN_".$selechain.
	      "_OVERLAP> $rkid :: $$conn{$rkid} <<>> $ac :: $$ac_ligand_entries{$ac} ";
	    printf ">> %5.2f  %5.2f %5.2f %4d %4d %4d\n",@{$mdpa{$rkid}{'OVERLAP'}{$ac}};
        }
    }

      print "_GSITE_$pdbid","_ACLIGAND__ALLLIST>> @{$ac_list}\n";
      foreach my $ac (sort keys %$ac_ligand_entries){
          print "_GSITE_$pdbid","_ACLIGAND_____LIST>> $ac :: $$ac_ligand_entries{$ac}\n";
      }
      print "_GSITE_$pdbid","_ACNONLIGAND__LIST>> @$ac_nonligand\n" if(@$ac_nonligand);
      print "_GSITE_$pdbid","_ACNONLIGAND__LIST>> NULL\n" if(! @$ac_nonligand);
      foreach my $ac (keys %$ac_ligand_atoms) {
        foreach my $ligand (keys %{$$ac_ligand_atoms{$ac}}) {
            foreach my $ii (@{$$ac_ligand_atoms{$ac}{$ligand}}){
#                print "_GSITE_$pdbid","_ACLIGAND__ATOMCRD>> $ac :: $ligand >>> $ii\n";
            }
        }
      }
      foreach my $ac (keys %$ac_protein_atoms) {  #LBS: ligand-binding sites
          print "_GSITE_$pdbid"."_ACLIGAND_LBS-LIST >> $ac :: @{$$ac_protein_resi{$ac}}\n";
          foreach my $resi (keys %{$$ac_protein_atoms{$ac}}) {
             foreach my $ii (@{$$ac_protein_atoms{$ac}{$resi}}){
#                 print "_GSITE_$pdbid","_ACLIG_LBS_ATOMCRD >> $ac :: $resi >>> $ii\n";
              }
          }   
      }
      foreach my $ac (keys %$ac_nonligand_protein_resi) {
          print "\n\n\n_GSITE_$pdbid","_ACNONLIG_LBS-LIST > $ac :: @{$$ac_nonligand_protein_resi{$ac}}\n";
      }

    foreach my $rkid (sort {$a<=>$b} keys %{$conn}){
        printf "REMARK  %3d :: %s\n",$rkid,$$conn{$rkid};}
    my $iatom=0; #my @chainid_list=[O,P,Q,R,S,T,U,V,W,X,Y,Z]; 
    my $chainid='';my @chainid_list=("O".."Z"); my $resname='';
    foreach my $cid (sort {$a<=>$b} keys %$clusters) {
        $chainid=@chainid_list[($cid-1)%scalar(@chainid_list)] if($cid >=1);
        $resname='TOP';$resname='NOS' if($cid eq 0);
        $chainid = "x" if($cid eq 0);
        foreach my $i (sort {$a<=>$b} keys %{$$clusters{$cid}}){
            $iatom++; my  $resseq=$i;
            my $x = @{$$clusters{$cid}{$i}}[0];
            my $y = @{$$clusters{$cid}{$i}}[1];
            my $z = @{$$clusters{$cid}{$i}}[2];
            my $pid =@{$$clusters{$cid}{$i}}[-1]; $pid =9999 if $pid >= 9999;
           #printf  "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%4d\n", 'HETATM',
           #$iatom,' CA ',' ',$resname,$chainid,$resseq,' ',$x,$y,$z,1.0,10.0,'','MING',$pid;
           printf "TT>%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%4d\n", 'HETATM',
           $iatom,' CA ',' ',$resname,$chainid,$resseq,' ',$x,$y,$z,1.0,10.0,'','MING',$pid;
         }
    }


    exit();

}

sub overlap_two_point_sets {
    my ($nlga,$np,$lig,$clup,$bcut) = @_;
    my ($precision,$recall); 
    my $ncnt = 0; my $distsq=9999.; my $bcutsq=$bcut*$bcut;
    my ($x,$y,$z,$dx,$dy,$dz);
    foreach my $i (keys %$clup) {
        $x = $$clup{$i}[0];$y = $$clup{$i}[1];$z = $$clup{$i}[2];
        foreach my $j (keys %$lig) {
            $dx=$$lig{$j}[0]-$x;$dy=$$lig{$j}[1]-$y;$dz=$$lig{$j}[2]-$z;
            $distsq =  $dx*$dx + $dy*$dy + $dz*$dz;
            if($distsq <= $bcutsq) {
                $ncnt++;last;
	    }
        }
    }
    $precision = $ncnt / $np;
    my $nfind=0;
    foreach my $i (keys %$lig) {
        $x = $$lig{$i}[0];$y = $$lig{$i}[1];$z = $$lig{$i}[2];
        foreach my $j (keys %$clup) {
            $dx=$$clup{$j}[0]-$x;$dy=$$clup{$j}[1]-$y;$dz=$$clup{$j}[2]-$z;
            $distsq =  $dx*$dx + $dy*$dy + $dz*$dz;
            if($distsq <= $bcutsq) {
                $nfind++;last;
            }
        }
    }
    $recall = $nfind/$nlga;
    return $precision,$recall,$ncnt,$nfind; 
}

sub overlap_two_sets{
    my ($predict,$expri,$nca)=@_;
    my ($precision,$recall,$np,$nexpr);
    $np = scalar @{$predict}; #total prediction-data
    $nexpr = scalar @{$expri}; #total experiment-data
    my $tp = 0; #true-positive
    foreach my $i (@{$predict}){
	foreach my $j (@{$expri}) {
            $tp++ if ($i == $j);
	}
    }
    $precision = -1; $recall = -1;
    $precision = $tp / $np if($np >= 1);
    $recall = $tp / $nexpr if($nexpr >= 1);

    my $fp = $np - $tp;  #false-positive
    my $fn = $nexpr - $tp; #false-nagtive
    my $tn = $nca - $nexpr - $np + $tp; #true-nagtive
    my $mcc = $np*$nexpr*($tn+$fp)*($tn+$fn);
    if($mcc <= 0) {
	$mcc = -1}
    else {
 	$mcc = ($tp * $tn - $fp * $fn)/sqrt($mcc)
    }

    return $mcc,$precision,$recall,$tp,$np,$nexpr;
}

sub get_pdbsite{
    my ($selechain,$structuredatadir,$workdir,$clu2ligcutoff,$flag_wsite,$flag_check)
	 =@_;

    my $pdbfile=$structuredatadir."/".$pdbid.".pdb";
    #print "\n\nGPSITE>  $pdbfile -- $selechain\n";

    my ($finfo,$file_data)= &get_pdb_file_data($pdbid,$pdbfile); #fdata:including  "REMARK", "ATOM  /HETATM", etc.

    my $cnt_patt_outnamef=$pdbid.".patt";
    my $outname = $pdbid.".site";

#get ATOM/HETATM records, for NMR structure ONLY the first MODEL is chosen.
    my ($atoms,$hetatms,$atomhs,$hetatmhs,$atom_hetatm_data) = &choose_atom_and_hetatm( $file_data );

    if($flag_check eq 'true') {
        foreach my $i (@$atoms){    print "_GSITECRD_$pdbid","_ATOM__> $i\n";}
        foreach my $i (@$hetatms){  print "_GSITECRD_$pdbid","_HETATM> $i\n";}
        foreach my $i (@$atomhs){   print "_GSITECRD_$pdbid","_ATOM_H> $i\n";}
        foreach my $i (@$hetatmhs){ print "_GSITECRD_$pdbid","_HETATH> $i\n";}
    }
#get the active center, associated ligands, and other-type site, such as phosporylation-site, etc.
    my ($ac_list,$ac_ligand_entries,$ac_nonligand)=&get_ligand_for_each_AC($file_data);
    if($flag_check eq 'true') {
        print "\n";
        print "_GSITE_$pdbid","_ACLIGAND__ALLLIST>> @{$ac_list}\n";
        foreach my $ac (sort keys %$ac_ligand_entries){ 
	    print "_GSITE_$pdbid","_ACLIGAND_____LIST>> $ac :: $$ac_ligand_entries{$ac}\n";
        }
        print "_GSITE_$pdbid","_ACNONLIGAND__LIST>> @$ac_nonligand\n" if(@$ac_nonligand);
        print "_GSITE_$pdbid","_ACNONLIGAND__LIST>> NULL\n" if(! @$ac_nonligand);
    }

#get the bound ligand atom data from the HETATM records
    my ($ac_ligand_atoms) = &get_ligand_atoms($ac_ligand_entries,$hetatms,$hetatmhs);
    if($flag_check eq 'true') {
      print "\n";
      foreach my $ac (keys %$ac_ligand_atoms) {
        foreach my $ligand (keys %{$$ac_ligand_atoms{$ac}}) {
            foreach my $ii (@{$$ac_ligand_atoms{$ac}{$ligand}}){
                print "_GSITE_$pdbid","_ACLIGAND__ATOMCRD>> $ac :: $ligand >>> $ii\n";
            }
        }
      }
    }

#extract atoms coordinates of the AAs in the SITE
    my ($ac_protein_atoms,$ac_protein_resi,$ac_nonligand_protein_resi) =
        &get_protein_atoms($file_data,$ac_list,$ac_nonligand,$atoms,$hetatms,$atomhs,$hetatmhs,$atom_hetatm_data);

    if($flag_check eq 'true') {
        foreach my $ac (keys %$ac_protein_atoms) {  #LBS: ligand-binding sites
            print "\n_GSITE_$pdbid"."_ACLIGAND_LBS-LIST >> $ac :: @{$$ac_protein_resi{$ac}}\n";
            foreach my $resi (keys %{$$ac_protein_atoms{$ac}}) {
                foreach my $ii (@{$$ac_protein_atoms{$ac}{$resi}}){
                    print "_GSITE_$pdbid","_ACLIG_LBS_ATOMCRD >> $ac :: $resi >>> $ii\n";
                }
            }    
        }
        foreach my $ac (keys %$ac_nonligand_protein_resi) {
            print "\n\n\n_GSITE_$pdbid","_ACNONLIG_LBS-LIST > $ac :: @{$$ac_nonligand_protein_resi{$ac}}\n";
        }
    }
#
#--

   return   $ac_ligand_atoms,$ac_protein_atoms,$ac_protein_resi,$ac_nonligand_protein_resi,
		$ac_list,$ac_ligand_entries,$ac_nonligand;


#
#-- END OF PDB-SITE-INFORMATION-PROCESSING------


#
#
print "STARTING CONTACT-PATTERN--ANALYSIS.....\n";
#
#-- START CALCULATE ATOM-CONTACT PATTERNS
#

#---
#if there are no binding ligand in this protein ,
#print the residue in the SITE then exit
#open( ALL,">>","all_sites" ) || die "Cannot open the file all_sites:$!\n";
#print ALL "_GPSITE_$pdbid",">>> BEGIN FUNCTION-SITES ANNOTATION\n"; 
    if($flag_wsite eq 'true') {
        open( SITE,">","$workdir/$outname" )
            or die "CANNOT OPEN WRITE site file $workdir/$outname:$!\n";
        printf SITE ( "#SIT RES A ###  <=> AC1 LGD A ###  >>> %3s %3s %3s %3s %3s %3s %4s %s\n",
                      'COV','COO','ELE','HDD','HDA',' PI','VAN','CAT');

        if( %$ac_nonligand_protein_resi ){
            my $null = "*** * *** ";
            foreach my $ac ( keys %$ac_nonligand_protein_resi){
                foreach my $a ( @{$$ac_nonligand_protein_resi{$ac}} ){
                    printf SITE ( "SITE %s <:> %s %s >>> %3d %3d %3d %3d %3d %3d %4d %d\n",
                                  $a,$ac,$null,0,0,0,0,0,0,0,0);
#           printf ALL  ( "SITE %s <:> %s %s >>> %3d %3d %3d %3d %3d %3d %4d %d\n",
#                        $a,$ac,$null,0,0,0,0,0,0,0,0);
                }    
            }
        }
        if( @$ac_list == 0 ){
            print "_GPSITE_$pdbid> ","NO BINDING LIGANDS RECORDED FOR THIS ENTRY\n";
            print SITE "END\n";
            close SITE;
            #system("rm -f $pdbid.pdb");
            #exit;
        }
    }

#
# VDW interactions
#get the atom pairs within 5A(defined by $ligand_contact_cutoff) and label the VDW contacts by "8" end-tag
#
my($atom_cnt_list,$ac_cnt_resnames,$ac_cnt_res_entries,$cnt_dist)=
    &get_atom_contact_pattern($ac_ligand_atoms,$ac_protein_atoms,
       $ac_list,$file_data,$ac_ligand_entries,$atom_vdws);
if($flag_check eq 'true') {
    print "\n";
    print "_GPSITE_$pdbid","_AC_CNT_RESIDUES>>> @$ac_cnt_resnames\n";
    print "_GPSITE_$pdbid","_AC_CNT_RESENTRY>>> @$ac_cnt_res_entries\n";
    foreach my $i (sort {$a<=>$b} keys %$atom_cnt_list){
        print "_GPSITE_$pdbid","_ATOM_CNT_PATT> $$atom_cnt_list{$i} \n";
    }
}


#
# H bonds
#find the hydrogen bond between the AC residues and the ligands using the software HBPLUS
# attach "6" FOR HM & HS, or Heteroatom with Mainchain/Sidechain Atoms, A.A. as receptor
#        "5" FOR MH & SH, or Mainchain/Sidechain Atoms with Heteroatom, A.A. as donor
#
my ($atom_cnt_list2,$nhb) = &get_hydrogen($atom_cnt_list,$pdbfile);
if($flag_check eq 'true') {
    print "\n";
    foreach my $i (sort {$a<=>$b} keys %$atom_cnt_list2){
        print "_GPSITE_$pdbid","_HBOND_ATOM_CNT_PATT> $$atom_cnt_list2{$i} \n" if($$atom_cnt_list2{$i} =~ /\|\| .*(5|6)$/);
    }
    print "_GPSITE_$pdbid","_NHBOND> $nhb HYDROGEN BOND FOUND \n";
}

#
#-- metal coordination bonds 
# find the coordinations bonds between the AAs and the metal ligands
#
my $link_ref;
my ($atom_cnt_list3,$ncoord,$links) = &get_coordination($atom_cnt_list2,$file_data);
if($flag_check eq 'true') {
    print "\n";
    foreach my $i (sort {$a<=>$b} keys %$atom_cnt_list3){
        print "_GPSITE_$pdbid","_COORD_ATOM_CNT_PATT> $$atom_cnt_list3{$i} \n" if($$atom_cnt_list3{$i} =~ /\|\| .*2$/);
    }
    print "_GPSITE_$pdbid","_NCOORD> $ncoord METAL COORDINATION BONDs FOUND \n";
}
#
#-- find the covalent bonds between the AAs and ligands
#
my ($atom_cnt_list4,$ncov) =
    &get_covalent($links,$ac_cnt_res_entries,$atom_cnt_list3,$covalence_cutoff);

if($flag_check eq 'true') {
    print "\n";
    foreach my $i (sort {$a<=>$b} keys %$atom_cnt_list4){
        print "_GPSITE_$pdbid","_COVALENT_ATOM_CNT_PATT> $$atom_cnt_list4{$i} \n" if($$atom_cnt_list4{$i} =~ /\|\| .*1$/);
    }
    print "_GPSITE_$pdbid","_NCOVALENT> $ncov COVALENCE BOND FOUND \n";
}





return 0;

}

###################################################################
###################################################################
###################################################################
#choose ATOM/HETATM, only the first MODEL in NMR structure is used
sub get_pdb_file_data{  #get pdb from pdb-database
    my ($pdbid,$filename) = @_;
    my @filedata   = ();

    if(! -e $filename) {
        print "_GPSITE_$pdbid> ERROR: CANNOT OPEN READ PDB-FILE: $filename\n";
        return -9999,\@filedata;
    }
    open(PDBFILE,'<',$filename) || die "_GPSITE_$pdbid> CANNOT OPEN-READ PDB-FILE: $filename\n";
    @filedata = <PDBFILE>;
    close PDBFILE;
    if (!grep(/^SITE/,@filedata) ){
        print "_GPSITE_$pdbid> WARNING: NO SITE RECORDS IN PDB-FILE: $filename!\n";
    }
    return (-1,\@filedata );
}



sub choose_atom_and_hetatm{
    my ( $file_data) = @_;
    my $atom_hetatm;
    my $records=join('',grep( /^(ATOM\s\s|HETATM|ENDMDL)/,@$file_data)); #select ALL ATOM/HETATM

    my @copy=();
    if ( grep( /^ENDMDL/,@$file_data)){     #FOR NMR MODELs
	my @split_rec=split("ENDMDL",$records);
	@copy=split("\n",$split_rec[0]);}   #select ONLY the first NMR MODEL
    else {
	@copy=split("\n",$records);  
    }
    
    my @atoms; my @hetatms;
    my @atomhs; my @hetatmhs;
    foreach my $i (@copy){
	if($i =~ /^ATOM/) { 
	    if ($i =~ /^.{76}\sH/) {
		push (@atomhs,$i)}
	    else {
		push (@atoms,$i)
	    }}
	else {
	    if ($i =~ /^.{76}\sH/) {
		push (@hetatmhs,$i)}
	    else {
		push (@hetatms,$i)
	    }
	}
    }
    

    return \@atoms,\@hetatms,\@atomhs,\@hetatmhs,\@copy;
}
     

# subroutine of getting the active center AC and associated ligands
sub get_ligand_for_each_AC{
    my ($file_data) = @_;
    my @ac_list     = ();    
    my @ac_nonligand = ();  #NON-LIGAND LISTED FOR THIS AC, SUCH AS PHOSPHORYLATION-SITES
    my %ac_ligand_entries= (); undef %ac_ligand_entries;  #RECORDING LIGAND FOR EACH AC
                                                          #usually, ONE LIGAND FOR EACH AC ENTRY
    
    my @ac_id_line  = grep( /^REMARK 800.*SITE_IDENTIFIER/,@$file_data);
    my @ac_ligand_records = grep( /^REMARK 800.*SITE_DESCRIPTION:/,@$file_data );
    

    for( my $i = 0;$i <= $#ac_id_line;$i++){
	$ac_id_line[$i] =~ /SITE_IDENTIFIER: (...)/;
	my $ac_id=$1; 
	if ( $ac_ligand_records[$i] =~ /SITE_DESCRIPTION: BINDING SITE FOR RESIDUE ([ 0-9A-Z]*) ([A-Z]) ?([0-9]*)/ ){
	    push( @ac_list,$ac_id );
	    #reset ligand name based on PDB "ATOM" records
	    my $a = $1;my $b = $2; my $c = $3; my $ligand = "";
	    if (length $a == 2){ $a = " ".$a; } 
	    if (length $c == 1){ $ligand = $a." ".$b."   ".$c." ";}
	    elsif (length $c == 2){$ligand = $a." ".$b."  ".$c." ";}
	    elsif (length $c== 3){$ligand = $a." ".$b." ".$c." ";}
	    else{
		$ligand = $a." ".$b.$c." ";
	    }  	
	    $ac_ligand_entries{$ac_id}=$ligand;  #ligand residue entry
	}
	else {
	    push( @ac_nonligand,$ac_id);    #SUCH AS phosphorylation-site, polysarcchride-modification-
	}
    }
 
    return ( \@ac_list,\%ac_ligand_entries,\@ac_nonligand );
}

# subroutine of getting the binding-ligand atoms, CRD-entries: HETATM/ATOM entries
sub get_ligand_atoms{
    my ($ac_ligand_entries,$hetatms,$hetatmhs)=@_;

    my %ac_ligand_atoms; undef %ac_ligand_atoms;
    foreach my $ac ( keys %$ac_ligand_entries ){
	my $ligand=$$ac_ligand_entries{$ac}; my @lg_atoms=(); undef @lg_atoms;

	foreach my $het (@$hetatms) {
	    if($het =~ /$ligand/) { 
		my $a_atom=substr($het,12,4);
		next if(grep (/$a_atom/,@lg_atoms));
		push(@lg_atoms,$a_atom);             #remove alternative ATOMS	    
		push( @{$ac_ligand_atoms{$ac}{$ligand}},$het);
	    }	    
	}
    }


    return \%ac_ligand_atoms;
}


# get the activesite atoms from the ATOM in this pdbfile
sub get_protein_atoms{
    my ($file_data,$ac_list,$ac_nonligand,$atoms,$hetatms,$atomhs,$hetatmhs,$atom_hetatm_data) = @_;

    my %ac_protein_resi=(); undef %ac_protein_resi;
    my %ac_nonligand_protein_resi=(); undef %ac_nonligand_protein_resi;
    my %ac_protein_atoms=(); undef %ac_protein_atoms;


    my @allsites=grep(/^SITE/,@$file_data);
    
    #get the a.a. in the site with binding ligand
    foreach my $ac ( @$ac_list ){     
	my @ac_aa_resi=(); undef @ac_aa_resi;
	foreach my $line ( grep(/$ac/,@allsites) ){ 
	    next if(substr($line,11,3) !~ /$ac/);
	    for( my $a = 18; $a < 62; $a = $a + 11 ){
		my $aresi=substr($line,$a,10);
		next if $aresi =~ /^ {8}/; #remove empty entries
		push( @ac_aa_resi,$aresi);
	    }
	}
#	push( @{$ac_protein_resi{$ac}},@ac_aa_resi);
	foreach my $aresi (@ac_aa_resi) {
	    my $addresi=0; my @aa_atoms=(); undef @aa_atoms;
	    foreach my $atm (@$atoms) {
		if($atm =~ /$aresi/) {
		    my $a_atom=substr($atm,12,4);  
		    next if(grep (/$a_atom/,@aa_atoms));
		    push(@aa_atoms,$a_atom);             #remove alternative ATOMS		    
		    push( @{$ac_protein_atoms{$ac}{$aresi}},$atm);
		    if(! $addresi) {
			$addresi++;
			push( @{$ac_protein_resi{$ac}},$aresi);	
		    }
		}
	    }
	}
    }
    foreach my $ac ( @$ac_nonligand ){   
	my @ac_aa_resi=(); undef @ac_aa_resi;
	foreach my $line ( grep(/$ac/,@allsites) ){ 
	    next if(substr($line,11,3) !~ /$ac/);
	    for( my $a = 18; $a < 62; $a = $a + 11 ){
		my $aresi=substr($line,$a,10);
		next if $aresi =~ /^ {8}/; #remove empty entries
		push( @ac_aa_resi,$aresi);
	    }
	}

#	push( @{$ac_nonligand_protein_resi{$ac}},@ac_aa_resi);
	foreach my $aresi (@ac_aa_resi) {
	    my $addresi=0;
	    foreach my $atm (@$atoms) {
		if($atm =~ /$aresi/) {
		    push( @{$ac_protein_atoms{$ac}{$aresi}},$atm);
		    if(! $addresi) {
			$addresi++;
			push( @{$ac_nonligand_protein_resi{$ac}},$aresi);	
		    }
		}
	    }
	}
    }

    return \%ac_protein_atoms, \%ac_protein_resi, \%ac_nonligand_protein_resi;
}


sub constants {
    my %atom_vdws  = (# atoms and their van der waals radus
                      'C' => '1.68', 'N' => '1.55', 'O' => '1.50', 'I' => '2.05',
                      'F' => '1.55', 'NA'=> '2.30', 'MG'=> '1.70', 'AL'=> '1.25',
                      'SI'=> '2.10', 'P' => '1.85', 'S' => '1.80', 'CL'=> '1.80',
                      'K' => '2.80', 'CA'=> '1.74', 'MN'=> '1.17', 'FE'=> '1.17',
                      'CO'=> '1.16', 'CU'=> '1.40', 'ZN'=> '1.40', 'LI'=> '1.80',
                      'SM'=> '1.66', 'PB'=> '2.00', 'PT'=> '1.75', 'Y' => '1.62',
                      'LA'=> '1.69', 'DY'=> '1.59', 'BA'=> '1.98', 'HG'=> '1.50',
                      'SR'=> '1.92', 'GD'=> '1.61', 'CS'=> '2.35', 'NI'=> '1.60',
                      'PD'=> '1.63', 'IR'=> '1.26', 'CR'=> '1.17', 'GA'=> '1.90',
                      'PR'=> '1.65', 'LU'=> '1.56', 'CE'=> '1.65', 'EU'=> '1.85',
                      'CD'=> '1.60', 'W' => '1.30', 'YB'=> '1.70', 'RH'=> '1.25',
                      'TL'=> '2.00', 'AG'=> '1.72', 'RB'=> '2.16', 'AU'=> '1.70',
                      'BR'=> '1.90'
        );
    my %atom_covradii  = (# atoms and their covalent radii
                          # from: Beatriz Cordero etc., 2008, "Covalent radii revisited", Dalton Trans. (21): 2832-2838
                          #       Pekka Pyykko etc., 2009, "Molecular Double-Bond Covalent Radii for Elements Li-E112". 
                          #                    Chemistry: A European Journal 15 (46): 12770-12779. 
                          # unit: pm (10^-12), or 0.01A; only "Single Bonds" are recorded

                          'H' =>  '31', 'HE' =>  '28', 'LI' => '128', 'BE' =>  '96',  'B' =>  '84',  'C' => '76',
                          'N' =>  '71',  'O' =>  '66',  'F' =>  '57', 'NE' =>  '58', 'NA' => '166', 'MG' => '141',
                          'AL' => '121', 'SI' => '111',  'P' => '107',  'S' => '105', 'CI' => '102', 'AR' => '106',
                          'K' => '203', 'CA' => '176', 'SC' => '170', 'TI' => '160',  'V' => '153', 'CR' => '139',
                          'MN' => '150', 'FE' => '142', 'CO' => '138', 'NI' => '124', 'CU' => '132', 'ZN' => '122',
                          'GA' => '122', 'GE' => '120', 'AS' => '119', 'SE' => '120', 'BR' => '120', 'KR' => '116',
                          'RB' => '220', 'SR' => '195',  'Y' => '190', 'ZR' => '175', 'NB' => '164', 'MO' => '154',
                          'TC' => '147', 'RU' => '146', 'RH' => '142', 'PD' => '139', 'AG' => '145', 'CD' => '144',
                          'IN' => '142', 'SN' => '139', 'SB' => '139', 'TE' => '138',  'I' => '139', 'XE' => '140',
                          'CS' => '244', 'BA' => '215', 'LA' => '207', 'CE' => '204', 'PR' => '203', 'ND' => '201',
                          'PM' => '199', 'YB' => '187',  'W' => '162', 'PT' => '136', 'AU' => '136', 'HG' => '132',
                          'PB' => '146', 'TH' => '206'
        );

    my $pdb_positive_charge_ele;  #we remove (O 1+) in D3O, (N 1+) in ND4
    $pdb_positive_charge_ele=$pdb_positive_charge_ele."_CA_ZN_CU_MG_MN_K_NA_SR_CD_W_YB_FE_NI_HG_CO";
    $pdb_positive_charge_ele=$pdb_positive_charge_ele."_MO_GA_CE_PB_CS_LI_AU_BA_PT_Y_TL_RB_SM_OS_CR";
    $pdb_positive_charge_ele=$pdb_positive_charge_ele."_AG_AL_LU_PR_PD_TB_EU_LA_RH_ER_RU_D_IR_IN_SB";
    $pdb_positive_charge_ele=$pdb_positive_charge_ele."_GD_BI_HO_DY_";

    my %aromatic_ring   = (
        'PHE'=>'CG|CD1|CD2|CE1|CE2|CZ',
        'TYR'=>'CG|CD1|CD2|CE1|CE2|CZ',
        'TRP'=>'CG|CD1|NE1|CE2|CD2|CE3|CE3|CZ3|CH2|CZ2',
        'HIS'=>'CG|CD2|ND1|NE2|CE1'
        );

    return \%atom_vdws,\%atom_covradii,$pdb_positive_charge_ele,\%aromatic_ring;
}

#----------------------------------------------
#
# generate ligand-protein atom-contact pairs, and sort them by their distances
#
sub get_atom_contact_pattern{
    my ($ac_ligand_atoms,$ac_protein_atoms,$ac_list,$file_data,$ac_ligand_entries,
	$atom_vdws) = @_;
    my @ac_resnames = ();  
    my @ac_res_entries  = ();     
    my %atom_contact_pattern = ();  my %cnt_dist=();
    my %atom_contact_list=();   #this is the same as %atom_contact_pattern except reorder according to the distance increase

    #reading PDB-SITE annotations  #AC contact protein residue names and their entries
    foreach my $line ( grep(/^SITE/,@$file_data)){ 
	for( my $a = 18; $a < 62; $a = $a + 11 ){
	    my $tmpaa=substr($line,$a,3);
	    if($tmpaa !~ /HOH|\ /) {                    #remove water and blanks
		push( @ac_resnames,substr($line,$a,3) );   #get 3-letter a.a.name
		push( @ac_res_entries,substr($line,$a,10) ); #get all 10-letter a.a. name
	    }
	}
    }
       
    #clear the repeat a.a. 
    my %count     = ();
    @ac_resnames  = grep{ (++$count{$_}) < 2 }@ac_resnames; 
    my %seen2=(); 
    foreach (keys %$ac_ligand_entries){ 
	$seen2{substr($$ac_ligand_entries{$_},0,3)}=1; 
    }
    @ac_resnames = grep{ ! $seen2{$_} }@ac_resnames;

    
#   remove the ligand residues/entries recorded in SITE list: @ac_res_entries
#   thus ONLY protein a.a. left in SITE list record: @ac_res_entries
    undef %count; %count=();
    @ac_res_entries = grep{ (++$count{$_}) < 2 }@ac_res_entries;   
    my %seen=(); foreach (keys %$ac_ligand_entries){ $seen{$$ac_ligand_entries{$_}}=1; } 
    @ac_res_entries = grep{ ! $seen{$_} }@ac_res_entries;   
    

#
#-- NOW build ligand-protein atom-contact pairs using ligand_contact_cutoff measurement
#
    my $icontact=0;
    foreach my $ac(keys %$ac_ligand_atoms){
	foreach my $lg (keys %{$$ac_ligand_atoms{$ac}}) {
	    foreach my $lgatm (@{$$ac_ligand_atoms{$ac}{$lg}}){  #FOR each ligand atom
		
		my $lgres=substr($lgatm,17,3);
		my $lgx = substr($lgatm,30,8);        #
		my $lgy = substr($lgatm,38,8);        #ligand atom crd
		my $lgz = substr($lgatm,46,8);        # 
		my $lgatmname = substr($lgatm,76,2); $lgatmname =~ s/\s//g;     #ligand atom name
		if(! $lgatmname) {
		    print "_GPSITE_$pdbid","_LGPROT_ATOM_CNT> WARNING_5 LIGAND ATOMNAME (\"$lgres\",\"$lgatmname\") NOT DEFINED IN ELEMENT COLUMN\n";
		    $lgatmname=substr($lgatm,12,2); $lgatmname =~ s/\s//g; 
		}
		my $rlgatm=9999.; #radius
		if(exists $$atom_vdws{$lgatmname}) {
		    $rlgatm=$$atom_vdws{$lgatmname};}
		else{
		    print "_GPSITE_$pdbid","_LGPROT_ATOM_CNT> WARNING_50 NO VDW RADIUS FOR LIGAND ATOM (\"$lgres\",\"$lgatmname\"), 9999 Angstrom USED\n";
		} 
		
		foreach my $resi (keys %{$$ac_protein_atoms{$ac}}) {
		    foreach my $atm (@{$$ac_protein_atoms{$ac}{$resi}}){
			
			my $resname=substr($atm,17,3);
			my $atmx = substr($atm,30,8);        #
			my $atmy = substr($atm,38,8);        #a.a. atom crd
			my $atmz = substr($atm,46,8);        # 
			my $atmname = substr($atm,76,2); $atmname =~ s/\s//g; #a.a. atom name
			if(! $atmname) {
			    print "_GPSITE_$pdbid","_LGPROT_ATOM_CNT> WARNING_5 A.A. ATOMNAME (\"$resname\", $atmname) NOT DEFINED IN ELEMENT COLUMN\n";
			    $atmname=substr($atm,12,4); $atmname =~ s/\s//g; $atmname=substr($atmname,0,1);
			}

			my $ratm=9999.;
			if(exists $$atom_vdws{$atmname}) {
			    $ratm=$$atom_vdws{$atmname};}
			else{
			    print "_GPSITE_$pdbid","_LGPROT_ATOM_CNT> WARNING_50 NO VDW RADIUS FOR A.A. ATOM (\"$resname\",\"$atmname\"), 9999 Angstrom USED\n";
			} 

			my $d = sqrt(($lgx-$atmx)**2+($lgy-$atmy)**2+($lgz-$atmz)**2);  

			print "DDD> $d  :::: $ligand_contact_cutoff << $resi :: substr($atm,12,4) ---- $lg :: substr($lgatm,12,4)\n";

			if ( $d< $ligand_contact_cutoff ){

#			    my $rand1=rand()/9999.;
#			    $d = sprintf("%6.5f",$d+$rand1); 

			    $atom_contact_pattern{++$icontact}=
				substr($atm,12,4)." ".substr($atm,17,11).
				substr($lgatm,12,4)." ".substr($lgatm,17,11).
				$ac." ".$d." || ";
			    $cnt_dist{$icontact}=$d;

			    if ( $d < $rlgatm+$ratm+$vdm_gap) {
				$atom_contact_pattern{$icontact}=$atom_contact_pattern{$icontact}."8";
			    }
			    print "DDD>ATM_PATTERN> ($d < $rlgatm+$ratm+$vdm_gap) $atom_contact_pattern{$icontact} ||  $atom_contact_pattern{$d}\n";    exit;
			}

		    }
		}

	    }
	}	

    }

#SORT ACCORDING TO CONTACT DISTANCE
    my $i=0;
    foreach my $ii (sort {$cnt_dist{$a}<=>$cnt_dist{$b}} keys  %cnt_dist) {
	$atom_contact_list{++$i}=$atom_contact_pattern{$ii};
    }
                                     
    return (\%atom_contact_list,\@ac_resnames,\@ac_res_entries,\%cnt_dist);
}


#subroutine of getting hydrogen bonds
sub get_hydrogen{
    my ($atom_cnt_list1,$pdbfilename) = @_;
    system( "hbplus $pdbfilename > _tmp_hbplus_87601_pdbfilename" );
    system ("rm -f _tmp_hbplus_87601_pdbfilename");

    my $hbfile=$pdbfilename; 
    $hbfile =~ s/\.pdb$/\.hb2/;
    
    if(! -e $hbfile) {
	print "_GPSITE_$pdbid","_HBOUND_BUILD> WARNING_10 CANNOT GENERATE HBOND FILE WITH HBPLUS\n";
	return $atom_cnt_list1,0;
    }
    open(HBPLUS,"<","$hbfile");
    my @hbplus = <HBPLUS>;
    close HBPLUS;    
    system("rm -f $hbfile hbdebug.dat") if($flag_check ne 'true');

 
    my @hydrogen = grep{$_ !~ /HOH/}@hbplus;  #remove water-involved HYDROGEN BONDS
    @hydrogen = grep(/MH|HM|SH|HS/,@hydrogen); #select the records including the protein atom
                                               #M-mainchain, S-sidechain, H-heteroatoms
    if ( !@hydrogen ){
	print STDOUT "_GPSITE_$pdbid","_HBOUND_BUILD> NO LIGAND-INVOLVED HYDROGEN BOND FOUND\n";
	return $atom_cnt_list1;
    }

    my $nhb=0;
    foreach my $hbline ( @hydrogen ){
	my $resseq=substr($hbline,1,4);
	my $icode=substr($hbline,5,1); $icode=' ' if $icode eq "-";
	my $resname=substr($hbline,6,3);
	my $chain=substr($hbline,0,1);
	my $atmtype=substr($hbline,9,4);
	if($resseq =~ /^0*(.*)/) {
	    $resseq=$1; 
	    $resseq='    '.$resseq if(length $resseq ==0); 
	    $resseq='   '.$resseq if(length $resseq ==1); 
	    $resseq='  '.$resseq if(length $resseq ==2); 
	    $resseq=' '.$resseq if(length $resseq ==3); 
	}

	my $resseq2=substr($hbline,15,4);
	my $icode2=substr($hbline,19,1); $icode2=' ' if $icode2 eq "-";
	my $resname2=substr($hbline,20,3);
	my $chain2=substr($hbline,14,1);
	my $atmtype2=substr($hbline,23,4);
	if($resseq2 =~ /^0*(.*)/) {
	    $resseq2=$1; 
	    $resseq2='    '.$resseq2 if(length $resseq2 ==0); 
	    $resseq2='   '.$resseq2 if(length $resseq2 ==1); 
	    $resseq2='  '.$resseq2 if(length $resseq2 ==2); 
	    $resseq2=' '.$resseq2 if(length $resseq2 ==3); 
	}

	if ($hbline =~ /HM|HS/){  #Hetero-atom + a.a. atom
	    my $hcnt=$atmtype2." ".$resname2." ".$chain2.$resseq2.$icode2.' ';
	    $hcnt=$hcnt.$atmtype." ".$resname." ".$chain.$resseq.$icode.' ';
	    foreach my $i (keys %$atom_cnt_list1){
		if ($$atom_cnt_list1{$i} =~ /$hcnt/){
		    $$atom_cnt_list1{$i} .= "6";  #FOR HM & HS
		    $nhb++;
		    last;
		}
	    }
	}
	else{	   #a.a. atom + hetero-atom
	    my $hcnt=$atmtype." ".$resname." ".$chain.$resseq.$icode.' ';
	    $hcnt=$hcnt.$atmtype2." ".$resname2." ".$chain2.$resseq2.$icode2.' ';
	    foreach my $i (keys %$atom_cnt_list1){
		if ($$atom_cnt_list1{$i} =~ /$hcnt/){
		    $$atom_cnt_list1{$i} = $$atom_cnt_list1{$i}."5";  #FOR MH & SH
		    $nhb++;
		}
	    }

	}  
    }
    
    return $atom_cnt_list1,$nhb;
}


#subroutine of getting the coordination bonds
sub get_coordination{
    my ($atom_cnt_list1,$file_data) = @_;

    my ($metallist) = &check_metal_contained($file_data);  #get metal list, either from REMARK 620 OR HETATM METALS
    my @link = grep(/^LINK/,@$file_data);

    my @coordination = ();                 #record a.a. atom entry + metal atom entry, a atom_contact_pattern format
    foreach my $i ( @$metallist ){   # $i is metal-atom entry
	foreach my $lk (@link) {
	    next if($lk !~ /$i/);
	    my $l1=substr($lk,12,16);  
	    if($l1 =~ /$i/) {               #link line: "metal atom" list first 
		push( @coordination,substr($lk,42,16).$i) } #list a.a. atom entry + "metal atom" entry
	    else {                          #link line: "metal atom" list last
		push( @coordination,$l1.$i)                 #list a.a. atom entry + "metal atom" entry
	    }
	}
    }
    @coordination = grep(!/HOH/,@coordination);  #remove water

    my $nlink=0;
    foreach my $a_coord ( @coordination ){
        print "\n\nGPSITE_COVAL> $a_coord\n";
	my $match=0;
	foreach my $i (keys %$atom_cnt_list1) {
	    my $pattline=$$atom_cnt_list1{$i};
            print "GPSITE_COVAL>........PATTLINE: $pattline\n";
	    if ($pattline =~ /$a_coord/){
		$$atom_cnt_list1{$i} = $$atom_cnt_list1{$i}."2";
		$match=1; $nlink++;
		last;
	    }
	}
	if(! $match) {
	    print "_GPSITE_$pdbid"."_METAL_COORD> WARNING_60 COORDINATION:\"$a_coord\" NOT FOUND IN ATOMS_PATTERN\n";
	}
        print ">GPSITE_COCAL> $a_coord\n\n\n";
    }
    
    return ($atom_cnt_list1,$nlink,\@link);
}
 

#subroutine of getting covalent bonds 
sub get_covalent{
    my ($links,$cnt_res_entry1,$atom_cnt_list1,$covalence_cutoff) = @_;
    my $ncov=0;

    my @links1     = ();    
    foreach my $a_link ( @$links  ){
	if (substr($a_link,73,5) < $covalence_cutoff){  #check the distance
	    push( @links1,$a_link);
	}
    }
    my @covalent    = ();
    foreach my $a_res ( @$cnt_res_entry1 ){
	push( @covalent,grep( /$a_res/,@links1 ) );
    }
    if(! @covalent) {
	return $atom_cnt_list1,$ncov;
    }

    foreach my $i (keys %$atom_cnt_list1) {
	my $pattline=$$atom_cnt_list1{$i};
	
	next if ($pattline =~ /| .*2$/);  #NO COVALENT BOND FOR METAL-COORDINATE COMPLEX

	foreach my $a_cova( @covalent ){
	    my $part1=substr($a_cova,12,15);
	    my $part2=substr($a_cova,42,15);
	    if ( $pattline =~ /$part1/ && $pattline =~ /$part2/) {
		$$atom_cnt_list1{$i} = $$atom_cnt_list1{$i}."1";
		$ncov++;
		last;
	    }
	}
    }
    
    return $atom_cnt_list1,$ncov;
}

#
#-- subroutine of figuring out the atom patterns with electrical force
#
sub get_electro_interactions{
    my ($atom_cnt_list1,$ligand_atoms,$ac_cnt_res_entry1,$atom_vdws1) = @_;    
    my $regex;

    my $ne=0;
    foreach my $i (keys %$atom_cnt_list1) {  #atom-contact-pattern list
	my $a_atom_cnt=$$atom_cnt_list1{$i};
#	print "PATTEF>$a_atom_cnt\n"; 

#	my $bonds;
#	if($a_atom_cnt =~ /\|\|\s(.*)/) {
#	    $bonds=$1; $bonds =~ s/^8//g;
#	}
#	next if($bonds);   #exis one or a few of HYDROGEN-, COVELENCE-, COORD-BONDS

	my ($is_electro)=&check_electro_interaction( $a_atom_cnt,$ligand_atoms,$atom_vdws1 );

	if($is_electro) {
 	    $a_atom_cnt = $a_atom_cnt."4";
	    $$atom_cnt_list1{$i}=$$atom_cnt_list1{$i}."4";
	
	    $ne++;
	}

    }
    
    return $atom_cnt_list1, $ne;
}

sub check_electro_interaction {   #check electro-static interaction for each atom-atom pair
    my ( $a_atom_cnt,$ac_all_ligand_atoms,$atom_vdws) = @_;


    my $nelectro=0;
    my $resname=substr($a_atom_cnt,5,3);
    my $charge_atm='';
#
#-- charged a.a. side-chain groups
    if ($resname=~/ARG/){
	 $charge_atm = "NH..ARG";}   # "+" positve charge
    elsif ($resname =~ /LYS/){	
	$charge_atm = "NZ..LYS";}    # "+" positive charge   
    elsif ($resname =~ /HIS/){
	$charge_atm = "NE2.HIS";}    # "+" positive charge
    elsif ($resname =~ /GLU/){
	$charge_atm = "OE..GLU";}    # "-" negative charge
    elsif ($resname =~ /ASP/){
	$charge_atm = "OD..ASP";     # "-" negative charge
    }
    return 0 if(! $charge_atm);

    return 0 if ($a_atom_cnt !~ /$charge_atm/);

    if ( $charge_atm =~ /^N/ ){     #POSITIVELY CHARGED SIDE-CHAIN IN A.A.
	
	if ( $a_atom_cnt =~ /^(.{16})( F| I|CL|BR)/ ) { #electro-interaction between F/I/CL/BR in ligand and positive a.a. atoms
	    return 1;}	
	elsif ($a_atom_cnt =~ /^.{16}( S| O)/ ){  #for O/S in ligand and positive charged a.a. atom	    
	    my $lgatm=substr($a_atom_cnt,16,16);  # my $lgresname=substr($lgatm,5,3);
	    my ($multi_conn,$nconn)=&is_atom_multiple_connection($lgatm,$ac_all_ligand_atoms); #check if ligand atom in multiple-connection

#	    print "CHECK_ELE>>> $a_atom_cnt---> $multi_conn,$nconn\n";
#
# check if ligand is of particular compound
#	    my ($is_miscell_cmpnd)=&is_paticular_compund($lgatm);

#	    print ">>>$lgatm<<< ---->>>>>$multi_conn :: $nconn\n"; exit;
	    
	    return 1 if(! $multi_conn)  # || $is_miscell_cmpnd);
	}
	else {
	    
	    print "_GPSITE_$pdbid","_ELECTRO> WARNING_70 NO DEFINED ELEC_INTERACTION (+): $a_atom_cnt\n" if $a_atom_cnt !~ /^.{16} (C|N)/;	    
	    return 0;
	}
    }
    elsif( $charge_atm =~ /^O/ ){  
	my $lgatm=substr($a_atom_cnt,16,16); 
	my $lgatmname=substr($lgatm,0,2); $lgatmname=~ s/\s//g;

	if($pdb_positive_charge_ele =~/\_$lgatmname\_/      #for metal atoms
	   || $a_atom_cnt =~ /^.{16} N/ ){     	#for N in ligand and negative charged a.a. atom	    

	    return 2;
	}
	else {	    
	    print "_GPSITE_$pdbid","_ELECTRO> WARNING_70 NO DEFINED ELEC_INTERACTION (-): $a_atom_cnt\n" if $a_atom_cnt !~ /^.{16} (C|O)/;	    
	    return 0;
	}
    }
    else {
	return 0;
    }
    
}


#check if given ligand is particular PDB ligand cluster, in which
#      electro-interaction exist even when O-/S- are in multiple connections...
sub is_paticular_compund {
    my ($lgname)=@_;
    my @paticular_pdb_ligands=(
	
	
	);

    return 0;
}

#check how many atoms link to the given atom in the given ligand
sub is_atom_multiple_connection {
    my($atom,$compound)=@_;              #the given atom
    my $lgres_entry=substr($atom,5,10);  #the given ligand

    my $multi_conn=0;
    my $nconn=0;     #counting how many atoms in the given ligand linked to given atom
    
    my $ac0='###'; my $lg0='###';  

    foreach my $ac (keys %$compound) {         #find ligand in AC-SITE LIST
	foreach my $lg (keys %{$$compound{$ac}}){

	 #   print "\n\n>OO>>> $ac --->>$lg<<>>>$lgres_entry<<< ::>> @{$$compound{$ac}{$lg}}\n"; 	    

	    if ($lg =~ /$lgres_entry/ or $lgres_entry =~ /$lg/){
		$ac0=$ac; $lg0=$lg;
		last;
	    }
	}
	last if $ac0 ne '###';
    }
    if( $ac0 eq '###' || $lg0 eq '###') {
	print "_$pdbid> WARNING:: CANNOT ANALYZE LIGAND-ATOM: \"$atom\", LIGAND RESIDUE NOT DEFINED!\n";
    }

    
    my @crd0=grep (/$atom/,@{$$compound{$ac0}{$lg0}});
    if(! @crd0) {
	
	#$atom=substr($atom,0,4)."(A|B|C|D)".substr($atom,5,10);  #sometime, A/B/C/D/.. are inserted before HETATM resname
	$atom=substr($atom,0,4)."([a-zA-Z0-9]{1})".substr($atom,5,10);  #sometime, A/B/C/D/.. are inserted before HETATM resname
	@crd0=grep (/$atom/,@{$$compound{$ac0}{$lg0}});
    }
    

    my $atmcrd=$crd0[0]; 
    my $r_cov= &get_atom_covalence_radii($atmcrd);       #for given atom
    
    my @lgcrds=grep (! /$atom/,@{$$compound{$ac0}{$lg0}});
    foreach my $a_crd (@lgcrds) {                        #for other atoms in the given ligand
	
	my $r_cov1= &get_atom_covalence_radii($a_crd);
	
	my $dd=&atom_dist($atmcrd,$a_crd); 
	if($dd >= 9999) {
	    print "_$pdbid> UNDEFINED DISTANCE FOR UNDEFINED ATOMS IN LIGAND RESIDUE: \"$lg0\"\n";
	    next;
	}
	$nconn++ if($dd < $r_cov1+$r_cov+$covalence_bond_gap);
    }
    
    $multi_conn=1 if $nconn >= 2;             #finish calc. on given ligand
    return $multi_conn,$nconn;
}
sub atom_dist{
    my ($atm1,$atm2)=@_;
    my $dist=99999.0;
    return $dist if (! $atm1 || ! $atm2);

    my $x1=substr($atm1,30,8); my $y1=substr($atm1,38,8); my $z1=substr($atm1,46,8);
    my $x2=substr($atm2,30,8); my $y2=substr($atm2,38,8); my $z2=substr($atm2,46,8);

    return $dist if(! $x1 || ! $y1 || ! $z1 || ! $x2 || ! $y2 || ! $z2);
	
    $dist=sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);

    return $dist;
}

sub get_atom_name {
    my ($pdb_atom_hetatm_line)=@_;
    my $resname1=substr($pdb_atom_hetatm_line,17,3);
    my $atmname1 = substr($pdb_atom_hetatm_line,76,2); $atmname1 =~ s/\s//g;     #ligand atom name
    if(! $atmname1) {
	print "_$pdbid> WARNING:: ATOM OF RESIDUE \"$resname1\" NOT NAMED BY ELEMENT $atmname1\n";
	$atmname1=substr($_,12,2); $atmname1 =~ s/\s//g; 
    }
    return $atmname1;
}
sub get_atom_covalence_radii {
    my ($pdb_atom_hetatm_line)=@_;
    my ($atmname)=&get_atom_name($pdb_atom_hetatm_line);
    my $r_cov=-9999.;
    my $resname1=substr($pdb_atom_hetatm_line,17,3);
    if(! $atmname) {
	print "_$pdbid> WARNING:: COVALENCE RADII NOT DEFINED FOR UNDEFINED ATOM IN RESIDUE $resname1, NEGATIVE RADII USED!\n";
    }
    elsif(exists $$atom_covradii{$atmname}) {
	$r_cov=	$$atom_covradii{$atmname}/100.   #the unit is pm (0.01A) in the table
    }
    else {
	print "_$pdbid> WARNING:: COVALENCE RADII NOT DEFINED FOR ATOM: $atmname, NEGATIVE RADII USED!\n";
    }

    return $r_cov;
}
#subroutine of check if there are electrical atoms in the ligands
sub label_e_atom{
    my ( $atom_cnt_list1,$regex,$ligand_atoms,$atom_vdws1 ) = @_;
    
    if ( $regex =~ /^N/ ){
	foreach my $a_pat ( grep(/$regex/,@$atom_cnt_list1) ){
	   # my $dis=substr( $a_pat,36,7 ); ">>>$dis<<<\n"; exit;
	   # next if (  $dis > 4.0 );

	    if ( $a_pat =~ /.{16,17}(F|CL|BR|I)/ ){
		$a_pat = $a_pat."4";

	    }elsif ( $a_pat =~ /.{17}(O|S)/ ){
		my $index = substr( $a_pat,34,1 );
		my $lig_atom = substr( $a_pat,16,15 );
		my @matched = grep(/$lig_atom/,@{$ligand_atoms->[$index]});
		unless ( @matched ){
		    next;
		}	
		my $a_matched = shift @matched;	
		my $x = substr( $a_matched,30,8);
		my $y = substr( $a_matched,38,8);
		my $z = substr( $a_matched,46,8);
		my $count = 0;
		foreach my $line ( @{$ligand_atoms->[$index]} ){
		    my $a = substr($line,30,8);
		    my $b = substr($line,38,8);
		    my $c = substr($line,46,8);
		    my $d = sqrt( ($x-$a)**2 + ($y-$b)**2+ ($z-$c)**2 ); 
		    if ( $d < 1.8 && $d != 0 ){
			$count++;
		    }
		}
		if ( $count < 2){
		    $a_pat = $a_pat."4";
		}
	    }else {
		next;
	    }
	}
    }elsif ( $regex =~ /^O/ ) {
	foreach my $a_pat ( grep(/$regex/,@$atom_cnt_list1) ){
	    if ( substr( $a_pat,35,6 ) > 4.0 ){
		next;
	    }
	    my $ligand = substr( $a_pat,14,3 );
	    if ( grep(/$ligand/,grep(!/O|S|F|CL|I|C/,keys %$atom_vdws1) ) ){
		$a_pat = $a_pat."4";
	    }else{
		next;
	    }
	}
    }else {
	print "There are no electrical force between this protein and its ligand!\n";
    }
}
			
#subroutine of finding the pi-pi interaction 
sub get_pi_pi{

    my ($filedata,$ac_protein_resi,$ac_ligand_entries,$atom_cnt_list1,$ligand_atoms) = @_;
    
    #check wether there are organic compounds in the 
    #ligands  according to the record FORMUL 
    my @formul = grep(/^FORMUL/,@$filedata); 
    my @formula = ();
    my $lg_pi_resi   = '';
    
    foreach my $a_formul ( @formul ){
	chomp($a_formul); $a_formul =~ s/^\s+|\s+$//g;
	my $lgformul = substr( $a_formul,19);       # print "LGFORMUL>>> $lgformul\n";
	if ( $lgformul =~ /C([0-9]+)/ && $1 >= 5){   #select ligands that have >= 5 carbon atoms
	    $lg_pi_resi=$lg_pi_resi."_".substr($a_formul,12,3)."_";
	}
    }
    unless( $lg_pi_resi){
	print "_$pdbid> WARNING_81: NO CARBON-RING FOUND IN LIGAND RESIDUES\n";
	return $atom_cnt_list1;
    }

    #find the ring in the ligands 
    #check if there are aromatic AAs ,if YES then get the ring in the associated ligand 

    my $npi=0;

    my @aromatic_cnt_list=(); undef @aromatic_cnt_list;
    foreach my $i (keys %$atom_cnt_list1) {
	push(@aromatic_cnt_list,$$atom_cnt_list1{$i});
    }	
    
    my %match=(); undef %match;
    my %match0=(); undef %match0;
    foreach my $ac (keys %$ligand_atoms){  #for a given AC

	#get the aromatic a.a. in this AC 
	my @aromatic_aa = ();
	push(@aromatic_aa,grep(/TRP|PHE|TYR|HIS/,@{$$ac_protein_resi{$ac}}));
	next if (! @aromatic_aa);

	my @lgatoms=(); undef @lgatoms;  #ALL ATOMS OF THIS LIGAND RESIDUE
	foreach my $lgres ( keys %{$$ligand_atoms{$ac}}){

	    my $lgresname=substr($lgres,0,3);
	    next if($lg_pi_resi !~ /$lgresname/);
	    @lgatoms=@{$$ligand_atoms{$ac}{$lgres}};
	    next if(! @lgatoms);
	    
	    #get the ring atom list in the ligand residue
	    my (@lg_ring_atoms) = &get_ligand_ring_atoms( @lgatoms);  
	    next if(! @lg_ring_atoms);
	    
	#    print "\n\nRING>>$ac>>$lgres>>>\n  @lg_ring_atoms\n"; print "\n@aromatic_aa\n";


	    foreach my $an_aromatic (@aromatic_aa){
		my $aa  = substr( $an_aromatic,0,3 );	
		my @aa_ring  = grep( /^$$aromatic_ring{$aa}/,grep( /$an_aromatic/,@aromatic_cnt_list ) ); #a.a. aromatic atoms
		my $count    = 0; 
		foreach my $a_ring_atom ( @lg_ring_atoms ){
		    $count++ if grep( /$a_ring_atom/,@aa_ring);    #ligand ring atom   
		    #my $n=grep( /$a_ring_atom/,@aa_ring);
		    #$count=$count+$n;
		}
		if ($count >=3) {
		    $match0{$ac}{$lgres}{$an_aromatic}=$count;
		    $match{$lgres}{$an_aromatic}=$count;
		    $npi++;
		}
	    }
	}
    }

    if($npi >= 1) {
	foreach my $i (keys %$atom_cnt_list1) {	
	    my $a_atom_cnt=$$atom_cnt_list1{$i};
	    my $resentry=substr($a_atom_cnt,5,10);
	    my $lgresentry=substr($a_atom_cnt,21,10);
	    if(exists $match{$lgresentry}{$resentry}){
		$$atom_cnt_list1{$i}=$$atom_cnt_list1{$i}."7";
	    }	
	}
    }

    return $atom_cnt_list1,$npi;
    
}


#subroutine of getting the ligand's aromatic ring atoms
sub get_ligand_ring_atoms{
    my (@lig_ring) = @_;
    my $times    = 0;
    
    for( my $i = 0 ;$i <= $#lig_ring;$i++){
	my $a = substr( $lig_ring[$i],30,8 );
	my $b = substr( $lig_ring[$i],38,8 );
	my $c = substr( $lig_ring[$i],46,8 );
	my $count = 0;
	for( my $j = 0;$j <= $#lig_ring;$j++){
	    my $x = substr( $lig_ring[$j],30,8 ); 
	    my $y = substr( $lig_ring[$j],38,8 );
	    my $z = substr( $lig_ring[$j],46,8 );
	    my $d = sqrt( ($x-$a)**2+($y-$b)**2+($c-$z)**2 );
	    if ( $d<1.8 && $d != 0 ){
		$count++;
	    }
	}
	if ( $count == 1 ){
	    splice( @lig_ring,$i,1 );
	    $times++;
	}
	unless( @lig_ring ){ return; }
    }
    if ( @lig_ring ){
	my @ring_atoms = ();
	if ( $times == 0 ){
	    foreach my $an_atom ( @lig_ring ){
		my $lig_atom = substr( $an_atom,12,4 )." ".
		    substr( $an_atom,17,11 );
		push( @ring_atoms,$lig_atom ); 
	    }
	    return @ring_atoms;
	}else{
	    &get_ligand_ring_atoms( @lig_ring );
	}	
    }
}

#subroutine of labelling the active sites of the enzyme in the PDB 
sub label_catalytic_sites{
    
    my $csa_dbf="/databases/structures/csa/csa.1";
    
    my %cat_list=(); undef %cat_list;
    if(! -e $csa_dbf){
	print "_$pdbid> WARNING_10 CSA FILE NOT FOUND: $csa_dbf\n";
	return \%cat_list;
    }
    open( FILE,"<",$csa_dbf); 
    my @catasites = <FILE>;
    chomp @catasites;
    close FILE;
    foreach my $a_csite (@catasites) {
	my @stmp = split("\ ",$a_csite);
	if($stmp[0] =~ /$pdbid/i) {
	    my $cat_entry=uc($stmp[2])." ".uc($stmp[3]);
	    if($stmp[4] >=1000) {
		$cat_entry=$cat_entry.$stmp[4];
	    }
	    elsif($stmp[4] >=100) {
		$cat_entry=$cat_entry." ".$stmp[4];
	    }
	    elsif($stmp[4] >=10) {
		$cat_entry=$cat_entry."  ".$stmp[4];
	    }
	    elsif($stmp[4] >=0) {
		$cat_entry=$cat_entry."   ".$stmp[4];
	    }
	    elsif($stmp[4] >=-9) {
		$cat_entry=$cat_entry."  ".$stmp[4];
	    }
	    elsif($stmp[4] >=-10) {
		$cat_entry=$cat_entry." ".$stmp[4];
	    }
	    elsif($stmp[4] >=-100) {
		$cat_entry=$cat_entry.$stmp[4];
	    }
	    
	    $cat_list{$cat_entry." "}=1;  #no icode
	}
    }
    
    return \%cat_list;
}


sub check_metal_contained{
    my ($file_data)=@_;
#
#-- TO CHECK IF "XX" IS METAIL ION, CHANGE IT TO "_XX_" AND THEN USE STRING-MATCH
#
#

    my @metallist    = grep(/REMARK 620 {10,}\S+/,@$file_data);

    if (@metallist){
	@metallist = map{$_ =substr($_,50,5).substr($_,39,10);} @metallist;  #METALLIST: "ATOMNAME RESIDUE_IDENTITY"
	my %count = ();
	@metallist = grep{(++$count{$_}) < 2} @metallist;    #clear repeatition
    }

    my @hetatms = grep(/^HETATM/,@$file_data);
    foreach my $i (@hetatms) {
	my ($el)=&get_atom_name($i);
	if($pdb_positive_charge_ele =~/\_$el\_/) {
	    my $metal=substr($i,12,15);  
	    if( ! grep (/$metal/,@metallist)) { 
		push(@metallist,$metal)
	    }
	}	    
    }


#    my $pdb_nagive_charge_ele="_F_CL_BR_I_";

    return \@metallist;  #return metal entries
}

return 1;
