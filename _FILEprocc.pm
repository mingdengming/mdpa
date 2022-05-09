package _FILEprocc;


sub lineoffile{
    my ($dir,$file)=@_;
    my $nline=0;
    my $tmpf=$file.'nline_tmp0000';

    return -9999 if(! -e  "$dir/$file") ;
    `wc $dir/$file > $dir/$tmpf`;
    open(otmpf,"< $dir/$tmpf");
    while(<otmpf>){
	my @tmp=split /\s+/,$_; 
	$nline=@tmp[0]; $nline=@tmp[1] if $tmp[0] eq '';}
    close(tmpf);
    `rm -f $dir/$tmpf`;

    return $nline;
}

1;
