use strict;
use warnings;
use File::Spec;
use FindBin;
use Getopt::Long;
use Pod::Usage;
use Math::CDF qw(pbinom);
use Data::Dumper;
use Storable;
use File::Basename;
use Bio::DB::Sam;

# manual and version
our	( $help, $man, $usage, $version );
#: input and output
our $SoftClipSort;
our	$bam;
our	$outdir = "./";
our	$clusters_file = "out.txt";
our $structure_file="struc.txt";
# running parameters
our $max_gap_of_cluster = "10";
our $min_supported_reads = "3";
our $min_supported_PE="3";
our	$bino_pvalue = "0.05";
our $max_cluster_length="5";
our $heter_ratio=0.2;
our $insert_size=300;
our $insert_size_vari=30;
our $bwa="bwa-0.7.8";
our $reference="/PROJ/CR/share/pipeline/cancerPipelineV2/DB/humanRef/hg19_24chr/human_g1k_v37c_decoy.fasta";
our $cap3 = "/PROJ/CR/share/pipeline/cancerPipeline/Bin/var_sv_CREST/cap3";
our $get_sc=0;
our $samtools='samtools';
#get options
my $optionOK = GetOptions(
	# common help parameters
	'h|help|?'		=> \$help,
	# input and output
	'bam:s'			=>	\$bam,
	'outdir:s'		=>	\$outdir,
	'clusters_file:s'		=>	\$clusters_file,
	'structure_file:s'	=>\	$structure_file,
	'reference:s'		=>	\$reference,
	#tools
	'bwa:s'			=>	\$bwa,
	'cap3:s'			=>	\$cap3,
	'samtools'		=>	\$samtools,
	# running parameters
	'insert_size:i'	=>	\$insert_size,
	'insert_size_vari:i'	=>	\$insert_size_vari,
	'max_gap_of_cluster:i'	=>	\$max_gap_of_cluster,
	'min_supported_reads:i'	=>	\$min_supported_reads,
	'min_supported_PE:i'	=>	\$min_supported_PE,
	'heter_ratio:f'	=> \$heter_ratio,
	'get_sc'	=>	\$get_sc,
);

die(usage()) if($help or $version);

my $sam = Bio::DB::Sam->new(-bam  =>"$bam",
                             -fasta=>"$reference",
                             );

mkdir ("$outdir") unless (-d $outdir);
#global variable declare
our $clusters; #used to store cluster results
our $cur_cluster;#current cluster, cur stand for current, the same below.
our $clu_number = 1;#cluster number


##SoftClipSort
our $fh;
#open ($fh,"$samtools view -h $bam 17:19294680-27240265 | ");
open ($fh,"$samtools view -h $bam | ");

$SoftClipSort=basename($bam);
$SoftClipSort=~s/\.bam$/\.sorftclip\.bam/;
$SoftClipSort=File::Spec->catdir($outdir,"$SoftClipSort");
our $sort_bam=$SoftClipSort;
$sort_bam=~s/\.bam$/\.sorted\.bam/;
SoftClipSort($fh,$SoftClipSort) if ($get_sc);

cluster();
store[$clusters],"$outdir/temp.txt";
#my ($asdf)=retrieve "$outdir/temp.txt";$clusters=$asdf->[0];
assemble_checking();
#
PE_checking();
#
print_cluster ();
#open (DUMP,">$outdir/data_dumper.txt");
#print DUMP Dumper($clusters);
#
our $pair={
	"r1[r2["=>"]r2]r1",
	"r1]r2]"=>"r1]r2]",
	"[r2[r1"=>"[r2[r1",
 	"]r2]r1"=>"r1[r2[",
	};

our (%site,%sv_cluster,%sv);
cluster_merge();
###############################################
############# s u b r o u t i n e #############
###############################################


 
sub SoftClipSort{
	my $fh=$_[0];
	my $SoftClipSort=$_[1];
	open OUT," | $samtools view -bS - >$SoftClipSort" or die;
	my $chr="initial";
	CYC:while (<$fh>){
		chomp;
		if($_=~m/^@/){
			print OUT $_,"\n";
			next;
			}
		next if ($_ =~ /GL\d+\.\d/);
		my @info = split /\s+/,;
		next if ($info[5] !~ /S/);
		next if ($info[5]!~/^\d+\w\d+\w$/);
		#next if ($info[11] !~ /^SA/);
		my ($unalign_side)=$info[5]=~/(\d+)S/;
		if ($unalign_side<=5){
			next CYC;
			}
		my $soft_clip_length=0;
		($soft_clip_length)	=	$info[5]	=~	/^(\d+)M\d+S/;
		my $site	=	$info[3];
		$site	+=	$soft_clip_length if ($soft_clip_length);
		$site--;
		my $site_adjust=0;
		{
			$info[5]=~/(\d+)M/;
			my $ad1=$1;
			my $temp=join"\t",@info;
			my $region2="";
			if ($temp=~/\s(SA:Z:.*?\,.*?\,.*?\,.*?\,.*?\,.)/){
				$region2=$1;
				}
			$region2=~/(\d+)S/ if ($region2);
			my $ad2=$1;
			$site_adjust=$ad1-$ad2 if ($ad1>$ad2);
			if ($ad1<$ad2){
				next CYC;
				} 
			$site_adjust*=(-1) if ($info[5]=~/\d+S\d+M/);
			my @temp=$info[5]=~/(\D+)/g;
			if ($#temp>1){
				next CYC;
				}
			}
		$site-=$site_adjust;
		$info[3]=$site;
		my $info=join"\t",@info;
		print OUT "$info\n";
		}
	close (OUT);
	my $sort_bam=$SoftClipSort;
	$sort_bam=~s/\.bam$/\.sorted/;
	system("$samtools sort -@ 2  $SoftClipSort $sort_bam");
	$sort_bam.="\.bam";
	system("$samtools index $sort_bam");
	return ($sort_bam);
	}

sub cluster{
	open IN,"$samtools view $sort_bam |" or die ("Can not open input file\n");
	my $temp=0;
	while ( <IN> ){
		chomp;
		#$temp++;last if ($temp>20000);
		next if ($_ =~ /GL\d+\.\d/);
		my @temp = split /\s+/, $_;
		my @info=@temp[0..5];
		push @info,$temp[3];
		push @info,@temp[6..$#temp];
		next if ($info[5]!~/^\d+\w\d+\w$/);
		my $cluster_flag = if_exist_cur_clusters($cur_cluster);
		if ($cluster_flag == 0){
			$cur_cluster = create_new_cur_cluster(\@info);
			next;
			}
		my $belong_flag = if_belong_cur_cluster ( \@info, $cur_cluster ) ;
		if ( $belong_flag == 0 ) {
			handle_cur_cluster($cur_cluster);
			$cur_cluster = create_new_cur_cluster(\@info);
			} else {
			add_info_to_cur_cluster ( \@info , $cur_cluster );
			}
		}
	handle_cur_cluster($cur_cluster)	if (if_exist_cur_clusters($cur_cluster));
	}

sub if_exist_cur_clusters{
	return 0 if (scalar(keys(%{$_[0]})) == 0);
	return 1 if (scalar(keys(%{$_[0]})) > 0);
	}

sub create_new_cur_cluster{
	my @info	=	@{$_[0]};
	my %cur_cluster;
	$cur_cluster{r1}->{start}	=	"$info[2]\:$info[6]";
	$cur_cluster{r1}->{end}	=	$cur_cluster{r1}->{start};
	push @{$cur_cluster{align}}	,	\@info;
	return \%cur_cluster;
	}

sub if_belong_cur_cluster{
	my @info	=	@{$_[0]};
	my %cur_cluster	=	%{$_[1]};
	my $distance	=	calculate_distance("$info[2]\:$info[6]",$cur_cluster{r1}->{end});
	if ($distance	>	$max_gap_of_cluster){
		return 0
		}
	return 1
	}

sub calculate_distance{
	my @pos1=split/\:/,$_[0];
	my @pos2=split/\:/,$_[1];
	if ($pos1[0] ne $pos2[0]){
		# different chromosome, just return a very large number to make check false
		return 1000
		}
	my $distance = abs($pos1[1] - $pos2[1]);
	return $distance;
	}

sub add_info_to_cur_cluster{
	my @info = @{$_[0]};
	$cur_cluster->{r1}->{end}="$info[2]\:$info[6]";
	push @{$cur_cluster->{align}},\@info;
	}

sub print_cluster{
	#while (my ($key,$value)=each %{$clusters}){
	my $out_put	=	File::Spec->catdir($outdir,$clusters_file);
	open (OUT , ">$out_put") or die ("Can not create out put file\n");
	for my $key(sort{$a<=>$b}keys %{$clusters}){
		my $value=$clusters->{$key};
		if ($value->{better} eq "original"){
			print OUT ">$key\t$value->{chr1}\t$value->{start1}\t$value->{end1}\t$value->{support_reads}\t$value->{local_align}\t$value->{pvalue}\t$value->{chr2}\t$value->{start2}\t$value->{end2}\t$value->{connection_type}\t$value->{su_pe}\n";
			}else{
			my $a=$value->{alternative_connection}->{$value->{better}};
			my $b=$a->{type};# `b` stand for better;
			print OUT ">$key\t$b->[0]\t$b->[1]\t$b->[1]\t$value->{support_reads}\t$value->{local_align}\t$value->{pvalue}\t$b->[2]\t$b->[3]\t$b->[3]\t$b->[4]\t$a->{su_pe}\n";
			print OUT "original\t$key\t$value->{chr1}\t$value->{start1}\t$value->{end1}\t$value->{support_reads}\t$value->{local_align}\t$value->{pvalue}\t$value->{chr2}\t$value->{start2}\t$value->{end2}\t$value->{connection_type}\t$value->{su_pe}\n";
			}
#		print OUT "longest_assemble\t$value->{longest}\n";
		for (@{$value->{align}}){
			print OUT "@{$_}\n";			
			}
		for (@{$value->{align_non_sc}}){
			print OUT "@{$_}\n";
			}
		print OUT "\nPE\n";
		for (@{$value->{PE_align}}){
			#print OUT  "@{$_}\n"
			}
		}
	close(OUT);
	}

sub handle_cur_cluster{
	#check the status of cur_cluster, add cur_cluster to cluster if it is eligible.
	#cluster if global variation.
	my %cur_cluster=%{$_[0]};
	return 0 if (scalar(@{$cur_cluster->{align}})	<	$min_supported_reads);
	#calculate reads depth
	my ($chr,$start)	=	$cur_cluster->{r1}->{start}	=~	/(\w*)\:(\d+)/;
	my ($end)	=	$cur_cluster->{r1}->{end}	=~	/\w*\:(\d+)/;
	my $region2	=	r2_chromosome_divide($cur_cluster->{align});
	my $sub_cur_cluster	=	get_sub_cur_cluster($region2);
	for my $chr_r2(sort{sort_site($a,$b)}keys %{$sub_cur_cluster}){
#	while (my ($chr_r2,$value1)	=	each %{$sub_cur_cluster}){
		my $value1=$sub_cur_cluster->{$chr_r2};
		while (my ($number,$value2)	=	each %{$value1}){
			cluster_filter($value2);# push $value2 info $cluster if it eligible
			}
		}
	}

sub r2_chromosome_divide{
	my @align=@{$_[0]};
	my %region2;
	CYC:for my $info(@align){
		my $temp_join=join"\t",@{$info};
		if ($temp_join=~/\t(SA\:Z.*?)/){
			my ($r2_chr,$r2_site,$chain,$flag)	=	$temp_join	=~	/SA\:Z\:(\w+)\,(\w+)\,(\-|\+),(\w+)\,/;
			my @temp=$flag=~/(\D+)/g;
			if ($#temp>1){
				next CYC;
				}
			my $soft_length=0;
			($soft_length)=$flag=~/^(\d+)M\d+S/;
			$r2_site+=$soft_length if ($soft_length);
			$r2_site--;
			my $site_adjust=0;
			{
				$flag=~/(\d+)M/;
				my $ad1=$1;
				$info->[5]=~/(\d+)S/;
				my $ad2=$1;
				$site_adjust=$ad1-$ad2 if ($ad1>$ad2);
				if ($ad1<$ad2){
					next CYC;
					}
				$site_adjust*=(-1) if ($flag=~/\d+S\d+M/);
				}
			$r2_site-=$site_adjust;
			push @{$region2{$r2_chr}->{$r2_site}}	,	$info;
			}else{
			push @{$region2{unalign}}	,	$info;
			}
		}
	return \%region2;
	}

sub get_sub_cur_cluster{
	my %region2=%{$_[0]};
	my $sub_cur_cluster={};
	while (my ($chr,$value)	=	each %region2){
		next if ($chr=~/unalign/);
		my %site	=	%{$value};
		for my $site(sort{	$a	<=>	$b	} keys %site){
			$sub_cur_cluster	=	check_belonging($sub_cur_cluster,$chr,$site,$site{$site});
			}
		}
	for my $unalign(@{$region2{unalign}}){
		my $soft_clip=get_soft_clip($unalign->[5],$unalign->[10]);
		next if (length($soft_clip)<5);
		while(my ($chr,$value)=each %{$sub_cur_cluster}){
			while (my ($number,$value1)=each %{$value}){
				for (@{$value1->{align_sc}}){
					my $soft_clip_mapped=get_soft_clip($_->[5],$_->[10]);
					my $align_stat=align_stat($soft_clip_mapped,$soft_clip);
					if ($align_stat){
						push @{$sub_cur_cluster->{$chr}->{$number}->{align_non_sc}},$unalign;#non_sc stand for non-soft-clip reads
						last;
						}
					}
				}
			}
		}
	return $sub_cur_cluster;
	}

sub align_stat{
	my $soft_clip1=shift;
	my $soft_clip2=shift;
	#soft_clip1 is short one
	($soft_clip1,$soft_clip2) = ($soft_clip2,$soft_clip1) if (length($soft_clip1) > length($soft_clip2));
	if ($soft_clip2=~/$soft_clip1/){
		return 1;
		}
	return 0;
	}

sub get_soft_clip{
	my $flag=$_[0];
	my $total_reads=$_[1];
	my $sub_reads;
	if ($flag=~/^(\d+)S(\d+)M/){
		$sub_reads=substr($total_reads,0,$1);
		}elsif($flag=~/^(\d+)M(\d+)S/){
		$sub_reads=substr($total_reads,$1,$2);
		}
	return $sub_reads;
	}

sub check_belonging{
	my $sub_cur_cluster=$_[0];
	my $chr=$_[1];
	my $site=$_[2];
	my @align=@{$_[3]};
	for my $align(@align){
		my $belong="false";
		my $connection_type=analysis_flag_soft_clip($align);
		while (my ($key,$value)=each %{$sub_cur_cluster->{$chr}}){
			if(($site - $value->{end} < $max_gap_of_cluster) && ($value->{connection_type} eq $connection_type)){
				$belong="true";
				$sub_cur_cluster->{$chr}->{$key}->{end}=$site;
				push @{$sub_cur_cluster->{$chr}->{$key}->{site}},$site;
				push @{$sub_cur_cluster->{$chr}->{$key}->{align_sc}},$align;
				}
			}
		if ($belong eq "false"){
			my $number = scalar(keys(%{$sub_cur_cluster->{$chr}}));
			push @{$sub_cur_cluster->{$chr}->{$number}->{site}},$site;
			$sub_cur_cluster->{$chr}->{$number}->{start}=$site;
			$sub_cur_cluster->{$chr}->{$number}->{end}=$site;
			$sub_cur_cluster->{$chr}->{$number}->{connection_type}=$connection_type;
			push @{$sub_cur_cluster->{$chr}->{$number}->{align_sc}},$align; #sc stand for soft_clip
			}
		}
	return $sub_cur_cluster;
	}

sub cluster_filter{
	my %sub_cluster	=	%{$_[0]};
	my $pvlaue;
	my $start1	=	0;
	my $end1	=	0;
	my $start2	=	0;
	my $end2	=	0;
	my $chr1;
	my $chr2;
	my $location;
	return 0 if (abs($sub_cluster{end} - $sub_cluster{start} > $max_cluster_length));
	my $support_reads=scalar(@{$sub_cluster{align_sc}});
	if ($sub_cluster{align_non_sc}){
		$support_reads+=scalar(@{$sub_cluster{align_non_sc}})
		};
	return 0 if ( $support_reads < $min_supported_reads);
	my @sites;
	for (@{$sub_cluster{align_sc}}){
		my @info	=	@{$_};
		for (@info){
			$chr1	=	$info[2];
			push @sites,$info[6];
			my $temp=join"\t",@info;
			($chr2,$location) = $temp=~/SA\:Z\:(\w+)\,(\w+)\,/;
			}
		}
	($start1,$end1)=get_start_end(\@sites);
	###bino filter
	my $start_tmp=$start1-2;
	my $end_tmp=$end1+2;
#	my @local_align =`samtools view $bam $chr1:$start_tmp-$end_tmp`;
	my @local_align = $sam->get_features_by_location(-seq_id => "$chr1",
                                                 -start  => "$start_tmp",
                                                 -end    => "$end_tmp");
	my $local_align=scalar(@local_align);
=pre
	for(@local_align){
		my @temp=split/\s+/,;
		$local_align=$temp[2] if ($temp[2]>$local_align);
		}
=cut
	#########CAUTION##########
	###problems can be derived here!!!!!!!
	my $pvalue= pbinom($support_reads , $local_align , 0.5);#0.5 for heterozygosis
	$pvalue=0 unless ($pvalue);
	#return 0 if ( $pvalue < $bino_pvalue);
	@sites=@{$sub_cluster{site}};
	($start2,$end2)=get_start_end(\@sites);
	my $connection_type	=	check_uniformity($sub_cluster{align_sc}); #if uniform get the connection_type, or 0;
	my ($strand1_side,$strand2_side) = get_side($sub_cluster{align_sc}->[0]);
	if ($strand1_side<0){
		$start1++;
		$end1++;
		}
	if ($strand2_side<0){
		$start2++;
		$end2++;
		}
	return 0 if ($connection_type eq "0");
	
	#push this cluster to $cluster;
	my $number=scalar(keys(%{$clusters}));
	print "$number\n";
	$clusters->{$number}->{chr1}		=	$chr1;
	$clusters->{$number}->{start1}	=	$start1;
	$clusters->{$number}->{end1}		=	$end1;
	$clusters->{$number}->{length1}	=	$end1	-	$start1;
	
	$clusters->{$number}->{chr2}		=	$chr2;
	$clusters->{$number}->{start2}	=	$start2;
	$clusters->{$number}->{end2}		=	$end2;
	$clusters->{$number}->{length2}	=	$end2	-	$start2;
	
	$clusters->{$number}->{support_reads}	=	$support_reads;
	$clusters->{$number}->{local_align}	=	$local_align;
	$clusters->{$number}->{pvalue}	=	sprintf("%0.2f",$pvalue);
	
	$clusters->{$number}->{connection_type}	=	$connection_type;
	$clusters->{$number}->{align}	=	$sub_cluster{align_sc};
	$clusters->{$number}->{align_non_sc}	=	$sub_cluster{align_non_sc};
	return 1;
	}

sub get_start_end{
	my @sites=@{$_[0]};
	my %sites;
	for (@sites){
		$sites{$_}++;
		}
	my @keys=sort{$sites{$b}<=>$sites{$a}}keys %sites;
	if ($#keys==0){
		return ($keys[0],$keys[0]);
		}
	if ($sites{$keys[0]}/$sites{$keys[1]}>=3){
		return ($keys[0],$keys[0]);
		}else{
		my @res=sort{$a<=>$b}($keys[0],$keys[1]);
		return @res;
		}
	}

sub get_side{
	my @info=@{$_[0]};
	my $strand1_side=1;
	my $strand2_side=-1;
	if ($info[5]=~/\d+(S|H)\d+M/){
		$strand1_side*=(-1);
		$strand2_side*=(-1);
		}
	
	my $temp=join"\t",@info;
	my ($region2)=$temp=~/\s(SA:Z:.*?\,.*?\,.*?\,.*?\,.*?\,.)/;

	if (($info[1]&0x10) and ($region2=~/\+/)){
		$strand2_side*=(-1);
		}
	
	if (!($info[1]&0x10) and ($region2=~/\-/)){
		$strand2_side*=(-1);
		}
	return ($strand1_side,$strand2_side);
	}

sub check_uniformity{
#return connection type if uniform, or return 0;
	my @align=@{$_[0]};
	my $connection_type_total="none";
	for (@align){
		my $connection_type=analysis_flag_soft_clip($_);
		if ($connection_type_total eq "none"){
			$connection_type_total = $connection_type;
			}else{
			return 0 if ($connection_type ne $connection_type_total);
			}
		}
	return $connection_type_total;
	}

sub analysis_flag_soft_clip{
	my @info=@{$_[0]};
	my $flag=1;
	my $sign="[";
	my $strand1="r1";
	my $strand2="r2";
	my @ALT=(\$strand1,\$sign,\$strand2,\$sign);
	if ($info[5]=~/\d+(S|H)\d+M/){
		@ALT = reverse (@ALT);
		$flag*=(-1);
		}

	my $temp=join"\t",@info;
	my ($region2)=$temp=~/\s(SA:Z:.*?\,.*?\,.*?\,.*?\,.*?\,.)/;
		
	if (($info[1]&0x10) and ($region2=~/\+/)){
		$flag*=(-1);
		}
	
	if (!($info[1]&0x10) and ($region2=~/\-/)){
		$flag*=(-1);
		}
	$sign="]" if ($flag==-1);
	my $ALT="";
	for (@ALT){
        $ALT.="${$_}";
		}
	return $ALT;
	}

sub add_seq{
	my $key=$_[0];
	my $info=$_[1];
	for (@{$info}){
		push @{$clusters->{$key}->{sclip}},[$_->[0],$_->[10],$_->[11]];
		}
	}

sub assemble_checking {
	my $cap3_options=" -h 70 -y 10 > /dev/null"; 
	my %sclip;
	my $flag=0;
	my $index;
	my $assemble_file=File::Spec->catdir($outdir,"assemble.fa");
#=asdf
	open (ASM,">$assemble_file");
	while(my ($key,$value1)=each %{$clusters}){
		add_seq($key,$value1->{align});
		add_seq($key,$value1->{align_non_sc});
		my $fa_file=prepare_reads($clusters->{$key}->{sclip});
		system(join(" ", ($cap3, $fa_file, $cap3_options)));
		my $longest_reads = get_longest_reads("$fa_file.cap.contigs");
		if (!$longest_reads){
			delete $clusters->{$key};
			next;
			};
		$clusters->{$key}->{longest}=$longest_reads;
		print ASM ">$key\n$longest_reads\n";
		}
	close (ASM);
#=cut
	my $bwa_out=File::Spec->catdir($outdir,"bwa1.out");
	system("$bwa mem -t 5 -a $reference $assemble_file> $bwa_out");
	my $fh; 
	open ($fh,"$bwa_out");
	analysis_bwa($fh);
	}

sub analysis_bwa{
	my $fh=$_[0];
	my $region2="";
	my $chr1="";
	my $start1="";
	my $chr2="";
	my $start2="";
	CYC:while (<$fh>){
		next if ($_=~/^\@/);
		my @temp=split/\s+/,;
		next if ($_ =~ /GL\d+\.\d/);
		next if (mismatch($temp[11])>1);
		next if ($temp[5]!~/^\d+\w\d+\w$/);
		my @info;
		my $temp=join"\t",@temp;
		if($temp=~/\s(SA:Z:.*?\,.*?\,.*?\,.*?\,.*?\,.)/){
			$region2=$1;
			}
		push @info,@temp[0..5];
		my $soft_clip_length=0;
		($soft_clip_length)	=	$info[5]	=~	/^(\d+)M\d+(S|H)$/;
		my @flag_number=$info[5]=~/(\D+)/g;
		if ($#flag_number>1){
			next CYC;
			}
		$info[6]=$temp[3]-1;
		$info[6]+=$soft_clip_length if ($soft_clip_length);
		push @info,@temp[6..11];
		$info[12]=$region2 if ($region2);
		$chr1=$info[2];
		$start1=$info[6];
		my ($chain,$flag);
		($chr2,$start2,$chain,$flag)=$info[12]=~/SA:Z:(.*?),(\d+),(.*?),(.*?),.*/;
		@flag_number=$flag=~/(\D+)/g;
		if ($#flag_number>1){
			next CYC;
			}
		my $site_adjust=0;
		{
			$info[5]=~/(\d+)M/;
			my $ad1=$1;
			$info[12]=~/(\d+)(S|H)/;
			my $ad2=$1;
			$site_adjust=$ad1-$ad2 if ($ad1>$ad2);
			if ($ad1<$ad2){
				next CYC;
				}
			$site_adjust*=(-1) if ($info[5]=~/\d+(S|H)\d+M/);
			}
		$start1-=$site_adjust;
		$soft_clip_length=0;
		($soft_clip_length)	=	$info[12]	=~	/\,(\d+)M\d+(S|H)\,/;
		$start2+=$soft_clip_length if ($soft_clip_length);
		$start2--;
		$site_adjust=0;
		{
			$info[12]=~/(\d+)M/;
			my $ad1=$1;
			$info[5]=~/(\d+)(S|H)/;
			my $ad2=$1;
			$site_adjust=$ad1-$ad2 if ($ad1>$ad2);
			$site_adjust*=(-1) if ($info[12]=~/\d+(S|H)\d+M/);
			}
		$start2-=$site_adjust;
		my $connection_type=analysis_flag_soft_clip(\@info);
		my ($strand1_side,$strand2_side) = get_side(\@info);
		if ($strand1_side<0){
			$start1++;
			}
		if ($strand2_side<0){
			$start2++;
			}
		my $alt_number=0;
		$alt_number=scalar(keys(%{$clusters->{$info[0]}->{alternative_connection}})) if ($clusters->{$info[0]}->{alternative_connection});
		$clusters->{$info[0]}->{alternative_connection}->{$alt_number}->{type}=[$chr1,$start1,$chr2,$start2,$connection_type];
		}
	}
	
sub mismatch{
	my @info=split/\:/,$_[0];
	return $info[2];
	}
	
sub prepare_reads{
	my@inf=@{$_[0]};
	my $fa_file=File::Spec->catdir($outdir,"temp.fa");
	open (my $fh,">$fa_file");
	my$qual_file=$fa_file.".qual";
	open my$qual_fh, ">$qual_file" || die $!;
 	foreach my$i(@inf){
		print $fh ">",$i->[0],"\n";
		print $fh $i->[1],"\n";
		print $qual_fh ">",$i->[0],"\n";
		my@qual=map {ord} split //,$i->[2];
		print $qual_fh join(" ",@qual[0..$#qual]),"\n";
		}
	close $qual_fh;
	return $fa_file;
	}

sub get_longest_reads {
	my $file = shift;
	my $reads;
	open (INP,"$file");
	<INP>;
	my $flag=1;
	while (<INP>){
		chomp;
		if ($_=~/^>/){
			$flag=0;
			}
		last if ($flag==0);
		$reads.=$_;
		}
	return $reads;
	}

sub PE_checking{
	for my $key(sort{$a<=>$b}keys %{$clusters}){
		my $max;
		my $value=$clusters->{$key};
		my $start1_tmp=$value->{start1}	-	$insert_size;
		my $end1_tmp=$value->{end1}	+	$insert_size;
		my $site_original_1="$value->{chr1}\:$value->{start1}";
		my $site_original_2="$value->{chr2}\:$value->{start2}";
		my @first_point = $sam->get_features_by_location(-seq_id => "$value->{chr1}",
                                                 -start  => $start1_tmp,
                                                 -end    => $end1_tmp);
		#my @first_point=`samtools view $bam $value->{chr1}:$start1_tmp-$end1_tmp`;
		my $first	=	PE(\@first_point,$value->{chr2},$value->{start2},$value->{end2},$value->{connection_type});
		$first->{true} =0 if (!$first->{true});
		$clusters->{$key}->{su_pe}=$first->{true};#su_pe stand for support Pair-End;
		$clusters->{$key}->{PE_align}=$first->{align};
		$clusters->{$key}->{better}="original";
		$max=$clusters->{$key}->{su_pe};
		my $fuzzy_number=scalar(keys(%{$value->{alternative_connection}}));
		#next if ($fuzzy_number==2);
		my $cal_number=0;
		while (my ($number,$value2)=each %{$value->{alternative_connection}}){
			$cal_number++;
			#last if ($cal_number>200);
			my $site1="$value2->{type}->[0]\:$value2->{type}->[1]";
			my $site2="$value2->{type}->[2]\:$value2->{type}->[3]";
			next if (calculate_distance($site1,$site_original_1)+calculate_distance($site2,$site_original_2) < $insert_size);
			next if (calculate_distance($site1,$site_original_2)+calculate_distance($site2,$site_original_1) < $insert_size);
			my $start1_tmp=$value2->{type}->[1]-$insert_size;
			my $end1_tmp=$value2->{type}->[1]+$insert_size;
			my @alt_first_point = $sam->get_features_by_location(-seq_id => "$value2->{type}->[0]",
                                                 -start  => $start1_tmp,
                                                 -end    => $end1_tmp);
			#my @alt_first_point=`samtools view $bam $value2->{type}->[0]:$start1_tmp-$end1_tmp `;
			my $first	=	PE(\@alt_first_point,$value2->{type}->[2],$value2->{type}->[3],$value2->{type}->[3],$value2->{type}->[4]);
			$first->{true} =0 if (!$first->{true});
			$value2->{su_pe}=$first->{true};#su_pe stand for support Pair-End;
			$value2->{PE_align}=$first->{align};
			if (choose_better($value2->{su_pe},$max)){
				$max=$value2->{su_pe};
				$clusters->{$key}->{better}=$number;
				}
			last if ($max>20);
			}
		}
	}

sub choose_better{
	my $select=$_[0];
	my $max=$_[1];
	if ($select==0){
		return 0;
		}
	if ($max==0&&$select>$max){
		return 1;
		}
	if ($max>=20){
		return 0;
		}
	if ($select/$max>3){
		return 1;
		}
	if ($select>$max){
		return 1;
		}
	return 0;
	}

sub PE{
	my $pe_reads=$_[0];
	my $chr2=$_[1];
	my $start2=$_[2];
	my $end2=$_[3];
	my $connection_type=$_[4];
	my $all=0;
	my $inter=0;
	my $second=0;
	my $eligible=0;
	my $first={};
	for my $align(@{$pe_reads}){
		my $chr_pe;
		$chr_pe=$align->mate_seq_id;
		next if (!$chr_pe);
		next if ($chr_pe!~/\w/);
		my $start=$align->mate_start;
		my $distance1=calculate_distance("$chr2\:$start2","$chr_pe\:$start");
		my $distance2=calculate_distance("$chr2\:$end2","$chr_pe\:$start");
		next if ($distance1 > $insert_size+3*$insert_size_vari);##########CAUTION##########
		next if ($distance2 > $insert_size+3*$insert_size_vari);
		my $flag=analysis_flag_PE($align->flag);
		next if ($flag ne $connection_type);
		$first->{true}++;
		push @{$first->{align}},$align;
#		print "$_\n";
		}
	return $first;
	}

sub analysis_flag_PE{
	my $symbol=$_[0];
	my $strand1="r1";
	my $strand2="r2";
	my $flag=1;
	my $sign="[";
	my @ALT=(\$strand1,\$sign,\$strand2,\$sign);
	if((($symbol&0x10)/0x10)!=(($symbol&0x20)/0x20)){
		if ($symbol&0x10){
			@ALT=reverse(@ALT);
			$flag*=(-1);
			}
		}
	if((($symbol&0x10)/0x10)==(($symbol&0x20)/0x20)){
		$flag*=(-1);
		if ($symbol&0x10){
			@ALT=reverse(@ALT);
			$flag*=(-1);
			}
		}
	$sign="]" if ($flag<0);
	my $ALT="";
	for (@ALT){
        $ALT.="${$_}";
		}
	return $ALT;	
	}

sub cluster_merge{
	my $cluster_files	=	File::Spec->catdir($outdir,$clusters_file);
	open (INP,$cluster_files);
	my %cluster;
	while (<INP>){
		chomp($_);
		next if ($_!~/^>/);
		my @temp=split/\s+/,;
		my $Untrusted=get_Untrusted(\@temp);
		next if ($Untrusted==0);
		$cluster{$temp[0]}->{chr1}=$temp[1];
		$cluster{$temp[0]}->{start1}=$temp[2];
		$cluster{$temp[0]}->{end1}=$temp[3];
		$cluster{$temp[0]}->{support_reads}=$temp[4];
		$cluster{$temp[0]}->{local_align}=$temp[5];
		$cluster{$temp[0]}->{pvalue}=$temp[6];
		$cluster{$temp[0]}->{chr2}=$temp[7];
		$cluster{$temp[0]}->{start2}=$temp[8];
		$cluster{$temp[0]}->{end2}=$temp[9];
		$cluster{$temp[0]}->{connection_type}=$temp[10];
		$cluster{$temp[0]}->{su_pe}=$temp[11];
		}
	while (my ($key,$value)=each %cluster){
		my $site1="$value->{chr1}\:$value->{start1}";
		my $site2="$value->{chr2}\:$value->{start2}";
		my $connection_type=$value->{connection_type};
		$site{$site1}->{$connection_type}->{$site2}->{info}=[$value->{support_reads},$value->{local_align},$value->{pvalue},$value->{su_pe}];
		$site{$site1}->{$connection_type}->{$site2}->{true}=1;
		if (!$site{$site2}->{$pair->{$connection_type}}->{$site1}->{true}){
			$site{$site2}->{$pair->{$connection_type}}->{$site1}->{info}=[$value->{support_reads},$value->{local_align},$value->{pvalue},$value->{su_pe}];
			$site{$site2}->{$pair->{$connection_type}}->{$site1}->{true}=0;
			}
		}

	my @sort_site=sort{sort_site($a,$b)}keys %site;
	for my $site1(@sort_site){
		while (my ($connection_type,$value1)=each %{$site{$site1}}){
			while (my ($site2,$value2)=each %{$value1}){
				next if ($value2->{info}->[4]);
				my $key=scalar(keys(%sv_cluster));
				push @{$sv_cluster{$key}->{info}},[$site1,$connection_type,$site2,@{$value2->{info}}];
				$sv_cluster{$key}->{number}++;
				$value2->{info}->[4]++;
				search_others($site1,$key);
				search_others($site2,$key);
				estimate_sv_structrue($key);
				}
			}
		}
	my $structure_file=File::Spec->catdir($outdir,$structure_file);
	open (OUT,">$structure_file");
	for my $key(sort{$a<=>$b}keys(%sv)){
		my $value1=$sv{$key};
		my $info;
		my $number=0;
		while (my ($site1,$value2)=each (%{$value1->{info}})){
			while (my ($type,$site2)=each (%{$value2})){
				if ($site{$site1}->{$type}->{$site2}->{true}){
					$number++;
					$info.="$site1\t$type\t$site2\t";
					my @temp=@{$site{$site1}->{$type}->{$site2}->{info}};
					for (my $n=0;$n<=3;$n++){
						$info.="$temp[$n]\t";
						}
					$info.="\n";
					}
				my $reverse_type=$pair->{$type};
				if ($site{$site2}->{$reverse_type}->{$site1}->{true}){
					$number++;
					$info.="$site2\t$reverse_type\t$site1\t";
					my @temp=@{$site{$site2}->{$reverse_type}->{$site1}->{info}};
					for (my $n=0;$n<=3;$n++){
						$info.="$temp[$n]\t";
						}
					$info.="\n";
					}
				}
			}
		print OUT ">$key\t$number\t$value1->{type}";
		if ($value1->{type}!~/unknown/){
			print OUT "\t$value1->{site}\t";
			for (@{$value1->{source}}){
				print OUT "$_\t";
				}
			print OUT "$value1->{connection_type}";
			}
		print OUT "\n";
		print OUT "$info";
		}
	close (OUT);
	}

sub get_Untrusted{
	my $temp=$_[0];
	return 0 if ($temp->[11])<7;
	return 0 if ($temp->[4])<2;
	return 0 if ($temp->[5])<30;
	return 0 if ($temp->[2]!=$temp->[3]);
	return 0 if ($temp->[8]!=$temp->[9]);
	return 1;
	}

sub search_others{
	my $site=$_[0];
	my $key=$_[1];
	my @temp=split/\:/,$site;
	for (my $n=-$max_gap_of_cluster;$n<=$max_gap_of_cluster;$n++){
		my $temp=$temp[1]+$n;
		my $site1="$temp[0]\:$temp";
		if ($site{$site1}){
			my @asdf=keys(%{$site{$site1}});
			while (my ($connection_type,$value1)=each %{$site{$site1}}){
				while (my ($site2,$value2)=each %{$value1}){
					next if ($value2->{info}->[4]);
					push @{$sv_cluster{$key}->{info}},[$site1,$connection_type,$site2,@{$value2->{info}}];
					$sv_cluster{$key}->{number}++;
					# delete ($site{$site1}->{$connection_type}->{$site2});
					$value2->{info}->[4]++;
					search_others($site2,$key);
					}
				}
			}
		}
	}

sub estimate_sv_structrue{
	my $sv_key=$_[0];
	my %type;
	for (@{$sv_cluster{$sv_key}->{info}}){
		$type{$_->[0]}->{$_->[2]}=$_->[1];
		my $reverse_conn=$pair->{$_->[1]};
		$type{$_->[2]}->{$_->[0]}=$reverse_conn;
		}
	my @sites=sort{sort_site($a,$b)}(keys %type);
	my $seeds=search_seeds(\@sites,\%type);
	if ($#{$seeds->{single}}==-1 and scalar(keys(%type))==2){ #potential deletion or duplication
		my $site1=$sites[0];
		my $site2=$sites[1];
		my $same_chrom=same_chrom($site1,$site2);
		my $flag=1;
		if (!$same_chrom){
			return 0;
			}
		unless ($type{$site1}->{$site2} && $type{$site2}->{$site1}){
			return 0;
			}
		if ($type{$site1}->{$site2} ne $pair->{$type{$site2}->{$site1}}){
			return 0;
			}
		if($type{$site1}->{$site2}=~/(\^\[\w+\[\w+$)|(\w+\]\w+\]$)/){
			return 0;
			}
		my $key=scalar(keys(%sv));
		$sv{$key}->{info}->{$site1}->{$type{$site1}->{$site2}}=$site2;
		$sv{$key}->{type}=judge_del_dup($site1,$site2,$type{$site1}->{$site2});
		my @temp=sort{sort_site($a,$b)}($site1,$site2);
		$sv{$key}->{site}=$temp[0];
		$temp[2]=calculate_distance($temp[0],$temp[1]);
		$sv{$key}->{source}=\@temp;
		$sv{$key}->{connection_type}="+/+";
		return 1;
		}
	my @dirty;
	return 0 if ($#{$seeds->{single}}<0);
	my $type_cis=seed_extend_process(\%type);
	my $type_choice=type_choice($type_cis,$seeds);
	for my $seed(@{$seeds->{single}}){
		my $key=scalar(keys(%sv));
		my $site1=$seed->[0];
		my $site2=$seed->[1];
		next if (grep/$site1/,@dirty);
		next if (grep/$site2/,@dirty);
		push @dirty,($seed->[0],$seed->[1]);
		my $flag=judge_flag($type_choice,$seed,$seeds,\%type);
		next if (!$flag);
		my @site3=keys($type{$site1});
		my @site4=keys($type{$site2});
		my ($site3,$site4)=get_site3_4(\@site3,\@site4);
		if ($site3=~/false/){
			next;
			}
		if (calculate_distance($site3,$site4)<$max_gap_of_cluster){#potential inversion
			if ($type{$site1}->{$site3} ne $pair->{$type{$site3}->{$site1}}){
				next
				}
			if ($type{$site2}->{$site4} ne $pair->{$type{$site4}->{$site2}}){
				next;
				}
			if ($type{$site1}->{$site3}=~/\^\]\w+\]\w+$|\w+\[\w+\[$/){
				next;
				}
			if ($type{$site2}->{$site4}=~/\^\]\w+\]\w+$|\w+\[\w+\[$/){
				next;
				}
			push @dirty,($site3,$site4);
			$sv{$key}->{info}->{$site1}->{$type{$site1}->{$site3}}=$site3;
			$sv{$key}->{info}->{$site2}->{$type{$site2}->{$site4}}=$site4;
			$sv{$key}->{type}="inversion";
			$sv{$key}->{site}=$site1;
			my @temp=sort{sort_site($a,$b)}($site1,$site3);
			$temp[2]=calculate_distance($temp[0],$temp[1]);
			$sv{$key}->{source}=\@temp;
			$sv{$key}->{connection_type}="+/-";
			next;
			}
		if ($type{$site1}->{$site3} ne $pair->{$type{$site2}->{$site4}}){
			next;
			}
		$sv{$key}->{info}->{$site1}->{$type{$site1}->{$site3}}=$site3;
		$sv{$key}->{info}->{$site2}->{$type{$site2}->{$site4}}=$site4;
		$sv{$key}->{type}="translocation_copy";
		$flag=copy_cut($site3,$site4,$type_choice);#copy or cut; 1 for cut, 0 for copy
		if ($flag){
			$sv{$key}->{type}="translocation_cut";
			$sv{$key}->{info}->{$type_choice->[0]}->{$type_choice->[2]}=$type_choice->[1];
			}
		$sv{$key}->{site}=$site1;

		my @temp=sort{sort_site($a,$b)}($site3,$site4);
		$temp[2]=calculate_distance($temp[0],$temp[1]);
		$sv{$key}->{source}=\@temp;
		$sv{$key}->{connection_type}="+/+";
		$sv{$key}->{connection_type}="+/-" if ($type{$site1}->{$site3}=~/\^\[\w+\[\w+$|\w+\]\w+\]$/);
		}
	return 0;
	}

sub get_site3_4{
	my @site3=@{$_[0]};
	my @site4=@{$_[1]};
	my $length=0;
	my @choice=("non","non");
	for (my $n=0;$n<=$#site3;$n++){
		for(my $m=0;$m<=$#site4;$m++){
			next if (!same_chrom($site3[$n],$site4[$m]));
			if (!$length){
				$length=calculate_distance($site3[$n],$site4[$m]);
				@choice=($n,$m);
				next;
				}
			my $temp=calculate_distance($site3[$n],$site4[$m]);
			if ($temp<$length){
				$length=$temp;
				@choice=($n,$m);
				}
			}
		}
	if ($choice[0]=~/non/){
		return "false";
		}
	return ($site3[$choice[0]],$site4[$choice[1]]);
	}

sub copy_cut{
	my ($site3,$site4,$type_choice)=@_;
	if(calculate_distance($site3,$type_choice->[0])<=$max_cluster_length+5 && calculate_distance($site4,$type_choice->[1])<=$max_cluster_length+5){
		return 1;
		}
	if(calculate_distance($site4,$type_choice->[0])<=$max_cluster_length+5 && calculate_distance($site3,$type_choice->[1])<=$max_cluster_length+5){
		return 1;
		}
	return 0;
	}

sub judge_flag{
	my ($type_choice,$seed_refind,$seeds,$type)=@_;
	my $flag1=0;
	if (!$type_choice->[0]){
		return 1;
		}
	for my $type_choice_s(@{$type_choice}){
		next if ($type_choice_s=~/\[|\]/);
		for (@{$seed_refind}){	
			if(calculate_distance($_,$type_choice_s)<$max_cluster_length+5){
				$flag1=1
				}
			}
		}
	my $flag2=0;
	for my  $type_choice_s(@{$type_choice}){
		next if ($type_choice_s=~/\[|\]/);
		my ($seed_p)=keys %{$type->{$type_choice_s}};
		if ($seeds->{hash}->{$seed_p}){
			$flag2=1;
			}
		}
	if ($flag1*$flag2>0){
		return 0;
		}else{
		return 1
		}
	}

sub seed_extend_process{
	my $seed_extend=$_[0];
	my $seed_result={};
	while (my ($start,$value)=each %{$seed_extend}){
		while (my ($end,$type)=each %{$value}){
			if ($type=~/(\[r2\[r1)|(\]r2\]r1)/){
				next;
				}else{
				$seed_result->{$start}->{$end}=$type;
				}
			}
		}
	return $seed_result;
	}

sub type_choice{
	my $seed_extend=$_[0];
	my $seeds=$_[1];
	my $type_choice=[0,0,0];
	while (my ($start,$value)=each %{$seed_extend}){
		while (my ($stop,$type)=each %{$value}){
			my $flag1=0;
			my $flag2=0;
			for (@{$seeds->{all}}){
				if(calculate_distance($_,$start)<$max_cluster_length){
					$flag1=1
					}
				if(calculate_distance($_,$stop)<$max_cluster_length){
					$flag2=1
					}
				}
			next if ($type!~/r1\[r2\[/);
			next unless ($flag1*$flag2);
			my $flag=sort_site($start,$stop);
			next if ($flag>0);
			$type_choice=mk_choice($type_choice,$start,$stop,$type);
			}
		}
	return $type_choice;
	}

sub mk_choice{
	my $type_choice=$_[0];
	my $start=$_[1];
	my $stop=$_[2];
	my $type=$_[3];
	my @start=split/\:/,$start;
	my @stop=split/\:/,$stop;
	if ($start[0] ne $stop[0]){
		return $type_choice;
		}
	if (!$type_choice->[0]){
		return [$start,$stop,$type];
		}
	my @start1=split/\:/,$type_choice->[0];
	my @stop1=split/\:/,$type_choice->[1];
	my $len1=$start1[$#start1]-$stop1[$#stop1];
	$len1=abs($len1);
	my $len2=$start[1]-$stop[1];
	$len2=abs($len2);
	if ($len1>=$len2){
		return [$start,$stop,$type];
		}else{
		return $type_choice;
		}
	}

sub judge_del_dup{
	my ($r1,$r2,$type)=@_;
	if($type!~/r1\[r2\[/){
		($r1,$r2)=($r2,$r1);
		}
	my @r1=split/\:/,$r1;
	my @r2=split/\:/,$r2;
	my $sv="deletion";
	$sv="duplication" if ($r1[1]>$r2[1]);
	return $sv;
	}

sub search_seeds{
	my @sites=@{$_[0]};
	my $type=$_[1];
	my %seed;
	for (my $n=0;$n<$#sites;$n++){
		my $distance=calculate_distance($sites[$n],$sites[$n+1]);
		my $flag=seed_check($sites[$n],$sites[$n+1],$type);
		next if (!$flag);
		next if ($distance>$max_gap_of_cluster);
		push @{$seed{single}},[$sites[$n],$sites[$n+1]];
		push @{$seed{all}},($sites[$n],$sites[$n+1]);
		$seed{hash}->{$sites[$n]}=$sites[$n+1];
		$seed{hash}->{$sites[$n+1]}=$sites[$n];
		$n++;
		}
	return \%seed;
	}

sub seed_check{
	my $seed1=$_[0];
	my $seed2=$_[1];
	my %type=%{$_[2]};
	my $flag=0;
	my @type1;
	my @type2;
	while (my ($end,$type)=each %{$type{$seed1}}){
		push @type1,$type;
		}
	while (my ($end,$type)=each %{$type{$seed2}}){
		push @type2,$type;
		}
	for my $type1(@type1){
		for my $type2 (@type2){
			my ($site1)=$type1=~/(r\d)/;
			my ($site2)=$type2=~/(r\d)/;
			next if ($site1 eq $site2);
			$flag=1;
			return $flag;
			}
		}
	return 0;
	}

sub same_chrom{
	my @pos1=split/\:/,$_[0];
	my @pos2=split/\:/,$_[1];
	if ($pos1[0] ne $pos2[0]){
		return 0
		}
	return 1;
	}

sub sort_site{
	my @a=split/\:/,$_[0];
	my @b=split/\:/,$_[1];
	return 1 if ($a[0] gt $b[0]);
	return -1 if ($a[0] lt $b[0]);
	if($a[0] eq $b[0]){
		return 1 if ($a[1] > $b[1]);
		return -1 if ($a[1] < $b[1]);
		return 0;
		}
	}

sub print_log{
	my ($log_message, $log_flag) = @_;
	chomp(my $time=`date`);
	print STDERR "$time:\t$log_message\n";
	if ($log_flag=~/E/){
		die ("Sorry program has been suicide...\n");
		}
	}
	
sub usage{
		# common help parameters
my $help=<<""
	
	Version 1.0
	
	'h|help|?'		=> print thist message
	# input and output
	'SoftClipSort:s'	=>	soft clip file sort by site
	'bam:s'			=>	input bam file
	'outdir:s'		=>	output directory
	'clusters_file:s'		=>	output cluster file
	'structure_file'	=>	output structure file
	'reference'		=>	reference in fasta format
	#tools
	'bwa'			=>	path to bwa
	'cap3'			=>	path to cap3
	# running parameters
	'insert_size'	=>	estimated insert size
	'insert_size_vari'	=>	estimated insert size variance
	'max_gap_of_cluster:i'	=>	maximum gap reads in an cluster
	'min_supported_reads:i'	=>	minimum supported reads
	'min_supported_PE'	=>	minimum supported PE
	'heter_ratio:f'	=> heterozygosis ratio

	}

