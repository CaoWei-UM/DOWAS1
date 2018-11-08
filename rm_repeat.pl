#!usr/bin/perl
$project_name=$ARGV[0];
$offspring=$ARGV[1];
$reference=$ARGV[2];
$file1="./$project_name/${offspring}_merge_offtarget_sites.txt";
open(FD,"$file1")||die("Can not open the file!$!n");
open FILE1, ">./$project_name/${offspring}_final_result.txt";

@base=("A","T","C","G");
$i=1;
#$cut=<FD>;
while($line=<FD>){
    @info=split(/\s+/,$line);
    @info1=split(/_/,$info[0]);
    
    #------------------deletion or insertion---------------------------
    if(length($info[1]) > length($info[2])){
	$length=length($info[1])-1;
	$find=substr($info[1],1,$length);
	$position=$length+$info1[1];
    }
    else{
	$length=length($info[2])-1;
	$find=substr($info[2],1,$length);
	$position=$info1[1];
    }
    #----------------get sequence from ref--------------------------------------------------------------
   
	$temp=$position+1;
	$temp1=$position+5;
	$seq_repeat_after=`samtools faidx $reference $info1[0]:$temp-$temp1 | sed '1d' | tr -d '\n' | tr a-z A-Z`;
	$temp=$info1[1]-5;
	$temp1=$info1[1]-1;
	$seq_repeat_before=`samtools faidx $reference $info1[0]:$temp-$temp1 | sed '1d' | tr -d '\n' | tr a-z A-Z`;
    #----------------judge the repeate---------------------------------------------------------------------

    $seq_repeat_after1=$seq_repeat_after;
    $seq_repeat_before1=$seq_repeat_before;
    $find1=$find;
    $j=0;
    while($j<4){
	$seq_repeat_after1=~s/$base[$j]//g;
	$after_length=length($seq_repeat_after1);
	$seq_repeat_before1=~s/$base[$j]//g;
	$before_length=length($seq_repeat_before1);
	$find1=~s/$base[$j]//g;
	$find_length=length($find1);


	if( $after_length == 0 && $find_length == 0 ){
	    last;
	}
	if( $before_length == 0 && $find_length == 0 ){
	    last;
	}
	if( $before_length == 0 && $after_length < 3 ){
	    last;
	}
	if( $after_length == 0 && $before_length < 3 ){
	    last;
	}


	$seq_repeat_after1=$seq_repeat_after;
	$seq_repeat_before1=$seq_repeat_before;
	$find1=$find;
	$j+=1;

    }

    #-------------print results---------------------------------------------
    if($j == 4){
	print FILE1 "$line";
    }

    $i+=1;
}
