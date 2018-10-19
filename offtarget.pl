#!usr/bin/perl                                          
$project_name=$ARGV[0];
$offspring=$ARGV[1];
$guide_RNA_seq=$ARGV[2];  
$ref=$ARGV[3];                              

$PAM_GG=13;
$PAM_AG=15;


$on=$guide_RNA_seq;
$on=reverse $on;
$re_on=$on;
$re_on=reverse $on;
$re_on=~s/A/M/g;
$re_on=~s/T/A/g;
$re_on=~s/M/T/g;
$re_on=~s/C/N/g;
$re_on=~s/G/C/g;
$re_on=~s/N/G/g;
open(FD,"./$project_name/father/father.vcf")||die("Can not open the file!$!n");
open(FD1,"./$project_name/mother/mother.vcf")||die("Can not open the file!$!n");
open(FD2,"./$project_name/$offspring/$offspring.vcf")||die("Can not open the file!$!n");

#-----------------input father------------------------------------------------


while($line=<FD>){
    @info=split(/\s+/,$line);
    $key=$info[0]."_".$info[1]."_".$info[3]."_".$info[4];
    $storefather{$key}=1;
}
#----------------------input mother-------------------------------------------------

while($line1=<FD1>){
    @info1=split(/\s+/,$line1);
    $key=$info1[0]."_".$info1[1]."_".$info1[3]."_".$info1[4];
    $storemother{$key}=1;
      
}


#-------------------input offspring-------------------------------------------------------------------------------------------------------------------------
mkdir("./$project_name/$offspring") unless(-d "./$project_name/$offspring");
open FILE1, ">./$project_name/$offspring/4_node_GG.txt";
open FILE2, ">./$project_name/$offspring/4_noin_GG.txt";
open FILE3, ">./$project_name/$offspring/4_node_AG.txt";
open FILE4, ">./$project_name/$offspring/4_noin_AG.txt";
print FILE1 "site\tref\tALT\tfather\tmother\n";
print FILE2 "site\tref\tALT\tfather\tmother\n";
print FILE3 "site\tref\tALT\tfather\tmother\n";
print FILE4 "site\tref\tALT\tfather\tmother\n";


while($line2=<FD2>){
    @info2=split(/\s+/,$line2);  
    $site=$info2[0]."_".$info2[1]."_".$info2[3]."_".$info2[4];
    if(exists $storefather{$site} || exists $storemother{$site}){
    }
    else{      
#	$count_num+=1;


	
	$loce = "";
	$loce1 = "";
	$loc = "";
	$loc1 = "";
	$depth="";
	$ratio="";
	$define_sum="";
	@infor2=split(/:/,$info2[9]);           

	$refl=length($info2[3]);
	$ALTl=length($info2[4]);
    

#	print "2222222\t$infor2[2]\n";
 #   if($infor2[3] > 5 && $depth > 2 && $ratio > 0.2){
#--------------------------------------insertion-------------------------------------------------------------------------------  
	if($infor2[2] > 8 && $ALTl> $refl){
	    #$filename=uc($info2[0]);
	    #open (FILE,"$ref_path/$filename")||die "cannot open word_file: $!";   
       	#    $seq=<FILE>;
	    $seq_ss=$info2[1]-23;
	    #$seqs=substr($seq,$seq_ss,23);
	    #$seqe=substr($seq,$info2[1],23);
		$temp=$seq_ss+1;
		$temp1=$temp+22;
		$seqs=`samtools faidx $ref $info2[0]:$temp-$temp1 | sed '1d' | tr -d '\n' |tr a-z A-Z `;
		$temp=$info2[1]+1;
		$temp1=$temp+22;
		$seqe=`samtools faidx $ref $info2[0]:$temp-$temp1 | sed '1d' | tr -d '\n' |tr a-z A-Z `;
	    if($seqs=~m/CC/ || $seqs=~m/CT/){
		$loc=0;
		$loc1=0;
		$number_C=0;
		$number_C1=0;
		@number_CC="";
		@number_CC1="";
		@change_C="";
		$seqs1=$seqs;
		while($loc>-1){
		    $loc = index($seqs,"CC");
		    @change_C=split(//,$seqs);
		    $change_C[$loc]="O";
		    $seqs=join "",@change_C;
		    #		    print "aaaa $number_C\t $loc\n";
		    $number_CC[$number_C]=$loc;
		    $number_C+=1;
		}
		$number_C1=0;
		while($loc1>-1){
		    $loc1 = index($seqs1,"CT");
		    $seqs1=~s/CT/OO/;
                    #    print "aaaa $number_C1\t $loc1\n";
		    $number_CC1[$number_C1]=$loc1;
		    $number_C1+=1;
		}
		                
		$loop=0;
		while ($number_CC[$loop]>-1){
		    $new_ss=$seq_ss+$number_CC[$loop]+3;       
		    #$new_seqs=substr($seq,$new_ss,20);
			$temp=$new_ss+1;
			$temp1=$temp+19;
			$new_seqs=`samtools faidx $ref $info2[0]:$temp-$temp1 | sed '1d' | tr -d '\n'|tr a-z A-Z`;
		                   
		    $num=0;
		    $yes=0;
		    while($num<20){
			$c1=substr($re_on,$num,1);
			$c2=substr($new_seqs,$num,1);
			if($c1 eq $c2){
			    $yes+=1;
			}
			$num+=1;
		    }
		    if($yes > $PAM_GG){
			$distant=23-$number_CC[$loop]-3;
			print FILE2 "$info2[0]_$info2[1]\t$info2[3]\t$info2[4]\t$distant | \t$storefather{$site}\t$storemother{$site}\n";
#                                print "result1  $site\t$info2[3]\t$info2[4]\t$distant | \t$storefather{$site}\t$storemother{$site}\n";
		    }   
		    $loop+=1;
		}
		                
		                
		$loop1=0;
		while ($number_CC1[$loop1]>-1){
		    $new_ss1=$seq_ss+$number_CC1[$loop1]+3;
		    #$new_seqs1=substr($seq,$new_ss1,20);
			$temp=$new_ss1+1;
			$temp1=$temp+19;
			$new_seqs1=`samtools faidx $ref $info2[0]:$temp-$temp1 | sed '1d' | tr -d '\n'|tr a-z A-Z`;
			$num=0;
		    $yes=0;
		    while($num<20){
			$c1=substr($re_on,$num,1);
			$c2=substr($new_seqs1,$num,1);
			if($c1 eq $c2){
			    $yes+=1;
			}
			$num+=1;
		    }
		    if($yes > $PAM_AG){
			$distant=23-$number_CC1[$loop1]-3;
			print FILE4 "$info2[0]_$info2[1]\t$info2[3]\t$info2[4]\t$distant | \t$storefather{$site}\t$storemother{$site}\n";
		    }    
		    $loop1+=1;     
		}
		                    
		                                      
	    }
	                
	                
	                
	    if($seqe=~m/GG/ || $seqe=~m/AG/){
		$loce=0;
		$loce1=0;
		$number_G=0;
		$number_G1=0;
		@number_GG="";
		@number_GG1="";
		@change_G="";
		$seqe1=$seqe;
		while($loce>-1){
#print "aaaa      $loce\t$number_G\n";
		    $loce = index($seqe,"GG");
		    @change_G=split(//,$seqe);
		    $change_G[$loce]="O";
		    $seqe=join "",@change_G;
		    $number_GG[$number_G]=$loce;
		    $number_G+=1;
		}
		while($loce1>-1){

		    $loce1 = index($seqe1,"AG");
		    $seqe1=~s/AG/OO/;
		    $number_GG1[$number_G1]=$loce1;
		    $number_G1+=1;
		}
		$loope=0;
		while ($number_GG[$loope]>-1){
		    $new_se=$info2[1]+$number_GG[$loope]-21;
		    #$new_seqe=substr($seq,$new_se,20);
			$temp=$new_se+1;
			$temp1=$temp+19;
			$new_seqe=`samtools faidx $ref $info2[0]:$temp-$temp1 | sed '1d' | tr -d '\n'|tr a-z A-Z`;
#		    print "I need seq $new_seqe\tsite$new_se\n";
		    $num=0;
		    $yes=0;
		    while($num<20){
			$c1=substr($on,$num,1);
			$c2=substr($new_seqe,$num,1);
			if($c1 eq $c2){
			    $yes+=1;
			}
			$num+=1;
		    }
		    if($yes > $PAM_GG){
			$distant=$number_GG[$loope]-1;
			print FILE2 "$info2[0]_$info2[1]\t$info2[3]\t$info2[4]\t$distant | \t$storefather{$site}\t$storemother{$site}\n";
		    }   
		    $loope+=1;
		}
		                
		                
		$loope1=0;
		while ($number_GG1[$loope1]>-1){
		    $new_se1=$$info2[1]+$number_GG1[$loope1]-21;
		   #$new_seqe1=substr($seq,$new_se1,20);
			$temp=$new_se1+1;
			$temp1=$temp+19;
			$new_seqe1=`samtools faidx $ref $info2[0]:$temp-$temp1 | sed '1d' | tr -d '\n'|tr a-z A-Z`;
		    $num=0;
		    $yes=0;
		    while($num<20){
			$c1=substr($on,$num,1);
			$c2=substr($new_seqe1,$num,1);
			if($c1 eq $c2){
			    $yes+=1;
			}
			$num+=1;
		    }
		    if($yes > $PAM_AG){
			$distant=$number_GG1[$loope1]-1;
			print FILE4 "$info2[0]_$info2[1]\t$info2[3]\t$info2[4]\t$distant | \t$storefather{$site}\t$storemother{$site}\n";
		    }    
		    $loope1+=1;       
		}
		                
		                           
		                                                
	    }
                         
	}
       
	                                   
	                                   
	                                   
#------------------------------deletion---------------------------------------------------------------------------------------------   
	                                   
	                                   
	if($infor2[2] > 8 && $refl > $ALTl){
	        
 #   print "deletion\n";
	        
	    #$filename=uc($info2[0]);
	    #open (FILE,"$ref_path/$filename")||die "cannot open word_file: $!";
	    #$seq=<FILE>;
	 #   $seq=uc($seq);
	    $l=length($info2[3])-1;
	    $seq_ss=$info2[1]-23;
	    $seq_se=$info2[1]+$l;
	    #$seqs=substr($seq,$seq_ss,23);
		$temp=$seq_ss+1;
		$temp1=$temp+22;
		$seqs=`samtools faidx $ref $info2[0]:$temp-$temp1 | sed '1d' | tr -d '\n'|tr a-z A-Z`;
		
	    #$seqe=substr($seq,$seq_se,23);
		$temp=$seq_se+1;
		$temp1=$temp+22;
		$seqe=`samtools faidx $ref $info2[0]:$temp-$temp1 | sed '1d' | tr -d '\n'|tr a-z A-Z`;
#	      print "seqs\t$seqs\n";
#    print "seqe\t$seqe\n";
  
	    $stop=-5;
	    while($stop<0){
		
		if($seqs=~m/CC/ || $seqs=~m/CT/){
		    $loc=0;
		    $loc1=0;
		    $number_C=0;
		    $number_C1=0;
		    @number_CC="";
		    @number_CC1="";
		    @change_C="";
		    $seqs1=$seqs;
		    while($loc>-1){
			$loc = index($seqs,"CC");
			@change_C=split(//,$seqs);
			$change_C[$loc]="O";
			$seqs=join "",@change_C;
			$number_CC[$number_C]=$loc;
			$number_C+=1;
		    }
		            
		    while($loc1>-1){
			$loc1 = index($seqs1,"CT");
			$seqs1=~s/CT/OO/;
			$number_CC1[$number_C1]=$loc1;
			$number_C1+=1;
		    }
		    

                    
		    $loop=0;
		    while ($number_CC[$loop]>-1){
			$new_ss=$seq_ss+$number_CC[$loop]+3;
			#$new_seqs=substr($seq,$new_ss,20);
			$temp=$new_ss+1;
			$temp1=$temp+19;
			$new_seqs=`samtools faidx $ref $info2[0]:$temp-$temp1 | sed '1d' | tr -d '\n'|tr a-z A-Z`;
			$num=0;
			$yes=0;

#print "CCA      $new_seqs\n";
			while($num<20){
			    $c1=substr($re_on,$num,1);
			    $c2=substr($new_seqs,$num,1);
			    if($c1 eq $c2){
				$yes+=1;
			    }
			    $num+=1;
			}
			if($yes > $PAM_GG){
			    $distant=23-$number_CC[$loop]-3;
			    print FILE1 "$info2[0]_$info2[1]\t$info2[3]\t$info2[4]\t$distant | \t$storefather{$site}\t$storemother{$site}\n";
			               # print       "$info2[0]_$info2[1]\t$info2[3]\t$info2[4]\t$distant | \t$storefather{$site}\t$storemother{$site}\n";
			}   
			$loop+=1;
		    }
	

		    $loop1=0;
		    while ($number_CC1[$loop1]>-1){
			$new_ss1=$seq_ss+$number_CC1[$loop1]+3;
			#$new_seqs1=substr($seq,$new_ss1,20);
			$temp=$new_ss1+1;
			$temp1=$temp+19;
			$new_seqs1=`samtools faidx $ref $info2[0]:$temp-$temp1 | sed '1d' | tr -d '\n'|tr a-z A-Z`;
			$num=0;
			$yes=0;
			
#print "CTA $new_seqs1\n";
			
			while($num<20){
			    $c1=substr($re_on,$num,1);
			    $c2=substr($new_seqs1,$num,1);
			    if($c1 eq $c2){
				$yes+=1;
			    }
			    $num+=1;
			}
			if($yes > $PAM_AG){
			    $distant=23-$number_CC1[$loop1]-3;
			    print FILE3 "$info2[0]_$info2[1]\t$info2[3]\t$info2[4]\t$distant | \t$storefather{$site}\t$storemother{$site}\n";
			}    
			$loop1+=1;       
		    }
		                             
		                                  
		}
		
		if($seqe=~m/GG/ || $seqe=~m/AG/){
#		    print "GGGG\n";
		    $loce=0;
		    $loce1=0;
		    $number_G=0;
		    $number_G1=0;
		    @number_GG="";
		    @number_GG1="";
		    @change_G="";
		    $seqe1=$seqe;
		    while($loce>-1){
			$loce = index($seqe,"GG");
			@change_G=split(//,$seqe);
			$change_G[$loce]="O";
			$seqe=join "",@change_G;
			                      
			$number_GG[$number_G]=$loce;
			$number_G+=1;
		    }
		    $number_G1=0;
		    while($loce1>-1){
			$loce1 = index($seqe1,"AG");
			$seqe1=~s/AG/OO/;
			$number_GG1[$number_G1]=$loce1;
			$number_G1+=1;
		    }
		    $loope=0;
		         #   print "aaaaa:$number_GG[0]\t";
		    while ($number_GG[$loope]>-1){
			$new_se=$seq_se+$number_GG[$loope]-21;
			#$new_seqe=substr($seq,$new_se,20);
			$temp=$new_se+1;
			$temp1=$temp+19;
			$new_seqe=`samtools faidx $ref $info2[0]:$temp-$temp1 | sed '1d' | tr -d '\n'|tr a-z A-Z`;
			$num=0;
			$yes=0;


#print "NGG      $new_seqe\non\t$on\n";

			while($num<20){
			    $c1=substr($on,$num,1);
			    $c2=substr($new_seqe,$num,1);
			    if($c1 eq $c2){
				$yes+=1;
			    }
			    $num+=1;
			}
			if($yes > $PAM_GG){
			    $distant=$number_GG[$loope]-1;
			    print FILE1 "$info2[0]_$info2[1]\t$info2[3]\t$info2[4]\t$distant | \t$storefather{$site}\t$storemother{$site}\n";
			}   
			$loope+=1;
		    }
		                        
		                        
		    $loope1=0;
		    while ($number_GG1[$loope1]>-1){
			$new_se1=$info2[1]+$number_GG1[$loope1]-21;
			#$new_seqe1=substr($seq,$new_se1,20);
			$temp=$new_se1+1;
			$temp1=$temp+19;
			$new_seqe1=`samtools faidx $ref $info2[0]:$temp-$temp1 | sed '1d' | tr -d '\n'|tr a-z A-Z`;
			$num=0;
			$yes=0;
			
#print "NAG $new_seqe1\n";
			
			while($num<20){
			    $c1=substr($on,$num,1);
			    $c2=substr($new_seqe1,$num,1);
			    if($c1 eq $c2){
				$yes+=1;
			    }
			    $num+=1;
			}
			if($yes > $PAM_AG){
			    $distant=$number_GG1[$loope1]-1;
			    print FILE3 "$info2[0]_$info2[1]\t$info2[3]\t$info2[4]\t$distant | \t$storefather{$site}\t$storemother{$site}\n";
			}    
			$loope1+=1;        
		    }
		                        
		                                            
		}
		
		
		
		
		$seq_ss=$info2[1]-1;
		$seq_se=$info2[1]-1;
		$seqs=$info2[3];
		$seqe=$info2[3];
		
		$stop+=4;
#print "aaaaaaaaa\n";
	    }
	        
	        
	        
	        
	        
	}
	                                   
	                                   
	                                   
	            
 #   }                
            
  
    
    }
}
