use XML::Simple;
$rawdata_list=$ARGV[0];
$thread_num=$ARGV[1];
#get project name
@temp=split(/\//,$rawdata_list);
@temp=split(/\./,$temp[-1]);
$project_name=$temp[0];
undef @temp;
mkdir("./$project_name") unless(-d "./$project_name");
mkdir("./$project_name/metadata") unless(-d "./$project_name/metadata");

#get reference and data location 
$xml = XMLin($rawdata_list,ForceArray => 1);
%member_tree=%{${${$xml}{'rawdata_list'}}[0]};
$reference=${${${${$xml}{'refdata_list'}}[0]}{'reference'}}[0];
$guide_RNA_seq=${${${${$xml}{'refdata_list'}}[0]}{'guide_RNA_seq'}}[0];
$dbsnp=${${${${$xml}{'refdata_list'}}[0]}{'dbsnp'}}[0];
$OneKG_indel=${${${${$xml}{'refdata_list'}}[0]}{'OneKG_indel'}}[0];
$Mills_indel=${${${${$xml}{'refdata_list'}}[0]}{'Mills_indel'}}[0];
$gatk=${${${${$xml}{'refdata_list'}}[0]}{'gatk3_8'}}[0];
$gatk4=${${${${$xml}{'refdata_list'}}[0]}{'gatk4'}}[0];
$picard=${${${${$xml}{'refdata_list'}}[0]}{'picard'}}[0];
$samtools=${${${${$xml}{'refdata_list'}}[0]}{'samtools'}}[0];
$bwa=${${${${$xml}{'refdata_list'}}[0]}{'bwa'}}[0];
$java=${${${${$xml}{'refdata_list'}}[0]}{'java'}}[0];
$fastp=${${${${$xml}{'refdata_list'}}[0]}{'fastp'}}[0];
#make bwa reference and single chromosome reference from reference file
mkdir("./$project_name/reference") unless(-d "./$project_name/reference");
$reference_dict= $reference;
$reference_dict=~s/.fa/.dict/;
if( -e $reference ){
	if(( -e $reference.amb ) && ( -e $reference.ann ) && ( -e $reference.bwt ) && ( -e $reference_dict ) && ( -e $reference.fai ) && ( -e $reference.pac ) && ( -e $reference.sa ) ){

	}else{
		@temp=split(/\//,$reference);
		@temp=split(/\./,$temp[-1]);
		$reference_prefix=$temp[0];
		undef @temp;
		if(! ((-e "./$project_name/reference/$reference_prefix.fa" )&&(-e "./$project_name/reference/$reference_prefix.dict" )&&(-e "./$project_name/reference/$reference_prefix.fa.fai" ))){
			print $reference_prefix."\n"; 
			`cp $reference ./$project_name/reference`;
			$reference="./$project_name/reference/$reference_prefix.fa";
			`$samtools faidx $reference`;
			`$bwa index $reference`;
			`$java -Xmx20g -jar $gatk4 CreateSequenceDictionary -R ./$project_name/reference/$reference_prefix.fa -O ./$project_name/reference/$reference_prefix.dict`;
		}
	}
}else{
	print "reference not found!\n";
}

if( -e $dbsnp ){
	if(! -e $dbsnp.idx ){
		@temp=split(/\//,$dbsnp);
		$new_dbsnp="./$project_name/reference/$temp[-1]";
		undef @temp;
		if(! (-e $new_dbsnp )&& !(-e $new_dbsnp.idx )){
			`cp $dbsnp ./$project_name/reference`;
			`$java -jar $gatk4 IndexFeatureFile -F $new_dbsnp`;
		}
		$dbsnp=$new_dbsnp;
	}
}else{
	print "dbsnp not found!\n";
}

if( -e $OneKG_indel ){
	if(! -e $OneKG_indel.idx ){
		@temp=split(/\//,$OneKG_indel);
		$new_OneKG_indel="./$project_name/reference/$temp[-1]";
		undef @temp;
		if(! (-e $new_OneKG_indel )&& !(-e $new_OneKG_indel.idx )){
			`cp $OneKG_indel ./$project_name/reference`;
			`$java -jar $gatk4 IndexFeatureFile -F $new_OneKG_indel`;
		}
		$OneKG_indel=$new_OneKG_indel;
	}
}else{
	print "OneKG_indel not found!\n";
}
if( -e $Mills_indel ){
	if(! -e $Mills_indel.idx ){
		@temp=split(/\//,$Mills_indel);
		$new_Mills_indel="./$project_name/reference/$temp[-1]";
		undef @temp;
		if(! (-e $new_Mills_indel )&& !(-e $new_Mills_indel.idx )){
			`cp $Mills_indel ./$project_name/reference`;
			`$java -jar $gatk4 IndexFeatureFile -F $new_Mills_indel`;
		}
		$Mills_indel=$new_Mills_indel;
	}
}else{
	print "Mills_indel not found!\n";
}

if($thread_num>16){
	$fastp_thread_num=16;
}else{
	$fastp_thread_num=$thread_num;
}
$samtools_memory=int(100/$thread_num).'G';
#qualify, mapping and sort data
foreach $member(keys %member_tree){
	mkdir("./$project_name/$member") unless(-d "./$project_name/$member");
	$merge_count=0;
	$temp_file_list='';
	foreach $lane (@{$member_tree{$member}}){
		$merge_count++;
		$names=${\split(/\//,${$lane}{R1})};
		`$fastp -w $fastp_thread_num -c -h "./$project_name/metadata/$member.$names.qc.html" -i ${$lane}{R1} -o "./$project_name/$member/temp.r1.qc.fq.gz" -I ${$lane}{R2} -O "./$project_name/$member/temp.r2.qc.fq.gz"`;
		`$bwa mem -t $thread_num -M $reference ./$project_name/$member/temp.r1.qc.fq.gz ./$project_name/$member/temp.r2.qc.fq.gz | $samtools view -@ $thread_num -bS - > "./$project_name/$member/temp.$merge_count.bam"`;
		`$samtools sort -m $samtools_memory -@ $thread_num "./$project_name/$member/temp.$merge_count.bam" -o "./$project_name/$member/temp.$merge_count.sorted.bam"`;
		$temp_file_list.="./$project_name/$member/temp.$merge_count.sorted.bam\n";
	}
	open (FILE, ">./$project_name/$member/temp_file_list.txt");
	print FILE $temp_file_list;
	close FILE;
	`$samtools merge -b ./$project_name/$member/temp_file_list.txt "./$project_name/$member/$member.sorted.bam"`;
	`rm ./$project_name/$member/temp.*`;
}

#call vcf from bam file (mutiple process)
%fork_child_list;
foreach $member(keys %member_tree){
	$pid = fork();
	if (!defined($pid)) {
		print "Error in fork: $!";
		exit 1;
	}
	if ($pid == 0) {
		#2. AddOrReplaceReadGroups
		`$java -Xmx20g -jar $picard AddOrReplaceReadGroups I="./$project_name/$member/$member.sorted.bam" O="./$project_name/$member/$member.sorted.RG.bam" RGID="1" RGLB=$member"L" RGPL=illimina RGPU=run RGSM=$member`;
		#3. MarkDuplicates
		`$java -Xmx20g -jar $picard MarkDuplicates I="./$project_name/$member/$member.sorted.RG.bam" O="./$project_name/$member/$member.sorted.RG.markdup.bam" M="./$project_name/$member/$member.markdup.metrics" CREATE_INDEX=true REMOVE_DUPLICATES=true`;	
		#4. BaseRecalibration Calculate
		`$java -Xmx20g -jar $gatk4 BaseRecalibrator -R $reference -I "./$project_name/$member/$member.sorted.RG.markdup.bam" -O "./$project_name/$member/$member.BQSR.metrics" --known-sites $OneKG_indel --known-sites $Mills_indel --known-sites $dbsnp`; 
		#5. Print Recalibrated Reads
		`$java -Xmx20g -jar $gatk4 ApplyBQSR -bqsr "./$project_name/$member/$member.BQSR.metrics" -I "./$project_name/$member/$member.sorted.RG.markdup.bam" -O "./$project_name/$member/$member.sorted.RG.markdup.BQSR.bam" `;
		#6. Genotype Caller 
		`$java -Xmx20g -jar $gatk -glm INDEL -R $reference -T UnifiedGenotyper -I "./$project_name/$member/$member.sorted.RG.markdup.BQSR.bam" -o "./$project_name/$member/$member.vcf" -metrics "./$project_name/$member/$member.metrics" -stand_call_conf 20 -minIndelFrac 0.03 -minIndelCnt 2`; 
		# 7. remove temp files
		`rm ./$project_name/$member/$member.sorted.RG.bam`;
		`rm ./$project_name/$member/$member.sorted.RG.markdup.bam`;
		print "$member call vcf end\n";
		exit 0;
	}else{
		$fork_child_list{$pid}=1;
	}	
}
sleep(3);
while( scalar(keys(%fork_child_list))>0 ){
   while (($exitPid = waitpid(-1, WNOHANG)) > 0){
      delete($fork_child_list{$exitPid});
   }
}
undef %fork_child_list;

#find offtarget position from vcf file (mutiple process)
%fork_child_list;
foreach $member(keys %member_tree){
	if($member ne 'father'&&$member ne 'mother'){
		$pid = fork();
		if (!defined($pid)) {
			print "Error in fork: $!";
			exit 1;
		}
		if ($pid == 0) {
			`perl ./offtarget.pl $project_name $member $guide_RNA_seq $reference`;
			`cat ./${project_name}/$member/*.txt >  ./${project_name}/$member/merge_offtarget_sites.txt`;
			`perl ./rm_repeat.pl ${project_name} $member`;
			print "$member find offtaret sites end\n";
			exit 0;
		}else{
			$fork_child_list{$pid}=1;
		}	
	}
	
}
sleep(3);
while( scalar(keys(%fork_child_list))>0 ){
   while (($exitPid = waitpid(-1, WNOHANG)) > 0){
      delete($fork_child_list{$exitPid});
   }
}
print "all over\n"; 
