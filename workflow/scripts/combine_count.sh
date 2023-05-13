END=$(ls -1 $2*/abundance.tsv | wc -l)

cnt=0; iter=0; str="1,2,"
for i in $(seq 4 $((END+4)))
    do
    # echo $((i+cnt-iter))
    iter=$((i+1))
    cnt=$((cnt+5))
    str=$str$((i+cnt-iter))","
    done
str=${str%?}

ls -1 $2*/abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\ttest_counts$_\n"' | sed "s/results\/$1\/quantification\/kallisto\///g" > $3
paste $2*/abundance.tsv | cut -f $str | grep -v "est_counts" >> $3