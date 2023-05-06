while getopts i:o: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        c) out_tmp=${OPTARG};;
        t) out_count=${OPTARG};;
    esac
done

ls -1 $input*/abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' | sed 's/results\/test_paried_end\/quantification\/kallisto\///g' > $out_tmp
ls -1 $input*/abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\test_counts$_\n"' | sed 's/results\/test_paried_end\/quantification\/kallisto\///g' > $out_count
paste $input*/abundance.tsv | cut -f 1,2,5,10,15,20,25,30 | grep -v "tpm" >> $out_tmp
paste $input*/abundance.tsv | cut -f 1,2,4,9,14,19,24,29 | grep -v "tpm" >> $out_count