while getopts i:o: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        o1) out1=${OPTARG};;
        o2) out2=${OPTARG};;
    esac
done

ls -1 $input*/abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' | sed 's/results\/test_paried_end\/quantification\/kallisto\///g' > $out1
paste $input*/abundance.tsv | cut -f 1,2,5,10,15,20,25,30 | grep -v "tpm" >> $out1

ls -1 $input*/abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\test_counts$_\n"' | sed 's/results\/test_paried_end\/quantification\/kallisto\///g' > $out2
paste $input*/abundance.tsv | cut -f 1,2,4,9,14,19,24,29 | grep -v "tpm" >> $out2