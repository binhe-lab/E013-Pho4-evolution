# Use EBNOSS/merger to merge Sanger sequencing results
# 7 oct 2018
# Bin He

file=../data/sequencing/2018-10-07-all-sanger-trimmed-50bp-40pct.txt

while read n
do
	echo "Processing $n..."
	# set variables
	fwd="${n}-p355"
	rev1="${n}-p400"
	rev2="${n}-p077|${n}-p032"

	# extract forward sequence
	bioawk -c fastx -v pat="$fwd" '$name ~ pat {print ">"$name;print $seq}' $file > ./tmp/fwd.fa
	echo "${n}-p355"
	# extract reverse sequence
	m=`grep $rev1 $file`
	if [[ $m != "" ]] # if the strain has been sequenced with p400
	then              # merge p355 with p400, output as new fwd.fa, then merge with p077 or p032
		echo "${n}-p400"
		bioawk -c fastx -v pat="$rev1" '$name ~ pat {print ">"$name;print revcomp($seq)}' $file > ./tmp/rev.fa
		merger ./tmp/fwd.fa ./tmp/rev.fa -outseq ./tmp/fwd.fa -outfile ./tmp/${n}-p400.merger
	fi
	echo "${n}-p077 or ${n}-p032"
	bioawk -c fastx -v pat="$rev2" '$name ~ pat {print ">"$name;print revcomp($seq)}' $file > ./tmp/rev.fa
	# assemble sequence
	merger ./tmp/fwd.fa ./tmp/rev.fa -outseq ./tmp/${n}.fa -outfile ./tmp/${n}.merger
done < strain_names.txt

cat ./tmp/XZ*.fa > ../data/sequencing/2018-10-07-assembled-seq.fa
echo "done"
