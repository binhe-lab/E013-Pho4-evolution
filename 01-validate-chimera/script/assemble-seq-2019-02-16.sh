# Use EBNOSS/merger to merge Sanger sequencing results
# 7 oct 2018
# Bin He
# modified 16 fev 2019
# sh assemble-seq.sh INPUT
# INPUT should contain:
# - line 1, fasta input file name
# - line 2, output file name
# - line 3+, strain names to be merged

IN=$(sed '1q;d' $1);        echo "Input:  " $IN
OUT=$(sed '2q;d' $1);       echo "Output: " $OUT
NAMES=$(sed -n '3,$ p' $1); echo "Seq Names: " $NAMES

for n in $NAMES
do
	echo "Processing $n..."
	# set variables
	fwd="${n}-p355"
	rev1="${n}-p400"
	rev2="${n}-p077|${n}-p032"

	# extract forward sequence
	bioawk -c fastx -v pat="$fwd" '$name ~ pat {print ">"$name;print $seq}' $IN > ./tmp/fwd.fa
	echo "${n}-p355"
	# extract reverse sequence
	m=$(grep $rev1 $IN); echo $m;
	if [[ $m != "" ]] # if the strain has been sequenced with p400
	then              # merge p355 with p400, output as new fwd.fa, then merge with p077 or p032
		echo "${n}-p400"
		bioawk -c fastx -v pat="$rev1" '$name ~ pat {print ">"$name;print revcomp($seq)}' $IN > ./tmp/rev.fa
		merger ./tmp/fwd.fa ./tmp/rev.fa -outseq ./tmp/fwd.fa -outfile ./tmp/${n}-p400.merger
	fi
	echo "${n}-p077 or ${n}-p032"
	bioawk -c fastx -v pat="$rev2" '$name ~ pat {print ">"$name;print revcomp($seq)}' $IN > ./tmp/rev.fa
	# assemble sequence
	merger ./tmp/fwd.fa ./tmp/rev.fa -outseq ./tmp/${n}.fa -outfile ./tmp/${n}.merger
	cat ./tmp/${n}.fa >> $OUT
done

#cat ./tmp/XZ*.fa > ../data/sequencing/2018-10-07-assembled-seq.fa
echo "done"
