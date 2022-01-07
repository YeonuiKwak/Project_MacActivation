#hojoong code for AU content in the 3UTR table.
#cat figS5D.sh 
cat $1 |
awk '{	seq = "";
	l = length($2);
	if(l < 5000) for(i = 1; i <= 5000-l; ++i) seq = seq"N";
	if(l > 5000) seq = seq toupper(substr($2,l-4999,5000));
	else seq = seq toupper($2);
	split("", au);
	for(i = 1; i <= 5000; ++i) {
		c = substr(seq, i, 1);
		if(c == "A" || c == "T") au[int((i-1)/20)] += 0.05;
		if(c == "N") au[int((i-1)/20)] = -1;
	}
	printf $1"\t"l;
	for(i = 0; i < 250; ++i) printf "\t"0 + au[i];
	printf "\n";
}' > cutom3UTR.AUposition.txt 


