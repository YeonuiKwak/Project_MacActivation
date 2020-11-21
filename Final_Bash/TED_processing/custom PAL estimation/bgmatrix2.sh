#BGMatrix second
  #!/bin/bash

  if [ $# -lt 2 ]; then
  	echo -e "Usage:\t bgmatrix.sh [options] -a <bed> -b <bedgraph>"
  	echo -e "Options:"
  	echo -e "\t-bin\tbin size (default=10)"
  	echo -e "\t-split\tTreat BED12 entries as distinct BED intervals."
  	exit
  fi

  SPLIT=false
  BIN=10

  while [[ $# -ge 1 ]]; do
  	key="$1"
  	case $key in
  		-a)
  		BED="$2"
  		shift
  		;;
  		-b)
  		BG="$2"
  		shift
  		;;
  		-bin)
  		BIN="$2"
  		shift
  		;;
  		-split)
  		SPLIT=true
  		;;
  		--default)
  		;;
  		*)
  		;;
  	esac
  	shift
  done

  BEDCOL="$(head -1 $BED | awk '{print NF}')"

  if [ $BEDCOL -eq 12 ]; then
  	if $SPLIT; then
  		bedtools map -a $BED -b $BG -c 2,4 -o collapse -sorted -split -nobuf -delim " " > _tmp
  		cat _tmp | \
  			awk -v bin=$BIN ' 
  			BEGIN {
  				bin=0+bin;
  			}
  			{
  				start=0+$2;
  				split($12,exs,",");
  				split($11,exb,",");
  				nex=0+$10
  				exs[nex+1]=999999999;
  				n=(NF-12)/2;
  				prev=0;
  				v=0;
  				curexs=exs[1];
  				nextexs=exs[2];
  				curexi=1;
  				exsum=0;
  				if($13!=".") for(i=13;i<13+n;++i) {
  					rp=$i-start;
  					while(rp>=nextexs) { curexs=nextexs; exsum+=exb[curexi]; ++curexi; nextexs=exs[curexi+1]; }
  					idx=int((rp-curexs+exsum)/bin);
  					if(idx!=prev) {
  						printf v"\t";
  						if(idx>prev+1) for(j=prev;j<idx-1;++j) printf "0\t";
  						prev=idx;
  						v=0;
  					}
  					v+=$(i+n);
  				}
  				printf v;
  				for(j=curexi;j<=nex;++j) exsum+=exb[j];
  				N=int(exsum/bin)+1
  				for(j=prev+1;j<N;++j) printf "\t0";
  				printf "\n";
  			}'
  		rm _tmp
fi

fi

