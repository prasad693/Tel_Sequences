#!/usr/bin/env bash

threads=10
output="T2T_sequences"
show_help=false


# Function to display usage information
display_help() {
  echo -e "\nUsage: $0 [-a Asm_fasta] [-r Reference_fasta] [-o Output prefix] [-m motif] [-t threads] [-h]"
  echo -e "Options:"
  echo -e "\t  -a file         Assembly fasta file [required]"
  echo -e "\t  -r file         Reference fasta file [optional]"
  echo -e "\t  -o file         Output file prefix [default: T2T_sequences]"
  echo -e "\t  -m string       Telomere motif to search [default: TTAGGG]"
  echo -e "\t  -t int          Number of threads used for minigraph, seqkit [default: 10]"
  echo -e "\t  -h,             Show this help message\n"
  echo -e "\nmailto: sarashettp@gis.a-star.edu.sg\n"
  exit 0
}

# Process command-line options using getopts
while getopts "a:r:o:t:m:h" opt; do
  case "$opt" in
    a) assembly=$OPTARG ;;
    r) reference=$OPTARG ;;
    o) output=$OPTARG ;;
	m) tmotif=$OPTARG ;;
    t) threads=$OPTARG ;;
    h) display_help ;;
    *) display_help ;;

  esac
done


# Check if required options are missing
if [ -z "$assembly" ]; then
  echo "Error: -a option is required."
  display_help
fi

if [[ ! -e $assembly || ! -s $assembly ]]; then
	echo -e "\n$assembly does not exists or an empty file. Please check the input file\n\n";
	exit 0
else
	echo -e "Proceeding analysis with"
	echo -e "\tInput fasta:             $assembly"
fi

if [ -n "$reference" ]; then
  if [[ ! -e $reference || ! -s $reference ]]; then
    echo -e "\n$reference does not exist or is an empty file. Please check the input file\n\n"
    exit 0  # Exit with a non-zero status to indicate an error
  else
	echo -e "\tReference fasta:         $reference"
	mapping=1
  fi
fi

if [ -n "$threads" ]; then
  if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo -e "\nA non-integer value is given for threads option. Switching back to the default value of 10\n";
  fi
fi

if [ -n "$tmotif" ]; then
  if ! [[ "$tmotif" =~ ^[A-Za-z]+$ ]] || [[ ${#tmotif} -lt 4 ]]; then
    echo -e "\nInvalid motif. Motif must have atleast 4 alphabetic characters. Switching back to the default value\n";
	tmotif="TTAGGG"
	echo -e "\tTelomere motif:          $tmotif"
  else
	tmotif=$(echo "$tmotif" | tr '[:lower:]' '[:upper:]')
	echo -e "\tTelomere motif:          $tmotif"
  fi
else
	tmotif="TTAGGG"
	echo -e "\tTelomere motif:          $tmotif"
fi

#echo -e "Formating the input file\n"
#echo -e "seqkit seq -w 0 -j $threads -o $output.format.fa $assembly\n\n"

#seqkit seq -w 0 -j $threads -o $output.format.fa $assembly

echo -e "\n\nExtracting sequence information\n"
echo -e "seqkit fx2tab -j $threads -o $output.seqinfo.txt $assembly\n\n"
seqkit fx2tab -j $threads -C N -l -n -o $output.seqinfo.txt $assembly

echo -e "Extracting 1000bps from start and end of sequences for telomere motif identification\n"
echo -e "seqkit subseq -w 0 -j $threads -r 1:1000/-1000:-1 $assembly >$output.teloinput.fa\n\n"

seqkit subseq -w 0 -j $threads -r 1:1000 $assembly | awk '/^>/ {print $0 "_Start"; next} {print}' >$output.teloinput.fa
seqkit subseq -w 0 -j $threads -r -1000:-1 $assembly | awk '/^>/ {print $0 "_End"; next} {print}' >>$output.teloinput.fa

echo -e "searching telomere motif\n"
echo -e "tidk search -d . -f $output.teloinput.fa -o $output -s $tmotif -w 1000\n\n"

tidk search -d . -f $output.teloinput.fa -o $output -s $tmotif -w 1000

tail -n +2 $output"_telomeric_repeat_windows.csv" | awk -F "," '$3 + $4 >=15{print $1}' | while read line; do contig=`echo $line | sed 's/_.*//'`; egrep --color ${contig} $output".seqinfo.txt" | tr -s "[:blank:]" "\t" ; done | cut -f1,2 | sort | uniq -c | sort -rgk1 | tr -s "[:blank:]" "\t" | sed 's/^\t//' | grep "^2" | cut -f2- >$assembly"_motif_"$output".txt"

if [[ $mapping -eq 1 ]]; then
	echo -e "\n\nReference sequence provided, proceeding with mapping bases T2T filtering\n"
	minigraph -x asm --show-unmap=yes -t $threads -K1.9g $reference $assembly >$output".paf"
	cov_cal -T $output".paf" 2>$output"_alignment_T2T.txt"
	tail -n +2 $output"_telomeric_repeat_windows.csv" | awk -F "," '$3 + $4 >= 15{print $1}' | while read line; do contig=$(echo "$line" | sed 's/_.*//'); grep $contig $output".seqinfo.txt" | tr -d '\n'; echo -ne '\t'; grep  $contig $output"_alignment_T2T.txt" | tr -d '\n'; echo -ne '\n' ; done | cut -f 1-3,9,10 | awk 'NF>4{print}' | sort | uniq -c | tr -s "[:blank:]" "\t" | sed 's/^\t//' | grep "^2" | cut -f2- >$assembly"_alignment_"$output".txt"
fi	

echo -e "\n\nMotif based T2T sequences are                 : "$assembly"_motif_"$output".txt\n"	
echo -e "Alignment and motif based T2T sequences       : "$assembly"_alignment_"$output".txt\n"
	