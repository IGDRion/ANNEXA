#!/bin/bash


# Format a gtf file with only 12 fields if more (to be used by overlap)

#!/bin/bash


# Format a bed file to an tabulated bed

# USAGE

usage() {
    echo "#" >&2
    echo "# Format a gtf 12+ in cleaned gtf 12 to be used by overlap...">&2
    echo "# USAGE: `basename $0` <file.gtf>">&2
    
}
################################################################################
infile=$1

if [ -p /dev/stdin ];then

    #from STDIN
    cat /dev/stdin | awk '{for (i=1;i<9;i++) {printf("%s\t",$i) }; for (i=9;i<=NF;i++) {printf("%s ",$i)}; print ""}' 

elif [ $# -eq  1 ];then

    awk '{for (i=1;i<9;i++) {printf("%s\t",$i) }; for (i=9;i<=NF;i++) {printf("%s ",$i)}; print ""}' $infile

else
    echo "#Error! no argument  file or file empty !"
    usage;
    exit 0;
fi
