#!/bin/bash
if [ $# -lt 1 ]; then
    INFILE=lr.in
else
    INFILE=$1
fi

if [ ! -e $INFILE ]; then
    echo "Error: Input file \"$INFILE\" does not exist."
    exit 1
fi

FOLDER=`grep 'out_folder' $INFILE`
FOLDER=${FOLDER##* }
FOLDER=${FOLDER#\'}
FOLDER=${FOLDER%\'}

from=`grep 'freq_min' $INFILE`
from=${from##* }

to=`grep 'freq_max' $INFILE`
to=${to##* }

step=`grep 'step ' $INFILE`
step=${step##* }

cat $INFILE > $INFILE.tmp1
echo "read_col  3" >> $INFILE.tmp1
echo "max_order  1" >> $INFILE.tmp1
echo "calc_diffs  .TRUE." >> $INFILE.tmp1

FILE=$FOLDER/chi1.dat
rm -f $FILE

i=0
for freq in `seq $from $step $to`; do
    i=$((i+1))

    cat $INFILE.tmp1 > $INFILE.tmp
    echo "omega_ev  $freq" >> $INFILE.tmp
    echo "name  '$freq'" >> $INFILE.tmp

    echo "Processing frequency $freq ..."
    ./bin/F90/rho_prop.exe $INFILE.tmp

    nf=1
# Uncomment this section to gather the energy_differences from the parameters
# file. These frequencies will also be included in the process. However, if
# the laser frequency is very close to an energy gap, then a large error
# occurs. In this case, try increasing the 'tol' in the subroutine 
# 'remove_duplicates' which compares two frequencies and eliminates any that
# are 'too close' to each other.
#    if [ "$i" == 1 ]; then
#        from=`grep -nr 'Levels' $FOLDER/${freq}.params` 
#        from=${from%%:*}
#        from=$((from +3))
#        to=`grep -nr 'Director' $FOLDER/${freq}.params` 
#        to=${to%%:*}
#        to=$((to -3))
#        nf=$((to - from + 2))
#        ens=`sed -n ${from},${to}p $FOLDER/${freq}.params`
#    fi

    FFF=`grep 'Laser Freq' $FOLDER/${freq}.params`
    FFF=${FFF##* }
    echo "in_file  '$FOLDER/$freq.out'" >> $INFILE.tmp
    echo "probe_freq  $FFF" >> $INFILE.tmp
    echo "num_freqs  $nf" >> $INFILE.tmp
    echo freqs  $FFF $ens >> $INFILE.tmp

    ./bin/F90/post_proc.exe $INFILE.tmp >> $FILE
done

E0=`grep 'Laser Ampli' $FOLDER/${freq}.params`
E0=${E0##* }

awk -F ' ' '{OFMT="%+.14e"; print $1/1, $2/1, $3/'$E0', $4/'$E0', $5}' $FILE > $FILE.tmp
mv $FILE.tmp $FILE
rm -rf $FILE.tmp
rm -rf $INFILE.tmp1
rm -rf $INFILE.tmp

sed -i '1ifrequency (a.u.)      frequency (eV)        chi (real)            chi (imaginary)       estimated relative error (%)' $FILE
