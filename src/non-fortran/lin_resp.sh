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

FOLDER=`grep '^out_folder\s' $INFILE`
FOLDER=${FOLDER##* }
FOLDER=${FOLDER#\'}
FOLDER=${FOLDER%\'}

from=`grep '^freq_min\s' $INFILE`
from=${from#* }

to=`grep '^freq_max\s' $INFILE`
to=${to#* }

step=`grep '^step\s' $INFILE`
step=${step#* }

CHIIN=`grep '^chi_out\s' $INFILE`
CHIIN=${CHIIN##* }
CHIIN=${CHIIN#\'}
CHIIN=${CHIIN%\'}
CHIOUT=${CHIIN}.chi1

cat $INFILE > $INFILE.tmp1
#echo "max_order  1" >> $INFILE.tmp1
echo "calc_diffs  .FALSE." >> $INFILE.tmp1

FILE=$FOLDER/$CHIOUT
Q_fin=$FOLDER/Q.dat

if [ -e $FOLDER ]; then
    echo -n 'Do you want to generate new data? (Y/n) '
    read GD
fi

FREQS=`seq $from $step $to`

#echo -n 'Process all files in folder or just those generated? (a/G)'
#read PRO
#if [[ "$PRO" == 'a' ]]; then
#    FREQS=`ls $FOLDER/*.out`
#fi

if [[ "$GD" != 'n' ]]; then

    Q=$Q_fin

    if [ ! -e $FOLDER ]; then
        mkdir -p $FOLDER
    fi
    touch $Q

    Q_tmp=$Q
    i=0
    while [[ -e $Q_tmp ]]; do
        i=$((i + 1))
        Q_tmp=${Q}_$i
    done
    Q=$Q_tmp

#    rm -f $Q
    echo 'Generating data...'
    echo
    for freq in $FREQS; do

        freq=${freq##*$FOLDER}
        freq=${freq##*/}
        freq=${freq%*.out}
        cat $INFILE.tmp1 > $INFILE.tmp
        echo "omega_ev  $freq" >> $INFILE.tmp
        echo "name  '$freq'" >> $INFILE.tmp

        echo "Generating frequency $freq ..."
        if [[ "$EXISTS" != 'y' ]]; then
            PATHTOBIN/rho_prop.exe $INFILE.tmp >> $Q
        fi

    done

    if [ ! -e $Q_fin ]; then
        touch $Q_fin
        touch $Q_fin.tmp
    fi

fi

rm -f $Q_fin
touch $Q_fin
echo 'frequency (eV)        Q_MNP (W)             Q_SQD (W)             Q = Q_SQD+Q_MNP (W)   ' > $Q_fin
touch $Q_fin.tmp
for foo in `ls ${FOLDER}/Q.dat_*`; do
    cat $foo >> $Q_fin.tmp
done
sort -gk 1 $Q_fin.tmp >> $Q_fin

rm -f $Q_fin.tmp
echo
echo "*** Q_mnp written to file $Q_fin ***"

##############################
##############################
exit
##############################
##############################

echo
echo 'Post-processing...'
echo

#if [ ! -e $FILE ]; then
#    HEAD='yes'
#fi

#if [[ "$GD" == 'n' ]]; then
#    echo -n 'Process all files in folder or just those generated? (a/G)'
#    read PRO
#fi

i=0
for freq in $FREQS; do
    i=$((i+1))
    
    cat $INFILE.tmp1 > $INFILE.tmp
    freq=${freq##*$FOLDER}
    freq=${freq##*/}
    freq=${freq%*.out}
    echo "name  '${freq}-$CHIIN'" >> $INFILE.tmp

    FFF=`grep 'Laser Freq' $FOLDER/${freq}.params`
    FFF=${FFF##* }
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
    echo "in_file  '$FOLDER/$freq.out'" >> $INFILE.tmp
    echo "probe_freq  $FFF" >> $INFILE.tmp
    echo "num_freqs  $nf" >> $INFILE.tmp
    echo freqs  $FFF $ens >> $INFILE.tmp

    echo "Processing frequency $freq ..."
    PATHTOBIN/post_proc.exe $INFILE.tmp > /dev/null
done


E0=`grep 'Laser Ampli' $FOLDER/${freq}.params`
E0=${E0##* }

eps_eff1=`grep 'eps_eff1' $FOLDER/${freq}.params`
eps_eff1=${eps_eff1##* }

grep ! ${FOLDER}/*${CHIIN}-pp.freqs > $FILE.new
awk -F ' ' '{OFMT="%+.14e"; print $2/1, $3/1, $4/'$E0'/'$eps_eff1', $5/'$E0'/'$eps_eff1', $6}' $FILE.new > $FILE.tmp

#MAX=`cat ./max_sqd`
#MAX=${MAX##* }
#MAX=${MAX##* }
#MAX=${MAX##* }
#MAX=${MAX%%????????}

#if [ ! -e $FOLDER/sqd.chi1 ];then
#    MAX=1.0
#fi

#####################~!!!!!!!!!!!!!!!!!!1
MAX=1.0
#####################~!!!!!!!!!!!!!!!!!!1

awk -F ' ' '{OFMT="%+.14e"; print $1/1, $2/1, $3/'$MAX', $4/'$MAX', $5}' $FILE.tmp > $FILE.new
cat $FILE.new > $FILE.tmp

#if [[ "$HEAD" == 'yes' ]]; then
    echo 'frequency (a.u.)      frequency (eV)        chi (real)            chi (imaginary)       estimated relative error (%)' > $FILE
#fi

cat $FILE.tmp >> $FILE
sort -gk 1 $FILE > $FILE.tmp
mv $FILE.tmp $FILE

echo
echo "Finished."
echo "*** chi data written to ${FILE}. ***"

rm -rf $FILE.new
rm -rf $FILE.tmp
rm -rf $INFILE.tmp1
rm -rf $INFILE.tmp
