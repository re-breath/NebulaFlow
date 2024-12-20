#!/bin/sh
#lipai@mail.ustc.edu.cn

begin=$1
if [ $begin =='' ]; then begin=0; fi    
awk -v begin=$begin '/E0/{if ( i<begin  ) i++;else print $0 }' OSZICAR >temp.e
    
awk '/POSITION/,/drift/{
    if(NF==6) print $4,$5,$6;
    else if($1=="total") print $1 }' OUTCAR >temp.f

awk '{if($4=="F"||$4=="T") print $4,$5,$6}' CONTCAR >temp.fix
flag=`wc temp.fix|awk '{print $1}'`
steps=`grep E0 OSZICAR |tail -1 |awk '{print $1}'`
if [ flag != '0' ] ; then
    if [ -f temp.fixx ] ; then rm temp.fixx ; fi
    for i in `seq $steps`;do
        cat temp.fix >>temp.fixx
        echo >>temp.fixx
    done
    paste temp.f temp.fixx >temp.ff
fi

awk  '{ if($1=="total") {print ++i,a;a=0}
        else {
            if($4=="F") x=0; else x=$1;
            if($5=="F") y=0; else y=$2;
            if($6=="F") z=0; else z=$3;
            force=sqrt(x^2+y^2+z^2);
            if(a<force) a=force} }' temp.ff >force.conv

gnuplot <<EOF 
set term dumb
set title 'Energy of each ionic step'
set xlabel 'Ionic steps'
set ylabel 'Energy(eV)'
plot 'temp.e' u 1:5 w l  t "Energy in eV"
set title 'Max Force of each ionic step'
set xlabel 'Ionic steps'
set ylabel 'Force (eV/Angstrom)'
plot 'force.conv' w l t "Force in eV/Angstrom"
EOF

if [ $steps -gt 8 ] ; then
tail -5 force.conv >temp.fff
tail -5 temp.e >temp.ee
gnuplot <<EOF 
set term dumb
set title 'Energy of each ionic step for the last few steps'
set xlabel 'Ionic steps'
set ylabel 'Energy(eV)'
plot 'temp.ee' u 1:5 w l  t "Energy in eV"
set title 'Max Force of each ionic step for the last few steps'
set xlabel 'Ionic steps'
set ylabel 'Force (eV/Angstrom)'
plot 'temp.fff' w l  t "Force in eV/Angstrom"
EOF
rm temp.fff temp.ee
fi

rm temp.e temp.f temp.ff temp.fix temp.fixx