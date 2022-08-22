set terminal pngcairo size 2700,1500 enhanced font 'arial, 30' lw 2

set palette defined (0  0.0 0.0 0.5, 1  0.0 0.0 1.0, 2  0.0 0.5 1.0, 3  0.0 1.0 1.0,\
     4  0.5 1.0 0.5, 5  1.0 1.0 0.0, 6  1.0 0.5 0.0, 7  1.0 0.0 0.0, 8  0.5 0.0 0.0)

set grid front lt 2 dt 2 lc rgb "white" lw 2
#set grid mxtics lt 2 dt 2 lc rgb "white" lw 2
set tics out nomirror
set yrange [90:600]
set key off
set xdata time
set timefmt "%Y%m%d\t%H:%M"
LB=sprintf("%s\t0:0", ARG1)
HB=sprintf("%s\t24:0", ARG2)
Step=3600*ARG3
print(LB)
print(HB)
set xrange [LB:HB]
set xtics LB, Step, HB
set xtics timedate
set mxtics 3

set xtics format " "
set ylabel "Height, km"

date=LB[1:6]
day(d)=d<10 ? sprintf("0%d", d) : sprintf("%d", d)
file(date, d, ch, fr, le)=sprintf("%s%s%s_temps_%d_%d_%d.dat", dir, date, d, ch, fr, le);
file1(date, d, fr)=sprintf("%s%s_temps_%d.dat", date, d, fr);
d1=LB[7:8]
d2=HB[7:8]
print(date)
print(d1)
print(d2)
ch=2
fr=ARG4+0
le=700

dir=sprintf("/mnt/data/users/tashlykov/%s/", date)
print(sprintf("%sISdata_%s-%s_%s_%d.png", dir, ARG1, ARG2, ARG4, le))
set output sprintf("%sISdata_%s-%s_%s_%d.png", dir, ARG1, ARG2, ARG4, le)

set multiplot layout 4,1\
    margins 0.1,0.85,0.15,0.95\
    spacing 0.03,0.035

#set cbrange [500:3000]

set ytics format "%g"
#set label "155.5 MHz | 700 us" at graph 0.35,1.075
set logscale cb
set label "S/N [rel.un.]" at graph 0.4,1.075
set cbrange [1e-3:3e0]
plot for [i=d1:d2] file(date, day(i), 0, fr, le) using 1:"Height":"S/N" w image,\

unset logscale cb
# unset label
# set label "P_c_o_n_v" at graph 0.4,1.075
# unset cbrange
# # set cbrange [1e3:3e3]
# plot for [i=d1:d2] file(date, day(i), 0, fr, le) using 1:"Height":"P_conv" w image

unset label
set label "Ne 10^3[cm^-^3]" at graph 0.4,1.075
set cbrange [1e0:1.3e1]
plot for [i=d1:d2] file(date, day(i), 0, fr, le) using 1:"Height":"Ne" w image


unset label
set label "Te [K]" at graph 0.4,1.075
set cbrange [100:3500]
plot for [i=d1:d2] file(date, day(i), ch, fr, le) using 1:"Height":"Te_acf" w image

# unset label
# set cbrange [1e5:3e5]
# plot for [i=d1:d2] file(date, day(i), ch, fr, le) using 1:"Height":"P_exp" w image

set xtics format "%d/%m"
set xlabel sprintf("%s", LB[1:4])
unset label
set label "Tr" at graph 0.4,1.075
set cbrange [1:3]
plot for [i=d1:d2] file(date, day(i), ch, fr, le) using 1:"Height":"Tr_acf" w image

# plot for [i=d1:d2] file(date, day(i), ch, fr, le) using 1:"Height":"P_model" w image
