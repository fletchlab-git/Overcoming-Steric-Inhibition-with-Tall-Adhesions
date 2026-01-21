set terminal jpeg font ",20";set output "Data.jpg";
set xlabel "MCSweeps"; set ylabel "Relative population [%]";
set xrange [100:];set xtics 10000;
set yrange [:1.1];
set format y "%.1f"
set key at screen 0.55, 0.9
unset arrow; set arrow from 9000, 0.0 to 9000,0.8 nohead lw 2;set arrow from 17000,0.0 to 17000,0.8 nohead lw 2;
p "/Data/FletcherMembrane/DataFiles/Data_Nodes_1600_NumPTypes_2_kc_100.00_A_25600.0_Run_3_.dat" u 1:( ($14 + $17)/($12 + 0.00000001)) title "Long NBP" w l lw 3, "" u 1:( ($15 + $18)/($12 + 0.0000001)) w l lw 3 title "Short BP"

set autoscale
set output "Data2.jpg";
set xlabel "MCSweeps"; set ylabel "z/L_{long}";set xrange [1000:];
unset arrow; set arrow from 1000, 2.0 to 40000,2.0 nohead lw 2;set arrow from 1000, 1.0 to 40000,1.0 nohead lw 2;
set arrow from 9000, 0.0 to 9000,2.0 dt 0 lw 2;set arrow from 17000, 0.0 to 17000,1.3 nohead dt 0 lw 2;
set arrow from 22000, 0.0 to 22000,1.0 nohead dt 0 lw 2;

p "/Data/FletcherMembrane/DataFiles/Data_Nodes_1600_NumPTypes_2_kc_100.00_A_25600.0_Run_3_.dat" u 1:( $11/16. ) title "z" w l lw 2

set output "Data3.jpg"
set sample 1000;
set xlabel "z/L"; set ylabel "{/Symbol bD}F";set xrange [0:3.0];set xtics autofreq;set yrange [:2.1]
unset ytics;

set style rect fc lt -1 fs solid 0.15 noborder
set obj rect from 2.0, graph 0 to 3.1, graph 1;
set label "Noninteracting" at 2.12, 1.8 font ",16" rotate by 0

set style rect fc lt 1 fs solid 0.3 noborder
set obj rect from 1.0, graph 0 to 2.0, graph 1;
set label "Interdigitating" at 1.15, 1.8 font ",16" rotate by 0

set style rect fc lt 2 fs solid 0.2 noborder
set obj rect from 0.0, graph 0 to 1.0, graph 1;
set label "Steric depletion" at 0.08, 1.8 font ",16" rotate by 0

set label "{/Symbol Dm}_1={/Symbol m}(2c_{NB}) - {/Symbol m}(c_{NB})" at 0.25, 0.3 font ",16" rotate by 0
set arrow from 1.5, 0.0 to 1.5, 2.0/3.0 heads lw 1.5

set label "{/Symbol Dm}_2={/Symbol m}(c_{NB}) - {/Symbol m}(0)" at 1.05, 1.1 font ",16" rotate by 0
set arrow from 0.8, 2.0/3.0 to 0.8, 1.0 + 2.0/3.0 heads lw 1.5

p erfc(15.0*(x-2.0))/3 + erfc(15.0*(x-1.0))/2 lw 3 title ""


