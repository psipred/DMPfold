# Note these files can also be found in your CNS install directory e.g. /usr/local/cns_solve_1.3
wget -O gseq.inp "http://cns-online.org/cgi-bin/cns_solve_1.3/cns_view.cgi?&file=inputs/general/generate_seq.inp"
wget -O extn.inp "http://cns-online.org/cgi-bin/cns_solve_1.3/cns_view.cgi?&file=inputs/nmr_calc/generate_extended.inp"
wget -O dgsa.inp "http://cns-online.org/cgi-bin/cns_solve_1.3/cns_view.cgi?&file=inputs/nmr_calc/dg_sa.inp"
patch dgsa.inp dgsa.diff
patch extn.inp extn.diff
patch gseq.inp gseq.diff
