#!/bin/bash
clear
# sed -ne '3p;5p' -s ./**/speed
# find ./**/speed -type f | xargs -i echo {} | sed -r 's#(.\/)(.*)#cat &\|sed  "s:^:file \2 :g"#ge'


path="/work/bjb2chen/740_project/claire_tiara_spectra/water_Claire_2/mctdh/op_nh36Q_5st/op_nh36Q_5st_PBF20_tf10.00/op_nh36Q_5st/output"
# path="/work/bjb2chen/740_project/claire_tiara_spectra/water_Claire_2/mctdh/**/**/**/output"
tail -n 7 ${path} | head -n 1
