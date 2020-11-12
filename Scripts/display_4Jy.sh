#!/bin/bash

src=$1

if [[ -e FIRST_${src}.fits ]]
then
    ds9 -tile yes -tile mode column -cmap BB -grid yes -grid load good.grd -colorbar vertical -lock frame wcs -scale limits -0.0005 0.01 VLASS_${src}.fits TGSS_${src}.fits -scale limits -0.1 1 FIRST_${src}.fits -log -scale limits -0.0001 0.1 -rgb -red  WISE12um_${src}.fits -linear -scale mode 99.5 -green WISE4.6um_${src}.fits -linear -scale mode 99.5 -blue WISE3.4um_${src}.fits -linear -scale mode 99.5 -zoom 2 -saveimage ${src}.png -region load all VLASS_${src}.reg

#-crosshair ${src:1:2}:${src:3:2}:${src:5:5} ${src:10:3}:${src:13:2}:${src:15:4} wcs fk5 -crosshair lock wcs
else
    ds9 -tile yes -tile mode column -cmap BB -grid yes -grid load good.grd -colorbar vertical -lock frame wcs -scale limits -0.0005 0.01 VLASS_${src}.fits TGSS_${src}.fits -scale limits -0.1 1 -rgb -red  WISE12um_${src}.fits -scale mode 99.5 -green WISE4.6um_${src}.fits -scale mode 99.5 -blue WISE3.4um_${src}.fits -scale mode 99.5 -zoom 2 -saveimage ${src}.png -region load all VLASS_${src}.reg
fi
