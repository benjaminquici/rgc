#!/usr/bin/bash

usage()
{
echo "Simple script to display .fits cutouts in DS9. Requires DS9 software. "
echo "convolve.sh [-n] [-s]
	-n source_name : Source name
	-s cutout_size : Cutout size in arcseconds (e.g. MIDAS_Jhhmmss-ddmmss_survey_cutoutsize.fits)" 1>&2;
exit 1;
}

# G23_pipeline_path="/mnt/e/Thesis/Paper_II/Pipeline/"
G23_pipeline_path="/mnt/e/Thesis/Paper_II/G23_sample_first_pass/"

while getopts ':n:s:' OPTION
do
	case "$OPTION" in
		n)
			source_name=${OPTARG}
			;;
		s)
			cutout_size=${OPTARG}
			;;
		? | : | h )
			usage
			;;
	esac
done

dir=$G23_pipeline_path"/${source_name}/"

ds9 -zscale -frame lock wcs -crosshair lock wcs ${dir}"/${source_name}_midas_${cutout_size}.fits" ${dir}"/${source_name}_gmrt_${cutout_size}.fits" ${dir}"/${source_name}_emu_${cutout_size}.fits" ${dir}"/${source_name}_glass_${cutout_size}.fits" ${dir}"/${source_name}_viking_${cutout_size}.fits" -cmap cubehelix0