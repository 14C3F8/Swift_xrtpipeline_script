#!/bin/sh

scriptdir=$(pwd)

if [ ! -d "flux_log" ]; then
	mkdir $scriptdir/flux_log
fi

for line in $(cat "$scriptdir/wtlist_new_2.txt")
do
	cd $scriptdir/$line

	if [ ! -e "po_fit.xcm" ]; then
		cp $scriptdir/po_fit.xcm $scriptdir/$line/
	fi

	xspec - po_fit.xcm

	if [ -e "./flux_out.txt" ] && [ -s "./flux_out.txt" ]; then
		cp ./flux_out.txt $scriptdir/flux_log/"$line".txt
		echo "$line" >> $scriptdir/flux_log/success_obsid_list.txt
	elif [ ! -e "./flux_out.txt" ] || [! -s "./flux_out.txt" ]; then
		echo "$line" >> $scriptdir/flux_log/failed_obsid_list.txt
	fi

	cd $scriptdir
done
