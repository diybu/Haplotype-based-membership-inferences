
for ((i=20000;i<=135700000-150000;i=i+150000));
do
	sfile="10_$(($i - 20000))-$(($i + 150000))"
	#echo $sfile
	java -jar Haploview.jar -nogui -memory 2000 -out /data2/diybu/beaconHaplopAttack/haploblock/hapOut/${sfile} -pedfile /data2/diybu/beaconHaplopAttack/haploblock/genPed/${sfile}.ped -info /data2/diybu/beaconHaplopAttack/haploblock/genPed/${sfile}.info -dprime -ldvalues RSQ -blockoutput GAB

done
