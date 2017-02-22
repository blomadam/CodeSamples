#!/bin/bash

# remove files from previous runs
rm temp* MAIDlookuptable.dat

#initialize variable k as float
k=0.0

# loop over Q2
# bash loops only work with integers, so define range there
# and divide later (into var k defined above)
for (( j=5; j <= 20; j=j+1 ))
do
# convert loop variable to Q2 in MeV^2 using bc
# bc does integer division, so multiply by a decimal instead
	k=`echo "$j * 0.01" | bc`


# loop over W.  Need to use trick above to get non-integer values
	for (( i=1150; i <= 1300; i=i+5 ))
	do

# use curl to DL html file with data for a single W,Q2 point
# wget may work too, but I don't know if the syntax changes
		curl -o tempfile.out "https://maid.kph.uni-mainz.de/cgi-bin/maid1?switch=211&param2=1&param50=3&value35=$k&value36=$i&value37=0&value41=10&value42=180&param99=0&param11=1&param12=1&param13=1&param14=1&param15=1&param16=1&param17=1&param18=1&param19=1&param20=1&param21=1&param22=1&param23=1&param24=1&param25=1&param26=1&value11=1.0&value12=1.0&value13=1.0&value51=1.0&value52=1.0&value53=1.0&value54=1.0&value55=1.0&value56=1.0&value57=1.0&value58=1.0&value59=1.0&value60=1.0&value61=1.0&value62=1.0&value63=1.0&value64=1.0&value65=1.0&value66=1.0&value67=1.0&value68=1.0&value69=1.0&value70=1.0&value71=1.0&value72=1.0&value73=1.0&value74=1.0&value75=1.0&value76=1.0&value77=1.0&value78=1.0&value79=1.0&value80=1.0&value81=1.0&value82=1.0&value83=1.0&value84=1.0"


#sed to copy out only the relevant lines 
		sed -n '/^  (deg)/,/^<\/FONT>/{/^<\/FONT>/d;/^  (deg)/d;p}' \
		   tempfile.out > tempfilemod.out
#cut to delete last four columns, keep only first 7
		cat tempfilemod.out | tr -s ' ' | cut -d ' ' -f 1-7 > tempfilemod2.out

#awk to convert theta to theta_proton = 180-theta_pion
		cat tempfilemod2.out |awk '{printf "%f %f %f %f %f %f\n" \
			,180-$1,$2,$3,$4,$5,$6}' >tempfilemod.out

#sed to reverse line order so theta is increasing again
		sed -n '1!G;h;$p' tempfilemod.out > tempfilemod2.out

#sed to instert Q2 and W at beginning of line
		sed "s/^/$k $i /" tempfilemod2.out > tempfile.out

# cat to append this info to complete table
		cat tempfile.out >> MAIDlookuptable.dat
		mv tempfile.out "temp-$k-$i.out"

	done

done

rm temp*
