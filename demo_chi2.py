# Demo file for "graphical fitting"
# this is used to check sanity of results
# obtained from Minuit fitting algorithm
# and also to explore the stability
# of the fitted result

# The graphs plotted include two cuts on the data
# lb is filtered: 0.005 < lb < 0.015
# and chi2 < 20 or chi2 < 5
# these can be edited in the lookup table loop
# tighter cuts lead to distribution approaching 
# gaussian, as seen with the green and red lines

# production version of code in C++ uses Monte Carlo
# selection of data opints with interpolation algorithm
# between lookup table points
# table has evenly spaced grid opints, so this throw
# will be flat and not introducing bias


from numpy import empty, trim_zeros, full, sqrt, mean, var 
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import sqlite3 as lite
import sys

def asym( sig1, sig2):
	"This takes two cross sections and returns the asymmetry"
	return 100*(sig1-sig2)/(sig1+sig2)
	
def asym_err(sig1,dsig1,sig2,dsig2):
	"This takes two cross sections and their uncertainty ( in pb units) and returns relative uncert for asym"
	return sqrt((pow(dsig1,2)+pow(dsig2,2))*(pow(sig2-sig1,-2)+pow(sig2+sig1,-2))) 


# connect to database of cross sections
con = lite.connect('demo.db')

with con:    

# select all the cross sections with "final" in the notes    
    cur = con.cursor()    
    cur.execute("SELECT Id FROM Notes WHERE Notes LIKE '%final%'")

    rows = cur.fetchone()
    
    ID = int(rows[0])

# Get each sigma, and identify which of the four kinematics it corresponds to
    cur.execute("SELECT * FROM Sigmas WHERE Id = " + str(ID))
    rows = cur.fetchall()

    for row in rows:
		if row[1] < 133:
			if row[2] < 90:
				dat_sigI0  = row[3]
				dat_sigI0_err  = row[4]
			else:
				dat_sigI180   = row[3]
				dat_sigI180_err = row[4]
		if row[1] > 133:
			if row[2] < 90:
				dat_sigII0 = row[3]
				dat_sigII0_err  = row[4]
			else:
				dat_sigII180  = row[3]
				dat_sigII180_err = row[4]
 

# calculate kinematic asymetries from cross sections
dat_asymI = asym(dat_sigI0,dat_sigI180)
dat_asymI_err = asym_err(dat_sigI0,dat_sigI0_err,dat_sigI180,dat_sigI180_err)
dat_asymII = asym(dat_sigII0,dat_sigII180)
dat_asymII_err = asym_err(dat_sigII0,dat_sigII0_err,dat_sigII180,dat_sigII180_err)

# print asymetries to screen for validation with model
print "kinI asym: ", dat_asymI
print "kinII asym: ", dat_asymII




# read the lookup table for the model and store in arrays
# generate arrays of zeros to save time appending to lists
# be careful with large n not to overrun available RAM
# the table this size requires roughly 8*9*n bytes = 4.5 MB
n = 63023

s1_grid =      empty( n, dtype=float)
la_grid =      empty( n, dtype=float)
chi2_grid =    empty( n, dtype=float)
s1_cut =       empty( n, dtype=float)
la_cut =       empty( n, dtype=float)
chi2_cut =     empty( n, dtype=float)
s1_cut2 =      empty( n, dtype=float)
la_cut2 =      empty( n, dtype=float)
chi2_cut2 =    empty( n, dtype=float)
i = 0
j = 0
k = 0

# open table, read parameters on each line and calculate chi2
# apply cuts on lb and on chi2 to make reasonable plots
# to judge location and stability of the minimum
with open('lookup.dat','r') as f:
	for r in f:
		r = r.strip()
		s1 ,la ,lb ,sigI0 ,sigI180 ,asymI ,sigII0 ,sigII180 ,asymII  = map(float, r.split(' '))
		chi2 = pow((sigI0 -dat_sigI0)/dat_sigI0_err,2) + pow((sigI180 -dat_sigI180)/dat_sigI180_err,2) + pow((sigII0 -dat_sigII0)/dat_sigII0_err,2) + pow((sigII180 -dat_sigII180)/dat_sigII180_err,2) + pow((dat_asymI-asymI )/dat_asymI_err,2) +  pow((dat_asymII-asymII )/dat_asymII_err,2) 
		if lb > 0.005 and lb < 0.015:
			la_grid[i] = la
			s1_grid[i] = s1
			chi2_grid[i] = chi2
			i=i+1
			if chi2 < 20 :
				la_cut[j] = la
				s1_cut[j] = s1
				chi2_cut[j] = chi2
				j=j+1
			if chi2 < 5 :
				la_cut2[k] = la
				s1_cut2[k] = s1
				chi2_cut2[k] = chi2
				k=k+1
		if not r: 
				continue


# clear the trailing zeroes
# have not tested to see if this 
# removes time savings from declaring numpy.empty arrays
s1_grid =     trim_zeros(s1_grid  ,'b')
la_grid =     trim_zeros(la_grid  ,'b')
chi2_grid =   trim_zeros(chi2_grid, 'b')
s1_cut =      trim_zeros(s1_cut  , 'b')
la_cut =      trim_zeros(la_cut  , 'b')
chi2_cut =    trim_zeros(chi2_cut , 'b')
s1_cut2 =     trim_zeros(s1_cut2  , 'b')
la_cut2 =     trim_zeros(la_cut2  , 'b')
chi2_cut2 =   trim_zeros(chi2_cut2 , 'b')



print "events passing cuts: ", i,j,k




#plot the data
# keep 0 out of range to keep axes sane  (not required if numpy.trim_zeroes is used!)
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')

ax1.set_ylabel(r"$\chi^2$")
ax1.plot(la_grid,chi2_grid,'k.')
ax1.plot([0, 2], [5, 5], 'r-')
ax1.axis([0, 2, 0, 50])


ax2.plot(s1_grid,chi2_grid,'k.')
ax2.plot([0, 1], [5, 5], 'r-')
ax2.axis([0, 1, 0, 50])


# add weighting scale factor for difference in bin width between S1 and la
ax3.set_ylabel("Counts")
ax3.set_xlabel(r"$\Lambda_\alpha$")
ax3.hist( la_grid, bins=101, range=[0.5,1.5], histtype='step', weights=full(i,8) )
ax3.hist( la_cut,  bins=101, range=[0.5,1.5], histtype='step', weights=full(j,8) )
ax3.hist( la_cut2, bins=101, range=[0.5,1.5], histtype='step', weights=full(k,8) )


ax4.set_xlabel(r"$^\frac{3}{2}S_1$")
ax4.hist( s1_grid, bins=25, range=[0,1], histtype='step' )
ax4.hist( s1_cut,  bins=25, range=[0,1], histtype='step' )
result = ax4.hist( s1_cut2, bins=25, range=[0,1], histtype='step' )
dx = result[1][1] - result[1][0]  # dx is the bin width
scale = len(s1_cut2)*dx  # scale the PDF to match events in histo
gaus = scale*mlab.normpdf(s1_grid,mean(s1_cut2),sqrt(var(s1_cut2)))
ax4.plot(s1_grid, gaus)

plt.suptitle('"Monte Carlo" fitting example')
f.subplots_adjust(hspace=0.05, wspace=0.1)
plt.show()

