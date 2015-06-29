# Solution to Project Euler problem 3: https://projecteuler.net/problem=3
#
# Completed by Adam Blomberg



# function roots finds all the roots of number passed to it
# it checks these roots for additional roots using recursion
# and prints the largest prime factor of the original number

# can be speed improved by dividing out known factors from the test pool.
# and only doing odd numbers (check 2 separately).
# ex:  84 --> 2 * 42
#     2*  42 --> 2* 21
#     2*2* 21 --> 3 * 7
#  84-->2*2*3*7

# store max value of while loop as a variable that shrinks as divide out previous factors
# move i=i+2 to happen if not a factor


def roots( num , maxprime = 1 ) :
	i=2   # test variable to find divisors
	count=0  # count number of roots... will equal zero for primes
	factors=[]  # list to store the factors
	
	# first find all the factors, sqrt(N) is largest factor
	while i <= num**0.5 :
		if num % i == 0 :
			factors.extend([i,num/i])
			count = count + 1
		i=i+1
		
	# if there is are factors, check each of them for factors
	if count > 0 :
		for trial in factors :      
			maxprime = roots(trial,maxprime)   # note recursion here!
	
	# if this factor has no others, then check if it is the
	# largest prime found so far and store it
	elif count == 0 :
		if num > maxprime:
			maxprime = num
	
	# return largest prime so far.... for printing or recursion
	return maxprime
	

#call roots for the test case and the question
var = 13195
print var,"largest prime factor is: ",roots(var)
var = 600851475143
print var,"largest prime factor is: ",roots(var)
