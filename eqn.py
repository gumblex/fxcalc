# Source code for the Quartic/Cubic/Quadratic Equation Calculator (quartic_cubic_quadratic.js)
# --> CURRENT VERSION: V3.5 <--
# IF YOU WANT UPDATES, DO NOT HOST THIS SOURCE FILE ON YOUR SERVER!
# Created By Brian Kieffer
# http://www.freewebs.com/brianjs

# Converted By Gumble
import math
from fparser import *
from decimal import *

D = Decimal
sqrt = lambda x : x ** 0.5

def cos(x):
	"""Return the cosine of x as measured in radians.

	The Taylor series approximation works best for a small value of x.
	For larger values, first compute x = x % (2 * pi).

	>>> print(cos(Decimal('0.5')))
	0.8775825618903727161162815826
	>>> print(cos(0.5))
	0.87758256189
	>>> print(cos(0.5+0j))
	(0.87758256189+0j)

	"""
	getcontext().prec += 2
	i, lasts, s, fact, num, sign = 0, 0, 1, 1, 1, 1
	while s != lasts:
		lasts = s
		i += 2
		fact *= i * (i-1)
		num *= x * x
		sign *= -1
		s += num / fact * sign
	getcontext().prec -= 2
	return +s

def sin(x):
	"""Return the sine of x as measured in radians.

	The Taylor series approximation works best for a small value of x.
	For larger values, first compute x = x % (2 * pi).

	>>> print(sin(Decimal('0.5')))
	0.4794255386042030002732879352
	>>> print(sin(0.5))
	0.479425538604
	>>> print(sin(0.5+0j))
	(0.479425538604+0j)

	"""
	getcontext().prec += 2
	#if abs(x) > 2 * pi:
		#x = x % (2 * pi)
	i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
	while s != lasts:
		lasts = s
		i += 2
		fact *= i * (i-1)
		num *= x * x
		sign *= -1
		s += num / fact * sign
	getcontext().prec -= 2
	return +s

tan = lambda x : sin(x) / cos(x)
cot = lambda x : cos(x) / sin(x)
sec = lambda x : 1 / cos(x)
csc = lambda x : 1 / sin(x)

toplace = 12
def convert(number):
	if (number.toString().indexOf("e") != -1 or abs(parseFloat(number)) > 1E+6):
		number = parseFloat(number)
		if (abs(number) < 1):
			place = number.toExponential(10)
			place1 = place
			length = place.length
			place3 = place.toString().indexOf("e-")
			place2 = place1.toString().substring(place3 + 1, length)
			place4 = place1.toString().substring(0, place3)
			str = "[" + eval(place4) + " x 10<sup>" + place2 + "</sup>] "
			return str
		elif (abs(number) > 10000):
			place = number.toExponential(6)
			place1 = place
			length = place.length
			place3 = place.toString().indexOf("e+")
			place2 = place1.toString().substring(place3 + 2, length)
			place4 = place1.toString().substring(0, place3)
			str = "[" + eval(place4) + " x 10<sup>" + place2 + "</sup>] "
			return str
		else:
			str = number
			return str
		
	else:
		return number
	

def calcmult(a2, b2, c2, d2, e2):
	real = a2 * c2 - b2 * d2
	img = b2 * c2 + a2 * d2
	if (e2 == 0):
		return real
	else:
		return img
	

def isquareroot(a1, b1, n1):
	y = sqrt((a1 * a1) + (b1 * b1))
	y1 = sqrt((y - a1) / 2)
	x1 = b1 / (2 * y1)
	if (n1 == 0):
		return x1
	else:
		return y1
	

def n_rt(x, y):
	if (x < 0):
		if (math.floor(y / 2) == (y / 2).toFixed(12)):
			rootone = "Input error detected. Using n_rt() for negative numbers with y being even is not allowed."
			roottwo = "Input error detected. Using n_rt() for negative numbers with y being even is not allowed."
			rootthree = "Input error detected. Using n_rt() for negative numbers with y being even is not allowed."
			rootfour = "Input error detected. Using n_rt() for negative numbers with y being even is not allowed."
			scriptexec = 0
			return "x"
		else:
			return -pow(-x, (1 / y))
		
	else:
		return pow(x, (1 / y))
	

def log(x, y):
	return math.log(x) / math.log(y)
	

def arcsin(x):
	if (x < -1 or x > 1):
		rootone = "Input error detected. Using arcsin() for numbers out of the range of -1 and 1 is not allowed."
		roottwo = "Input error detected. Using arcsin() for numbers out of the range of -1 and 1 is not allowed."
		rootthree = "Input error detected. Using arcsin() for numbers out of the range of -1 and 1 is not allowed."
		rootfour = "Input error detected. Using arcsin() for numbers out of the range of -1 and 1 is not allowed."
		scriptexec = 0
	else:
		return math.asin(x)
	

def arctan(x):
	return math.atan(x)

def eqn(self, a=0, b=0, c=0, d=0, e=0):
	'''Solve quartic/cubic/quadratic equations: ax^4 + bx^3 + cx^2 + dx + e = 0'''
	
	N = lambda x : self.ntype(x, 'frac')

	a3, b3, c3, d3, e3 = N(0), N(0), N(0), N(0), N(0)
	a = N(a)
	b = N(b)
	c = N(c)
	d = N(d)
	e = N(e)
	rootone, roottwo, rootthree, rootfour = (None,)*4
	scriptexec = 1
	specialexec = 0
	# Catch user errors
	if a == 0 and b == 0 and c != 0 and d == 0 and e == 0:
		rootone = 0
		roottwo = 0
		scriptexec = 0
		return (0,0)
	elif a == 0 and b == 0 and c == 0:
		if d != 0:
			rootone = -e / d
		else:
			raise MathERROR
		scriptexec = 0
		return (rootone,)
	
	if (e == 0 and a != 0):
		specialexec = 2
		scriptexec = 2
	
	if (d == 0 and e == 0):
		specialexec = 1
		scriptexec = 3
	
	if a == 0 and b == 0:
		scriptexec = 3 # Recongize as a Quadratic Equation
	elif a == 0:
		scriptexec = 2 # Recongize as a Cubic Equation
	
	if scriptexec:
		if (scriptexec == 1):
			discrim = (b * b) - (4 * a * c)
			x1 = (-1 * b + sqrt(discrim)) / (2 * a)
			x2 = (-1 * b - sqrt(discrim)) / (2 * a)
			alreadydone2 = 0
			divide = b / a
			alreadydone3 = 0
			a2, b2, c2, d2, e2 = N(0), N(0), N(0), N(0), N(0)
			dontalert = 0
			ripart = 0
			qipart = 0
			# Extract X^4 Coefficent
			aq = a
			aq2 = aq # Keeps Orignial AQ value
			# Extract X^3 Coefficent
			bq = b
			bq2 = bq # Keeps Orignial BQ Value
			# Extract X^2 Coefficent
			cq = c
			# Extract X Coefficent
			dq = d
			# Extract Constant
			eq = e
			# Define Perfect Quartic Vars
			perfect = 0
			perfectbiquadratic = 0
			# The Bi-Quadratic 2 Perfect Squares that are negative test
			if (cq * cq - 4 * aq * eq == 0 and (cq / aq) > 0 and bq == 0 and dq == 0):
				perfectbiquadratic = 1
			# Divide Equation by the X^4 Coefficent to make equation in the form of X^4 + AX^3 + BX^2 + CX + D
			bq /= aq
			cq /= aq
			dq /= aq
			eq /= aq
			aq = 1
			f2 = cq - (3 * bq * bq / 8)
			g2 = dq + (bq * bq * bq / 8) - (bq * cq / 2)
			h2 = eq - (3 * bq * bq * bq * bq / 256) + (bq * bq * (cq / 16)) - (bq * dq / 4)
			a = 1
			b = f2 / 2
			c = (f2 * f2 - (4 * h2)) / 16
			d = -1 * ((g2 * g2) / 64)
			if (b == 0 and c == 0 and d == 0):
				perfect = 1
			# Cubic routine starts here.....
			f = (((3 * c) / a) - ((b * b) / (a * a))) / 3
			g = (((2 * b * b * b) / (a * a * a)) - ((9 * b * c) / (a * a)) + ((27 * d) / a)) / 27
			h = ((g * g) / 4) + ((f * f * f) / 27)
			z = 1 / 3
			alreadydone2 = 0
			ipart = 0
			p, q, r, s = 0, 0, 0, 0
			if (h <= 0):
				_exec = 2
				i = sqrt(((g * g) / 4) - h)
				j = i ** z
				k = N(math.acos(-(g / (2 * i))))
				l = -j
				m = cos(k / 3)
				n = sqrt(N(3)) * N(sin(k / 3))
				p = -(b / (3 * a))
				xoneterm = (2 * j) * N(cos(k / 3)) - (b / (3 * a))
				xtwoterm = l * (m + n) + p
				xthreeterm = l * (m - n) + p
			if (h > 0):
				_exec = 1
				R = -(g / 2) + sqrt(h)
				if (R < 0):
					S = -((-1 * R) ** z)
				else:
					S = R ** z
				
				T = -(g / 2) - sqrt(h)
				if (T < 0):
					U = -((-T) ** z)
				else:
					U = T ** z
				
				xoneterm = (S + U) - (b / (3 * a))
				xtwoterm = (-(S + U) / 2) - (b / (3 * a))
				ipart = ((S - U) * sqrt(3)) / 2
				xthreeterm = xtwoterm
			
			if (f == 0 and g == 0 and h == 0):
				if ((d / a) < 0):
					xoneterm = (-(d / a)) ** z
					xtwoterm = xoneterm
					xthreeterm = xoneterm
				else:
					xoneterm = -((d / a) ** z)
					xtwoterm = xoneterm
					xthreeterm = xoneterm
			
			# ....and ends here.
			# if (abs(ipart) < 5E-7):
				# ipart = 0
			# if (abs(xoneterm) < 5E-7):
				# xoneterm = 0
			# if (abs(xtwoterm) < 5E-7):
				# xtwoterm = 0
			# if (abs(xthreeterm) < 5E-7):
				# xthreeterm = 0
			
			# Return to solving the Quartic.
			if (ipart == 0 and xoneterm == 0):
				alreadydone2 = 1
				p2 = sqrt(xtwoterm)
				q = sqrt(xthreeterm)
				r = -g2 / (8 * p2 * q)
				s = bq2 / (4 * aq2)
			
			if (ipart == 0 and xtwoterm == 0 and alreadydone2 == 0 and alreadydone2 != 1):
				alreadydone2 = 2
				p2 = sqrt(xoneterm)
				q = sqrt(xthreeterm)
				r = -g2 / (8 * p2 * q)
				s = bq2 / (4 * aq2)
			
			if (ipart == 0 and xthreeterm == 0 and alreadydone2 == 0 and alreadydone2 != 1 and alreadydone2 != 2):
				alreadydone2 = 3
				p2 = sqrt(xoneterm)
				q = sqrt(xtwoterm)
				r = -g2 / (8 * p2 * q)
				s = bq2 / (4 * aq2)
			
			if (alreadydone2 == 0 and ipart == 0):
				if (xthreeterm < 0):
					alreadydone2 = 4
					p2 = sqrt(xoneterm)
					q = sqrt(xtwoterm)
					r = -g2 / (8 * p2 * q)
					s = bq2 / (4 * aq2)
				else:
					alreadydone2 = 5
					p2 = sqrt(xoneterm)
					q = sqrt(xthreeterm)
					r = -g2 / (8 * p2 * q)
					s = bq2 / (4 * aq2)
				
			if (ipart != 0):
				if (xoneterm < 0):
					xoneterm /= -1
					p2 = 0
					p2ipart = sqrt(xoneterm)
					q = 0
					qipart = 0
					mult = 0
					r = 0
					s = bq2 / (4 * aq2)
					approx = 1
					xoneterm /= -1
				else:
					p2 = isquareroot(xtwoterm, ipart, 0)
					p2ipart = isquareroot(xtwoterm, ipart, 1)
					q = p2
					qipart = -p2ipart
					mult = calcmult(p2, p2ipart, q, qipart, 0)
					r = -g2 / (8 * mult)
					s = bq2 / (4 * aq2)
				
			if (ipart == 0 and xtwoterm < 0 and xthreeterm < 0):
				xtwoterm /= -1
				xthreeterm /= -1
				p2 = 0
				q = 0
				p2ipart = sqrt(xtwoterm)
				qipart = sqrt(xthreeterm)
				mult = calcmult(p2, p2ipart, q, qipart, 0)
				r = -g2 / (8 * mult)
				s = bq2 / (4 * aq2)
				ipart = 1
			
			if (xoneterm > 0 and xtwoterm < 0 and xthreeterm == 0 and ipart == 0):
				alreadydone5 = 1
				xtwoterm /= -1
				p2 = sqrt(xoneterm)
				q = 0
				p2ipart = 0
				qipart = sqrt(xtwoterm)
				mult = calcmult(p2, p2ipart, q, qipart, 0)
				mult2 = calcmult(p2, p2ipart, q, qipart, 1)
				r = -g2 / (8 * mult)
				if (mult2 != 0):
					ripart = g2 / (8 * mult2)
					r = 0
				
				s = bq2 / (4 * aq2)
				ipart = 1
			
			if (xoneterm > 0 and xtwoterm < 0 and xthreeterm == 0 and ipart == 0 and alreadydone5 != 1):
				alreadydone5 = 2
				xtwoterm /= -1
				p2 = sqrt(xoneterm)
				q = 0
				p2ipart = 0
				qipart = sqrt(xtwoterm)
				mult = calcmult(p2, p2ipart, q, qipart, 0)
				mult2 = calcmult(p2, p2ipart, q, qipart, 1)
				r = -g2 / (8 * mult)
				if (mult2 != 0):
					ripart = g2 / (8 * mult2)
					r = 0
				
				s = bq2 / (4 * aq2)
				ipart = 1
			
			if (xoneterm > 0 and xtwoterm < 0 and xthreeterm == 0 and ipart == 0 and alreadydone5 != 1 and alreadydone5 != 2):
				alreadydone5 = 3
				xtwoterm /= -1
				p2 = sqrt(xoneterm)
				q = 0
				p2ipart = 0
				qipart = sqrt(xtwoterm)
				mult = calcmult(p2, p2ipart, q, qipart, 0)
				mult2 = calcmult(p2, p2ipart, q, qipart, 1)
				r = -g2 / (8 * mult)
				if (mult2 != 0):
					ripart = g2 / (8 * mult2)
					r = 0
				
				s = bq2 / (4 * aq2)
				ipart = 1
			
			if (xtwoterm == 0 and xthreeterm == 0 and ipart == 0):
				if (xoneterm < 0):
					ipart = 1
					p2 = 0
					p2ipart = sqrt(-xoneterm)
					q = 0
					r = 0
					s = bq2 / (4 * aq2)
				else:
					p2 = sqrt(xoneterm)
					q = 0
					r = 0
					s = bq2 / (4 * aq2)
				
			
			if (xoneterm == 0 and xtwoterm == 0 and ipart == 0):
				if (xthreeterm < 0):
					ipart = 1
					p2 = 0
					p2ipart = sqrt(-xthreeterm)
					q = 0
					r = 0
					s = bq2 / (4 * aq2)
				else:
					p2 = sqrt(xthreeterm)
					q = 0
					r = 0
					s = bq2 / (4 * aq2)
				
			
			if (xoneterm.toFixed(8) == 0 and xthreeterm.toFixed(8) == 0 and ipart == 0):
				if (xtwoterm < 0):
					ipart = 1
					p2 = 0
					p2ipart = sqrt(-xtwoterm)
					q = 0
					r = 0
					s = bq2 / (4 * aq2)
				else:
					p2 = sqrt(xtwoterm)
					q = 0
					r = 0
					s = bq2 / (4 * aq2)
				
			
			sum = abs(xtwoterm) + abs(ipart)
			if (xtwoterm == 0 and ipart == 0):
				sum = abs(xthreeterm) + abs(ipart)
			
			if (sum <= 0 and approx != 1 and ipart != 0):
				if (xoneterm > 0):
					p2 = sqrt(xoneterm)
					q = 0
					p2ipart = 0
					qipart = 0
					mult = 0
					r = 0
					s = bq2 / (4 * aq2)
				else:
					p2 = 0
					q = 0
					p2ipart = sqrt(-xoneterm)
					qipart = 0
					mult = 0
					r = 0
					s = bq2 / (4 * aq2)
				
			# No need for float workaround.
			# if (p2 < 5.0E-8):
				# p2 = 0
			# if (q < 5.0E-8):
				# q = 0
			
			if (ipart == 0 and h > 0 and abs(xoneterm) > 0 and xtwoterm != 0):
				if (xoneterm > 0 and xtwoterm > 0):
					p2 = sqrt(xoneterm)
					p2ipart = 0
					q = sqrt(xtwoterm)
					qipart = 0
					r = -g2 / (8 * p2 * q)
					s = bq2 / (4 * aq2)
				elif (xoneterm < 0 and xtwoterm < 0):
					p2 = 0
					p2ipart = sqrt(-xoneterm)
					q = 0
					qipart = sqrt(-xtwoterm)
					mult = calcmult(p2, p2ipart, q, qipart, 0)
					mult2 = calcmult(p2, p2ipart, q, qipart, 1)
					r = -g2 / (8 * mult)
					if (mult2 != 0):
						ripart = g2 / (8 * mult2)
						r = 0
					
					s = bq2 / (4 * aq2)
				elif (xoneterm < 0):
					p2 = 0
					p2ipart = sqrt(-xoneterm)
					q = sqrt(xtwoterm)
					qipart = 0
					mult = calcmult(p2, p2ipart, q, qipart, 0)
					mult2 = calcmult(p2, p2ipart, q, qipart, 1)
					r = -g2 / (8 * mult)
					if (mult2 != 0):
						ripart = g2 / (8 * mult2)
						r = 0
					
					s = bq2 / (4 * aq2)
				else:
					p2 = sqrt(xoneterm)
					p2ipart = 0
					q = 0
					qipart = sqrt(-xtwoterm)
					mult = calcmult(p2, p2ipart, q, qipart, 0)
					mult2 = calcmult(p2, p2ipart, q, qipart, 1)
					r = -g2 / (8 * mult)
					if (mult2 != 0):
						ripart = g2 / (8 * mult2)
						r = 0
					
					s = bq2 / (4 * aq2)
				
			
			# Now output the answers!
			if (ipart == 0):
				rootone = N(p2 + q + r - s)
				roottwo = N(p2 - q - r - s)
				rootthree = N(-p2 + q - r - s)
				rootfour = N(-p2 - q + r - s)
			
			if (perfect == 1):
				rootone = N(-bq / 4)
				roottwo = N(-bq / 4)
				rootthree = N(-bq / 4)
				rootfour = N(-bq / 4)
			
			if (ipart != 0):
				x1imag = N(p2ipart + qipart + ripart)
				x2imag = N(p2ipart - qipart - ripart)
				x3imag = N(-p2ipart + qipart - ripart)
				x4imag = N(-p2ipart - qipart + ripart)
				x1real = N(p2 + q + r - s)
				x2real = N(p2 - q - r - s)
				x3real = N(-p2 + q - r - s)
				x4real = N(-p2 - q + r - s)
				x1imagc = p2ipart + qipart + ripart
				x2imagc = p2ipart - qipart - ripart
				x3imagc = -p2ipart + qipart - ripart
				x4imagc = -p2ipart - qipart + ripart
				x1realc = p2 + q + r - s
				x2realc = p2 - q - r - s
				x3realc = -p2 + q - r - s
				x4realc = -p2 - q + r - s
				# if (abs(x1imagc) < 4E-8):
					# x1imag = 0
				# if (abs(x2imagc) < 4E-8):
					# x2imag = 0
				# if (abs(x3imagc) < 4E-8):
					# x3imag = 0
				# if (abs(x4imagc) < 4E-8):
					# x4imag = 0
				
				rootone = cfrac(x1real, x1imag)
				roottwo = cfrac(x4real, x4imag)
				rootthree = cfrac(x2real, x2imag)
				rootfour = cfrac(x3real, x3imag)
				
			if (perfectbiquadratic == 1):
				perfect = sqrt(cq / 2)
				if (perfect == 1):
					rootone = cfrac(0, 1)
					roottwo = cfrac(0, 1)
					rootthree = cfrac(0, -1)
					rootfour = cfrac(0, -1)
				else:
					rootone = cfrac(0, perfect)
					roottwo = cfrac(0, perfect)
					rootthree = cfrac(0, -perfect)
					rootfour = cfrac(0, -perfect)
			
			if (cq == 0 and dq == 0 and eq == 0):
				rootone = -divide
				roottwo = 0
				rootthree = 0
				rootfour = 0
			
			# if (approx == 2):
				# rootone += " [Approx]"
				# roottwo += " [Approx]"
				# rootthree += " [Approx]"
				# rootfour += " [Approx]"
			
			# if (rootone == 0 and roottwo == 0 and rootthree == 0 and rootfour == 0):
				# rootone = 0
				# roottwo = 0
				# rootthree = 0
				# rootfour = 0
			
			# if (rootone == "NaN" and roottwo == "NaN" and rootthree == "NaN" and rootfour == "NaN"):
				# rootone = 0
				# roottwo = 0
				# rootthree = 0
				# rootfour = 0
		
		elif (scriptexec == 2):
			if (specialexec != 2):
				a = b
				b = c
				c = d
				d = e
			
			divide = -b / a
			a2 = a
			b2 = b
			c2 = c
			d2 = d
			if (b < 0):
				b2 /= -1
			if (c < 0):
				c2 /= -1
			if (d < 0):
				d2 /= -1
			if (b == 0):
				b2 = ""
				x2show = 0
			if (c == 0):
				c2 = ""
				xshow = 0
			if (d == 0):
				d2 = ""
			if (a == 1):
				a2 = ""
			if (a == -1):
				a2 = "-"
			if (abs(b) == 1):
				b2 = ""
			if (abs(c) == 1):
				c2 = ""
			
			
			f = (((3 * c) / a) - ((b * b) / (a * a))) / 3
			g = (((2 * b * b * b) / (a * a * a)) - ((9 * b * c) / (a * a)) + ((27 * d) / a)) / 27
			h = ((g * g) / 4) + ((f * f * f) / 27)
			z = 1 / 3
			if (discrim == 0 and d == 0):
				xoneterm = 0
				xtwoterm = ((-1 * b) / (2 * a))
			
			if (h <= 0):
				_exec = 2
				i = sqrt(((g * g) / 4) - h)
				j = pow(i, z)
				k = math.acos(-1 * (g / (2 * i)))
				l = -1 * j
				m = cos(k / 3)
				n = sqrt(3) * sin(k / 3)
				p = (b / (3 * a)) * -1
				xoneterm = (2 * j) * cos(k / 3) - (b / (3 * a))
				xtwoterm = l * (m + n) + p
				xthreeterm = l * (m - n) + p
			
			if (h > 0):
				_exec = 1
				R = (-1 * (g / 2)) + sqrt(h)
				if (R < 0):
					S = -1 * (pow((-1 * R), z))
				else:
					S = pow(R, z)
				
				T = (-1 * (g / 2)) - sqrt(h)
				if (T < 0):
					U = -1 * (pow((-1 * T), z))
				else:
					U = pow(T, z)
				
				xoneterm = (S + U) - (b / (3 * a))
				xtwoterm = (-1 * (S + U) / 2) - (b / (3 * a))
				ipart = ((S - U) * sqrt(3)) / 2
				xthreeterm = xtwoterm
			
			if (f == 0 and g == 0 and h == 0):
				if ((d / a) < 0):
					xoneterm = (pow((-1 * (d / a)), z))
				else:
					xoneterm = -1 * (pow((d / a), z))
				
			rootone = N(xoneterm)
			if (h > 0):
				if (xtwoterm == 0 and xthreeterm == 0):
					if (ipart == 1):
						roottwo = cfrac(0,1)
						rootthree = cfrac(0,-1)
					elif (ipart == 0):
						roottwo = 0
						rootthree = 0
					else:
						roottwo = cfrac(0,ipart)
						rootthree = cfrac(0,-ipart)
				else:
					if (ipart == 1):
						roottwo = cfrac(0,xtwoterm)
						rootthree = cfrac(0,-xthreeterm)
					else:
						roottwo = cfrac(xtwoterm,ipart)
						rootthree = cfrac(xthreeterm,-ipart)

			if (h < 0):
				roottwo = N(xtwoterm)
				rootthree = N(xthreeterm)
			
			if (discrim == 0 and d == 0):
				roottwo = N(xtwoterm)
				rootthree = roottwo
			
			if (a != 0 and b != 0 and c == 0 and d == 0):
				alreadydone = 1
				rootone = 0
				roottwo = -b / a
				rootthree = 0
			
			if (b == 0 and c == 0 and d == 0):
				alreadydone = 1
				roottwo = 0
				rootthree = 0
			
			if (h == 0 and alreadydone != 1):
				roottwo = N(xtwoterm)
				rootthree = N(xthreeterm)
			
			if (f == 0 and g == 0 and h == 0):
				alreadydone2 = 1
				roottwo = xoneterm
				rootthree = xoneterm
			
			if (c == 0 and d == 0 and alreadydone != 1):
				alreadydone = 2
				if (x2 == 0):
					rootone = N(x1)
					roottwo = 0
					rootthree = 0
				else:
					alreadydone3 = 1
					rootone = 0
					roottwo = N(x2)
					rootthree = 0
			if (specialexec == 2):
				rootfour = 0
		
		elif (scriptexec == 3):
			a2 = a
			b2 = b
			c2 = c
			d2 = d
			if (c < 0):
				c2 /= -1
			if (d < 0):
				d2 /= -1
			if (c == 0):
				c2 = ""
				xshow = 0
			if (d == 0):
				d2 = ""
			if (a == 1):
				a2 = ""
			if (a == -1):
				a2 = "-"
			if (abs(b) == 1):
				b2 = ""
			if (abs(c) == 1):
				c2 = ""
			
			if (specialexec == 1):
				xoneterm1 = (-1 * b + sqrt(b * b - 4 * a * c)) / (2 * a)
				xtwoterm1 = (-1 * b - sqrt(b * b - 4 * a * c)) / (2 * a)
				sqrtterm = (b * b) - (4 * a * c)
				realpart = (-1 * b / (2 * a))
				imag = (sqrt(-sqrtterm) / (2 * a))
			else:
				xoneterm1 = (-1 * d + sqrt(d * d - 4 * c * e)) / (2 * c)
				xtwoterm1 = (-1 * d - sqrt(d * d - 4 * c * e)) / (2 * c)
				sqrtterm = (d * d) - (4 * c * e)
				realpart = -1 * d / (2 * c)
				imag = sqrt(-sqrtterm) / (2 * c)
			
			rootone = xoneterm1
			roottwo = xtwoterm1
			if specialexec:
				rootthree = 0
				rootfour = 0
			
			if (specialexec == 0):
				if (sqrtterm < 0):
					if (d == 0 or realpart == 0):
						sqrtterm = -1 * sqrtterm
						rootone = cfrac(0, abs(imag))
						roottwo = cfrac(0, -abs(imag))
					else:
						sqrtterm = -1 * sqrtterm
						rootone = cfrac(realpart, abs(imag))
						roottwo = cfrac(realpart, -abs(imag))
						if (((sqrt(sqrtterm)) / (2 * c)) == 1):
							rootone = cfrac(realpart, 1)
							roottwo = cfrac(realpart, -1)
			else:
				if (sqrtterm < 0):
					if (b == 0):
						sqrtterm = -1 * sqrtterm
						rootone = cfrac(0, imag)
						roottwo = cfrac(0, -imag)
					else:
						sqrtterm = -1 * sqrtterm
						rootone = cfrac(realpart, abs(imag))
						roottwo = cfrac(realpart, -abs(imag))
						if (((sqrt(sqrtterm)) / (2 * a)) == 1):
							rootone = cfrac(realpart, 1)
							roottwo = cfrac(realpart, -1)
			if (specialexec == 0):
				if (xoneterm1 == xtwoterm1):
					rootone = cfrac(xoneterm1)
					roottwo = cfrac(xoneterm1)
			if (realpart == 0 and imag == 1):
				rootone = cfrac(0,1)
				roottwo = cfrac(0,-1)
	return tuple(filter(lambda x: x != None, (rootone, roottwo, rootthree, rootfour)))

print(eqn(Parser(),1,2,3,4,5))