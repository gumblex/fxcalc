#Source code for the Quartic/Cubic/Quadratic Equation Calculator (quartic_cubic_quadratic.js)
#--> CURRENT VERSION: V3.5 <--
#IF YOU WANT UPDATES, DO NOT HOST THIS SOURCE FILE ON YOUR SERVER!
#Created By Brian Kieffer
#http://www.freewebs.com/brianjs
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
		y = math.sqrt((a1 * a1) + (b1 * b1))
		y1 = math.sqrt((y - a1) / 2)
		x1 = b1 / (2 * y1)
		if (n1 == 0):
			return x1
		else:
			return y1
		
	
	def extractcoefficents(execexample):
		scriptexec = 1
		specialexec = 0
		def sqrt(x):
			if (parseFloat(x) < 0):
				document.getElementById('rootone').innerHTML = "Input error detected. Using sqrt() for negative numbers is not allowed."
				document.getElementById('roottwo').innerHTML = "Input error detected. Using sqrt() for negative numbers is not allowed."
				document.getElementById('rootthree').innerHTML = "Input error detected. Using sqrt() for negative numbers is not allowed."
				document.getElementById('rootfour').innerHTML = "Input error detected. Using sqrt() for negative numbers is not allowed."
				scriptexec = 0
				return "x"
			else:
				return math.sqrt(x)
			
		
		def n_rt(x, y):
			if (x < 0):
				if (math.floor(y / 2) == (y / 2).toFixed(12)):
					document.getElementById('rootone').innerHTML = "Input error detected. Using n_rt() for negative numbers with y being even is not allowed."
					document.getElementById('roottwo').innerHTML = "Input error detected. Using n_rt() for negative numbers with y being even is not allowed."
					document.getElementById('rootthree').innerHTML = "Input error detected. Using n_rt() for negative numbers with y being even is not allowed."
					document.getElementById('rootfour').innerHTML = "Input error detected. Using n_rt() for negative numbers with y being even is not allowed."
					scriptexec = 0
					return "x"
				else:
					return -pow(-x, (1 / y))
				
			else:
				return pow(x, (1 / y))
			
		
		def log(x, y):
			return math.log(x) / math.log(y)
		
		def PI():
			return math.PI
		
		def E():
			return math.E
		
		def sin(x):
			return math.sin(x)
		
		def cosine(x):
			return math.cos(x)
		
		def tan(x):
			return math.tan(x)
		
		def arccos(x):
			if (x < -1 or x > 1):
				document.getElementById('rootone').innerHTML = "Input error detected. Using arccos() for numbers out of the range of -1 and 1 is not allowed."
				document.getElementById('roottwo').innerHTML = "Input error detected. Using arccos() for numbers out of the range of -1 and 1 is not allowed."
				document.getElementById('rootthree').innerHTML = "Input error detected. Using arccos() for numbers out of the range of -1 and 1 is not allowed."
				document.getElementById('rootfour').innerHTML = "Input error detected. Using arccos() for numbers out of the range of -1 and 1 is not allowed."
				scriptexec = 0
			else:
				return math.acos(x)
			
		
		def arcsin(x):
			if (x < -1 or x > 1):
				document.getElementById('rootone').innerHTML = "Input error detected. Using arcsin() for numbers out of the range of -1 and 1 is not allowed."
				document.getElementById('roottwo').innerHTML = "Input error detected. Using arcsin() for numbers out of the range of -1 and 1 is not allowed."
				document.getElementById('rootthree').innerHTML = "Input error detected. Using arcsin() for numbers out of the range of -1 and 1 is not allowed."
				document.getElementById('rootfour').innerHTML = "Input error detected. Using arcsin() for numbers out of the range of -1 and 1 is not allowed."
				scriptexec = 0
			else:
				return math.asin(x)
			
		
		def arctan(x):
			return math.atan(x)
		
		a3 = 0
		b3 = 0
		c3 = 0
		d3 = 0
		e3 = 0
		a_m = document.numbers2.a.value
		b_m = document.numbers2.b.value
		c_m = document.numbers2.c.value
		d_m = document.numbers2.d.value
		e_m = document.numbers2.e.value
		if (execexample == 1):
			a = 15
			b = -58
			c = 8
			d = -80
			e = 64
			document.numbers2.a.value = 15
			document.numbers2.b.value = -58
			document.numbers2.c.value = 8
			document.numbers2.d.value = -80
			document.numbers2.e.value = 64
		elif (execexample == 2):
			a = 0
			b = 3
			c = -10
			d = 14
			e = 27
			document.numbers2.a.value = 0
			document.numbers2.b.value = 3
			document.numbers2.c.value = -10
			document.numbers2.d.value = 14
			document.numbers2.e.value = 27
		elif (execexample == 3):
			a = 0
			b = 0
			c = 10
			d = -28
			e = 16
			document.numbers2.a.value = 0
			document.numbers2.b.value = 0
			document.numbers2.c.value = 10
			document.numbers2.d.value = -28
			document.numbers2.e.value = 16
		else:
			a = eval(document.numbers2.a.value)
			b = eval(document.numbers2.b.value)
			c = eval(document.numbers2.c.value)
			d = eval(document.numbers2.d.value)
			e = eval(document.numbers2.e.value)
		
		if (e == 0 and a != 0):
			specialexec = 2
			scriptexec = 2
		
		if (d == 0 and e == 0):
			specialexec = 1
			scriptexec = 3
		
		discrim = (b * b) - (4 * a * c)
		x1 = (-1 * b + math.sqrt(discrim)) / (2 * a)
		x2 = (-1 * b - math.sqrt(discrim)) / (2 * a)
		alreadydone2 = 0
		divide = b / a
		alreadydone3 = 0
		a2 = 0
		b2 = 0
		c2 = 0
		d2 = 0
		e2 = 0
		dontalert = 0
		ripart = 0
		qipart = 0
		if (a == 0):
			scriptexec = 2 # Recongize as a Cubic Equation
		
		if (a == 0 and b == 0):
			scriptexec = 3 # Recongize as a Quadratic Equation
		
		# Catch user errors
		if (scriptexec != 0):
			if (isNaN(a)):
				scriptexec = 0
				specialexec = 0
				document.getElementById('rootone').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('roottwo').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('rootthree').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('rootfour').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			
			if (isNaN(b)):
				scriptexec = 0
				specialexec = 0
				document.getElementById('rootone').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('roottwo').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('rootthree').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('rootfour').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			
			if (isNaN(c)):
				scriptexec = 0
				specialexec = 0
				document.getElementById('rootone').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('roottwo').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('rootthree').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('rootfour').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			
			if (isNaN(d)):
				scriptexec = 0
				specialexec = 0
				document.getElementById('rootone').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('roottwo').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('rootthree').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('rootfour').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			
			if (isNaN(e)):
				scriptexec = 0
				specialexec = 0
				document.getElementById('rootone').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('roottwo').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('rootthree').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
				document.getElementById('rootfour').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			
		if (a == 0 and b == 0 and c != 0 and d == 0 and e == 0):
			document.getElementById('rootone').innerHTML = 0
			document.getElementById('roottwo').innerHTML = 0
			document.getElementById('rootthree').innerHTML = "Equation is Quadratic, no 3rd solution exists."
			document.getElementById('rootfour').innerHTML = "Equation is Quadratic, no 4th solution exists."
			scriptexec = 0
		
		if (a == 0 and b == 0 and c == 0):
			document.getElementById('rootone').innerHTML = "Input error detected. The Equation is neither Quartic, Cubic, or Quadratic. Please enter a non-zero number into at least 1 of the first 3 boxes."
			document.getElementById('roottwo').innerHTML = "Input error detected. The Equation is neither Quartic, Cubic, or Quadratic. Please enter a non-zero number into at least 1 of the first 3 boxes."
			document.getElementById('rootthree').innerHTML = "Input error detected. The Equation is neither Quartic, Cubic, or Quadratic. Please enter a non-zero number into at least 1 of the first 3 boxes."
			document.getElementById('rootfour').innerHTML = "Input error detected. The Equation is neither Quartic, Cubic, or Quadratic. Please enter a non-zero number into at least 1 of the first 3 boxes."
			scriptexec = 0
		
		if (document.numbers2.a.value == "" or document.numbers2.b.value == "" or document.numbers2.c.value == "" or document.numbers2.d.value == "" or document.numbers2.e.value == ""):
			document.getElementById('rootone').innerHTML = "Input error detected. Please enter all fields."
			document.getElementById('roottwo').innerHTML = "Input error detected. Please enter all fields."
			document.getElementById('rootthree').innerHTML = "Input error detected. Please enter all fields."
			document.getElementById('rootfour').innerHTML = "Input error detected. Please enter all fields."
			scriptexec = 0
		
		if (a_m.indexOf('Math') != -1):
			scriptexec = 0
			specialexec = 0
			document.getElementById('rootone').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('roottwo').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('rootthree').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('rootfour').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
		
		if (b_m.indexOf('Math') != -1):
			scriptexec = 0
			specialexec = 0
			document.getElementById('rootone').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('roottwo').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('rootthree').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('rootfour').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
		
		if (c_m.indexOf('Math') != -1):
			scriptexec = 0
			specialexec = 0
			document.getElementById('rootone').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('roottwo').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('rootthree').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('rootfour').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
		
		if (d_m.indexOf('Math') != -1):
			scriptexec = 0
			specialexec = 0
			document.getElementById('rootone').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('roottwo').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('rootthree').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('rootfour').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
		
		if (e_m.indexOf('Math') != -1):
			scriptexec = 0
			specialexec = 0
			document.getElementById('rootone').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('roottwo').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('rootthree').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
			document.getElementById('rootfour').innerHTML = "Input error detected. There is a non-numerical value entered in one of the fields."
		
		if (scriptexec != 0):
			document.numbers2.a.value = a
			document.numbers2.b.value = b
			document.numbers2.c.value = c
			document.numbers2.d.value = d
			document.numbers2.e.value = e
			if (scriptexec == 1):
				# Extract X^4 Coefficent
				aq = document.numbers2.a.value
				aq2 = aq # Keeps Orignial AQ value
				# Extract X^3 Coefficent
				bq = document.numbers2.b.value
				bq2 = bq # Keeps Orignial BQ Value
				# Extract X^2 Coefficent
				cq = document.numbers2.c.value
				# Extract X Coefficent
				dq = document.numbers2.d.value
				# Extract Constant
				eq = document.numbers2.e.value
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
				h = eval(((g * g) / 4) + ((f * f * f) / 27))
				z = 1 / 3
				i
				j
				k
				l
				m
				n
				p
				xoneterm
				xtwoterm
				xthreeterm
				alreadydone
				alreadydone2 = 0
				ipart = 0
				p = 0
				q = 0
				r = 0
				s = 0
				if (h <= 0):
					_exec = 2
					i = math.sqrt(((g * g) / 4) - h)
					j = pow(i, z)
					k = math.acos(-1 * (g / (2 * i)))
					l = -1 * j
					m = math.cos(k / 3)
					n = math.sqrt(3) * math.sin(k / 3)
					p = (b / (3 * a)) * -1
					xoneterm = (2 * j) * math.cos(k / 3) - (b / (3 * a))
					xtwoterm = l * (m + n) + p
					xthreeterm = l * (m - n) + p
				
				if (h > 0):
					_exec = 1
					R = (-1 * (g / 2)) + math.sqrt(h)
					if (R < 0):
						S = -1 * (pow((-1 * R), z))
					else:
						S = pow(R, z)
					
					T = (-1 * (g / 2)) - math.sqrt(h)
					if (T < 0):
						U = -1 * (pow((-1 * T), z))
					else:
						U = pow(T, z)
					
					xoneterm = (S + U) - (b / (3 * a))
					xtwoterm = (-1 * (S + U) / 2) - (b / (3 * a))
					ipart = ((S - U) * math.sqrt(3)) / 2
					xthreeterm = xtwoterm
				
				if (f == 0 and g == 0 and h == 0):
					if ((d / a) < 0):
						xoneterm = (pow((-1 * (d / a)), z))
						xtwoterm = xoneterm
						xthreeterm = xoneterm
					else:
						xoneterm = -1 * (pow((d / a), z))
						xtwoterm = xoneterm
						xthreeterm = xoneterm
					
				
				# ....and ends here.
				if (abs(ipart) < 5E-7):
					ipart = 0
				
				if (abs(xoneterm) < 5E-7):
					xoneterm = 0
				
				if (abs(xtwoterm) < 5E-7):
					xtwoterm = 0
				
				if (abs(xthreeterm) < 5E-7):
					xthreeterm = 0
				
				# Return to solving the Quartic.
				if (ipart == 0 and xoneterm.toFixed(toplace) == 0):
					alreadydone2 = 1
					p2 = math.sqrt(xtwoterm)
					q = math.sqrt(xthreeterm)
					r = -g2 / (8 * p2 * q)
					s = bq2 / (4 * aq2)
				
				if (ipart == 0 and xtwoterm.toFixed(toplace) == 0 and alreadydone2 == 0 and alreadydone2 != 1):
					alreadydone2 = 2
					p2 = math.sqrt(xoneterm)
					q = math.sqrt(xthreeterm)
					r = -g2 / (8 * p2 * q)
					s = bq2 / (4 * aq2)
				
				if (ipart == 0 and xthreeterm.toFixed(toplace) == 0 and alreadydone2 == 0 and alreadydone2 != 1 and alreadydone2 != 2):
					alreadydone2 = 3
					p2 = math.sqrt(xoneterm)
					q = math.sqrt(xtwoterm)
					r = -g2 / (8 * p2 * q)
					s = bq2 / (4 * aq2)
				
				if (alreadydone2 == 0 and ipart == 0):
					if (xthreeterm.toFixed(toplace) < 0):
						alreadydone2 = 4
						p2 = math.sqrt(xoneterm)
						q = math.sqrt(xtwoterm)
						r = -g2 / (8 * p2 * q)
						s = bq2 / (4 * aq2)
					else:
						alreadydone2 = 5
						p2 = math.sqrt(xoneterm)
						q = math.sqrt(xthreeterm)
						r = -g2 / (8 * p2 * q)
						s = bq2 / (4 * aq2)
					
				
				if (ipart != 0):
					if (xoneterm.toFixed(toplace) < 0):
						xoneterm /= -1
						p2 = 0
						p2ipart = math.sqrt(xoneterm)
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
					
				
				if (ipart == 0 and xtwoterm.toFixed(toplace) < 0 and xthreeterm.toFixed(toplace) < 0):
					xtwoterm /= -1
					xthreeterm /= -1
					p2 = 0
					q = 0
					p2ipart = math.sqrt(xtwoterm)
					qipart = math.sqrt(xthreeterm)
					mult = calcmult(p2, p2ipart, q, qipart, 0)
					r = -g2 / (8 * mult)
					s = bq2 / (4 * aq2)
					ipart = 1
				
				if (xoneterm.toFixed(toplace) > 0 and xtwoterm.toFixed(toplace) < 0 and xthreeterm.toFixed(toplace) == 0 and ipart == 0):
					alreadydone5 = 1
					xtwoterm /= -1
					p2 = math.sqrt(xoneterm)
					q = 0
					p2ipart = 0
					qipart = math.sqrt(xtwoterm)
					mult = calcmult(p2, p2ipart, q, qipart, 0)
					mult2 = calcmult(p2, p2ipart, q, qipart, 1)
					r = -g2 / (8 * mult)
					if (mult2 != 0):
						ripart = g2 / (8 * mult2)
						r = 0
					
					s = bq2 / (4 * aq2)
					ipart = 1
				
				if (xoneterm.toFixed(15) > 0 and xtwoterm.toFixed(15) < 0 and xthreeterm.toFixed(15) == 0 and ipart == 0 and alreadydone5 != 1):
					alreadydone5 = 2
					xtwoterm /= -1
					p2 = math.sqrt(xoneterm)
					q = 0
					p2ipart = 0
					qipart = math.sqrt(xtwoterm)
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
					p2 = math.sqrt(xoneterm)
					q = 0
					p2ipart = 0
					qipart = math.sqrt(xtwoterm)
					mult = calcmult(p2, p2ipart, q, qipart, 0)
					mult2 = calcmult(p2, p2ipart, q, qipart, 1)
					r = -g2 / (8 * mult)
					if (mult2 != 0):
						ripart = g2 / (8 * mult2)
						r = 0
					
					s = bq2 / (4 * aq2)
					ipart = 1
				
				if (xtwoterm.toFixed(10) == 0 and xthreeterm.toFixed(10) == 0 and ipart == 0):
					if (xoneterm.toFixed(10) < 0):
						ipart = 1
						p2 = 0
						p2ipart = math.sqrt(-xoneterm)
						q = 0
						r = 0
						s = bq2 / (4 * aq2)
					else:
						p2 = math.sqrt(xoneterm)
						q = 0
						r = 0
						s = bq2 / (4 * aq2)
					
				
				if (xoneterm.toFixed(10) == 0 and xtwoterm.toFixed(10) == 0 and ipart == 0):
					if (xthreeterm.toFixed(10) < 0):
						ipart = 1
						p2 = 0
						p2ipart = math.sqrt(-xthreeterm)
						q = 0
						r = 0
						s = bq2 / (4 * aq2)
					else:
						p2 = math.sqrt(xthreeterm)
						q = 0
						r = 0
						s = bq2 / (4 * aq2)
					
				
				if (xoneterm.toFixed(8) == 0 and xthreeterm.toFixed(8) == 0 and ipart == 0):
					if (xtwoterm.toFixed(10) < 0):
						ipart = 1
						p2 = 0
						p2ipart = math.sqrt(-xtwoterm)
						q = 0
						r = 0
						s = bq2 / (4 * aq2)
					else:
						p2 = math.sqrt(xtwoterm)
						q = 0
						r = 0
						s = bq2 / (4 * aq2)
					
				
				sum = abs(xtwoterm) + abs(ipart)
				if (xtwoterm == 0 and ipart == 0):
					sum = abs(xthreeterm) + abs(ipart)
				
				if (sum < 5E-8 and approx != 1 and ipart != 0):
					if (xoneterm > 0):
						p2 = math.sqrt(xoneterm)
						q = 0
						p2ipart = 0
						qipart = 0
						mult = 0
						r = 0
						s = bq2 / (4 * aq2)
					else:
						p2 = 0
						q = 0
						p2ipart = math.sqrt(-xoneterm)
						qipart = 0
						mult = 0
						r = 0
						s = bq2 / (4 * aq2)
					
				
				if (p2 < 5.0E-8):
					p2 = 0
				
				if (q < 5.0E-8):
					q = 0
				
				if (abs(ipart) < 5.0E-8 and h > 0 and abs(xoneterm) > 1E-8 and xtwoterm.toFixed(15) != 0):
					if (xoneterm > 0 and xtwoterm > 0):
						p2 = math.sqrt(xoneterm)
						p2ipart = 0
						q = math.sqrt(xtwoterm)
						qipart = 0
						r = -g2 / (8 * p2 * q)
						s = bq2 / (4 * aq2)
					elif (xoneterm < 0 and xtwoterm < 0):
						p2 = 0
						p2ipart = math.sqrt(-xoneterm)
						q = 0
						qipart = math.sqrt(-xtwoterm)
						mult = calcmult(p2, p2ipart, q, qipart, 0)
						mult2 = calcmult(p2, p2ipart, q, qipart, 1)
						r = -g2 / (8 * mult)
						if (mult2 != 0):
							ripart = g2 / (8 * mult2)
							r = 0
						
						s = bq2 / (4 * aq2)
					elif (xoneterm < 0):
						p2 = 0
						p2ipart = math.sqrt(-xoneterm)
						q = math.sqrt(xtwoterm)
						qipart = 0
						mult = calcmult(p2, p2ipart, q, qipart, 0)
						mult2 = calcmult(p2, p2ipart, q, qipart, 1)
						r = -g2 / (8 * mult)
						if (mult2 != 0):
							ripart = g2 / (8 * mult2)
							r = 0
						
						s = bq2 / (4 * aq2)
					else:
						p2 = math.sqrt(xoneterm)
						p2ipart = 0
						q = 0
						qipart = math.sqrt(-xtwoterm)
						mult = calcmult(p2, p2ipart, q, qipart, 0)
						mult2 = calcmult(p2, p2ipart, q, qipart, 1)
						r = -g2 / (8 * mult)
						if (mult2 != 0):
							ripart = g2 / (8 * mult2)
							r = 0
						
						s = bq2 / (4 * aq2)
					
				
				# Now output them answers!
				if (ipart == 0):
					document.getElementById('rootone').innerHTML = convert(eval((p2 + q + r - s).toFixed(toplace)))
					document.getElementById('roottwo').innerHTML = convert(eval((p2 - q - r - s).toFixed(toplace)))
					document.getElementById('rootthree').innerHTML = convert(eval((-p2 + q - r - s).toFixed(toplace)))
					document.getElementById('rootfour').innerHTML = convert(eval((-p2 - q + r - s).toFixed(toplace)))
				
				if (perfect == 1):
					document.getElementById('rootone').innerHTML = convert(eval((-bq / 4).toFixed(toplace)))
					document.getElementById('roottwo').innerHTML = convert(eval((-bq / 4).toFixed(toplace)))
					document.getElementById('rootthree').innerHTML = convert(eval((-bq / 4).toFixed(toplace)))
					document.getElementById('rootfour').innerHTML = convert(eval((-bq / 4).toFixed(toplace)))
				
				if (ipart != 0):
					x1imag = convert(eval((p2ipart + qipart + ripart).toFixed(toplace)))
					x2imag = convert(eval((p2ipart - qipart - ripart).toFixed(toplace)))
					x3imag = convert(eval((-p2ipart + qipart - ripart).toFixed(toplace)))
					x4imag = convert(eval((-p2ipart - qipart + ripart).toFixed(toplace)))
					x1real = convert(eval((p2 + q + r - s).toFixed(toplace)))
					x2real = convert(eval((p2 - q - r - s).toFixed(toplace)))
					x3real = convert(eval((-p2 + q - r - s).toFixed(toplace)))
					x4real = convert(eval((-p2 - q + r - s).toFixed(toplace)))
					x1imagc = eval((p2ipart + qipart + ripart).toFixed(toplace))
					x2imagc = eval((p2ipart - qipart - ripart).toFixed(toplace))
					x3imagc = eval((-p2ipart + qipart - ripart).toFixed(toplace))
					x4imagc = eval((-p2ipart - qipart + ripart).toFixed(toplace))
					x1realc = eval((p2 + q + r - s).toFixed(toplace))
					x2realc = eval((p2 - q - r - s).toFixed(toplace))
					x3realc = eval((-p2 + q - r - s).toFixed(toplace))
					x4realc = eval((-p2 - q + r - s).toFixed(toplace))
					if (abs(x1imagc) < 4E-8):
						x1imag = 0
					
					if (abs(x2imagc) < 4E-8):
						x2imag = 0
					
					if (abs(x3imagc) < 4E-8):
						x3imag = 0
					
					if (abs(x4imagc) < 4E-8):
						x4imag = 0
					
					if (x1real == 0 and x1imag != 0):
						document.getElementById('rootone').innerHTML = ""
					else:
						document.getElementById('rootone').innerHTML = x1real
					
					if (x1imag == -1):
						document.getElementById('rootone').innerHTML += " - i"
					elif (x1imag < 0):
						x1imag /= -1
						document.getElementById('rootone').innerHTML += " - " + x1imag + "i"
					elif (x1imag == 0):
						pass # Do nothing
					elif (x1imag == 1):
						document.getElementById('rootone').innerHTML += " + i"
					else:
						document.getElementById('rootone').innerHTML += " + " + x1imag + "i"
					
					if (x2real == 0 and x2imag != 0):
						document.getElementById('rootthree').innerHTML = ""
					else:
						document.getElementById('rootthree').innerHTML = x2real
					
					if (x2imag == -1):
						document.getElementById('rootthree').innerHTML += " - i"
					elif (x2imag < 0):
						x2imag /= -1
						document.getElementById('rootthree').innerHTML += " - " + x2imag + "i"
					elif (x2imag == 0):
						pass # Do nothing
					elif (x2imag == 1):
						document.getElementById('rootthree').innerHTML += " + i"
					else:
						document.getElementById('rootthree').innerHTML += " + " + x2imag + "i"
					
					if (x3real == 0 and x3imag != 0):
						document.getElementById('rootfour').innerHTML = ""
					else:
						document.getElementById('rootfour').innerHTML = x3real
					
					if (x3imag == -1):
						document.getElementById('rootfour').innerHTML += " - i"
					elif (x3imag < 0):
						x3imag /= -1
						document.getElementById('rootfour').innerHTML += " - " + x3imag + "i"
					elif (x3imag == 0):
						pass # Do nothing
					elif (x3imag == 1):
						document.getElementById('rootfour').innerHTML += " + i"
					else:
						document.getElementById('rootfour').innerHTML += " + " + x3imag + "i"
					
					if (x4real == 0 and x4imag != 0):
						document.getElementById('roottwo').innerHTML = ""
					else:
						document.getElementById('roottwo').innerHTML = x4real
					
					if (x4imag == -1):
						document.getElementById('roottwo').innerHTML += " - i"
					elif (x4imag < 0):
						x4imag /= -1
						document.getElementById('roottwo').innerHTML += " - " + x4imag + "i"
					elif (x4imag == 0):
						pass # Do nothing
					elif (x4imag == 1):
						document.getElementById('roottwo').innerHTML += " + i"
					else:
						document.getElementById('roottwo').innerHTML += " + " + x4imag + "i"
					
				
				if (perfectbiquadratic == 1):
					perfect= eval(math.sqrt(cq / 2).toFixed(toplace))
					if (perfect== 1):
						document.getElementById('rootone').innerHTML = " + i"
						document.getElementById('roottwo').innerHTML = " + i"
						document.getElementById('rootthree').innerHTML = " - i"
						document.getElementById('rootfour').innerHTML = " - i"
					else:
						document.getElementById('rootone').innerHTML = " + " + perfect+ "i"
						document.getElementById('roottwo').innerHTML = " + " + perfect+ "i"
						document.getElementById('rootthree').innerHTML = " - " + perfect+ "i"
						document.getElementById('rootfour').innerHTML = " - " + perfect+ "i"
					
				
				if (cq == 0 and dq == 0 and eq == 0):
					document.getElementById('rootone').innerHTML = eval(-divide.toFixed(10))
					document.getElementById('roottwo').innerHTML = 0
					document.getElementById('rootthree').innerHTML = 0
					document.getElementById('rootfour').innerHTML = 0
				
				if (approx == 2):
					document.getElementById('rootone').innerHTML += " [Approx]"
					document.getElementById('roottwo').innerHTML += " [Approx]"
					document.getElementById('rootthree').innerHTML += " [Approx]"
					document.getElementById('rootfour').innerHTML += " [Approx]"
				
				if (document.getElementById('rootone').innerHTML == "0" and document.getElementById('roottwo').innerHTML == "0" and document.getElementById('rootthree').innerHTML == "0" and document.getElementById('rootfour').innerHTML == "0"):
					document.getElementById('rootone').innerHTML = 0
					document.getElementById('roottwo').innerHTML = 0
					document.getElementById('rootthree').innerHTML = 0
					document.getElementById('rootfour').innerHTML = 0
				
				if (document.getElementById('rootone').innerHTML == "NaN" and document.getElementById('roottwo').innerHTML == "NaN" and document.getElementById('rootthree').innerHTML == "NaN" and document.getElementById('rootfour').innerHTML == "NaN"):
					document.getElementById('rootone').innerHTML = 0
					document.getElementById('roottwo').innerHTML = 0
					document.getElementById('rootthree').innerHTML = 0
					document.getElementById('rootfour').innerHTML = 0
				
			
			if (scriptexec == 2):
				if (specialexec != 2):
					a = b
					b = c
					c = d
					d = e
				
				divide = -b / a
				asign = "+"
				bsign = "+"
				csign = "+"
				dsign = "+"
				a2 = a
				b2 = b
				c2 = c
				d2 = d
				if (b < 0):
					bsign = "-"
					b2 /= -1
				
				if (c < 0):
					csign = "-"
					c2 /= -1
				
				if (d < 0):
					dsign = "-"
					d2 /= -1
				
				if (b == 0):
					bsign = ""
					b2 = ""
					x2show = 0
				
				if (c == 0):
					csign = ""
					c2 = ""
					xshow = 0
				
				if (d == 0):
					dsign = ""
					d2 = ""
				
				if (a == 1):
					a2 = ""
				
				if (a == -1):
					a2 = "-"
				
				if (abs(b) == 1):
					b2 = ""
				
				if (abs(c) == 1):
					c2 = ""
				
				if (xshow == 0):
					xshow = ""
				else:
					xshow = "x"
				
				if (x2show == 0):
					x2show = ""
				else:
					x2show = "x<sup>2</sup>"
				
				f = (((3 * c) / a) - ((b * b) / (a * a))) / 3
				g = (((2 * b * b * b) / (a * a * a)) - ((9 * b * c) / (a * a)) + ((27 * d) / a)) / 27
				h = eval(((g * g) / 4) + ((f * f * f) / 27))
				z = 1 / 3
				i
				j
				k
				l
				m
				n
				p
				xoneterm
				xtwoterm
				xthreeterm
				alreadydone
				if (discrim == 0 and d == 0):
					xoneterm = 0
					xtwoterm = ((-1 * b) / (2 * a))
				
				if (h <= 0):
					_exec = 2
					i = math.sqrt(((g * g) / 4) - h)
					j = pow(i, z)
					k = math.acos(-1 * (g / (2 * i)))
					l = -1 * j
					m = math.cos(k / 3)
					n = math.sqrt(3) * math.sin(k / 3)
					p = (b / (3 * a)) * -1
					xoneterm = (2 * j) * math.cos(k / 3) - (b / (3 * a))
					xtwoterm = l * (m + n) + p
					xthreeterm = l * (m - n) + p
				
				if (h > 0):
					_exec = 1
					R = (-1 * (g / 2)) + math.sqrt(h)
					if (R < 0):
						S = -1 * (pow((-1 * R), z))
					else:
						S = pow(R, z)
					
					T = (-1 * (g / 2)) - math.sqrt(h)
					if (T < 0):
						U = -1 * (pow((-1 * T), z))
					else:
						U = pow(T, z)
					
					xoneterm = (S + U) - (b / (3 * a))
					xtwoterm = (-1 * (S + U) / 2) - (b / (3 * a))
					ipart = ((S - U) * math.sqrt(3)) / 2
					xthreeterm = xtwoterm
				
				if (f == 0 and g == 0 and h == 0):
					if ((d / a) < 0):
						xoneterm = (pow((-1 * (d / a)), z))
					else:
						xoneterm = -1 * (pow((d / a), z))
					
				
				document.getElementById('rootone').innerHTML = convert(parseFloat(xoneterm.toFixed(toplace)))
				if (h > 0):
					if (xtwoterm.toFixed(toplace) == 0 and xthreeterm.toFixed(toplace) == 0):
						if (ipart.toFixed(toplace) == 1):
							document.getElementById('roottwo').innerHTML = "i"
							document.getElementById('rootthree').innerHTML = "- i"
						elif (ipart.toFixed(toplace) == 0):
							document.getElementById('roottwo').innerHTML = "0"
							document.getElementById('rootthree').innerHTML = "0"
						else:
							document.getElementById('roottwo').innerHTML = convert(parseFloat(ipart.toFixed(toplace))) + "i"
							document.getElementById('rootthree').innerHTML = "-" + convert(parseFloat(ipart.toFixed(toplace))) + "i"
						
					else:
						if (ipart.toFixed(toplace) == 1):
							document.getElementById('roottwo').innerHTML = convert(parseFloat(xtwoterm.toFixed(toplace))) + " + i"
							document.getElementById('rootthree').innerHTML = convert(parseFloat(xthreeterm.toFixed(toplace))) + " - i"
						else:
							document.getElementById('roottwo').innerHTML = convert(parseFloat(xtwoterm.toFixed(toplace))) + " + " + convert(parseFloat(ipart.toFixed(toplace))) + "i"
							document.getElementById('rootthree').innerHTML = convert(parseFloat(xthreeterm.toFixed(toplace))) + " - " + convert(parseFloat(ipart.toFixed(toplace))) + "i"
						
					
				
				if (h < 0):
					document.getElementById('roottwo').innerHTML = convert(parseFloat(xtwoterm.toFixed(toplace)))
					document.getElementById('rootthree').innerHTML = convert(parseFloat(xthreeterm.toFixed(toplace)))
				
				if (discrim == 0 and d == 0):
					document.getElementById('roottwo').innerHTML = convert(parseFloat(xtwoterm.toFixed(toplace)))
					document.getElementById('rootthree').innerHTML = document.getElementById('roottwo').innerHTML
				
				if (a != 0 and b != 0 and c == 0 and d == 0):
					alreadydone = 1
					document.getElementById('rootone').innerHTML = 0
					document.getElementById('roottwo').innerHTML = eval((-b / a).toFixed(toplace))
					document.getElementById('rootthree').innerHTML = 0
				
				if (b == 0 and c == 0 and d == 0):
					alreadydone = 1
					document.getElementById('roottwo').innerHTML = "0"
					document.getElementById('rootthree').innerHTML = "0"
				
				if (h == 0 and alreadydone != 1):
					document.getElementById('roottwo').innerHTML = convert(parseFloat(xtwoterm.toFixed(toplace)))
					document.getElementById('rootthree').innerHTML = convert(parseFloat(xthreeterm.toFixed(toplace)))
				
				/* if (h == 0):
 if (f < 0 and g < 0):
document.getElementById('roottwo').innerHTML = parseFloat(xtwoterm.toFixed(toplace))
document.getElementById('rootthree').innerHTML = parseFloat(xthreeterm.toFixed(toplace))

alreadydone2 = 2
document.getElementById('roottwo').innerHTML = parseFloat(xtwoterm.toFixed(toplace))
document.getElementById('rootthree').innerHTML = parseFloat(xoneterm.toFixed(toplace))
 */
				if (f == 0 and g == 0 and h == 0):
					alreadydone2 = 1
					document.getElementById('roottwo').innerHTML = parseFloat(xoneterm.toFixed(toplace))
					document.getElementById('rootthree').innerHTML = parseFloat(xoneterm.toFixed(toplace))
				
				if (c == 0 and d == 0 and alreadydone != 1):
					alreadydone = 2
					if (x2.toFixed(toplace) == 0):
						document.getElementById('rootone').innerHTML = convert(parseFloat(x1.toFixed(toplace)))
						document.getElementById('roottwo').innerHTML = "0"
						document.getElementById('rootthree').innerHTML = "0"
					else:
						alreadydone3 = 1
						document.getElementById('rootone').innerHTML = "0"
						document.getElementById('roottwo').innerHTML = convert(parseFloat(x2.toFixed(toplace)))
						document.getElementById('rootthree').innerHTML = "0"
					
				
				if (specialexec != 2):
					document.getElementById('rootfour').innerHTML = "Equation is Cubic, no 4th solution exists."
				else:
					document.getElementById('rootfour').innerHTML = 0
				
			
			if (scriptexec == 3):
				asign = "+"
				bsign = "+"
				csign = "+"
				dsign = "+"
				a2 = a
				b2 = b
				c2 = c
				d2 = d
				if (c < 0):
					csign = "-"
					c2 /= -1
				
				if (d < 0):
					dsign = "-"
					d2 /= -1
				
				if (c == 0):
					csign = ""
					c2 = ""
					xshow = 0
				
				if (d == 0):
					dsign = ""
					d2 = ""
				
				if (a == 1):
					a2 = ""
				
				if (a == -1):
					a2 = "-"
				
				if (abs(b) == 1):
					b2 = ""
				
				if (abs(c) == 1):
					c2 = ""
				
				if (xshow == 0):
					xshow = ""
				else:
					xshow = "x"
				
				if (x2show == 0):
					x2show = ""
				else:
					x2show = "x<sup>2</sup>"
				
				if (specialexec == 1):
					xoneterm1 = (-1 * b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
					xtwoterm1 = (-1 * b - math.sqrt(b * b - 4 * a * c)) / (2 * a)
					sqrtterm = (b * b) - (4 * a * c)
					acterm = ((4 * a * c)).toString()
					realpart = (-1 * b / (2 * a)).toFixed(toplace)
					imag = (math.sqrt(-sqrtterm) / (2 * a)).toFixed(toplace)
				else:
					xoneterm1 = (-1 * d + math.sqrt(d * d - 4 * c * e)) / (2 * c)
					xtwoterm1 = (-1 * d - math.sqrt(d * d - 4 * c * e)) / (2 * c)
					sqrtterm = (d * d) - (4 * c * e)
					acterm = ((4 * c * e)).toString()
					realpart = (-1 * d / (2 * c)).toFixed(toplace)
					imag = (math.sqrt(-sqrtterm) / (2 * c)).toFixed(toplace)
				
				document.getElementById('rootone').innerHTML = convert(eval(xoneterm1.toFixed(toplace)))
				document.getElementById('roottwo').innerHTML = convert(eval(xtwoterm1.toFixed(toplace)))
				if (specialexec == 0):
					document.getElementById('rootthree').innerHTML = "Equation is Quadratic, no 3rd solution exists."
					document.getElementById('rootfour').innerHTML = "Equation is Quadratic, no 4th solution exists."
				else:
					document.getElementById('rootthree').innerHTML = 0
					document.getElementById('rootfour').innerHTML = 0
				
				if (specialexec == 0):
					if (sqrtterm < 0):
						if (d == 0 or realpart == 0):
							sqrtterm = -1 * sqrtterm
							document.getElementById('rootone').innerHTML = convert(eval(imag)) + "i"
							document.getElementById('roottwo').innerHTML = convert(eval(-1 * imag)) + "i"
						else:
							sqrtterm = -1 * sqrtterm
							document.getElementById('rootone').innerHTML = convert(eval(realpart)) + " + " + convert(eval(abs(imag))) + "i"
							document.getElementById('roottwo').innerHTML = convert(eval(realpart)) + " - " + convert(eval(abs(imag))) + "i"
							if (((math.sqrt(sqrtterm).toFixed(toplace)) / (2 * c).toFixed(toplace)) == 1):
								document.getElementById('rootone').innerHTML = convert(eval(realpart)) + " + i"
								document.getElementById('roottwo').innerHTML = convert(eval(realpart)) + " - i"
							
						
					
				else:
					if (sqrtterm < 0):
						if (b == 0):
							sqrtterm = -1 * sqrtterm
							document.getElementById('rootone').innerHTML = convert(eval(imag)) + "i"
							document.getElementById('roottwo').innerHTML = convert(eval(-1 * imag)) + "i"
						else:
							sqrtterm = -1 * sqrtterm
							document.getElementById('rootone').innerHTML = convert(eval(realpart)) + " + " + convert(eval(abs(imag))) + "i"
							document.getElementById('roottwo').innerHTML = convert(eval(realpart)) + " - " + convert(eval(abs(imag))) + "i"
							if (((math.sqrt(sqrtterm).toFixed(toplace)) / (2 * a).toFixed(toplace)) == 1):
								document.getElementById('rootone').innerHTML = convert(eval(realpart)) + " + i"
								document.getElementById('roottwo').innerHTML = convert(eval(realpart)) + " - i"
							
						
					
				
				if (specialexec == 0):
					if (eval(xoneterm1.toFixed(toplace)) == eval(xtwoterm1.toFixed(toplace))):
						document.getElementById('rootone').innerHTML = convert(eval(xoneterm1.toFixed(toplace)))
						document.getElementById('roottwo').innerHTML = convert(eval(xoneterm1.toFixed(toplace)))
						document.getElementById('rootthree').innerHTML = "Equation is Quadratic, no 3rd solution exists."
					
				
				if (realpart == 0 and imag == 1):
					document.getElementById('rootone').innerHTML = "i"
					document.getElementById('roottwo').innerHTML = "- i"
				
			
		
	
