#!/Applications/VisIt.app/Contents/Resources/bin/visit visit
class HaltException(Exception): pass
import math 


#### Constants: Make sure they are correct!!! ###
pi = math.pi 
omega=0.09793983350312327
numr=258
dr=1/(130-3.0) 
#factor=((numr-3)*dr)/(2*pi)  #0.315054
factor=(dr)/(2*pi/256) 
dz=dr
dphi=2*pi/256.0
rhoth = 1e-5


def findBounds():
# Find upper and lower bounds of the grid
	global xmin, xmax, zmin, zmax

        AddPlot("Pseudocolor", "X")
	DrawPlots()
	Query("MinMax")
	xmmstr = string.split(GetQueryOutputString(), " ")
	xmin = xmmstr[4]
	xmin = ''.join(xmin)
        xmin = float(xmin)	
	xmax = xmmstr[15]
	xmax = ''.join(xmax)
        xmax = float(xmax)	
	
	DeleteAllPlots()
#	print "xmin",xmin
#	print "xmax",xmax

        AddPlot("Pseudocolor", "Z")
	DrawPlots()
	Query("MinMax")
	zmmstr = string.split(GetQueryOutputString(), " ")
	zmin = zmmstr[4]
	zmin = ''.join(zmin)
        zmin = float(zmin)	
	zmax = zmmstr[15]
	zmax = ''.join(zmax)
        zmax = float(zmax)	
	
	DeleteAllPlots()

#	print "zmin",zmin
#	print "zmax",zmax


	

def rBounds():
# Reject the boundary points on the grid
	global xmin, xmax, zmin, zmax

	AddOperator("Isovolume")
	isoz = IsovolumeAttributes()
	isoz.variable = "Z"
	isoz.ubound = zmax-dr
	isoz.lbound = zmin+dr
	SetOperatorOptions(isoz)

	isor = IsovolumeAttributes()
	isor.variable = "r"
	isor.ubound = xmax-10*dr
	isor.lbound = xmin+dr
	SetOperatorOptions(isor)





def getCoords():
# Get the basic coordinates
	DefineVectorExpression("coords", "coords(mesh)")
	DefineScalarExpression("X", "coords[0]")
	DefineScalarExpression("Y", "coords[1]")
	DefineScalarExpression("Z", "coords[2]")	
	DefineScalarExpression("r",  "sqrt(X*X+Y*Y)" )


def totalMass():
# Find total mass on the grid
	global mt
	DefineScalarExpression("mass",  "<mesh_quality/volume>*den")
        AddPlot("Pseudocolor", "mass")
	rBounds()
	DrawPlots()
	Query("Variable Sum")
	mtstr = string.split(GetQueryOutputString(), " ")
	mt = mtstr[4]
	mt = float(mt)*factor**3 
	DeleteAllPlots()


def totalAM():
# Find total AM on the grid
	global jt
	DefineScalarExpression("jtF",  "<mesh_quality/volume>*den*(velphi+r*r*(%f)*(%f)*(%f))"%(factor,factor,omega))
        AddPlot("Pseudocolor", "jtF")
	rBounds()
	DrawPlots()
	Query("Variable Sum")
	jtstr = string.split(GetQueryOutputString(), " ")
	jt = jtstr[4]
	jt = float(jt)*factor**3 
	DeleteAllPlots()


def findL1():
# 1. Find the angle through which the primary is rotated wrt x-axis
# 2. Find location of the minimum potential for both stars
# 3. Find location and potential of L1 point  BUG!!
# 4. Find the plane dividing two stars 
	global theta, maxx, maxy, maxz, maxx2, maxy2, maxz2
	global l1x, l1y, l1z, l1_val, div

# find location of primary as deepent potential
       	AddPlot("Pseudocolor", "pot")
       	DrawPlots()

       	Query("Min", use_actual_data=1) 
	q = string.split(GetQueryOutputString(), " ")
	maxx = list(q[9])
       	maxy = list(q[10])
	maxz = list(q[11])

	maxx[0] = " "
        maxx[-1] = " "
	maxy[-1] = " "
	maxz[-6:-1] = " "

	maxx = ''.join(maxx)
        maxx = float(maxx)
        maxy = ''.join(maxy)
        maxy = float(maxy)
	maxz = ''.join(maxz)
        maxz = float(maxz)


#	print "Deepest potential --->", maxx, maxy, maxz


# Find the angle through which the primary is rotated 
	if maxx>0:
		if maxy>0: 
			theta = -math.atan(maxy/maxx)
		else:
			theta = -math.atan(maxy/maxx)
	if maxx<0:
		if maxy>0: 
			theta = -math.atan(maxy/maxx)-pi
		else:
			theta = -math.atan(maxy/maxx)-pi

		
#	print "theta = ", theta	
	

	DefineScalarExpression("rotX", "X*cos(%f)-Y*sin(%f)" % (theta, theta))
	DefineScalarExpression("rotY", "X*sin(%f)+Y*cos(%f)" % (theta, theta))	

# Find location of the other star with lower potential well
	AddOperator("Isovolume")
	iso1 = IsovolumeAttributes()
	iso1.variable = "rotX"
	iso1.ubound = 0
	SetOperatorOptions(iso1)

        DrawPlots()

	Query("Min", use_actual_data=1)
	q2 = string.split(GetQueryOutputString(), " ")
	maxx2 = list(q2[9])
        maxy2 = list(q2[10])
	maxz2 = list(q[11])


	maxx2[0] = " "
        maxx2[-1] = " "
	maxy2[-1] = " "
	maxz2[-6:-1] = " "

	maxx2 = ''.join(maxx2)
        maxx2 = float(maxx2)
        maxy2 = ''.join(maxy2)
        maxy2 = float(maxy2)
	maxz2 = ''.join(maxz2)
        maxz2 = float(maxz2)

#	print "Second deepest potential --->", maxx2, maxy2, maxz2
		
# Find the L1 point and L1 potential
	DeleteAllPlots()
        AddPlot("Pseudocolor", "rotX")
	DrawPlots()
	DeleteAllPlots()

        AddPlot("Pseudocolor", "pot")
	DrawPlots()


	Lineout(( maxx, maxy, maxz),( maxx2, maxy2, maxz2))
# Lineout plots are replotted with the origin at the beginning of the plot,
# hence the min/max locations will have an offset!!

	SetActiveWindow(2)

	Query("Max", use_actual_data=1)
 	q3 = string.split(GetQueryOutputString(), " ")
#        print q3

	DeleteAllPlots()
	DeleteWindow()

	l1_val = q3[4]
	l1x = list(q3[9])
        l1y = list(q3[10])
	l1z = list(q3[11])

	l1x[0] = " "
        l1x[-1] = " "
	l1y[-1] = " "
	l1z[-4:-1] = " "
	l1_val = ''.join(l1_val)
	l1_val = ''.join(l1_val)
        l1_val = float(l1_val)
	l1x = ''.join(l1x)
        l1x = float(l1x)
        l1y = ''.join(l1y)
        l1y = float(l1y)
        l1z = ''.join(l1z)
	l1z = float(l1z)

		
#	print "l1x = ", l1x,"l1y = ", l1y,"l1z = ",l1z, "l1val", l1_val

#  Find the dividing plane
	div = -(l1x - math.sqrt(maxx*maxx+maxy*maxy))
	corr =  math.sqrt(maxx*maxx+maxy*maxy)
#	print "corr ", corr, maxx*math.cos(theta)-maxy*math.sin(theta)
#	print "div = " , div

	DeleteAllPlots()
	



def getMasses():
# Get masses of individual stars
# Find mass of one star 
	global div, m1, m2

	DeleteAllPlots()

#	DefineScalarExpression("r",  "sqrt(X*X+Y*Y)" )
#	DefineScalarExpression("mass",  "den*r" )
	DefineScalarExpression("mass",  "<mesh_quality/volume>*den")

        AddPlot("Pseudocolor", "mass")

	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.lbound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "den"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)


	DrawPlots()
		
	Query("Variable Sum")
	m1str = string.split(GetQueryOutputString(), " ")
#	print m1str
	m1 = m1str[4]
	m1 = float(m1)*factor**3 

	DeleteAllPlots()

# Find mass of the other star
       	AddPlot("Pseudocolor", "mass")

	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.ubound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "den"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m2str = string.split(GetQueryOutputString(), " ")
#	print m2str
	m2 = m2str[4]
	m2 = float(m2)*factor**3 


	q = m2/m1
#	print m1, m2 , " q = ", m2/m1

        with open('qtest','a') as fil:

                fil.write(str(q))
        	fil.write("\n")
        fil.close()

	DeleteAllPlots()


def getCom():
# Actually get SUM m*r, i.e. divide by M to get CoM
# Find CoM of one star
	global div
	global rhoth
	global m1, m2, com1x, com1y, com1z, com2x,com2y, com2z, sep 
	print "Finding CoMs"

	DefineScalarExpression("mass",  "<mesh_quality/volume>*den")

	DefineScalarExpression("comx",  "mass*X" )
	DefineScalarExpression("comy",  "mass*Y" )
	DefineScalarExpression("comz",  "mass*Z" )

       	AddPlot("Pseudocolor", "comx")
	rBounds()
	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.lbound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "den"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m1xstr = string.split(GetQueryOutputString(), " ")
#	print m1xstr
	m1x = m1xstr[4]
	m1x = float(m1x)*factor**4


	DeleteAllPlots()
		

       	AddPlot("Pseudocolor", "comy")
	rBounds()
	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.lbound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "den"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m1ystr = string.split(GetQueryOutputString(), " ")
#	print m1ystr
	m1y = m1ystr[4]
	m1y = float(m1y)*factor**4

	DeleteAllPlots()


       	AddPlot("Pseudocolor", "comz")
	rBounds()
	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.lbound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "den"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m1zstr = string.split(GetQueryOutputString(), " ")
#	print m1zstr
	m1z = m1zstr[4]
	m1z = float(m1z)*factor**4

	DeleteAllPlots()


# Find CoM of the other star

       	AddPlot("Pseudocolor", "comx")
	rBounds()
	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.ubound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "den"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m2xstr = string.split(GetQueryOutputString(), " ")
#	print m2xstr
	m2x = m2xstr[4]
	m2x = float(m2x)*factor**4

	DeleteAllPlots()
		

       	AddPlot("Pseudocolor", "comy")
	rBounds()
	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.ubound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "den"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m2ystr = string.split(GetQueryOutputString(), " ")
#	print m2ystr
	m2y = m2ystr[4]
	m2y = float(m2y)*factor**4

	DeleteAllPlots()


       	AddPlot("Pseudocolor", "comz")
	rBounds()
	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.ubound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "den"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m2zstr = string.split(GetQueryOutputString(), " ")
#	print m2zstr
	m2z = m2zstr[4]
	m2z = float(m2z)*factor**4

	DeleteAllPlots()


	com1x = m1x/m1
	com1y = m1y/m1
	com1z =  m1z/m1

	com2x =  m2x/m2
	com2y =  m2y/m2
	com2z =  m2z/m2	
	
	sep = math.sqrt((m1x/m1-m2x/m2)**2+(m1y/m1-m2y/m2)**2+(m1z/m1-m2z/m2)**2) 

#	print "sep: ", sep
#	print com1x, com1y, com2x, com2y, div, rhoth


def findSpins(rhothv):
	global com1x, com1y, com2x, com2y, div, rhoth, j1, j2, theta

	print "Finding spins.."
#	print com1x, com1y, com2x, com2y, div, rhoth
#	com1x = -0.467819661435 
#	com1y = -0.0057412893278 
#	com2x = 0.666178785115 
#	com2y = 0.00817566537476 
#	div = -0.4417500465 
#	rhoth = 1e-05
#	theta = -3.1538645064

	print "com1x",com1x,"com1y", com1y,"com2x", com2x,"com2y", com2y

	DefineScalarExpression("rotX", "X*cos(%f)-Y*sin(%f)" % (theta, theta))
	DefineScalarExpression("rotY", "X*sin(%f)+Y*cos(%f)" % (theta, theta))	



	DefineScalarExpression("vrot", "r*(%f)*(%f)" % (omega,factor))
	DefineScalarExpression("vphinew", "vrot+velphi/r/(%f)" % (factor))
	DefineScalarExpression("sinF", "Y/r")
	DefineScalarExpression("cosF", "X/r")
	
#	DefineScalarExpression("dmF", "<mesh_quality/volume>*den*(%f)*(%f)*(%f)" % (factor,factor,factor))
	DefineScalarExpression("velrnew", "velr")#/r/(%f)" % (factor) )

#	DefineScalarExpression("jtF",  "<mesh_quality/volume>*den*(velphi+r*r*(%f)*(%f)*(%f))"%(factor,factor,omega))
#	jt = float(jt)*factor**3 



#	DefineScalarExpression("jF",  "<mesh_quality/volume>*den*(velphi+r*r*(%f)*(%f)*(%f))"%(factor,factor,omega))
#	DefineScalarExpression("j1F", \
#	 "<mesh_quality/volume>*den* (  (X*(%f)-(%f))*(sinF*velrnew+cosF*vphinew)   +   (Y*(%f)-(%f))*(cosF*velrnew-sinF*vphinew)   )" \
#	 % (factor, com1x, factor, com1y)) 

	DefineScalarExpression("j1F", \
	"<mesh_quality/volume>*den* (  (X*(%f)-(%f))*(sinF*velrnew+cosF*vphinew)   -   (Y*(%f)-(%f))*(cosF*velrnew-sinF*vphinew)   )" \
	 % (factor, com1x, factor, com1y)) 


	DefineScalarExpression("j2F",  \
	 "<mesh_quality/volume>*den* (  (X*(%f)-(%f))*(sinF*velrnew+cosF*vphinew)   -   (Y*(%f)-(%f))*(cosF*velrnew-sinF*vphinew)   )" \
	  % (factor,com2x, factor, com2y)) 
#	print "omega, factor", omega, factor

        AddPlot("Pseudocolor", "j1F")

	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.lbound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "den"
	iso5.lbound = rhothv
	SetOperatorOptions(iso5)


	DrawPlots()
		
	Query("Variable Sum")
	j1str = string.split(GetQueryOutputString(), " ")
#	print j1str
	j1 = j1str[4]
	j1 = float(j1)*factor**3


	DeleteAllPlots()

	print "j1",j1

        AddPlot("Pseudocolor", "j2F")

	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.ubound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "den"
	iso5.lbound = rhothv
	SetOperatorOptions(iso5)


	DrawPlots()
		
	Query("Variable Sum")
	j2str = string.split(GetQueryOutputString(), " ")
#	print j2str
	j2 = j2str[4]
	j2 = float(j2)*factor**3

	DeleteAllPlots()

	print "j2",j2

	s1 = j1
	s2 = j2


#J_point = SUM rho * r * dr dz dphi [(X-X_point) (sin(theta)*velr+cos(theta)*vphinew)
#+ (Y-Y_point) (cos(theta)*velr-sin(theta)*vphinew)]
#where vphinew = velphi + r*omega.



############### Program starts here ###############
try:
        for i in range(1000,1001,30):
                db = str(i)+".silo"
                print db

# Open the file
                OpenDatabase(db)

# Coordinates
                getCoords()

                findBounds()

# Find L1 point, angle of the binary, location of dividing plane
		findL1()

# Find individual masses and volumes
		getMasses()

# Find total mass
                totalMass()
#                print "mt", mt

# Find total angular momentum
                totalAM()
                print "jt", jt




# Find center of masses
		getCom()

# Find spin angular momenta of the stars

		findSpins(0.0)
		s1=j1
		s2=j2
                print "s12", s1,s2

		findSpins(rhoth) #1e-5
		s3=j1
		s4=j2
		print "s34", s3,s4

		findSpins(0.0001)
                s5=j1
                s6=j2
                print "s56", s5,s6

                findSpins(0.001)
                s7=j1
                s8=j2
                print "s78", s7,s8

                findSpins(0.01)
                s9=j1
                s10=j2
                print "s910", s9,s10

                findSpins(0.1)
		s11=j1
		s12=j2
                print "s1112", s11,s12



#		pirint "m1, m2", m1, m2
#		print  maxx, maxy, maxz,maxx2, maxy2, maxz2
		a = math.sqrt((maxx-maxx2)**2+(maxy-maxy2)**2+(maxz-maxz2)**2)
#		print "a", a

                CloseDatabase(db)
	        with open('spins','a') as fil1:
			fil1.write(str(db)+' '), fil1.write(str(i)+' '), \
			fil1.write(str(s1)+' '), fil1.write(str(s2)+' '),fil1.write(str(s3)+' '),\
			fil1.write(str(s4)+' '), fil1.write(str(s5)+' '), fil1.write(str(s6)+' '),\
			fil1.write(str(j1)+' '), fil1.write(str(j2)+' ')  
               		fil1.write("\n")
		fil1.close()

		with open('diags','a') as fil2:
			fil2.write(str(db)+' '), fil2.write(str(i)+' '), \
	                fil2.write(str(mt)+' '), fil2.write(str(m1)+' '), \
			fil2.write(str(m2)+' '), fil2.write(str(a)+' '), fil2.write(str(sep)+' '), \
			fil2.write(str(jt)+' '), fil2.write(str(s1)+' '), fil2.write(str(s2)+' ')
	        	fil2.write("\n")
		fil2.close()


        raise HaltException("Done!")


except HaltException as h:
        print(h)

############### Program ends here ###############

