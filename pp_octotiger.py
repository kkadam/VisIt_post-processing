#!/Applications/VisIt.app/Contents/Resources/bin/visit visit
class HaltException(Exception): pass
import math 

pi = math.pi 
rhoth = 1e-5


############## Define functions ################

def getCoords():
# Get the basic coordinates
	DefineVectorExpression("coords", "coords(mesh)")
	DefineScalarExpression("X", "coords[0]")
	DefineScalarExpression("Y", "coords[1]")
	DefineScalarExpression("Z", "coords[2]")	
	DefineScalarExpression("r",  "sqrt(X*X+Y*Y)" )



def findL1():
# 1. Find the angle through which the primary is rotated wrt x-axis
# 2. Find location of the minimum potential for both stars
# 3. Find location and potential of L1 point  BUG!!
# 4. Find the plane dividing two stars 
	global theta, maxx, maxy, maxz, maxx2, maxy2, maxz2
	global l1x, l1y, l1z, l1_val, div

# find location of primary as deepent potential
       	AddPlot("Pseudocolor", "phi")
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

        AddPlot("Pseudocolor", "phi")
	DrawPlots()


	Lineout(( maxx, maxy, maxz),( maxx2, maxy2, maxz2))
# Lineout plots are replotted with the origin at the beginning of the plot,
# hence the min/max locations will have an offset!!

	SetActiveWindow(2)

	Query("Max", use_actual_data=1)
 	q3 = string.split(GetQueryOutputString(), " ")

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

	DefineScalarExpression("mass",  "<mesh_quality/volume>*rho")

        AddPlot("Pseudocolor", "mass")

	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.lbound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "rho"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)


	DrawPlots()
		
	Query("Variable Sum")
	m1str = string.split(GetQueryOutputString(), " ")
	m1 = m1str[4]
	m1 = float(m1)

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
	iso5.variable = "rho"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m2str = string.split(GetQueryOutputString(), " ")
	m2 = m2str[4]
	m2 = float(m2)


	q = m2/m1
#	print m1, m2 , " q = ", m2/m1


	DeleteAllPlots()



def totalMass():
# Find total mass on the grid
	global mt
	DefineScalarExpression("mass",  "<mesh_quality/volume>*rho")
        AddPlot("Pseudocolor", "mass")
	DrawPlots()
	Query("Variable Sum")
	mtstr = string.split(GetQueryOutputString(), " ")
	mt = mtstr[4]
	mt = float(mt)
	DeleteAllPlots()


def totalAM():
# Find total AM on the grid
	global jt

	DefineScalarExpression("sinF", "Y/r")
	DefineScalarExpression("cosF", "X/r")

#	DefineScalarExpression("velphi","cosF + sinF")

	DefineScalarExpression("jtF",  "<mesh_quality/volume>* r * ( sy*cosF - sx*sinF )")
        AddPlot("Pseudocolor", "jtF")
	DrawPlots()
	Query("Variable Sum")
	jtstr = string.split(GetQueryOutputString(), " ")
	jt = jtstr[4]
	jt = float(jt)
	DeleteAllPlots()


def getCom():
# Actually get SUM m*r, i.e. divide by M to get CoM
# Find CoM of one star
	global div
	global rhoth
	global m1, m2, com1x, com1y, com1z, com2x,com2y, com2z, sep 
#	print "Finding CoMs"

	DefineScalarExpression("mass",  "<mesh_quality/volume>*rho")

	DefineScalarExpression("comx",  "mass*X" )
	DefineScalarExpression("comy",  "mass*Y" )
	DefineScalarExpression("comz",  "mass*Z" )

       	AddPlot("Pseudocolor", "comx")
	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.lbound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "rho"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m1xstr = string.split(GetQueryOutputString(), " ")
	m1x = m1xstr[4]
	m1x = float(m1x)


	DeleteAllPlots()
		

       	AddPlot("Pseudocolor", "comy")
	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.lbound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "rho"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m1ystr = string.split(GetQueryOutputString(), " ")
	m1y = m1ystr[4]
	m1y = float(m1y)

	DeleteAllPlots()


       	AddPlot("Pseudocolor", "comz")
	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.lbound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "rho"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m1zstr = string.split(GetQueryOutputString(), " ")
	m1z = m1zstr[4]
	m1z = float(m1z)

	DeleteAllPlots()


# Find CoM of the other star

       	AddPlot("Pseudocolor", "comx")
	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.ubound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "rho"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m2xstr = string.split(GetQueryOutputString(), " ")
	m2x = m2xstr[4]
	m2x = float(m2x)

	DeleteAllPlots()
		

       	AddPlot("Pseudocolor", "comy")
	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.ubound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "rho"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m2ystr = string.split(GetQueryOutputString(), " ")
	m2y = m2ystr[4]
	m2y = float(m2y)

	DeleteAllPlots()


       	AddPlot("Pseudocolor", "comz")
	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.ubound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "rho"
	iso5.lbound = rhoth
	SetOperatorOptions(iso5)

	DrawPlots()

	Query("Variable Sum")
	m2zstr = string.split(GetQueryOutputString(), " ")
	m2z = m2zstr[4]
	m2z = float(m2z)

	DeleteAllPlots()


	com1x = m1x/m1
	com1y = m1y/m1
	com1z =  m1z/m1

	com2x =  m2x/m2
	com2y =  m2y/m2
	com2z =  m2z/m2	
	
	sep = math.sqrt((m1x/m1-m2x/m2)**2+(m1y/m1-m2y/m2)**2+(m1z/m1-m2z/m2)**2) 

#	print com1x, com1y, com2x, com2y, div, rhoth, theta





def findSpins(rhothv):
	global com1x, com1y, com2x, com2y, div, rhoth, j1, j2, theta

#	print "Finding spins.."

	DefineScalarExpression("rotX", "X*cos(%f)-Y*sin(%f)" % (theta, theta))
	DefineScalarExpression("rotY", "X*sin(%f)+Y*cos(%f)" % (theta, theta))		


	DefineScalarExpression("j1F",  \
	"<mesh_quality/volume>* (  (X-(%f))*sy   -   (Y-(%f))*sx   )" \
	  % (com1x, com1y)) 


	DefineScalarExpression("j2F",  \
	"<mesh_quality/volume>* (  (X-(%f))*sy   -   (Y-(%f))*sx   )" \
	  % (com2x, com2y)) 



        AddPlot("Pseudocolor", "j1F")

	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.lbound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "rho"
	iso5.lbound = rhothv
	SetOperatorOptions(iso5)

	DrawPlots()
		
	Query("Variable Sum")
	j1str = string.split(GetQueryOutputString(), " ")
	j1 = j1str[4]
	j1 = float(j1)

	DeleteAllPlots()

#	print "j1",j1


        AddPlot("Pseudocolor", "j2F")

	AddOperator("Isovolume")
	iso4 = IsovolumeAttributes()
	iso4.variable = "rotX"
	iso4.ubound = div
	SetOperatorOptions(iso4)

	AddOperator("Isovolume")
	iso5 = IsovolumeAttributes()
	iso5.variable = "rho"
	iso5.lbound = rhothv
	SetOperatorOptions(iso5)

	DrawPlots()
		
	Query("Variable Sum")
	j2str = string.split(GetQueryOutputString(), " ")
	j2 = j2str[4]
	j2 = float(j2)

	DeleteAllPlots()

#	print "j2",j2
#	print com1x, com1y, com2x, com2y, div, rhoth, theta, rhothv




############### Program starts here ###############
try:
        for i in range(1000,1300,50):
                db = "X."+str(i)+".chk.silo"
                print db

# Open the file
                OpenDatabase(db)

# Coordinates
                getCoords()


# Find L1 point, angle of the binary, location of dividing plane
		findL1()

# Find individual masses and volumes
		getMasses()
#		print "m1", m1, "m2", m2
		

# Find total mass
                totalMass()
#                print "mt", mt

# Find total angular momentum
                totalAM()
#                print "jt", jt




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
		

#	print "s1234", s1, s2, s3,s4

        raise HaltException("Done!")


except HaltException as h:
        print(h)

############### Program ends here ###############
