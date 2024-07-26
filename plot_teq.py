import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as c
from astropy import units as u
from astropy.io import ascii 

from matplotlib.patches import Rectangle

# these are only used for fancy comet marker! replace it later if you want.
# conda install svgpathtools
# conda install svgpath2mpl

from svgpath2mpl import parse_path

# here is the comet marker

### https://stackoverflow.com/questions/14324270/matplotlib-custom-marker-symbol
# def generate_marker_from_svg(svg_path):
# 	from svgpathtools import svg2paths

# 	import matplotlib as mpl
# 	image_path, attributes = svg2paths(svg_path)

# 	image_marker = parse_path(attributes[0]['d'])

# 	image_marker.vertices -= image_marker.vertices.mean(axis=0)

# 	image_marker = image_marker.transformed(mpl.transforms.Affine2D().rotate_deg(180))
# 	image_marker = image_marker.transformed(mpl.transforms.Affine2D().scale(-1,1))

# 	return image_marker

#marker1 = generate_marker_from_svg("comet2.svg")



def generate_comet_marker():
	import matplotlib as mpl

	# the big string is the comet marker but you still need the svgpath2mpl module

	image_marker = parse_path("M22.707,10.658a1,1,0,0,1,0,1.414l-8.7,8.7A7.622,7.622,0,0,1,3.232,9.989L9.611,3.611a1,1,0,1,1,1.414,1.414L4.646,11.4A5.622,5.622,0,0,0,12.6,19.354l8.7-8.7A1,1,0,0,1,22.707,10.658ZM11.335,12.664a3.838,3.838,0,1,1-5.427,0A3.842,3.842,0,0,1,11.335,12.664ZM9.921,14.078a1.838,1.838,0,1,0-2.6,2.6,1.88,1.88,0,0,0,2.6,0A1.84,1.84,0,0,0,9.921,14.078ZM18.293,1.293l-7.7,7.7A1,1,0,1,0,12.011,10.4l7.7-7.7a1,1,0,0,0-1.414-1.414Zm-2,9.414a1,1,0,0,0,1.414,0l4-4a1,1,0,1,0-1.414-1.414l-4,4A1,1,0,0,0,16.293,10.707ZM4,5A1,1,0,1,0,3,6,1,1,0,0,0,4,5ZM21,19.5A1.5,1.5,0,1,0,19.5,21,1.5,1.5,0,0,0,21,19.5Z")

	image_marker.vertices -= image_marker.vertices.mean(axis=0)

	image_marker = image_marker.transformed(mpl.transforms.Affine2D().rotate_deg(180))
	image_marker = image_marker.transformed(mpl.transforms.Affine2D().scale(-1,1))

	return image_marker

marker1 = generate_comet_marker()

# we want to use LaTeX on the figure

plt.rcParams.update({
    "text.usetex": True,
 #   "font.family": "Helvetica",
    "font.size": 12
})


# def Teq_units(D,T=5700.*u.K,R=1*u.Rsun,albedo=0.5):
# 	return (T*np.power(1-albedo,0.25)*np.sqrt(R/(2*D))).to(u.K)

# def Deq_units(Teq,T=5700.*u.K,R=1*u.Rsun,albedo=0.5):
# 	alb = np.power(1-albedo,0.25)
# 	denom = 2*np.power((Teq/T)*(1./alb),2)
# 	return (R/denom).to(u.au)


# def Teq(D,T=5700.,R=1.,albedo=0.5):
# 	return T*np.power(1-albedo,0.25)*np.sqrt(R/(2*D))

# def Deq(Teq,T=5700.,R=1.,albedo=0.5):
# 	alb = np.power(1-albedo,0.25)
# 	denom = 2*np.power((Teq/T)*(1./alb),2)
# 	return R/denom


# ░█▀▀░█░█░█▀█
# ░▀▀█░█░█░█░█
# ░▀▀▀░▀▀▀░▀░▀


# in order to plot the different Axes scales on the panels, you need to define
# two functions (a forward and inverse function) so that matplotlib can make
# the scales.
#
# ALSO: you can't pass more than one variable in these functions, so you need a duplicate
# of the function for EACH different scale. Grah!

def TeqsolRstar(D):
	T=5700.
	R=1.
	albedo=0.0
	return T*np.power(1-albedo,0.25)*np.sqrt(R/(2*D))

def DeqsolRstar(Teq):
	T=5700.
	R=1.
	albedo=0.0
	alb = np.power(1-albedo,0.25)
	denom = 2*np.power((Teq/T)*(1./alb),2)
	return R/denom

def Teqsolau(D):
	T=5700.
	R=(1*u.Rsun/(1*u.au)).decompose()
	albedo=0.0
	return T*np.power(1-albedo,0.25)*np.sqrt(R/(2*D))

def Deqsolau(Teq):
	T=5700.
	R=(1*u.Rsun/(1*u.au)).decompose()
	albedo=0.0
	alb = np.power(1-albedo,0.25)
	denom = 2*np.power((Teq/T)*(1./alb),2)
	return R/denom


# GENERAL PARAMETERS FOR THE PLOT
# GENERAL PARAMETERS FOR THE PLOT
# GENERAL PARAMETERS FOR THE PLOT

Teq = (1e4,30) # temperature range for all the plots!

# width of the temperature scales 
left = 0.1
right = 0.9
wid = right-left



fig = plt.figure(figsize=(6,8))

axsun = fig.add_axes((left,0.57,wid,0.30))

#axsun.set_xlim(Deqsol(np.array(Teq)))
axsun.set_xlim(Teq)
axsun.set_ylim(0,1)
axsun.set_xscale('log')
#axsun.set_xlabel(r'$T\ [^oK]$')
axsun.text(50,-0.08,r'$T\ [K]$',fontsize=14)

axsun.tick_params(left=False,right=False)
#axteq.tick_params('x', top=True, labeltop=True)
axsun.set(yticklabels=[])

axsun_x = axsun.secondary_xaxis(1.2, functions=(DeqsolRstar, DeqsolRstar))
axsun_x.set_xlabel(r'$D\ [R_\odot]$')

axsun_x2 = axsun.secondary_xaxis('top', functions=(Deqsolau, Teqsolau))
axsun_x2.set_xlabel(r'$D\ [au]$')



solplanets='''
Mercury	0.4
Venus	0.7
Earth	1
Mars	1.5
Jupiter	5.2
Saturn	9.5
Uranus	19
Neptune	30
'''
for l in ascii.read(solplanets):
	(name, dist) = l
	axsun.text(Teqsolau(dist),0.5+0.03,name[0:2],
		rotation=90,verticalalignment='bottom',horizontalalignment='center')
	axsun.scatter(Teqsolau(dist),0.5)

solcomets='''
Kreutz sungrazers,0.018568,0.13926
Jupiter Family Comets,0.3,12
Main Asteroid Belt,2.1,3.3
Kuiper Belt,30,50
'''
comets=(0.1,0.4) # comet patches are put between these two limits on the y-axis
solcom = ascii.read(solcomets)
deltay = np.abs(comets[1]-comets[0]) / len(solcom)
hei = deltay*0.9

# calculate lower left corner for each Rect
for i, (l) in enumerate(solcom):
	(name, d1, d2) = l
	d1x = Teqsolau(d1)
	d2x = Teqsolau(d2)

	lower = comets[0] + i * deltay
	wi = d2x-d1x

	# draw the rectangle
	t = Rectangle((d1x,lower), wi.value, hei)
	axsun.add_patch(t)

	axsun.text(d1x,lower+hei,name+'  ',
		rotation=0,verticalalignment='bottom',horizontalalignment='right',fontsize=8)

solpoi='''
Halley perihelion,0.6
Halley aphelion,35
Rosetta operations,1.2
Most distant activity observed,25
Hale Bopp Na tail,1
'''
for l in ascii.read(solpoi):
	(name, dist) = l
	axsun.text(Teqsolau(dist),0.77,name,
		rotation=-35,verticalalignment='bottom',horizontalalignment='right',
		fontsize=5)
	axsun.scatter(Teqsolau(dist),0.75,marker=marker1,linewidth=0,color='blue')

# put in the solar surface
axsun.text(8000,0.8,"Sun",fontsize=20)

# add star surface (CAUTION this doesn't work because we use Teq as the scale and 
# scale to that, but this formula doesn't work for close to the star.)
# t2 = Rectangle((Teq[0],0.), 5700-Teq[0], 1.,facecolor='yellow')
# axsun.add_patch(t2)


# ░█▀▄░█▀▀░▀█▀░█▀█░░░█▀█░▀█▀░█▀▀░▀█▀░█▀█░█▀▄░▀█▀░█▀▀
# ░█▀▄░█▀▀░░█░░█▀█░░░█▀▀░░█░░█░░░░█░░█░█░█▀▄░░█░░▀▀█
# ░▀▀░░▀▀▀░░▀░░▀░▀░░░▀░░░▀▀▀░▀▀▀░░▀░░▀▀▀░▀░▀░▀▀▀░▀▀▀

def TeqbpicRstar(D):
	T=8052.
	R=1.
	albedo=0.0
	return T*np.power(1-albedo,0.25)*np.sqrt(R/(2*D))

def DeqbpicRstar(Teq):
	T=8052.
	R=1.
	albedo=0.0
	alb = np.power(1-albedo,0.25)
	denom = 2*np.power((Teq/T)*(1./alb),2)
	return R/denom

def Teqbpicau(D):
	T=8052.
	R=(1.8*u.Rsun/(1*u.au)).decompose()
	albedo=0.0
	return T*np.power(1-albedo,0.25)*np.sqrt(R/(2*D))

def Deqbpicau(Teq):
	T=8052.
	R=(1.8*u.Rsun/(1*u.au)).decompose()
	albedo=0.0
	alb = np.power(1-albedo,0.25)
	denom = 2*np.power((Teq/T)*(1./alb),2)
	return R/denom

axbpic = fig.add_axes((left,0.18,wid,0.345))

axbpic.set_xlim(Teq)
axbpic.set_ylim(0,1)
axbpic.xaxis.tick_top()

axbpic.set_xscale('log')
axbpic.set_xlabel('')
#axbpic.tick_params(top=False)
axbpic.axes.xaxis.set_ticklabels([])

axbpic.tick_params(left=False,right=False)
axbpic.set(yticklabels=[])

axbpic_x = axbpic.secondary_xaxis(-0.2, functions=(DeqbpicRstar, DeqbpicRstar))
axbpic_x.set_xlabel(r'$D\ [R_*]$')

axbpic_x2 = axbpic.secondary_xaxis('bottom', functions=(Deqbpicau, Teqbpicau))
axbpic_x2.set_xlabel(r'$D\ [au]$')

axbpic.text(5700,0.8,r"$\beta$ Pic",fontsize=20)

bpicplanets='''
b	10
c	2.68
'''
for l in ascii.read(bpicplanets):
	(name, dist) = l
	axbpic.text(Teqbpicau(dist),0.5+0.03,name[0:2],
		rotation=90,verticalalignment='bottom',horizontalalignment='center')
	axbpic.scatter(Teqbpicau(dist),0.5)


bpiccomets='''
inner belt,6.2,6.4	
main disc,16,1500
'''
comets=(0.1,0.4) # comet patches are put between these two limits on the y-axis
bpiccom = ascii.read(bpiccomets,format='fast_no_header')
print(bpiccom)
deltay = np.abs(comets[1]-comets[0]) / len(bpiccom)
hei = deltay*0.9

# calculate lower left corner for each Rect
for i, (l) in enumerate(bpiccom):
	(name, d1, d2) = l
	d1x = Teqbpicau(d1)
	d2x = Teqbpicau(d2)

	lower = comets[0] + i * deltay
	wi = d2x-d1x

	# draw the rectangle
	t = Rectangle((d1x,lower), wi.value, hei)
	axbpic.add_patch(t)

	axbpic.text(d1x,lower+hei,name+'  ',
		rotation=0,verticalalignment='bottom',horizontalalignment='right',fontsize=8)




sublim='''
$H_2O$	180
$CO$	25
$CO_2$	80
silicates	1400
$HCN$	95
'''
for l in ascii.read(sublim):
	(name, dist) = l
	axsun.text(dist,0.85,name,
		rotation=90,verticalalignment='top',
		horizontalalignment='right',alpha=0.5,zorder=-10)
#	axsun.scatter(dist,0.8,color='black',alpha=0.5,zorder=-10)
	axsun.vlines(dist,0,1,colors='black',alpha=0.5,zorder=-10)

	axbpic.text(dist,0.85,name,
		rotation=90,verticalalignment='top',
		horizontalalignment='right',alpha=0.5,zorder=-10)
#	axsun.scatter(dist,0.8,color='black',alpha=0.5,zorder=-10)
	axbpic.vlines(dist,0,1,colors='black',alpha=0.5,zorder=-10)

 


plt.draw()
plt.savefig("teq_exocomets.pdf")
plt.show()
