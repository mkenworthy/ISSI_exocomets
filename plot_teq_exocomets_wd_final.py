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



fig = plt.figure(figsize=(8,8))

axsun = fig.add_axes((left,0.57,wid,0.30))

axsun.set_xlim(Teq)
axsun.set_ylim(0,1)
axsun.set_xscale('log')
axsun.xaxis.set_visible(False)
axsun.yaxis.set_visible(False)

axsun_x2 = axsun.secondary_xaxis(1.2)
axsun_x2.set_xlabel(r'$\mathrm{T_{eq}\ [K]}$')

axsun_x3 = axsun.secondary_xaxis('top', functions=(DeqsolRstar, DeqsolRstar))
#axsun_x3.set_xlabel(r'$D\ [R_\odot]$')

axsun_x4 = axsun.secondary_xaxis('bottom', functions=(Deqsolau, Teqsolau))

axsun.text(0.05, 0.95, "Sun", transform=axsun.transAxes,
            fontsize=16, fontweight='bold', va='top')


solplanets='''
Mercury	0.4 blue	5
Venus	0.7	blue	10
Earth	1	blue	10
Mars	1.5	red	5
Jupiter	5.2	orange	40
Saturn	9.5	orange	20
Uranus	19	blue 	15
Neptune	30	blue	15
'''
for l in ascii.read(solplanets):
	(name, dist,col, siz) = l
	axsun.text(Teqsolau(dist),0.5+0.04,name[0:2],
		rotation=90,verticalalignment='bottom',horizontalalignment='center',
		bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.7))
	axsun.scatter(Teqsolau(dist),0.5,c=col,s=siz)

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
		rotation=0,verticalalignment='bottom',horizontalalignment='right',
		fontsize=8,bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.7))

solpoi='''
Halley perihelion,0.6
Halley aphelion,35
Rosetta operations,1.2
Most distant activity seen,25
Hale Bopp Na tail,1
'''
for l in ascii.read(solpoi):
	(name, dist) = l
	axsun.text(Teqsolau(dist),0.77,name,
		rotation=-35,verticalalignment='bottom',horizontalalignment='right',
		fontsize=5,
		bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.7))
	axsun.scatter(Teqsolau(dist),0.75,marker=marker1,linewidth=0,color='blue')

# put in the solar surface
#axsun.text(8000,0.8,"Sun",fontsize=20)

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

axbpic = fig.add_axes((left,0.34,wid,0.145))

axbpic.set_xlim(Teq)
axbpic.set_ylim(0,1)

axbpic.set_xscale('log')
axbpic.xaxis.set_visible(False)
axbpic.yaxis.set_visible(False)

axbpic_x2 = axbpic.secondary_xaxis('top', functions=(DeqbpicRstar, DeqbpicRstar))
#axbpic_x2.set_xlabel(r'$D\ [R_{*}]$')

axbpic_x3 = axbpic.secondary_xaxis('bottom', functions=(Deqbpicau, Teqbpicau))
#axbpic_x3.set_xlabel(r'$D\ [au]$')

axbpic.text(0.05, 0.95, r"$\beta$ Pic", transform=axbpic.transAxes,
            fontsize=16, fontweight='bold', va='top')

bpicplanets='''
b	10	orange	50
c	2.68	orange	70
'''
for l in ascii.read(bpicplanets):
	(name, dist,col,sci) = l
	axbpic.text(Teqbpicau(dist),0.5+0.07,name[0:2],
		rotation=90,verticalalignment='bottom',horizontalalignment='center' ,
		bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.7))
	axbpic.scatter(Teqbpicau(dist),0.5,c=col,s=sci)


bpiccomets='''
inner belt,6.2,6.4	
main disc,16,1500
exocomets,0.0248,0.215
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
		rotation=0,verticalalignment='bottom',horizontalalignment='right',fontsize=8 ,
		bbox=dict(facecolor='white', edgecolor='none', pad=1.0,alpha=0.7))

# WD
# WD
# WD

def TeqwdRstar(D):
	T=15500.
	R=1.
	albedo=0.0
	return T*np.power(1-albedo,0.25)*np.sqrt(R/(2*D))

def DeqwdRstar(Teq):
	T=15500.
	R=1.
	albedo=0.0
	alb = np.power(1-albedo,0.25)
	denom = 2*np.power((Teq/T)*(1./alb),2)
	return R/denom

def Teqwdau(D):
	T=15500.
	R=(0.0112*u.Rsun/(1*u.au)).decompose()
	albedo=0.0
	return T*np.power(1-albedo,0.25)*np.sqrt(R/(2*D))

def Deqwdau(Teq):
	T=15500.
	R=(0.0112*u.Rsun/(1*u.au)).decompose()
	albedo=0.0
	alb = np.power(1-albedo,0.25)
	denom = 2*np.power((Teq/T)*(1./alb),2)
	return R/denom



axwd = fig.add_axes((left,0.08,wid,0.145))

axwd.set_xlim(Teq)
axwd.set_ylim(0,1)
axwd.set_xscale('log')
axwd.xaxis.set_visible(False)
axwd.yaxis.set_visible(False)


axwd_x2 = axwd.secondary_xaxis('top', functions=(DeqwdRstar, DeqwdRstar))
#axwd_x2.set_xlabel(r'$D\ [R_{WD}]$')


axwd_x3 = axwd.secondary_xaxis('bottom', functions=(Deqwdau, Teqwdau))
#axwd_x3.set_xlabel(r'$D\ [au]$')

#axwd.text(9000,0.8,r"WD 1145+017",fontsize=20)

axwd.text(0.05, 0.95, r"WD 1145+017", transform=axwd.transAxes,
            fontsize=16, fontweight='bold', va='top')

# putting on D labels on all three

fs = 14
for (ax,r) in zip((axsun,axbpic,axwd),(r'$\mathrm{R_\odot}$',r'$\mathrm{R_*}$',r'$\mathrm{R_{WD}}$')):
#	ax.text(-0.01, -0.02, r"$au$", transform=ax.transAxes,
#            fontsize=fs, fontweight='bold', va='top',ha='right')

	ax.text(1.012, -0.02, r"$\mathrm{au}$", transform=ax.transAxes,
            fontsize=fs, fontweight='bold', va='top',ha='left')

#	ax.text(-0.01, 1.02, r, transform=ax.transAxes,
#            fontsize=fs, fontweight='bold', va='bottom',ha='right')

	ax.text(1.012, 1.02, r, transform=ax.transAxes,
            fontsize=fs, fontweight='bold', va='bottom',ha='left')




wdper='''
A, 4.49888
B, 4.60530
C, 4.78283
D, 4.55000
E, 4.82336
F, 4.85848
'''

@u.quantity_input
def Ptoa(P:u.year, m1:u.M_sun, m2:u.M_jup)->u.au:
    """calculate orbital radius from period

    Args:
        P: orbital period
        m1, m2: Primary and secondary masses

    Returns:
        a: semi-major axis

    >>> import astropy.units as u
    >>> Ptoa(11.86*u.year, 1.0*u.M_sun, 1.0*u.M_jup)
    <Quantity 5.20222482 AU>
    """

    # a^3/P^2 = (G/4pipi) (m1 + m2)
    const = c.G / (4.*np.pi*np.pi)
    mu = m1 + m2
    a3 = P*P*const*mu
    aa = np.power(a3, 1./3.)

    return aa

axzoom = fig.add_axes((0.4,0.175,0.2,0.037))

axzoom.set_xlim(1050,1015)
axzoom.set_ylim(0,1)

axzoom.set_xscale('log')
axzoom.xaxis.set_visible(False)
axzoom.yaxis.set_visible(False)
axzoom.set_ylim(0,1)

for pp in ascii.read(wdper):
	(name, per) = pp
	dist = Ptoa(per*u.hour, 0.707*u.Msun, 0.01*u.Mjup)
	#axwd.text(Teqwdau(dist.value),0.5+0.03,name[0:1],
	#	rotation=90,verticalalignment='bottom',horizontalalignment='center')
	axwd.scatter(Teqwdau(dist.value),0.5,s=20,marker=marker1,linewidth=0,color='blue')


	axzoom.text(Teqwdau(dist.value),0.37+0.05,name[0:1],
		rotation=0,verticalalignment='bottom',horizontalalignment='center')
	axzoom.scatter(Teqwdau(dist.value),0.36,s=50,marker=marker1,linewidth=0,color='blue')

# zoom lines
axwd.plot((950,270),(0.50,0.65),color='black',linewidth=1)
axwd.plot((1100,1121),(0.50,0.65),color='black',linewidth=1)

sublim='''
$\mathrm{H_2O}$	180
$\mathrm{CO}$	25
$\mathrm{CO_2}$	80
silicates	1400
$\mathrm{HCN}$	95
'''

toff = 0.90
for l in ascii.read(sublim):
	(name, dist) = l
	for ax in (axsun,axbpic,axwd):
		ax.text(dist+1.5,toff,name,
			rotation=90,verticalalignment='top',
			horizontalalignment='right',alpha=0.5,zorder=-10)
		ax.vlines(dist,0,1,colors='black',alpha=0.5,zorder=-10)

wdcomets='''
debris disk,1570,300
'''
comets=(0.1,0.4) # comet patches are put between these two limits on the y-axis
bpiccom = ascii.read(wdcomets,format='fast_no_header')
print(bpiccom)
deltay = np.abs(comets[1]-comets[0]) / len(bpiccom)
hei = deltay*0.9

# calculate lower left corner for each Rect
for i, (l) in enumerate(bpiccom):
	(name, d1, d2) = l
	d1x = d1 # the table values are in Teq so no conversion necessary
	d2x = d2 # the table values are in Teq so no conversion necessary

	lower = comets[0] + i * deltay
	wi = d2x-d1x

	# draw the rectangle
	t = Rectangle((d1x,lower), wi, hei)
	axwd.add_patch(t)

	axwd.text(d1x,lower+hei,name+'  ',
		rotation=0,verticalalignment='bottom',horizontalalignment='right',fontsize=8)


plt.draw()
plt.savefig("teq_exocomets_wd_final.pdf")

plt.show()
