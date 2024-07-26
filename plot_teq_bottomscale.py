import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as c
from astropy import units as u
from astropy.io import ascii 

from matplotlib.patches import Rectangle

# these are only used for fancy comet marker! replace it later if you want.
from svgpathtools import svg2paths
from svgpath2mpl import parse_path


### https://stackoverflow.com/questions/14324270/matplotlib-custom-marker-symbol
def generate_marker_from_svg(svg_path):
	import matplotlib as mpl
	image_path, attributes = svg2paths(svg_path)

	image_marker = parse_path(attributes[0]['d'])

	image_marker.vertices -= image_marker.vertices.mean(axis=0)

	image_marker = image_marker.transformed(mpl.transforms.Affine2D().rotate_deg(180))
	image_marker = image_marker.transformed(mpl.transforms.Affine2D().scale(-1,1))

	return image_marker

marker1 = generate_marker_from_svg("comet2.svg")



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

fig = plt.figure(figsize=(6,8))

# are the two functions actual inverses of each other?
#print(Deq(Teq(np.array([1,2,5,10,100]))))

Teq = (1e4,30) # temperature range for all the plots!

# width of the temperature scales 
left = 0.1
right = 0.9
wid = right-left

# # the equilibruim temperature scale
# axteq = fig.add_axes((left,0.9,wid,0.05))

# axteq.tick_params(left=False,right=False)
# axteq.tick_params('x', top=True, labeltop=True)
# axteq.set(yticklabels=[])
# axteq.set_xscale('log')
# axteq.set_xticks([1e4,1e3,1e2,1e1])
# axteq.set_xlim(Teq)

# Each star

# def plot_star(ax,Teq,Rstar=1.0*u.Rsun,Tstar=5700*u.K):
# 	'given teq ranges, plot out the upper and lower scales for Rstar and au'

axsun = fig.add_axes((left,0.57,wid,0.30))

#axsun.set_xlim(Deqsol(np.array(Teq)))
axsun.set_xlim(Teq)
axsun.set_ylim(0,1)
axsun.set_xscale('log')
axsun.set_xlabel(r'$T\ [^oK]$')

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
	axsun.scatter(Teqsolau(dist),0.75,marker=marker1)



# put in the solar surface
axsun.text(8000,0.8,"Sun",fontsize=20)

#\Pisymbol{astrosym}{134}


t2 = Rectangle((Teq[0],0.), TeqsolRstar(1.0)-Teq[0], 1.,facecolor='yellow')
axsun.add_patch(t2)


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

axbpic = fig.add_axes((left,0.08,wid,0.3))

axbpic.set_xlim(Teq)
axbpic.set_ylim(0,1)

axbpic.set_xscale('log')
axbpic.set_xlabel(r'$T\ [^oK]$')

axbpic.tick_params(left=False,right=False)
axbpic.set(yticklabels=[])

axbpic_x = axbpic.secondary_xaxis(1.2, functions=(DeqbpicRstar, DeqbpicRstar))
axbpic_x.set_xlabel(r'$D\ [R_*]$')

axbpic_x2 = axbpic.secondary_xaxis('top', functions=(Deqbpicau, Teqbpicau))
axbpic_x2.set_xlabel(r'$D\ [au]$')

axbpic.text(5700,0.8,r"$\beta$ Pic",fontsize=20)

t2 = Rectangle((Teq[0],0.), TeqbpicRstar(1.0)-Teq[0], 1.,facecolor='yellow')
axbpic.add_patch(t2)



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
plt.show()
