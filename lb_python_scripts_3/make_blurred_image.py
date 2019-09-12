import Image, ImageColor, pyfits, os.path, types, pdb
from make_color import *
from numpy import *
from scipy import ndimage
import matplotlib.pyplot as plt
import scipy.signal as signal

def gauss_kern(size, sizey=None):
	""" Returns a normalized 2D gauss kernel array for convolutions """
	size = int(size)
	if not sizey:
		sizey = size
	else:
		sizey = int(sizey)
	#print size, sizey    
	x, y = mgrid[-size:size+1, -sizey:sizey+1]
	g = exp(-(x**2/float(size)+y**2/float(sizey)))
	return g / g.sum()

def blur_image(im, n, ny=None) :
	""" blurs the image by convolving with a gaussian kernel of FWHM
		n. The optional keyword argument ny allows for a different
		size in the y direction.
	"""
	#g = gauss_kern(n, sizey=ny)
	#improc = signal.convolve(im,g)
	improc = ndimage.filters.gaussian_filter(im,n/(2*sqrt(2*log(2))),mode='constant')
	return(improc)


def make_blurred_image(image,return_jpeg=True,**kwargs):

	"""
    Blurs image with a gaussian and converts a data array to a color
    image using the Lupton scheme.

	Almost identical to Patrik's make_color.

	A typical example (with beam FWHM 5 pix) would be:

	imshow(make_blurred_image.make_blurred_image(im, return_jpeg=False, band=(6,4,2), \
           scale="autolum",autopercentile=0.1,beamsize=5.))
	"""

	if kwargs.has_key("band"):
		band=kwargs["band"]
	else:
		band = (6,4,2)

	image = extract_bands(image, band)

	if kwargs.has_key("alpha"):
		alpha=float(kwargs["alpha"])
	else:
		alpha=.3

	if kwargs.has_key("Q"):
		Q=float(kwargs["Q"])
	else:
		Q=9

	if kwargs.has_key("autopercentile"):
		autopercentile=kwargs["autopercentile"]
	else:
		autopercentile=0.0005

	if kwargs.has_key("beamsize"):
		beamsize=float(kwargs["beamsize"])
	else:
		beamsize=1.

	print "Convolving image with Gaussian kernel of FWHM = "+str(beamsize)

	m=0

	sz = image.shape

	slices = sz[0]
	xs = sz [1]
	ys = sz [2]

	# deal with nan
	image = where(image != image, 0, image)

	if kwargs.has_key("scale"):
		if kwargs["scale"]=="auto":
			# set scale
			scale=find_autoscale(image, autopercentile, lum=False)
			print "autoscaling:"
			#raise RuntimeError, `scale `
		elif kwargs["scale"]=="autolum":
			scale=find_autoscale(image, autopercentile, lum=True)
			print "autoscaling:"
		elif isinstance(kwargs["scale"], str):
			sc = kwargs["scale"].split(',',2)
			scale=(float(sc[0]),float(sc[1]),float(sc[2]))
		else:
			scale=kwargs["scale"]
	else:
		scale=(1.,1.,1.)

	print "scale",scale

	image*=array(scale)

	r=image[:,:,0]
	g=image[:,:,1]
	b=image[:,:,2]

	r=blur_image(r,beamsize)
	g=blur_image(g,beamsize)
	b=blur_image(b,beamsize)

	i = (r+g+b)/3+1e-20
	r[:,:] = r*arcsinh (alpha*Q*(i-m))/(Q*i)
	g[:,:] = g*arcsinh (alpha*Q*(i-m))/(Q*i)
	b[:,:] = b*arcsinh (alpha*Q*(i-m))/(Q*i)

	image=transpose(array([r,g,b]),axes=(1,2,0))

	image=clip(image,0.0,1.0)

	if not return_jpeg:
		# imshow wants 0-1 range numbers, so that's what we return
		return image

	# the jpeg encoding wants 0-255 uchars
	image = image*255;
	jpgimg = Image.fromarray(image.astype(uint8))
	jpg = jpgimg.tostring("jpeg","RGB", 95)

	return jpg


