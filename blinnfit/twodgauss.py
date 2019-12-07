import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt

#define model function and pass independant variables x and y as a list
def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x=xy[0]
    y=xy[1]
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g=offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g.ravel()

def gen_twoD_Gaussian(xs,ys, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    # Create x and y indices
    x = np.linspace(0, xs-1, xs)
    y = np.linspace(0, ys-1, ys)
    x, y = np.meshgrid(x, y)
    data=twoD_Gaussian((x,y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset)
    return data.reshape((ys,xs))

def fit_twoD_Gaussian(img,g_amplitude,g_xo,g_yo,g_sigma_x,g_sigma_y,g_theta,g_offset):
    data=img.ravel()
    x = np.linspace(0, img.shape[1]-1,img.shape[1])
    y = np.linspace(0, img.shape[0]-1,img.shape[0])
    x, y = np.meshgrid(x, y)
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x,y), data,
                               p0=(g_amplitude,g_xo,g_yo,g_sigma_x,g_sigma_y,g_theta,g_offset))
    data_fitted = twoD_Gaussian((x, y), *popt).reshape(img.shape)
    return tuple(popt)+(data_fitted,)

if __name__=="__main__":
    #create data
    xs=201
    ys=301
    data = gen_twoD_Gaussian(xs,ys, 3, 100, 150, 20, 40, np.radians(15), 10)

    # plot twoD_Gaussian data generated above
    plt.figure()
    plt.imshow(data,origin='bottom')
    plt.colorbar()

    # add some noise to the data and try to fit the data generated beforehand

    data_noisy = data + 0.2*np.random.normal(size=data.shape)

    (amp,xo,yo,sigx,sigy,theta,ofs,data_fitted)=fit_twoD_Gaussian(data_noisy,3,100,100,20,40,0,10)

    fig, ax = plt.subplots(1, 1)
    ax.imshow(data_noisy, origin='bottom', extent=(0,xs-1,0,ys-1))
    x = np.linspace(0, xs - 1, xs)
    y = np.linspace(0, ys - 1, ys)
    x, y = np.meshgrid(x, y)
    ax.contour(x, y, data_fitted, 8, colors='w')
    plt.show()