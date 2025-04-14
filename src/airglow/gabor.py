from numpy import *
import numpy as np
import multiprocessing as mp

NUM_PROCESSORS = 4

n_stds = None
Bf     = None
Btheat = None
DC     = None
x0_max = None
y0_max = None

def Gabor_Kernel(n_stds, theta_k, Bf, Btheta, frequency, DC = False):

    u     = frequency

    lam = 1./u

    sigma_x   = lam*sqrt(2*log(2))/(2*pi*tanh(Bf/2.*log(2)))
    sigma_y   = lam*sqrt(2*log(2))/(2*pi)*1./tan(Btheta/2.)

    x0 = np.ceil(max(np.abs(n_stds * sigma_x * np.cos(theta_k)),
                     np.abs(n_stds * sigma_y * np.sin(theta_k)), 1))
    y0 = np.ceil(max(np.abs(n_stds * sigma_y * np.cos(theta_k)),
                     np.abs(n_stds * sigma_x * np.sin(theta_k)), 1))

    Y, X = np.mgrid[-y0:y0 + 1, -x0:x0 + 1]
    g     = zeros([X.shape[0], X.shape[1]], dtype = 'complex')

    g = 1./(2*pi*sigma_x*sigma_y)* \
        exp(-1/2.*(X*cos(theta_k) + Y*sin(theta_k))**2/sigma_x**2-1/2.*(-X*sin(theta_k) + Y*cos(theta_k))**2/sigma_y**2)* \
        exp(1j*2*pi*1./lam*(X*cos(theta_k) + Y*sin(theta_k)))

    DC_Term = 1./(2*pi*sigma_x*sigma_y)* \
            exp(-1/2.*(X**2/sigma_x**2 + Y**2/sigma_y**2) - 2*pi**2*sigma_x**2/lam**2)

    if DC == False:
        return g - DC_Term
    else:
        return g

def Gabor_Kernel_PLL(args):

#    print(".", end = '')

    global n_stds
    global Bf
    global Btheta
    global DC
    global x0_max
    global y0_max

    frequency, theta_k = args

    u     = frequency

    lam = 1./u

    sigma_x   = lam*sqrt(2*log(2))/(2*pi*tanh(Bf/2.*log(2)))
    sigma_y   = lam*sqrt(2*log(2))/(2*pi)*1./tan(Btheta/2.)

    x0 = np.ceil(max(np.abs(n_stds * sigma_x * np.cos(theta_k)),
                     np.abs(n_stds * sigma_y * np.sin(theta_k)), 1))
    y0 = np.ceil(max(np.abs(n_stds * sigma_y * np.cos(theta_k)),
                     np.abs(n_stds * sigma_x * np.sin(theta_k)), 1))

    # make all of the kernels the same size
    x0 = x0_max
    y0 = y0_max

    #print 'x0, y0:', x0, y0

    Y, X = np.mgrid[-y0:y0 + 1, -x0:x0 + 1]
    g     = zeros([X.shape[0], X.shape[1]], dtype = 'complex')

    g = 1./(2*pi*sigma_x*sigma_y)* \
        exp(-1/2.*(X*cos(theta_k) + Y*sin(theta_k))**2/sigma_x**2-1/2.*(-X*sin(theta_k) + Y*cos(theta_k))**2/sigma_y**2)* \
        exp(1j*2*pi*1./lam*(X*cos(theta_k) + Y*sin(theta_k)))

    DC_Term = 1./(2*pi*sigma_x*sigma_y)* \
            exp(-1/2.*(X**2/sigma_x**2 + Y**2/sigma_y**2) - 2*pi**2*sigma_x**2/lam**2)

    if DC == False:
        return g - DC_Term
    else:
        return g

def Gabor_Kernel_Sigma(n_stds, theta_k, sigma_x, sigma_y, frequency, DC = False):

    u     = frequency
    lam = 1./u

    x0 = np.ceil(max(np.abs(n_stds * sigma_x * np.cos(theta_k)),
                     np.abs(n_stds * sigma_y * np.sin(theta_k)), 1))
    y0 = np.ceil(max(np.abs(n_stds * sigma_y * np.cos(theta_k)),
                     np.abs(n_stds * sigma_x * np.sin(theta_k)), 1))

    Y, X = np.mgrid[-y0:y0 + 1, -x0:x0 + 1]

    g     = zeros([X.shape[0], X.shape[1]], dtype = 'complex')


    g = 1./(2*pi*sigma_x*sigma_y)* \
        exp(-1/2.*(X*cos(theta_k) + Y*sin(theta_k))**2/sigma_x**2-1/2.*(-X*sin(theta_k) + Y*cos(theta_k))**2/sigma_y**2)* \
        exp(1j*2*pi*1./lam*(X*cos(theta_k) + Y*sin(theta_k)))

    DC_Term = 1./(2*pi*sigma_x*sigma_y)* \
            exp(-1/2.*(X**2/sigma_x**2 + Y**2/sigma_y**2) - 2*pi**2*sigma_x**2/lam**2)

    if DC == False:
        return g - DC_Term
    else:
        return g

def Plot_Gabor_Kernel(g):
    pcolormesh(abs(g), cmap = 'inferno'); cb=colorbar();
#     xlim(np.min(X), np.max(X)); ylim(np.min(Y), np.max(Y));
    xlabel('$x$'); ylabel('$y$'); cb.set_label('$|g|$', fontsize=10)

def Plot_Gabor_Kernel_FFT(g):
    pcolormesh(abs(fftshift(fftpack.fft2(g))), cmap = 'inferno'); cb=colorbar();
#     xlim(np.min(U), np.max(V)); ylim(np.min(U), np.max(V));
    xlabel('$u$'); ylabel('$v$'); cb.set_label('$|G|$', fontsize=10)

def Gabor_Kernel_FT(n_stds, theta_k, Bf, Btheta, frequency, DC = False):

    u     = frequency
    lam = 1./u

    sigma_x   = lam*sqrt(2*log(2))/(2*pi*tanh(Bf/2.*log(2)))
    sigma_y   = lam*sqrt(2*log(2))/(2*pi)*1./tan(Btheta/2.)

    x0 = np.ceil(max(np.abs(n_stds * sigma_x * np.cos(theta_k)),
                     np.abs(n_stds * sigma_y * np.sin(theta_k)), 1))
    y0 = np.ceil(max(np.abs(n_stds * sigma_y * np.cos(theta_k)),
                     np.abs(n_stds * sigma_x * np.sin(theta_k)), 1))

    Y, X = np.mgrid[-y0:y0 + 1, -x0:x0 + 1]
    V, U = np.mgrid[-1/2.:1/2.:1j*X.shape[0], -1/2.:1/2.:1j*X.shape[1]]

    G     = zeros([X.shape[0], X.shape[1]], dtype = 'complex')

    G = 1*exp(- 2*pi**2*sigma_x**2*((U - 1./lam*cos(theta_k))*cos(theta_k) + (V - 1./lam*sin(theta_k))*sin(theta_k))**2  \
        -2*pi**2*sigma_y**2*(-(U - 1./lam*cos(theta_k))*sin(theta_k) + (V - 1./lam*sin(theta_k))*cos(theta_k))**2)

    DC_Term = exp(-2*pi**2*(sigma_x**2*(U**2 + 1./lam**2) + V**2*sigma_y**2))

    if DC == False:
        return G - DC_Term
    else:
        return G

def Gabor_Kernel_FT_Sigma(n_stds, theta_k, sigma_x, sigma_y, frequency, DC = False):

    u     = frequency
    lam = 1./u

    x0 = np.ceil(max(np.abs(n_stds * sigma_x * np.cos(theta_k)),
                     np.abs(n_stds * sigma_y * np.sin(theta_k)), 1))
    y0 = np.ceil(max(np.abs(n_stds * sigma_y * np.cos(theta_k)),
                     np.abs(n_stds * sigma_x * np.sin(theta_k)), 1))

    Y, X = np.mgrid[-y0:y0 + 1, -x0:x0 + 1]
    V, U = np.mgrid[-1/2.:1/2.:1j*X.shape[0], -1/2.:1/2.:1j*X.shape[1]]

    G     = zeros([X.shape[0], X.shape[1]], dtype = 'complex')

    G = 1*exp(- 2*pi**2*sigma_x**2*((U - 1./lam*cos(theta_k))*cos(theta_k) + (V - 1./lam*sin(theta_k))*sin(theta_k))**2  \
        -2*pi**2*sigma_y**2*(-(U - 1./lam*cos(theta_k))*sin(theta_k) + (V - 1./lam*sin(theta_k))*cos(theta_k))**2)

    DC_Term = exp(-2*pi**2*(sigma_x**2*(U**2 + 1./lam**2) + V**2*sigma_y**2))

    if DC == False:
        return G - DC_Term
    else:
        return G

def Plot_Gabor_Kernel_Frequency(G):
    fig = plt.figure(0)
    ax = fig.add_subplot(111, aspect='equal')
    pcolormesh(abs(G), cmap = 'inferno'); cb=colorbar();
#     xlim(np.min(U), np.max(V)); ylim(np.min(U), np.max(V));
    xlabel('$u$'); ylabel('$v$'); cb.set_label('$|G|$', fontsize=10)


def Create_FilterBank_Gabor(n_stdsarg, carriers, orientations, Bfarg, Bthetaarg, DCarg = False):

#    print("Constructing filter bank...", end = '')
    global n_stds
    global Bf
    global Btheta
    global DC
    global x0_max
    global y0_max

    n_stds       = n_stdsarg
    Bf           = Bfarg
    Btheta       = Bthetaarg
    DC           = DCarg

    g     = []

    u_min     = np.min(carriers)
    cost_max  = np.max(np.cos(orientations))
    sint_max  = np.max(np.sin(orientations))
    max_lam = 1./u_min

    # determine the largest kernel size:
    sigma_x_max   = max_lam*sqrt(2*log(2))/(2*pi*tanh(Bf/2.*log(2)))
    sigma_y_max   = max_lam*sqrt(2*log(2))/(2*pi)*1./tan(Btheta/2.)

    xa1 = np.abs(n_stds * sigma_x_max * np.cos(orientations))
    xa2 = np.abs(n_stds * sigma_y_max * np.sin(orientations))
    xa  = np.concatenate([xa1, xa2, [1]])

    ya1 = np.abs(n_stds * sigma_y_max * np.cos(orientations))
    ya2 = np.abs(n_stds * sigma_x_max * np.sin(orientations))
    ya  = np.concatenate([ya1, ya2, [1]])

    x0_max = np.ceil(np.max(xa))
    y0_max = np.ceil(np.max(ya))

    print('x0_max: %d' % x0_max)
    print('y0_max: %d' % y0_max)

    for l, u_l in enumerate(carriers):
        g.append([])
        for k, theta_k in enumerate(orientations):
            g[-1].append([])

    args = []
    for l, u_l in enumerate(carriers):
        for k, theta_k in enumerate(orientations):
            args.append((u_l, theta_k))

#    pool = mp.Pool(processes = NUM_PROCESSORS)
#    print(args)
#    g_LIST = pool.map(Gabor_Kernel_PLL, args)
#    pool.close()
#    pool.join()

    g_LIST = []
    for a in args:
        g_LIST.append(Gabor_Kernel_PLL(a))

    p = 0
    for l, u_l in enumerate(carriers):
        for k, theta_k in enumerate(orientations):
            g[l][k] = g_LIST[p]
            p = p + 1

    print("finished.")

    return g

def Create_FilterBank_Gabor_Sigma(n_stds, carriers, orientations, sigma_x, sigma_y, DC = False):

    g     = []

    for l, u_l in enumerate(carriers):
        g.append([])
        for k, theta_k in enumerate(orientations):
            g[-1].append([])

    for l, u_l in enumerate(carriers):
        for k, theta_k in enumerate(orientations):

            g_kl = Gabor_Kernel_Sigma(n_stds, theta_k, sigma_x, sigma_y, u_l, DC = DC)

            g[l][k] = g_kl

    return g

def Create_FilterBank_Gabor_Frequency(n_stds, carriers, orientations, Bf, Btheta, DC = False):

    G     = []

    for l, u_l in enumerate(carriers):
        G.append([])
        for k, theta_k in enumerate(orientations):
            G[-1].append([])

    for l, u_l in enumerate(carriers):
        for k, theta_k in enumerate(orientations):

            G_kl = Gabor_Kernel_FT(n_stds, theta_k, Bf, Btheta, u_l, DC = DC)

            G[l][k] = G_kl

    return G

def Create_FilterBank_Gabor_Frequency_Sigma(n_std, carriers, orientations, sigma_x, sigma_y, DC = False):

    G     = []

    for l, u_l in enumerate(carriers):
        G.append([])
        for k, theta_k in enumerate(orientations):
            G[-1].append([])

    for l, u_l in enumerate(carriers):
        for k, theta_k in enumerate(orientations):

            G_kl = Gabor_Kernel_FT_Sigma(n_stds, theta_k, Bf, Btheta, u_l, DC = DC)

            G[l][k] = G_kl

    return G

from scipy import signal

def Convolve_FilterBank(image, g):
    # if x has length N and h has length K, then x convolved with h has length N+K-1, thus the output image
    # will have size N+K-1 x N+K-1

    output     = []

    for l in range(0, len(g)):
        output.append([])
        for k in range(0, len(g[l][:])):
            output[-1].append([])

    for l in range(0, len(g)):
        print('%d/%d |,' % (l, len(g)))
        for k in range(0, len(g[l][:])):
            output[l][k] = signal.fftconvolve(image, g[l][k].real, mode='full')
#             output[l][k] = signal.convolve2d(image, g[l][k].real, mode='full', boundary='fill', fillvalue=0)

    return output

def Energies(convolved_images):

    energies     = []

    for l in range(0, len(convolved_images)):
        energies.append([])
        for k in range(0, len(convolved_images[l][:])):
            energies[-1].append([])

    for l in range(0, len(convolved_images)):
        for k in range(0, len(convolved_images[l][:])):
            energies[l][k] = sqrt(sum(abs(convolved_images[l][k])**2))

    return energies

def Energies_Direct(args):

    # this version convolves and calculates the energy at the same time -- then throws away the convolution output
    # this is more space efficient; the convolutions can take up a ton of space

    image, g = args

#    print(".", end = '')

    energies     = []

    for l in range(0, len(g)):
        energies.append([])
        for k in range(0, len(g[l][:])):
            energies[-1].append([])

    for l in range(0, len(g)):
        for k in range(0, len(g[l][:])):
            cnv = signal.fftconvolve(image, g[l][k].real, mode='full')
            N = cnv.shape[0]
            M = cnv.shape[1]
            energies[l][k] = sqrt(sum(abs(cnv)**2)/(M*N))

    return energies

data_frames = None
g_kernels   = None

def spatial_cnv(args):

    global data_frames
    global g_kernels

    imgi, gl, gk = args

    image = data_frames[imgi]
    glk   = g_kernels[gl][gk].real

    cnv    = signal.fftconvolve(image, glk, mode='full')

    N      = cnv.shape[0]
    M      = cnv.shape[1]

    energy = sum(abs(cnv)**2)/(M*N)

    return energy

def Energies_Direct_PLL(i):

    global data_frames
    global g_kernels

    energies  = []

    # this version convolves and calculates the energy at the same time -- then throws away the convolution output
    # this is more space efficient; the convolutions can take up a ton of space

    for l in range(0, len(g_kernels)):
        energies.append([])
        for k in range(0, len(g_kernels[l][:])):
            energies[-1].append([])

    args = []
    for l in range(0, len(g_kernels)):
        for k in range(0, len(g_kernels[l][:])):
            args.append((i,l,k))

#    pool = mp.Pool(processes = NUM_PROCESSORS)
#    energies_LIST = pool.map(spatial_cnv, args)
#    pool.close()
#    pool.join()

    energies_LIST = []
    for a in args:
        energies_LIST.append(spatial_cnv(a))

    p = 0
    for l in range(0, len(g_kernels)):
        for k in range(0, len(g_kernels[l][:])):
            energies[l][k] = energies_LIST[p]
            p = p + 1

    return energies

def Energies_Direct_Calibrated(image, edge_image, g):

    # this version convolves and calculates the energy at the same time -- then throws away the convolution output
    # this is more space efficient; the convolutions can take up a ton of space

    energies     = []

    for l in range(0, len(g)):
        energies.append([])
        for k in range(0, len(g[l][:])):
            energies[-1].append([])

    for l in range(0, len(g)):
        print('%d/%d' % (l, len(g)))
        for k in range(0, len(g[l][:])):
            cnv  = signal.fftconvolve(image, g[l][k].real, mode='full')
            edge = signal.fftconvolve(edge_image, g[l][k].real, mode='full')
            energies[l][k] = sqrt(sum(abs(cnv-edge)**2))

    return energies

import matplotlib.gridspec as gridspec
def multiresolution_plot(image, image_convolved, g, u, theta, Bf, Btheta):

    gs = gridspec.GridSpec(len(theta), len(u)*2, wspace = 0.0, hspace = 0.0)

    for k, theta_k in enumerate(theta):

        offset = 0

        for l, u_l in enumerate(u):

            lam       = 1./u_l
            sigma_x   = lam*sqrt(2*log(2))/(2*pi*tanh(Bf/2.*log(2)))
            sigma_y   = lam*sqrt(2*log(2))/(2*pi)*1./tan(Btheta/2.)
            hp_ellipse = Ellipse([1./lam*cos(theta_k),1./lam*sin(theta_k)], width=2/(2.*pi*sigma_x)*sqrt(2*log(2)),\
                         height=2/(2.*pi*sigma_y)*sqrt(2*log(2)), angle=degrees(theta_k))

            plt.subplot(gs[k, offset]); offset += 1
            imshow(abs(image_convolved[l][k]), cmap = 'inferno')
            gca().axes.xaxis.set_ticklabels([])
            gca().axes.yaxis.set_ticklabels([])
            plt.subplot(gs[k, offset]); offset +=1
            imshow(abs(fftshift(fftpack.fft2(g[l][k]))), cmap = 'inferno')
            grid(color='w')
            gca().axes.xaxis.set_ticklabels([])
            gca().axes.yaxis.set_ticklabels([])

    show()

def multiresolution_plot_spatialkernel(image, image_convolved, g, u, theta, Bf, Btheta):

    gs = gridspec.GridSpec(len(theta), len(u)*2, wspace = 0.0, hspace = 0.0)

    for k, theta_k in enumerate(theta):

        offset = 0

        for l, u_l in enumerate(u):

            lam       = 1./u_l
            sigma_x   = lam*sqrt(2*log(2))/(2*pi*tanh(Bf/2.*log(2)))
            sigma_y   = lam*sqrt(2*log(2))/(2*pi)*1./tan(Btheta/2.)
            hp_ellipse = Ellipse([1./lam*cos(theta_k),1./lam*sin(theta_k)], width=2/(2.*pi*sigma_x)*sqrt(2*log(2)),\
                         height=2/(2.*pi*sigma_y)*sqrt(2*log(2)), angle=degrees(theta_k))

            plt.subplot(gs[k, offset]); offset += 1
            imshow(abs(image_convolved[l][k]), cmap = 'inferno')
            gca().axes.xaxis.set_ticklabels([])
            gca().axes.yaxis.set_ticklabels([])
            plt.subplot(gs[k, offset]); offset +=1
            imshow(g[l][k].real, cmap = 'inferno')
            grid(color='w')
            gca().axes.xaxis.set_ticklabels([])
            gca().axes.yaxis.set_ticklabels([])

    show()

def filterbank_plot(g, u, theta, Bf, Btheta):

    gs = gridspec.GridSpec(len(theta), len(u), wspace = 0.0, hspace = 0.0)

    for k, theta_k in enumerate(theta):

        offset = 0

        for l, u_l in enumerate(u):

            plt.subplot(gs[k,l])
            imshow(g[l][k].real, cmap = 'inferno')
            grid(color='w')
            gca().axes.xaxis.set_ticklabels([])
            gca().axes.yaxis.set_ticklabels([])

    show()

def flower_plot(u, theta, Bf, Btheta):

    fig = plt.figure(0)
    ax = fig.add_subplot(111, aspect='equal')

    for k, theta_k in enumerate(theta):
        for l, u_l in enumerate(u):

            lam       = 1./u_l
            sigma_x   = lam*sqrt(2*log(2))/(2*pi*tanh(Bf/2.*log(2)))
            sigma_y   = lam*sqrt(2*log(2))/(2*pi)*1./tan(Btheta/2.)
            hp_ellipse = Ellipse([1./lam*cos(theta_k),1./lam*sin(theta_k)], width=2/(2.*pi*sigma_x)*sqrt(2*log(2)),\
                         height=2/(2.*pi*sigma_y)*sqrt(2*log(2)), angle=degrees(theta_k))

            ax.add_artist(hp_ellipse)
            hp_ellipse.set_clip_box(ax.bbox)
            hp_ellipse.set_alpha(1)
            hp_ellipse.set_facecolor('none')

    xlim(-1,1); ylim(-1,1);
    title('Filter Bank ($B_f = %0.2f$, $B_{\\theta} = %0.2f^{\circ}$)' % (Bf, degrees(Btheta)));
    xlabel('$u$'); ylabel('$v$');
    axvline(0.5, color='k', linestyle='dashed', linewidth=2)
    axvline(-0.5, color='k', linestyle='dashed', linewidth=2)
    axhline(0.5, color='k', linestyle='dashed', linewidth=2)
    axhline(-0.5, color='k', linestyle='dashed', linewidth=2)
    tight_layout();

def flower_plot_sigma(u, theta, sigma_x, sigma_y):

    fig = plt.figure(0)
    ax = fig.add_subplot(111, aspect='equal')

    for k, theta_k in enumerate(theta):
        for l, u_l in enumerate(u):

            lam       = 1./u_l
            hp_ellipse = Ellipse([1./lam*cos(theta_k),1./lam*sin(theta_k)], width=2/(2.*pi*sigma_x)*sqrt(2*log(2)),\
                         height=2/(2.*pi*sigma_y)*sqrt(2*log(2)), angle=degrees(theta_k))

            ax.add_artist(hp_ellipse)
            hp_ellipse.set_clip_box(ax.bbox)
            hp_ellipse.set_alpha(1)
            hp_ellipse.set_facecolor('none')

    xlim(-1,1); ylim(-1,1);
    title('Filter Bank ($B_f = %0.2f$, $B_{\\theta} = %0.2f^{\circ}$)' % (Bf, degrees(Btheta)));
    xlabel('$u$'); ylabel('$v$');
    axvline(0.5, color='k', linestyle='dashed', linewidth=2)
    axvline(-0.5, color='k', linestyle='dashed', linewidth=2)
    axhline(0.5, color='k', linestyle='dashed', linewidth=2)
    axhline(-0.5, color='k', linestyle='dashed', linewidth=2)
    tight_layout(); show()

def calibrate(convolved, convolved_edges):

    output = []

    for l in range(0, len(convolved)):
        output.append([])
        for k in range(0, len(convolved[l][:])):
            output[-1].append([])

    for l in range(0, len(convolved)):
        for k in range(0, len(convolved[l][:])):
            output[l][k] = convolved[l][k] - convolved_edges[l][k]

    return output

def calibrate_energies(energies, energies_edges):

    output = []

    for l in range(0, len(energies)):
        output.append([])
        for k in range(0, len(energies[l][:])):
            output[-1].append([])

    for l in range(0, len(energies)):
        for k in range(0, len(energies[l][:])):
            output[l][k] = energies[l][k] - energies_edges[l][k]

    return output
