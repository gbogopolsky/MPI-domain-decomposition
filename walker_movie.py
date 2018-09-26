import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure()

data = np.loadtxt('position_movie.res')
print(data.shape)
domain_x = 20
domain_y = 50
nframes = 100

frames = np.zeros((nframes, domain_y, domain_x))
ims = []
for k in range(nframes):
	istart = k * domain_y
	iend = istart + domain_y
	frames[k,:,:] = data[istart:iend,:]
	im = plt.imshow(frames[k,:,:], animated=True)
	ims.append([im])

interval = 50
repeat_delay = interval * nframes
ani = animation.ArtistAnimation(fig, ims, interval=interval, blit=True, repeat_delay=repeat_delay)
plt.show()