import numpy as np
import utils

t = np.linspace(0,1,51)
vec = np.ones((51,5))
vec[:, 0] = t
vec[:, 1] = t*t
vec[:, 2] = np.exp(t)
vec[:, 3] = np.log(t+1)
vec[:, 4] = 1/(t+2)

dvec = utils.moving_derivative(vec, 4, 3, t)

vec_smooth = utils.savgol_filter_nonuniform(vec, 5, 3, t)

print(max(abs(dvec[:,0] - 1)))

print(max(abs(dvec[:,1] - 2*t)))

print(max(abs(dvec[:,2] - vec[:, 2])))

print(max(abs(dvec[:,3] - 1/(1+t))))


print(np.max(abs(vec_smooth - vec)))
