#Define forsberg_g function for the gravitational attraction for a single prism and gravimeter location

from math import sqrt, log, atan

def pointmass_g(sysc, Hint, Hfin, posX, posY, posZ, xc, yc, delr, delc):
      ugal = 1e8
      gma = 6.673e-11
      rho = 1000.0
      x = posX - xc
      y = posY - yc
      z = posZ - (Hfin + Hint)/2
      dh = Hfin - Hint
      rf = sqrt((x**2 + y**2 + z**2)**3)
      grav = ugal*gma*rho*sysc*dh*delc*delr*(z/rf)
      return grav
      
def forsberg_g(sysc, drMid, dcMid, Hint, Hfin, posX, posY, posZ, dr, dc):
    # sysc: density change in g/cm^3
    # drMid: x-coordinate of cell midpoint in the along-row (x) direction
    # dcMid: y-coordinate of cell midpoint in the along-column (y) direction
    # Hint: z-coordinate of the bottom of the cell (initial head)
    # Hfin: z-coordinate of the top of the cell (final head)
    # posX, posY, posZ: x-, y-, and z-coordinates of the gravimeter position
    # dr: dimension of the cell in the along-row (x) direction
    # dc: dimension of the cell in the along-column (y) direction

    gamma = 6.67e-8
    rho = 1
    x = [(dcMid - dc / 2) - posX, (dcMid + dc / 2) - posX]
    y = [(drMid - dr / 2) - posY, (drMid + dr / 2) - posY]
    z = [Hfin - posZ, Hint - posZ]
    # print(dcMid)
    # print(drMid)
    # print(dc)
    # print(dr)
    prism_sum = 0
    corner_idx = [0, 1]
    for x_idx in corner_idx:
        for y_idx in corner_idx:
            for z_idx in corner_idx:
                rf = sqrt(x[x_idx] ** 2 + y[y_idx] ** 2 + z[z_idx] ** 2)
                prism_sum = prism_sum + (-1) ** (x_idx + y_idx + z_idx) * (x[x_idx] * log(y[y_idx] + rf)
                                                               + y[y_idx] * log(x[x_idx] + rf)
                                                               - z[z_idx] * atan(x[x_idx] * y[y_idx] / z[z_idx] / rf))

    forsberg = gamma * sysc * rho * 1e8 * prism_sum  # gravity in microGal

    return forsberg


def timestep_gravity(modelgrid, sy, z_initial, z_final, x, y, z):
    ctr = 0
    g = 0
    for m in range(modelgrid.nrow):
        for k in range(modelgrid.ncol):
            if z_initial[0][m][k] > -998 and abs(z_initial[0][m][k] - z_final[0][m][k]) > 0.01:  # NaNs are -999
                # in fhd file (indicates cells outside model boundary)

                r2 = (x-modelgrid.xcellcenters[m][k])**2 + \
                     (y-modelgrid.ycellcenters[m][k])**2 + \
                     (z_initial[0][m][k] - z_final[0][m][k])**2
                dr2 = modelgrid.delr[k]**2 +  modelgrid.delc[m]**2 + (z_initial[0][m][k]-z_final[0][m][k])**2
                if ( r2/dr2 > 81):
                    g += pointmass_g(sy[0][m][k],
                                z_initial[0][m][k],
                                z_final[0][m][k],
                                x, y, z,
                                modelgrid.xcellcenters[m][k],
                                modelgrid.ycellcenters[m][k],
                                modelgrid.delr[k], modelgrid.delc[m])
                else:
                    g += forsberg_g(sy[0][m][k],
                                modelgrid.ycellcenters[m][k],
                                modelgrid.xcellcenters[m][k],
                                z_initial[0][m][k],
                                z_final[0][m][k],
                                x, y, z,
                                modelgrid.delc[m], modelgrid.delr[k])
                ctr += 1
    return g
