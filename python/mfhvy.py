import numpy as np


class HeavyObs(object):
    def __init__(self, klay, zloc, xloc, yloc, hvylbl):
        self.klay = klay
        self.xloc = xloc
        self.yloc = yloc
        self.zloc = zloc
        self.hvylbl = hvylbl

    def __repr__(self):
        return "{:6d} {:12.2f} {:12.2f} {:12.2f} {:>12s}\n".format(self.klay,
                                                                   self.zloc,
                                                                   self.xloc,
                                                                   self.yloc,
                                                                   self.hvylbl)


class Heavy(object):
    """
    Python class for reading, creating, editing, and writing
    heavy gravity observation files

    Parameters
    ----------
    probs : list
        list of HeavyObs objects for the forsberg equation
    ifun : int
        heavy output file unit
    pmobs : list
        list of HeavyObs objects for the point mass equation
    ipmun : int
        heavy point mass output unit and flag
    hvynoh : float
        default value to write if cell goes dry
    sy : float
        optional array of specific yield or a scalar
    """
    def __init__(self, probs, ifun=61, pmobs=(), ipmun=0, hvynoh=-999.,
                 sy=None):

        tobs = []
        for ob in probs:
            if isinstance(ob, HeavyObs):
                tobs.append(ob)
            else:
                tobs.append(HeavyObs(*ob))

        self.obs = tobs
        self.ifun = ifun
        self.ipmun = ipmun
        self.hvynoh = hvynoh
        self.pmobs = pmobs
        self.sy = sy

        if pmobs is not None:
            tobs = []
            for ob in pmobs:
                if isinstance(ob, HeavyObs):
                    tobs.append(ob)
                else:
                    tobs.append(HeavyObs(*ob))

            self.pmobs = tobs

    def write(self, fname):
        """
        Method to write a Heavy output file

        Parameters
        ----------
        fname : str
            file name path

        """
        with open(fname, "w") as foo:
            foo.write("# File created with Heavy Python\n")
            foo.write("{}  {}  {}  {}  {:.2f}\n".format(len(self.obs),
                                                        self.ifun,
                                                        len(self.pmobs),
                                                        self.ipmun,
                                                        self.hvynoh))
            for ob in self.obs:
                foo.write(str(ob))

            for ob in self.pmobs:
                foo.write(str(ob))

        if self.sy is not None:
            if isinstance(self.sy, (float, int)):
                foo.write("CONSTANT  {}".format(self.sy))

            elif isinstance(self.sy, (list, np.ndarray)):
                sy = np.array(self.sy)
                if len(sy.shape) == 2:
                    foo.write("INTERNAL 1.0 (FREE) -1\n")
                    np.savetxt(foo, sy, fmt='%.4f', delimiter='  ')

                elif len(sy.shape) == 3:
                    for ix, arr in sy:
                        foo.write("INTERNAL 1.0 (FREE) -1  # layer {}\n"
                                  .format(ix + 1))
                        np.savetxt(foo, sy, fmt='%.4f', delimiter='  ')

                else:
                    raise AssertionError("SY must be 2 or 3 dimensional")

            else:
                raise TypeError("unrecognized type for SY. Must be "
                                "int, float, list, or np.array")

    @staticmethod
    def load(f):
        """
        Method to load an existing Heavy input file
        for editing or appending

        Parameters
        ----------
        f : str
            filename of the Heavy file

        Returns
        -------
            Heavy object

        """
        with open(f) as f:
            # remove comment line
            while True:
                line = f.readline()
                if line.startswith("#"):
                    continue
                elif line.startswith("!"):
                    continue
                else:
                    break

            # read dataset 1
            t = line.strip().split()
            nhvy = int(t[0])
            ihvyun = int(t[1])
            npmhvy = int(t[2])
            ipmun = int(t[3])
            hvynoh = float(t[4])

            obs = []
            for _ in range(nhvy):
                line = f.readline()
                t = line.strip().split()
                klay = int(t[0])
                zloc = float(t[1])
                xloc = float(t[2])
                yloc = float(t[3])
                hvylbl = t[4]
                obs.append(HeavyObs(klay, zloc, xloc, yloc, hvylbl))

            pmobs = []
            if ipmun > 0 and npmhvy > 0:
                for _ in range(npmhvy):
                    line = f.readline()
                    t = line.strip().split()
                    klay = int(t[0])
                    zloc = float(t[1])
                    xloc = float(t[2])
                    yloc = float(t[3])
                    hvylbl = t[4]
                    obs.append(HeavyObs(klay, zloc, xloc, yloc, hvylbl))

        # todo: update load for optional specific yield for SS models
        return Heavy(obs, ihvyun, pmobs, ipmun, hvynoh)


if __name__ == "__main__":
    import numpy as np

    n = 43218
    zlocs = np.ones((n,)) * 200
    xlocs = np.random.random((n,)) * 100
    ylocs = np.random.random((n,)) * 200
    label = ["GRAV_{}".format(i) for i in range(1, n+1)]

    observations = [HeavyObs(1, zlocs[i], xlocs[i], ylocs[i], label[i]) for i
                    in range(n)]
    hvy = Heavy(observations)
    hvy.write("test2.hvy")

    # hy = HeavyObs(1, 1.5, 2.0, 'test')
