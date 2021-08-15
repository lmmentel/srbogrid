import numpy as np
from scipy.optimize import ridder


class SRBO:
    """
    Space-Reduced Bond-Order Grid generator for diatomics

    Based on:
        Rampino, S. (2016). Configuration-Space Sampling in Potential Energy
        Surface Fitting: A Space-Reduced Bond-Order Grid Approach.
        The Journal of Physical Chemistry A, 120(27), 4683–4692.
        http://doi.org/10.1021/acs.jpca.5b10018

    Args:
        Re : float
            Equlibrium bond length in bohr
        De : float
            Bond dissociation energy in hartree
        ke : float
            Force constant in hartree/bohr^2
        nrep : int
            Number of points in the repulsive part of the PES
        natt : int
            Number of points in the attractive part of the PES
        Vfact : float
            Parameter for getting `rmin`, Vfact = V(rmin)/De,
            default=1.5
        Vthrs : float
            Parameter for getting `rmax`, Vthrs = [De - V(rmax)]/De,
            default=0.001
        rmin : float
            Smallest bond length in bohr
        rmax : float
            Largest bond length in bohr

    .. note::

       `Vfact` and `Vthrs` are used to estimate the `rmin` and `rmax`
       parameters based on the Morse potential model, but `rmin` and
       `rmax` can be also passed explicitly and `Vfact`, `Vthrs` will
       be skipped
    """

    def __init__(
        self,
        Re: float,
        De: float,
        ke: float,
        nrep: int = 5,
        natt: int = 10,
        Vfact: float = 1.5,
        Vthrs: float = 0.001,
        rmin: float = None,
        rmax: float = None,
    ):
        self.Re = Re
        self.De = De
        self.ke = ke
        self.Vfact = Vfact
        self.Vthrs = Vthrs
        self.nrep = nrep
        self.natt = natt

        # calculate the remaining parameters
        self.npoints = self.nrep + self.natt + 1
        # Morse potential `alpha` parameter
        self.alpha = np.sqrt(self.ke / (2.0 * self.De))

        # `rmin` and `rmax` need to be initialized after `alpha`
        self.rmin = rmin
        self.rmax = rmax

        self.f = float(self.natt) / float(self.nrep)
        self.set_beta()

    @property
    def rmin(self):
        return self._rmin

    @rmin.setter
    def rmin(self, value):
        "Calculate the `rmin` value and update the value `Vfact`"

        if value is None:
            self._rmin = self.Re - np.log(1.0 + np.sqrt(self.Vfact)) / self.alpha
        else:
            self._rmin = value
            self.Vfact = self.morse(np.array([self._rmin]))[0] / self.De

    @property
    def rmax(self):
        return self._rmax

    @rmax.setter
    def rmax(self, value):
        "Calculate the `rmax` value and update the value `Vthrs`"
        if value is None:
            self._rmax = self.Re - np.log(1.0 - np.sqrt(1.0 - self.Vthrs)) / self.alpha
        else:
            self._rmax = value
            self.Vthrs = (self.De - self.morse(np.array([self._rmax]))[0]) / self.De

    @property
    def grid(self):
        "Alias for `get_bl_grid` method"
        return self.get_bl_grid()

    def set_beta(self):
        "Calculate the value of `beta` by finding a zero of `fbeta` function"

        def fbeta(beta, Re, rmin, rmax, f):

            nume = 1.0 - np.exp(-beta * (rmax - Re))
            deno = np.exp(-beta * (rmin - Re)) - 1.0
            return nume / deno - f

        args = tuple([self.Re, self.rmin, self.rmax, self.f])
        self.beta = ridder(fbeta, a=0.001, b=2.0, args=args)

    def evaluate(self, r):
        "Evaluate the BO coord. for a given BL coord. `r`"

        return np.exp(-self.beta * (r - self.Re))

    def evaluate_bl(self, n):
        "Evaluate the BL coord. for a gieven BO coord `n`"

        return self.Re - np.log(n) / self.beta

    def get_bo_grid(self):
        "Calculate the grid in BO coordinates"

        emin = np.exp(-self.beta * (self.rmin - self.Re))
        emax = np.exp(-self.beta * (self.rmax - self.Re))
        dn = (emin - emax) / (self.npoints - 1)

        return (emax + np.arange(self.npoints) * dn)[::-1]

    def get_bl_grid(self):
        "Calculate the grid in BL coordinates"

        return np.array([self.evaluate_bl(x) for x in self.get_bo_grid()])

    def summary(self):
        "Print internals, mainly for debugging"
        out = "System info:\n"
        for attr in ["Re", "De", "ke", "alpha"]:
            out += "\t{0:10s}: {1:10.6f}\n".format(attr, getattr(self, attr))

        out += "\nBoundaries:\n"
        for attr in ["rmin", "rmax", "Vfact", "Vthrs"]:
            out += "\t{0:10s}: {1:10.6f}\n".format(attr, getattr(self, attr))

        out += "\n\t{0:10s}: {1:10.6f}\n".format("Beta", self.beta)

        out += "\nGrid:\n"
        for attr in ["nrep", "natt", "npoints"]:
            out += "\t{0:10s}: {1:10d}\n".format(attr, getattr(self, attr))
        out += "\t{0:10s}: {1:10.6f}\n".format("f", getattr(self, "f"))
        out += "\nGrid points:\n"
        out += str(self.get_bl_grid())
        return out

    def morse(self, grid):
        """
        Morse potential on a discrete grid

        Args:
            grid : array_like
                Grid of points to evaluate the potential on
        """

        return self.De * (1.0 - np.exp(-self.alpha * (grid - self.Re))) ** 2.0

    def plot_morse(self, figsize=(16, 9), markersize: int = 50):
        "Plot the Morse potential with grid points indicated"

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        x = np.linspace(self.rmin, self.rmax, 100)
        ax.plot(x, self.morse(x), "k-")
        grid = self.get_bl_grid()

        ax.scatter(
            grid[self.nrep : self.nrep + 1],
            self.morse(grid[self.nrep : self.nrep + 1]),
            color="lightgreen",
            s=markersize,
        )
        ax.scatter(
            grid[: self.nrep], self.morse(grid[: self.nrep]), color="r", s=markersize
        )
        ax.scatter(
            grid[self.nrep + 1 :],
            self.morse(grid[self.nrep + 1 :]),
            color="b",
            s=markersize,
        )

        return fig, ax

    def save(self, filename):
        "Save the grid as numpy array"

        np.save(filename, self.get_bl_grid())
