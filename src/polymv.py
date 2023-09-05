import healpy as hp
import mpsolve as mp
import numpy as np

from functools import partial
from gmpy2 import bincoef
from gmpy2 import sqrt
from iteration_utilities import duplicates
from multiprocessing import Pool
from numba import jit

__author__ = (
    "Renan A. Oliveira, Thiago S. Pereira, Miguel Quartin"
)
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Renan A. Oliveira"
__email__ = "oliveirara@uel.br"

class mvs:
    @staticmethod
    def m_vectors(alms, l):
        """
        Get multipole vectors in spherical coordinates for a multipole given a map.

        Args:
            alms (complex array): multipole moments from a map (Healpy indexing).
            l (int): multipole.

        Returns:
            Float array [theta, phi] in radians.
        """

        # Context for MPSolve:
        ctx = mp.Context()
        poly = mp.MonomialPoly(ctx, 2 * l)

        # Get maximum multipole moment given an array.
        lmax = hp.sphtfunc.Alm.getlmax(len(alms))

        # Polynomial indexes:
        m_array = np.arange(-l, l + 1, 1)  # m values.
        neg_m = m_array[0:l]  # Negative m values.
        pos_m = m_array[l:]  # Positive m values.

        # Binomial terms:
        binom_term = [sqrt(bincoef(2 * l, l + i)) for i in range(-l, l + 1)]

        # Indexes for all alms list:
        index = hp.sphtfunc.Alm.getidx(lmax, l, abs(m_array))
        alm_ell_real = np.real(alms[index])
        alm_ell_imag = np.imag(alms[index])

        # Negative imaginary part for m:
        neg_rea_m_term = [(-1) ** int(i) for i in neg_m]
        # Positive imaginary part for m:
        neg_imag_m_term = [-((-1) ** int(i)) for i in neg_m]

        # Join all real signs:
        rea_sign = np.concatenate((neg_rea_m_term, np.ones(l + 1)))
        # Join all imaginary signs:
        imag_sign = np.concatenate((neg_imag_m_term, np.ones(l + 1)))

        rea = rea_sign * alm_ell_real * binom_term  # Real terms.
        ima = imag_sign * alm_ell_imag * binom_term  # Imaginary terms.

        # Set coefficients:
        [poly.set_coefficient(i, str(rea[i]), str(ima[i])) for i in range(2 * l + 1)]

        # Getting roots:
        ctx.solve(poly, algorithm=mp.Algorithm.SECULAR_GA)
        roots = np.array(ctx.get_roots())

        phi = np.angle(roots)
        r = np.absolute(roots)
        theta = 2 * np.arctan(1 / r)

        return np.vstack((theta, phi)).T

    @staticmethod
    def many_m_vectors(alm, multipoles, parallel=4):
        """
        Get several multipole vectors in spherical coordinates.

        Args:
            alms (complex array): multipolar coefficients from a map.
            multipoles (list, or range): list or range of multipoles.
            parallel (int): number of parallel process.

        Returns:
            Float array [multipole, [theta, phi]] in radians.
        """

        pool = Pool(processes=parallel)
        func = partial(mvs.m_vectors, alm)
        res = pool.map(func, multipoles)
        pool.close()
        pool.join()

        return np.vstack(res)


class otherfuncs:
    @staticmethod
    def get_multipole(alms, l):
        """
        Get specific multipolar coefficients.

        Args:
            alms (complex array): multipolar coefficients.
            l (int): multipole.

        Returns:
            Complex array.
        """
        ells, ms = hp.Alm.getlm(hp.Alm.getlmax(len(alms)))
        return alms[np.where(ells == l)]

    @staticmethod
    def to_cart(theta_phi_array):
        """
        Convert spherical to Cartesian coordinate. This function works only for a single multipole or stacked arrays, e.g. 2 and 3 MVs together.

        Args:
            theta_phi_array (float array): array containing in radians [theta, phi].

        Returns:
            Float array.
        """
        x = np.sin(theta_phi_array[:, 0]) * np.cos(theta_phi_array[:, 1])
        y = np.sin(theta_phi_array[:, 0]) * np.sin(theta_phi_array[:, 1])
        z = np.cos(theta_phi_array[:, 0])
        return np.vstack((x, y, z)).T

    @staticmethod
    def mvs_north(mvs):
        """
        Get multipole vectors on north hemisphere.

        Args:
            mvs (float array): array containing multipole vectors [theta, phi].

        Returns:
            Float array [theta, phi] in radians.
        """
        # First ordinary test:
        mvs_n = mvs[mvs[:, 0] < np.pi / 2]
        l = int(len(mvs) / 2)
        if len(mvs_n) == l:
            return mvs_n

        if len(mvs_n) == l + 1:

            # Find duplicates:
            duplicate_thetas = duplicates(mvs_n[:, 0])
            if len(duplicate_thetas) == 1:
                result = np.delete(
                    mvs_n,
                    np.random.choice(
                        np.where(mvs_n[:, 0] == duplicate_thetas)[0], size=1
                    )[0],
                    axis=0,
                )
                if len(result[:, 0]) == l:
                    return result

            # Delete highest value:
            else:
                result = np.delete(
                    mvs_n, np.where(mvs_n[:, 0] == np.max(mvs_n[:, 0]))[0], axis=0,
                )
                if len(result[:, 0]) == l:
                    return result

        if len(mvs_n) == l - 1:
            result = mvs[mvs[:, 0] >= np.pi / 2]
            sorted_result = result[result[:, 0].argsort()][0 : 1 + 1]

            duplicate_thetas = duplicates(sorted_result[:, 0])
            if len(duplicate_thetas) == 1:
                rand_theta = np.random.choice([0, 1])
                result1 = np.vstack((result, sorted_result[rand_theta]))

                if len(result1[:, 0]) == l:
                    return result1
            else:
                result2 = np.vstack((result, sorted_result[0]))
                if len(result2[:, 0]) == l:
                    return result2

        if (len(mvs_n) != l + 1) or (len(mvs_n) != l - 1):
            raise ValueError("Found a bug!")

    @staticmethod
    def many_mvs_north(mvs):
        """
        Get multipole vectors on north hemisphere for multiples multipoles.

        Args:
            mvs (float array): array containing multipole vectors [theta, phi].

        Returns:
            Float array [theta, phi] in radians.
        """
        mvs_north_full = []
        for i in range(len(mvs)):
            mvs_ell = mvs[i]
            mvs_north_full.append(otherfuncs.mvs_north(mvs_ell))
        return mvs_north_full


@jit(nopython=True)
def psi_first_estimation(mvs, theta_u, phi_u, x_u, y_u, z_u, type_="min"):
    """
    Fréchet Mean first estimator. Using Nside=64, this function gets the pixel which psi value is minimum/maximum.

    Args:
        mvs (float array): multipole vectors from an specific multipole.

    Returns:
        A vector.
    """
    f = np.sum(
        np.arccos(mvs[:, 0] * x_u[0] + mvs[:, 1] * y_u[0] + mvs[:, 2] * z_u[0]) ** 2
    )
    theta = theta_u[0]
    phi = phi_u[0]

    if type_ == "min":
        for i in range(1, len(theta_u)):
            frechet_mean = np.sum(
                np.arccos(mvs[:, 0] * x_u[i] + mvs[:, 1] * y_u[i] + mvs[:, 2] * z_u[i])
                ** 2
            )
            if frechet_mean < f:
                f = frechet_mean
                theta = theta_u[i]
                phi = phi_u[i]

        if theta > np.pi / 2:
            theta = np.pi - theta
            phi = np.pi + phi
            if phi > 2 * np.pi:
                phi = phi - 2 * np.pi
        return np.array([f, theta, phi])
    elif type_ == "max":
        for i in range(1, len(theta_u)):
            frechet_mean = np.sum(
                np.arccos(mvs[:, 0] * x_u[i] + mvs[:, 1] * y_u[i] + mvs[:, 2] * z_u[i])
                ** 2
            )
            if frechet_mean > f:
                f = frechet_mean
                theta = theta_u[i]
                phi = phi_u[i]
        if theta > np.pi / 2:
            theta = np.pi - theta
            phi = np.pi + phi
            if phi > 2 * np.pi:
                phi = phi - 2 * np.pi
        return np.array([f, theta, phi])
    else:
        return print("You must choose min or max.")


class fvs:
    # Index for first nside, which is 64.
    index_init = 6
    nside_first_estimator = 2 ** index_init

    # Number of pixels given nside.
    n_pixels = int(hp.nside2npix(nside_first_estimator) / 2)

    # Cartesian basis.
    x_u, y_u, z_u = hp.pix2vec(nside_first_estimator, np.arange(n_pixels))

    # Angles for each pixel.
    theta_u, phi_u = hp.pix2ang(nside_first_estimator, np.arange(n_pixels))

    @staticmethod
    def psi(mvs, type_="min"):
        """
        Fréchet Mean maximum or minimum.

        Args:
            mvs (float array): multipole vectors from an specific multipole.

        Returns:
            A vector.
        """
        # Vector guess:
        # Get the pixel center which is minimum given a map with Nside=64.
        mvs = otherfuncs.to_cart(mvs)
        first_estimation = psi_first_estimation(
            mvs,
            theta_u=fvs.theta_u,
            phi_u=fvs.phi_u,
            x_u=fvs.x_u,
            y_u=fvs.y_u,
            z_u=fvs.z_u,
            type_=type_,
        )[1:]
        x_vec = np.sin(first_estimation[0]) * np.cos(first_estimation[1])
        y_vec = np.sin(first_estimation[0]) * np.sin(first_estimation[1])
        z_vec = np.cos(first_estimation[0])

        # This part will get pixels inside the pixel guess that was calculated initially.
        pixels = hp.query_disc(
            nside=1024, vec=[x_vec, y_vec, z_vec], radius=np.deg2rad(2), inclusive=True
        )

        # We will minimize/maximize the Fréchet Mean with these pixels:
        x_u_pix, y_u_pix, z_u_pix = hp.pix2vec(1024, pixels)
        theta_u_pix, phi_u_pix = hp.pix2ang(1024, pixels)

        f = np.sum(
            np.arccos(
                mvs[:, 0] * x_u_pix[0] + mvs[:, 1] * y_u_pix[0] + mvs[:, 2] * z_u_pix[0]
            )
            ** 2
        )
        theta = theta_u_pix[0]
        phi = phi_u_pix[0]

        if type_ == "min":
            for i in range(1, len(x_u_pix)):
                frechet_mean_ = np.sum(
                    np.arccos(
                        mvs[:, 0] * x_u_pix[i]
                        + mvs[:, 1] * y_u_pix[i]
                        + mvs[:, 2] * z_u_pix[i]
                    )
                    ** 2
                )
                if frechet_mean_ < f:
                    f = frechet_mean_  # Fréchet-Mean Value
                    theta = theta_u_pix[i]
                    phi = phi_u_pix[i]
            if theta > np.pi / 2:
                theta = np.pi - theta
                phi = np.pi + phi
                if phi > 2 * np.pi:
                    phi = phi - 2 * np.pi
            return np.array([f, theta, phi])
        elif type_ == "max":
            for i in range(1, len(x_u_pix)):
                frechet_mean_ = np.sum(
                    np.arccos(
                        mvs[:, 0] * x_u_pix[i]
                        + mvs[:, 1] * y_u_pix[i]
                        + mvs[:, 2] * z_u_pix[i]
                    )
                    ** 2
                )
                if frechet_mean_ > f:
                    f = frechet_mean_  # Fréchet-Mean Value
                    theta = theta_u_pix[i]
                    phi = phi_u_pix[i]
            if theta > np.pi / 2:
                theta = np.pi - theta
                phi = np.pi + phi
                if phi > 2 * np.pi:
                    phi = phi - 2 * np.pi
            return np.array([f, theta, phi])
        else:
            return print("You must choose min or max.")
