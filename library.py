# ****************************************************************************
#       Copyright (C) 2022 Patricio Gallardo, Benjamin Schmidt
#       Contact: <pgallard@ucr.edu, schmbe@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ****************************************************************************

from sage.all import var, vector

x, y, z = var('x'), var('y'), var('z')

# This is the basis labelled e_0 to e_14 in the paper.
e = [vector([x**2 * y, -x**3, 0]),
     vector([x * y**2, -x**2 * y, 0]),
     vector([y**3, -x * y**2, 0]),
     vector([x**2 * z, 0, -x**3]),
     vector([x * z**2, 0, -x**2 * z]),
     vector([z**3, 0, -x * z**2]),
     vector([x * y * z, -x**2 * z, 0]),
     vector([x * y * z, 0, -x**2 * y]),
     vector([y**2 * z, -x * y * z, 0]),
     vector([y**2 * z, 0, -x * y**2]),
     vector([y * z**2, -x * z**2, 0]),
     vector([y * z**2, 0, -x * y * z]),
     vector([0, y**2 * z, -y**3]),
     vector([0, y * z**2, -y**2 * z]),
     vector([0, z**3, -y * z**2])]


def coeff_to_point(coeff):
    r"""Returns a vector in `math:H^0(\Omega(4))` defined through the
    coefficients `coeff` for our standard basis `math:e_1, \ldots, e_15`.
    """
    return sum(coeff[j] * e[j] for j in range(len(e)))


def group_action(coeff, matrix):
    r"""Computes the action of `matrix` on a point in the final model for
    `math: n = 7`. The point is determined through its coefficients in front
    of our standard basis (see paper) in the form of the list `coeff`.
    """
    # Check the remarks below Section 5.1 in the paper to understand the
    # group action.
    v = vector([x, y, z]) * matrix
    point = coeff_to_point(coeff)
    w = (matrix * point).subs(x=v[0], y=v[1], z=v[2])
    w = vector([entry.full_simplify().expand() for entry in w])

    # Determine the coefficients in front of the basis vectors after the
    # group action with matrix.
    return [w[0].coefficient(x, 2).coefficient(y, 1),  # e0
            w[0].coefficient(x, 1).coefficient(y, 2),  # e1
            w[0].coefficient(y, 3),  # e2
            w[0].coefficient(x, 2).coefficient(z, 1),  # e3
            w[0].coefficient(x, 1).coefficient(z, 2),  # e4
            w[0].coefficient(z, 3),  # e5
            -w[1].coefficient(x, 2).coefficient(z, 1),  # e6
            -w[2].coefficient(x, 2).coefficient(y, 1),  # e7
            -w[1].coefficient(x, 1).coefficient(y, 1).coefficient(z, 1),  # e8
            -w[2].coefficient(x, 1).coefficient(y, 2),  # e9
            -w[1].coefficient(x, 1).coefficient(z, 2),  # e10
            -w[2].coefficient(x, 1).coefficient(y, 1).coefficient(z, 1),  # e11
            w[1].coefficient(y, 2).coefficient(z, 1),  # e12
            w[1].coefficient(y, 1).coefficient(z, 2),  # e13
            w[1].coefficient(z, 3)]  # e14
