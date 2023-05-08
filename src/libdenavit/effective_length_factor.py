from math import inf, tan, pi, nextafter
from scipy.optimize import brentq
from typing import Literal

__all__ = [
    'effective_length_factor',
    'sidesway_inhibited_effective_length_factor',
    'sidesway_uninhibited_effective_length_factor',
]


def effective_length_factor(sidesway: Literal['inhibited', 'uninhibited'],
                            GA: float,
                            GB: float) -> float:
    """Calculate the effective length factor K.

    Parameters
    ----------
    sidesway : {'inhibited', 'uninhibited'}
        Whether sidesway is inhibited (e.g., braced frames) or uninhibited
        (e.g., moment frames).
    GA, GB : float
        Ratio of column flexural stiffness to girder flexural stiffness at
        ends A and B of the member, respectively. (See AISC Specification
        equation C-A-7-3.)

    Notes
    -----
    K is calculated using the AISC nomogram equations, which rely on the
    following assumptions (quoting the AISC commentary):

    1. Behavior is purely elastic.
    2. All members have constant cross section.
    3. All joints are rigid.
    4. For columns in frames with sidesway inhibited, rotations at opposite
       ends of the restraining beams are equal in magnitude and opposite in
       direction, producing single curvature bending.
    5. For columns in frames with sidesway uninhibited, rotations at opposite
       ends of the restraining beams are equal in magnitude and direction,
       producing reverse curvature bending.
    6. The stiffness parameter L√(P/EI) of all columns is equal.
    7. Joint restraint is distributed to the column above and below the joint
       in proportion to EI/L for the two columns.
    8. All columns buckle simultaneously.
    9. No significant axial compression force exists in the girders.
    10. Shear deformations are neglected.
    """
    # case_inf:     default value when GA = GB = inf
    # case_0_inf:   default value when GA = 0 and GB = inf, or vice-versa
    # case_0:       default value when GA = GB = 0
    # bracket:      bracket of values of K across which sign change occurs
    #               in the nomogram equation for the given case. Since there's
    #               only one solution, just use the endpoints of the nomogram.
    # fcn_one_zero: simplified nomogram equation as a function of K and G for
    #               when one G is 0.
    # fcn:          full nomogram equation as a function of K, GA, and GB.
    if sidesway == 'inhibited':
        case_inf = 1.0
        case_0_inf = 0.7
        case_0 = 0.5
        # Use the value just below 1 to make solver happy; evaluating the
        # nomogram equation with K = 1.0 causes divide by zero
        bracket = (0.5, nextafter(1.0, -inf))
        fcn_with_one_zero = _sidesway_inhibited_one_zero
        fcn = _sidesway_inhibited
    elif sidesway == 'uninhibited':
        case_inf = inf
        case_0_inf = 2.0
        case_0 = 1.0
        # Technically [1, inf), but that makes things unhappy
        bracket = (1.0, 1e10)
        fcn_with_one_zero = _sidesway_uninhibited_one_zero
        fcn = _sidesway_uninhibited
    else:
        raise ValueError(f'Unrecognized sidesway value {sidesway!r}')

    # Check basic cases for a quick return
    if GA == inf and GB == inf:
        return case_inf
    elif (GA == inf and GB == 0) or (GA == 0 and GB == inf):
        return case_0_inf
    elif GA == 0 and GB == 0:
        return case_0

    # Change infinity to large value to avoid trouble in solver
    if GA == inf:
        GA = 1e10
    if GB == inf:
        GB = 1e10

    # Solve for the effective length factor
    #
    # Since we have a well-defined bracket across which K is defined
    # in the nomogram equations, use Brent's method which is fast,
    # stable, and guarantees convergence.
    if GA == 0:
        K = brentq(fcn_with_one_zero, *bracket, args=GB)
    elif GB == 0:
        K = brentq(fcn_with_one_zero, *bracket, args=GA)
    else:
        K = brentq(fcn, *bracket, args=(GA, GB))

    return K


def sidesway_uninhibited_effective_length_factor(GA: float, GB: float) -> float:
    """Calculate the effective length factor for sidesway uninhibited frames
    (e.g., moment frames).

    Parameters
    ----------
    GA, GB : float
        Ratio of column flexural stiffness to girder flexural stiffness at
        ends A and B of the member, respectively. (See AISC Specification
        equation C-A-7-3.)
    """
    return effective_length_factor('uninhibited', GA, GB)


def sidesway_inhibited_effective_length_factor(GA: float, GB: float) -> float:
    """Calculate the effective length factor for sidesway inhibited frames
    (e.g., braced frames).

    Parameters
    ----------
    GA, GB : float
        Ratio of column flexural stiffness to girder flexural stiffness at
        ends A and B of the member, respectively. (See AISC Specification
        equation C-A-7-3.)
    """
    return effective_length_factor('inhibited', GA, GB)


#----------------------------------
# Sidesway uninhibited functions
#----------------------------------
def _sidesway_uninhibited(K, GA, GB):
    # Equation C-A-7-2 of the Commentary on the 2016 AISC Specification
    pK = pi / K
    return (GA * GB * pK ** 2 - 36) / (6 * (GA + GB)) - pK / tan(pK)


def _sidesway_uninhibited_one_zero(K, G):
    pK = pi / K
    return -6 / G - pK / tan(pK)


#----------------------------------
# Sidesway inhibited functions
#----------------------------------
def _sidesway_inhibited(K, GA, GB):
    # Equation C-A-7-1 of the Commentary on the 2016 AISC Specification
    pK = pi / K
    return 0.25*GA*GB*pK**2 + 0.5*(GA + GB)*(1 - pK/tan(pK)) + 2*tan(0.5*pK)/pK - 1


def _sidesway_inhibited_one_zero(K, G):
    pK = pi / K
    return 0.5*G*(1 - pK/tan(pK)) + 2*tan(0.5*pK)/pK - 1


def run_examples():
    target_K = 3.0
    GA = inf
    GB = 6/(pi/target_K)/tan(pi/target_K)
    K = sidesway_uninhibited_effective_length_factor(GA, GB)
    print(f'Target K = {target_K}')
    print(f'Calculated K = {K}')


if __name__ == "__main__":
    run_examples()
