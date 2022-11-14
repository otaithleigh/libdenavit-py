import openseespy.opensees as ops
import matplotlib.pyplot as plt
from math import ceil
from typing import Iterable, List, Literal, Tuple, Union


def uniaxial_material_analysis(
    definition: Iterable,
    peak_points,
    rate_type: Literal['None', 'StrainRate', 'Steps'] = 'None',
    rate_value: Union[float, int, None] = None,
    matTag: int = 1,
    plot_stress_strain: bool = False,
    compression_positive: bool = False,
) -> Tuple[List[float], List[float]]:
    """Analyze a uniaxial material.

    Parameters
    ----------
    definition : iterable
        Material definition as an iterable of arguments for
        `ops.uniaxialMaterial`.
    peak_points : array_like
        1-d array of imposed displacement targets.
    rate_type : {'None', 'StrainRate', 'Steps'}, optional
        How to interpolate between peaks in `peak_points`:
        - 'None' -- no interpolation is performed. (default)
        - 'StrainRate' -- interpolate between peaks at a fixed rate.
        - 'Steps' -- interpolate between peaks with a fixed number of steps.

    rate_value : float, optional
        Rate value, used when `rate_type` is 'StrainRate' or 'Steps'.
    matTag : int, optional
        Tag of the uniaxial material to analyze. (default: 1)
    plot_stress_strain : bool, optional
        If True, plot the recorded stress-strain response. (default: False)
    compression_positive : bool, optional
        If True, consider compression to be positive. (default: False)

    Returns
    -------
    (strain, stress) : (list[float], list[float])
        Recorded material response.
    """

    # ------------------------------
    # Build Model
    # ------------------------------
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)

    # Define Nodes
    ops.node(1,0.0)
    ops.node(2,1.0)

    # Define Boundary Conditions
    ops.fix(1, 1)

    # Define Elements
    ops.uniaxialMaterial(*definition)
    ops.element('Truss', 1, 1, 2, 1.0, matTag)

    # Define Loads
    ops.timeSeries('Linear', 854)
    ops.pattern('Plain', 979, 854)
    ops.load(2, 1.0)

    # ------------------------------
    # Build and Perform the Analysis
    # ------------------------------
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandSPD')
    ops.test('NormUnbalance', 1e-8, 10)
    ops.algorithm('Newton')
    ops.integrator('DisplacementControl', 2, 1, peak_points[0])
    ops.analysis('Static')
    ops.analyze(1)
    if compression_positive:
        stress = [-1*ops.getTime()]
        strain = [-1*ops.nodeDisp(2,1)]
    else:
        stress = [ops.getTime()]
        strain = [ops.nodeDisp(2,1)]

    for i in range(len(peak_points)-1):
        if rate_type == 'None':
            num_steps = 1
        elif rate_type == 'StrainRate':
            num_steps = ceil(abs(peak_points[i+1] - peak_points[i])/rate_value)
        elif rate_type == 'Steps':
            num_steps = rate_value
        else:
            raise Exception('Unknown rate_type: %s' % rate_type)

        incr = (peak_points[i+1] - peak_points[i])/num_steps

        for j in range(num_steps):
            ops.integrator('DisplacementControl', 2, 1, incr)
            ops.analyze(1)
            if compression_positive:
                stress.append(-1*ops.getTime())
                strain.append(-1*ops.nodeDisp(2,1))
            else:
                stress.append(ops.getTime())
                strain.append(ops.nodeDisp(2,1))

    if plot_stress_strain:
        fig = plt.figure()
        plt.plot(strain,stress)
        plt.xlabel('Strain')
        plt.ylabel('Stress')
        plt.show()

    return (strain,stress)
